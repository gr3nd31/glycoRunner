library(ggpubr)
library(tidyverse)

glycoRun <- function(
  # Whole protein gel profile ladder file
  ladderWP = "ladder_wp.csv",
  # Glycosylation protein gel profile ladder file
  ladderGP = "ladder_gp.csv",
  # Whole protein profile ladder values
  ladderValuesWP = c(250, 150, 100, 75, 50, 37, 25, 20, 15, 10),
  # Glycolsyation ladder values
  ladderValuesGP = c(75, 25),
  # Whole protein gel profile file
  wp = "wp.csv",
  # Glycosylation gel profile file
  gp = "gp.csv",
  # Name of the background correction csv for coomassie
  wp_background_file = "wp_background.csv",
  # Name of the background correction csv for glyco
  gp_background_file = "gp_background.csv",
  # Output file name
  outputCSV = "allData.csv",
  # Location of the peak data
  peakCSV = "csv/allData.csv",
  # Numerical threshold for peak calling
  givenThreshold = F,
  # The number of iterations attempted to get the ladder right before it gives up
  ladderAttempts = 10,
  # If the ladder attempts max out, should the peaks called be assumed as the upper or lower portion of the given ladder sizes?
  upperOrLowerLadder = "upper") {
  
  print("Running GlycoRunner...")
  cat("\n")
  #Lists all the files in the main folder
  currentSetup <- list.files()
  # If 'figures', 'bands', or 'csv' aren't present, these directories get made
  if(!"figures" %in% currentSetup){
    print("Setting up directories...")
    dir.create("figures")
    pngList <- list.files(pattern = ".png")
    for (a in pngList){
      file.rename(a, paste0("figures/", a))
    }
  }
  # if(!"bands" %in% currentSetup){
  #   dir.create("bands")
  # }
  if(!"csv" %in% currentSetup){
    dir.create("csv")
  }
  
  # If the corrections option is TRUE and the ladder file is present, the data gets corrected
  correctionsCheck <- list.files()
  # First the ladder gets corrected
  if (ladderWP %in% correctionsCheck){
    # The csv's are listed
    cList <- list.files(pattern = ".csv")
    # The csvs are read into a dataframe
    pxWP <- read.csv(ladderWP)
    wp <- read.csv(wp)
    # The dataframe is given a name and type on ladder
    # For non-ladder samples, these are taken from the file name: 'name'_'type'.csv
    # Dataframe names are changed
    names(pxWP)[1] <- "position"
    names(pxWP)[2] <- "Ladder_value"
    
    names(wp)[1] <- "position"
    names(wp)[2] <- "WP_value"
    
    if (wp_background_file %in% correctionsCheck){
      print("Applying background corrections to whole protein files...")
      background_data <- read.csv(wp_background_file)
      pxWP$Ladder_value <- pxWP$Ladder_value - background_data$Gray_Value
      wp$WP_value <- wp$WP_value - background_data$Gray_Value
    } else{
      print("Background file not found. Try again")
    }
    write.csv(wp, "csv/wp_corrected.csv", row.names = F)
    # The corrected WP points are then bound to the WP ladder dataframe
    if (nrow(wp) != nrow(pxWP)){
      print("The WP ladder is not the same length as the WP sample file. Script will now self-destruct...")
    }
    pxWP$WP_value <- wp$WP_value

    # Now we do the same for the GP ladder and sample
    gp <- read.csv(gp)
    names(gp)[1] <- "position"
    names(gp)[2] <- "GP_value"
    if (gp_background_file %in% correctionsCheck){
      print("Applying background corrections to glycoprotein files...")
      background_data <- read.csv(gp_background_file)
      gp$GP_value <- gp$GP_value - background_data$Gray_Value
    } else {
      print("Background file not found. Try again although it may not be necessary")
    }
    
    # If the ladders are different files, then the ladders and samples will be corrected
    if (ladderGP != ladderWP){
      # The GP ladder is read into a CSV
      pxGP <- read.csv(ladderGP)
      # Dataframe names are changed
      names(pxGP)[1] <- "position"
      names(pxGP)[2] <- "Ladder_value"
      if (gp_background_file %in% correctionsCheck){
        print("Applying background corrections to glycoprotein ladder...")
        background_data <- read.csv(gp_background_file)
        pxGP$Ladder_value <- pxGP$Ladder_value - background_data$Gray_Value
      } else{
        print("Background file not found. Try again although it may not be necessary")
      }
      print("Background correction is complete.")
      cat("\n")

      # The corrected GP values are saved
      write.csv(gp, "csv/gp_corrected.csv", row.names = F)
      if (nrow(gp) != nrow(pxGP)){
        print("The GP ladder is not the same length as the GP sample file. Script will now self-destruct...")
      }
      # The corrected GP values are saved to the GP ladder
      pxGP$GP_value <- gp$GP_value
      
      print("Attempting to equalize datasets...")

      #Now we have to correct everything
      # First, lets assign ladder peaks to their corresponding values
      print("Calling peaks on the whole protein ladder...")
      pxWP <- peakCaller(df = pxWP,
                 value = "Ladder_value",
                 threshold = F,
                 expectedPeaks = T,
                 expectedNumber = length(ladderValuesWP),
                 ladder = ladderValuesWP,
                 attemptsForPeaks = ladderAttempts, 
                 ladderUpOrDown = upperOrLowerLadder)
      print("Whole protein ladder peaks called...")
      cat("\n")
      # Now that ladder peaks are assigned, we'll fill out the rest of the positions by linear regression between peaks
      # First the pxWP values
      ladderValuesWP <- unique(pxWP[pxWP$size != 0,]$size)
      for (peakz in 1:length(ladderValuesWP)){
        x1 <- pxWP[pxWP$size == ladderValuesWP[peakz],]$position[1]
        y1 <- ladderValuesWP[peakz]
        
        if(peakz == 1){
          x2 <- pxWP[pxWP$size == ladderValuesWP[peakz+1],]$position[1]
          y2 <- ladderValuesWP[peakz+1]
          peakSlope <- (y2-y1)/(x2-x1)
          peakIntercept <- y1-x1*peakSlope
          pxWP[pxWP$size == 0 & pxWP$position < x2,]$size <- (pxWP[pxWP$size == 0 & pxWP$position < x2,]$position*peakSlope)+peakIntercept
        } else if (peakz == length(ladderValuesWP)){
          x2 <- max(pxWP$position)
          y2 <- 0.1
          peakSlope <- (y2-y1)/(x2-x1)
          peakIntercept <- y1-x1*peakSlope
          pxWP[pxWP$size == 0 & pxWP$position > x1,]$size <- (pxWP[pxWP$size == 0 & pxWP$position > x1,]$position*peakSlope)+peakIntercept
        } else {
          x2 <- pxWP[pxWP$size == ladderValuesWP[peakz+1],]$position[1]
          y2 <- ladderValuesWP[peakz+1]
          peakSlope <- (y2-y1)/(x2-x1)
          peakIntercept <- y1-x1*peakSlope
          pxWP[pxWP$size == 0 & pxWP$position > x1 & pxWP$position < x2,]$size <- (pxWP[pxWP$size == 0 & pxWP$position > x1 & pxWP$position < x2,]$position*peakSlope)+peakIntercept
        }
      }
      
      # Now to assign the ladder peaks to the pxGP dataset
      cat("\n")
      print("Calling peaks on the glycoprotein ladder...")
      pxGP <- peakCaller(df = pxGP,
                 value = "Ladder_value",
                 threshold = F,
                 expectedPeaks = T,
                 expectedNumber = length(ladderValuesGP),
                 ladder = ladderValuesGP,
                 attemptsForPeaks = ladderAttempts, 
                 ladderUpOrDown = upperOrLowerLadder)

      print("Glyco ladder peaks called...")
      # Now we adjust the positions so that the GP values and WP values match based on the ladder peaks
      # First, we'll align the GP position based on the ladder peaks 
      # Since not all ladder marks may have been found, we'll restrict our ladder bands to those that were assigned
      ladderValuesGP <- unique(pxGP[pxGP$size != 0,]$size)
      # Now we assign a placeholder value newPosition
      pxGP$newPosition <- 0.1
      # Then we give the max GP position a size equal to the smallest WP pixel size
      pxGP[pxGP$position == max(pxGP$position),]$size <- min(pxWP$size)
      # And we'll give the minimal GP position the maximal WP pixel size
      pxGP[pxGP$position == min(pxGP$position),]$size <- max(pxWP$size)
      # These will act as our upper and lower bounds or 'anchors' (see below)
      # While we could just use these anchors to bin or expand the GP data to fit WP positions, its better to do this on a stack-by-stack basis where
      # stacks are the pixels between known boundaries. We've got two boundaries, but the GP ladder peaks can provide us with other boundaries
      # For each band in the GP ladder we'll assign is a newPosition based on the position of the pixel in the pxWP with the closest size
      for (gpPeakNumber in length(ladderValuesGP)){
        gpPeak <- ladderValuesGP[gpPeakNumber]
        # Then we find the position of the WP pixel that is closest in size to the GP ladder peak and assign the GP ladder peak that position as a newPosition
        pxGP[pxGP$size == gpPeak,][1,]$newPosition <- pxWP[abs(pxWP$size-gpPeak) == min(abs(pxWP$size-gpPeak)),]$position[1]
      }
      # Now that we have some newposition anchors, we'll use math to calculate the newPosition of the interpeak pixels
      # first, we'll make an array of the anchor points
      anchors <- append(c(max(pxWP$size)), ladderValuesGP)
      anchors <- append(anchors, c(min(pxWP$size)))
      
      pxWP$GP_value <- 0
      runRegression <- F
      
      # Now we iterate through the stacks
      for (anchorNum in 2:length(anchors)){
        # First we generate wp and gp lists of the pixels that fit within the stack bounded by our known sizes
        # The WP ipList is easy, since sizes have already been determined for every position
        ipList_wp <- pxWP[pxWP$size <= anchors[anchorNum-1] & pxWP$size >= anchors[anchorNum],]
        # But since only certain positions have assigned sizes to the GP values, we take pixels with positions bound by the known points
        # So first we get the maximal position bound by the bigger size defining the stack
        upperGP <- pxGP[pxGP$size == anchors[anchorNum-1],]$position[1]
        # Then we get the minimal position bound by the smaller size defining the stack
        lowerGP <- pxGP[pxGP$size == anchors[anchorNum],]$position[1]
        # Then we just pull the GP pixels that have original position within the stack
        ipList_gp <- pxGP[pxGP$position >= upperGP & pxGP$position <= lowerGP,]
        # Now, there are three possibilities:
        # 1) The ipList_wp and the ipList_gp are the same length, wherein we can ignore increasing or decreasing pixels and just assign the GP_values directly
        if (nrow(ipList_gp) == nrow(ipList_wp)){
          pxWP[pxWP$position %in% ipList_wp$position,]$GP_value <- ipList_gp$value
          
          
        } else if (nrow(ipList_gp) < nrow(ipList_wp)){
        # 2) The ipList_wp could be larger than the correspond ipList_gp.
          # Here, we'll have to bin the wp pixels into a number of bins defined by the gp pixels
          gpBin <- split(ipList_gp, cut(ipList_gp$position, nrow(ipList_wp)))
          # Then we assign the average GP value of each bin to the corresponding wp pixel
          for (ip_wp in 1:nrow(ipList_wp)){
            pxWP[pxWP$position == ipList_wp[ip_wp,]$position,]$GP_value <- mean(gpBin[[ip_wp]]$GP_value)
          }
          
          # Finally, we'll make a linear equation of the GP_Values to fill in the gaps
          emptyNest <- pxWP[is.na(pxWP$GP_value),]
          for(nest in 1:nrow(emptyNest)){
            nestPosition <- emptyNest[nest,]$position
            upperValue <- pxWP[pxWP$position == nestPosition-1,]$GP_value
            lowerValue <- pxWP[pxWP$position == nestPosition+1,]$GP_value
            if(is.na(upperValue) | is.na(lowerValue)){
              print("Legitimate neighbors not for the missing pixel")
              runRegression <- T
            }else{
              pxWP[pxWP$position == emptyNest[nest,]$position,]$GP_value <- mean(upperValue, lowerValue)
            }
          }
        } else if (nrow(ipList_gp) > nrow(ipList_wp)){
        # 3) The ipList_wp is smaller than the ipList_gp.
          # In this final case, we'll bin the gp pixels in to a number of bins defined by the wp pixels in the size stack
          gpBin <- split(ipList_gp, cut(ipList_gp$position, nrow(ipList_wp)))
          # Then we assign the average GP value of each bin to the corresponding wp pixel
          for (ip_wp in 1:nrow(ipList_wp)){
            pxWP[pxWP$position == ipList_wp[ip_wp,]$position,]$GP_value <- mean(gpBin[[ip_wp]]$GP_value)
          }
        }
      }
      # So, there's a chance that, in binning the GP values across a larger WP pixel set, there will be gaps larger that one pixel
      #In these cases, NA's will still be present in the dataset, which will trigger the runRegression parameter
      # Here, we identify runs of NAs in the GP values, find the non-NA neighbors, generate a linear equation, and reassign pxWP GP_values based on the equation
      if (runRegression == T){
        # we collect all the pixels in pxWP that still have NA values in GP_value
        emptyNest <- pxWP[is.na(pxWP$GP_value),]
        # We iterate through them
        for (nest in 1:nrow(emptyNest)){
          # if a interim dataset 'stickList' doesn't exist, we make it exist
          if (!exists("stickList")){
            # Then we assign the pixel to the list
            stickList <- emptyNest[nest,]
            # we'll also grab the position to compare because complicated code hides mistakes
            stick <- emptyNest[nest,]$position
          } else {
            # if the next NA pixel is right next to the previous pixel, we add it to the list
            if (emptyNest[nest,]$position == stick+1){
              stickList <- rbind(stickList, emptyNest[nest,])
              # And we increment the positional comparison by one
              stick <- stick+1
            } else {
              # if the current pixel is *not* continuous with the current list of NAs, we calculate our new values and reset the process
              # First we'll grab the neighboring positions that aren't NA
              x1 <- min(stickList$position)-1
              x2 <- max(stickList$position)+1
              # Then we get their non-NA GP_value's
              y1 <- pxWP[pxWP$position == x1,]$GP_value
              y2 <- pxWP[pxWP$position == x2,]$GP_value
              # Then we calculate the slope and intercept
              stickSlope <- (y2-y1)/(x2-x1)
              stickIntercept <- y1-(stickSlope*x1)
              # Then we assign the calcuated values based on said slope and intercept
              pxWP[pxWP$position %in% stickList$position,]$GP_value <- stickList$position*stickSlope+stickIntercept
              
              # Now that the current set is done, we reset the NA set to the current pixel
              stickList <- emptyNest[nest,]
              stick <- emptyNest[nest,]$position
              }
            }
          }
        }
      # Finally, we can assign everything to px
      px <- pxWP
      
    } else {
      # if the same ladder is used for both GP and WP, then everything is simply bound together
      # Then the GP gel is read, corrected, and bound to the WP ladder
      if (nrow(gp) != nrow(pxWP)){
        print("The WP ladder and WP sample file is not the same length as the GP sample file. Script may self-destruct...")
      }
      pxWP$GP_value <- gp$GP_value
      #Assign the variable to px
      px <- pxWP

      #Run the peak-assigning script
      px <- peakCaller(df = px, 
                value = "Ladder_value",
                threshold = F,
                expectedPeaks = T,
                expectedNumber = length(ladderValuesWP),
                ladder = ladderValuesWP,
                attemptsForPeaks = ladderAttempts, 
                ladderUpOrDown = upperOrLowerLadder)
      
      ladderValuesWP <- unique(px[px$size != 0,]$size)
      for (peakz in 1:length(ladderValuesWP)){
        x1 <- px[px$size == ladderValuesWP[peakz],]$position[1]
        y1 <- ladderValuesWP[peakz]
        
       if(peakz == 1){
         x2 <- px[px$size == ladderValuesWP[peakz+1],]$position[1]
         y2 <- ladderValuesWP[peakz+1]
         peakSlope <- (y2-y1)/(x2-x1)
         peakIntercept <- y1-x1*peakSlope
         px[px$size == 0 & px$position < x2,]$size <- (px[px$size == 0 & px$position < x2,]$position*peakSlope)+peakIntercept
       } else if (peakz == length(ladderValuesWP)){
         x2 <- max(px$position)
         y2 <- 0.1
         peakSlope <- (y2-y1)/(x2-x1)
         peakIntercept <- y1-x1*peakSlope
         px[px$size == 0 & px$position > x1,]$size <- (px[px$size == 0 & px$position > x1,]$position*peakSlope)+peakIntercept
       } else {
         x2 <- px[px$size == ladderValuesWP[peakz+1],]$position[1]
         y2 <- ladderValuesWP[peakz+1]
         peakSlope <- (y2-y1)/(x2-x1)
         peakIntercept <- y1-x1*peakSlope
         px[px$size == 0 & px$position > x1 & px$position < x2,]$size <- (px[px$size == 0 & px$position > x1 & px$position < x2,]$position*peakSlope)+peakIntercept
       }
      }
    }
    
    #Now that the pixels are all aligned and such, lets calculate the relative glycosylation score
    # First, lets get rid of negative values
    write.csv(px, "csv/full_values.csv", row.names = F)
    if (TRUE %in% unique(is.na(px))){
      col_id <- names(px)[unique(which(is.na(px), arr.ind = T)[,2])]
      na_num <- round(100*nrow(px[is.na(px[col_id]),])/nrow(px), 2)
      print(paste0("Residual NA values(", na_num, ") found because Jason sucks. Running it anyway..."))
      px <- px[!is.na(px[col_id]),]
    }
    px$GP_value <- px$GP_value+abs(min(px$GP_value))
    px$WP_value <- px$WP_value+abs(min(px$WP_value))
    px$Ladder_value <- px$Ladder_value+abs(min(px$Ladder_value))
    # Now to subtract the GP_value from the WP_values and get rid of negatives
    px$relGylcoScore <- px$GP_value-px$WP_value
    px$relGylcoScore <- px$relGylcoScore +abs(min(px$relGylcoScore))
    #Now we'll lower quartile normalize because... Z scores don't look as good?
    px$relGylcoScore <- px$relGylcoScore/quantile(px$relGylcoScore)[2]
    
    # The dataframe is saved as the output file
    write.csv(px, "csv/full_values.csv", row.names = F)
    
    cat("\n")
    print("CSVs generated, creating graphs...")
    
    # First, a correlation scatter plot is generated and saved
    draft <- ggplot(data = px, aes(
      x = WP_value,
      y = GP_value,
    ))+
      geom_point(size = 4, aes(color = position))+
      geom_density_2d()+
      geom_abline(slope = 1, intercept = quantile(px$GP_value)[2], size = 2, linetype = 2)+
      theme_classic()+
      xlab("Whole protein signal")+
      ylab("Glycosylation signal")
    draft
    ggsave(paste0("figures/signalCorr_position.pdf"), width = 7.5, height = 7.5, units = "in")
    
    draft <- ggplot(data = px, aes(
      x = WP_value,
      y = GP_value,
    ))+
      geom_point(size = 4, aes(color = relGylcoScore))+
      geom_density_2d()+
      geom_abline(slope = 1, intercept = quantile(px$GP_value)[2], size = 2, linetype = 2)+
      theme_classic()+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = median(px$relGylcoScore))+
      xlab("Whole protein signal")+
      ylab("Glycosylation signal")
    draft
    ggsave(paste0("figures/signalCorr_relGlycoScore.pdf"), width = 7.5, height = 7.5, units = "in")
    
    #Now we make the super cool heat map histograms using the newly calculated diagrams
    draft <- ggplot(data = px, aes(
      x = position,
      y = WP_value,
      color = relGylcoScore,
      fill = relGylcoScore
    ))+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$relGylcoScore))+
      scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$relGylcoScore))+
      geom_bar(stat = "identity")+
      theme_classic2()+
      theme(legend.position = "top")+
      xlab("Relative position")+
      ylab("Whole protein gel intensity")
    draft
    ggsave("figures/relGlyco.pdf", width = 10, height = 4, units = "in")
    draft+geom_line(aes(y = Ladder_value), color = "#5c5c5cff")
    ggsave(paste0("figures/relGlyco_ladder.pdf"), width = 10, height = 4, units = "in")
    
    draft <- ggplot(data = px, aes(
      x = position,
      y = WP_value,
      color = WP_value,
      fill = WP_value
    ))+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$WP_value))+
      scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$WP_value))+
      geom_bar(stat = "identity")+
      theme_classic2()+
      theme(legend.position = "top")+
      xlab("Relative position")+
      ylab("Whole protein gel intensity")
    draft
    ggsave("figures/WP.pdf", width = 10, height = 4, units = "in")
    draft+geom_line(aes(y = Ladder_value), color = "#5c5c5cff")
    ggsave(paste0("figures/WP_ladder.pdf"), width = 10, height = 4, units = "in")
    
    draft <- ggplot(data = px, aes(
      x = position,
      y = GP_value,
      color = GP_value,
      fill = GP_value
    ))+
      scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$GP_value))+
      scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(px$GP_value))+
      geom_bar(stat = "identity")+
      theme_classic2()+
      theme(legend.position = "top")+
      xlab("Relative position")+
      ylab("Glyco protein gel intensity")
    draft
    ggsave("figures/GP.pdf", width = 10, height = 4, units = "in")
    draft+geom_line(aes(y = Ladder_value), color = "#5c5c5cff")
    ggsave(paste0("figures/GP_ladder.pdf"), width = 10, height = 4, units = "in")
    
    cat("\n")
    print("Glycorunner complete.")
  } else{
    print("Ladder not present. Try again; maybe better?")
  }
}

peakCaller <- function(df = px,
                       value = "Ladder_value",
                       threshold = 0,
                       expectedPeaks = F,
                       expectedNumber = 9,
                       ladder,
                       attemptsForPeaks=10,
                       ladderUpOrDown = "lower"){
  peaks <- df
  boomer <- df
  if(threshold != F){
    print(paste0("Threshold given of ", threshold, "."))
    peaks <- peaks[peaks[value] > threshold,]
    
    bandNumber <- 1
    peaks$band <- bandNumber
    for(a in 2:nrow(peaks)){
      if (peaks[a,]$position != peaks[a-1,]$position+1){
        bandNumber <- bandNumber+1
      }
      peaks[a,]$band <- bandNumber
    }
    print(paste0("Total peaks found: ", bandNumber, "."))
    
  } else {
    if(expectedPeaks == T){
      print(paste0("Peaks expected: ", expectedNumber))
      holder <- peaks
      threshold <- mean(unlist(peaks[value]))+sd(unlist(peaks[value]))
      peaks <- peaks[peaks[value] > threshold,]
      
      bandNumber <- 1
      peaks$band <- bandNumber
      for(a in 2:nrow(peaks)){
        if (peaks[a,]$position != peaks[a-1,]$position+1){
          bandNumber <- bandNumber+1
        }
        peaks[a,]$band <- bandNumber
      }
      
      if (max(peaks$band) != expectedNumber){
        counter <- 1
        print("Incorrect peaks found. Exploring other thresholds...")
        while (max(peaks$band) != expectedNumber){
          if (max(peaks$band) > expectedNumber){
            threshold <- threshold+1
          } else {
            threshold <- threshold-1
          }
          print(paste0("Setting threshold at ", round(threshold, 2)))
          peaks <- holder
          peaks <- peaks[peaks[value] > threshold,]
          bandNumber <- 1
          peaks$band <- bandNumber
          for(a in 2:nrow(peaks)){
            if (peaks[a,]$position != peaks[a-1,]$position+1){
              bandNumber <- bandNumber+1
            }
            peaks[a,]$band <- bandNumber
          }
          print(paste0("Found ", bandNumber, " bands"))
          print(paste0("Expected: ", expectedNumber))
          if(bandNumber==expectedNumber){
            break
          }
          counter <- counter+1
          if (counter > attemptsForPeaks){
            print("Max number of attempts reached...")
            print(paste0("Using ", bandNumber, " bands set for the ", ladderUpOrDown, " set of bands."))
            break
          }
        }
      }
      print(paste0("Peaks found: ", bandNumber))
    } else {
      print("Threshold not given.")
      threshold <- mean(unlist(peaks[value]))+sd(unlist(peaks[value]))
      print(paste0("Using mean+sd threshold of: ",threshold,"."))
      peaks <- peaks[peaks[value] > threshold,]
      
      bandNumber <- 1
      peaks$band <- bandNumber
      for(a in 2:nrow(peaks)){
        if (peaks[a,]$position != peaks[a-1,]$position+1){
          bandNumber <- bandNumber+1
        }
        peaks[a,]$band <- bandNumber
      }
      print(paste0("Peaks found: ", bandNumber))
    }
  }
  boomer$size <- 0
  
  if (value == "Ladder_value"){
    for (ladderBand in unique(peaks$band)){
      bandSize <- ladder[ladderBand]
      midPoint <- median(unlist(peaks[peaks$band == ladderBand,]["position"]))
      closestPoint <- min(abs(peaks[peaks$band == ladderBand,]["position"]-midPoint))
      peakCoord <- peaks[abs(peaks["position"]-midPoint) == closestPoint & peaks$band == ladderBand,]$position
      boomer[boomer$position %in% peakCoord,]$size <- bandSize
    }
  }
  
  # for (peakID in unique(peaks$band)){
  #   interim <- subset(peaks, band == peakID)
  #   if(!exists("bands")){
  #     bands <- data.frame(
  #       "name" = unique(interim$name),
  #       "type" = unique(interim$type),
  #       "bandID" = peakID,
  #       "position" = median(interim$position),
  #       "width" = nrow(interim),
  #       "intensity" = sum(interim$adjValue),
  #       "size" = 0
  #     )
  #   } else {
  #     turkey <- data.frame(
  #       "name" = unique(interim$name),
  #       "type" = unique(interim$type),
  #       "bandID" = peakID,
  #       "position" = median(interim$position),
  #       "width" = nrow(interim),
  #       "intensity" = sum(interim$adjValue),
  #       "size" = 0
  #     )
  #     bands <- rbind(bands, turkey)
  #   }
  # }
  # 
  # if (value == "Ladder_value"){
  #   bands$size <- 0
  #   if (bandNumber != length(ladder)){
  #     if (bandNumber > length(ladder)){
  #       bands <- bands[1:length(ladder),]
  #     }
  #     if (ladderUpOrDown == "lower"){
  #       ladder <- ladder[(1+length(ladder)-nrow(bands)):length(ladder)]
  #     } else {
  #       ladder <- ladder[1:length(bands)]
  #     }
  #   }
  #   for (a in 1:length(ladder)){
  #     bands[bands$bandID == a,]$size <- ladder[a]
  #   }
  # }
  # write.csv(bands, paste0("bands/", dfName, ".csv"), row.names = F)
  return(boomer)
}

banderSnatch <- function(){
  bList <- list.files(path = "bands", pattern = "csv")
  for (a in bList){
    if (!grepl("bands", a) & !grepl("ladder", a)){
      a <- paste0("bands/", a)
      if(!exists("allBands")){
        allBands <- read.csv(a)
      } else {
        someBands <- read.csv(a)
        allBands <- rbind(allBands, someBands)
      }
      for (gName in unique(allBands$name)){
        interim <- subset(allBands, name == gName)
        write.csv(interim, paste0("bands/",gName, "_bands.csv"), row.names = F)
        draft <- ggplot(data = interim, aes(
          x = size,
          y = intensity,
          color = type
        ))+
          geom_point(size = 4)+
          theme_classic()
        print(draft)
        ggsave(paste0("figures/",gName, "_bands.pdf"))
      }
    }
  }
}


dirGen <- function(workingDir="./"){
  if (workingDir != "./"){
    setwd(workingDir)
  }
  
  currentFiles <- list.files()
  if ("schema.csv" %in% currentFiles){
    print("Opening schema file")
    cat("\n")
    schema <- read_csv("schema.csv", progress = F, show_col_types = F)
    if(nrow(schema == 0)){
      print("No sample ID's detected. Please fill out the schema file")
      break
    }
    print("Creating directies...")
    for (sampleFile in schema$sample_id){
      if (!sampleFile %in% currentFiles){
        print(paste0("Creating directory for: ", sampleFile))
        dir.create(sampleFile)
      } else if(schema[schema$sample_id==sampleFile,]$Processed == "N"){
        print(paste0("Sample ID: ", sampleFile, " already exists and is unprocessed. Skipping directory generation."))
      } else {
        print(paste0("Sample ID: ", sampleFile, " already exists and is processed. Skipping directory generation."))
      }
    }
  } else {
    print("Schema file not detected.")
    cat("\n")
    break
  }
  cat("\n")
  print("Directories generated.")
  cat("\n")
}

glycoCat <- function(targets="all",
                     graphIt = T,
                     colorBy = "Metadata"){
  currentFiles <- list.files()
  if ("schema.csv" %in% currentFiles){
    print("Opening schema file")
    cat("\n")
    schema <- read_csv("schema.csv", progress = F, show_col_types = F)
    
    for(targetFile in schema$sample_id){
      if (schema[schema$sample_id==targetFile,]$Processed=="N"){
        print(paste0("Processing file: ", targetFile))
        setwd(targetFile)
        glycoRun()
        setwd("../")
        schema[schema$sample_id==targetFile,]$Processed <- "Y"
      }
    }
    write_csv(schema, "schema.csv")
    print("Processing of detected directories completed.")
    cat("\n")
    
    if (targets=="all"){
      targetList <- schema$sample_id
    } else {
      targetList <- schema[schema$Pair=="Y",]$sample_id
      schema <- schema[schema$Pair=="Y",]
    }
    schema$meanRatio <- 0
    schema$medianRatio <- 0
    schema$upperQuartile <- 0
    schema$lowerQuartile <- 0
    for (targetFile in targetList){
      print(paste0("Gathering data from: ", targetFile))
      if(!exists("gatheredData")){
        gatheredData <- read_csv(paste0(targetFile,"/csv/full_values.csv"), show_col_types = F)
        gatheredData$WP_value <- gatheredData$WP_value/quantile(gatheredData$WP_value)[2]
        gatheredData$GP_value <- gatheredData$GP_value/quantile(gatheredData$GP_value)[2]
        gatheredData$Ladder_value <- gatheredData$Ladder_value/quantile(gatheredData$Ladder_value)[2]
        gatheredData$ID <- targetFile
        gatheredData$Metadata <- schema[schema$sample_id == targetFile,]$Metadata
        gatheredData$Replicate <- schema[schema$sample_id == targetFile,]$Replicate
        
        schema[schema$sample_id == targetFile,]$meanRatio <- mean(gatheredData$relGylcoScore)
        schema[schema$sample_id == targetFile,]$medianRatio <- median(gatheredData$relGylcoScore)
        schema[schema$sample_id == targetFile,]$upperQuartile <- quantile(gatheredData$relGylcoScore)[4]
        schema[schema$sample_id == targetFile,]$lowerQuartile <- quantile(gatheredData$relGylcoScore)[2]
      } else {
        interimGather <- read_csv(paste0(targetFile,"/csv/full_values.csv"), show_col_types = F)
        interimGather$WP_value <- interimGather$WP_value/quantile(interimGather$WP_value)[2]
        interimGather$GP_value <- interimGather$GP_value/quantile(interimGather$GP_value)[2]
        interimGather$Ladder_value <- interimGather$Ladder_value/quantile(interimGather$Ladder_value)[2]
        interimGather$ID <- targetFile
        interimGather$Metadata <- schema[schema$sample_id == targetFile,]$Metadata
        interimGather$Replicate <- schema[schema$sample_id == targetFile,]$Replicate
        
        schema[schema$sample_id == targetFile,]$meanRatio <- mean(interimGather$relGylcoScore)
        schema[schema$sample_id == targetFile,]$medianRatio <- median(interimGather$relGylcoScore)
        schema[schema$sample_id == targetFile,]$upperQuartile <- quantile(interimGather$relGylcoScore)[4]
        schema[schema$sample_id == targetFile,]$lowerQuartile <- quantile(gatheredData$relGylcoScore)[2]
        gatheredData <- rbind(gatheredData, interimGather)
      }
    }
    write_csv(gatheredData, "glycoCated.csv")
    write_csv(schema, "annotated_glycoCated.csv")
    
    if(graphIt){
      print("Generating graphs...")
      middi <- median(gatheredData$relGylcoScore)
      draft <- ggplot(data = gatheredData, aes(
        x = position,
        y = WP_value,
        color = relGylcoScore,
        fill = relGylcoScore
      ))+
        scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(middi))+
        scale_fill_gradient2(low = "blue", mid = "green", high = "red", midpoint = mean(middi))+
        geom_bar(stat = "identity")+
        theme_classic2()+
        theme(legend.position = "top")+
        xlab("Relative position")+
        ylab("Relative whole protein gel intensity")+
        facet_wrap(~ID, ncol = 1)
      draft
      ggsave("relGlycoCated.pdf", width = 10, height = 4*nrow(schema), units = "in")
      draft+geom_line(aes(y = Ladder_value), color = "#5c5c5cff")
      ggsave(paste0("relGlycoCated_ladder.pdf"), width = 10, height = 4*nrow(schema), units = "in")
      
      draft <- ggplot(data = gatheredData, aes(
        x = WP_value,
        y = GP_value,
      ))+
        geom_point(size = 4, aes(color = position))+
        geom_density_2d()+
        geom_abline(slope = 1, intercept = 0, size = 2, linetype = 2)+
        theme_classic()+
        xlab("Relative whole protein signal")+
        ylab("Relative glycosylation signal")
      draft
      ggsave(paste0("signalCorr_position.pdf"), width = 7.5, height = 7.5, units = "in")
      
      draft <- ggplot(data = gatheredData, aes(
        x = WP_value,
        y = GP_value,
      ))+
        geom_point(size = 4, aes(color = relGylcoScore))+
        geom_density_2d()+
        geom_abline(slope = 1, intercept = 0, size = 2, linetype = 2)+
        theme_classic()+
        scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = median(gatheredData$relGylcoScore))+
        xlab("Relative whole protein signal")+
        ylab("Relative glycosylation signal")
      draft
      ggsave(paste0("signalCorr_relGlycoScore.pdf"), width = 7.5, height = 7.5, units = "in")
      
      draft <- ggplot(data = gatheredData, aes(
        x = WP_value,
        y = GP_value,
      ))+
        geom_point(size = 4)+
        geom_density_2d()+
        geom_abline(slope = 1, intercept = 0, size = 2, linetype = 2)+
        theme_classic()+
        xlab("Relative whole protein signal")+
        ylab("Relative glycosylation signal")
      if (colorBy == "Metadata"){
        draft <- draft+geom_point(size=4, aes(color = Metadata))
      } else {
        draft <- draft+geom_point(size=4, aes(color = ID))
      }
      draft
      ggsave(paste0("signalCorr_metadata.pdf"), width = 7.5, height = 7.5, units = "in")
    }
    
  } else {
    print("Schema file not detected.")
    cat("\n")
    break
  }
  print("Completed.")
}

runGlyco <- function(directoryLocation="./", 
                     graphing = T,
                     pairing = "all",
                     color_by = "Metadata"){
  genDirs <- readline(prompt = "Welcome to GlycoRunn. Would you like to generate new directories (Y/n)?")
  if (genDirs!="n"){
    dirGen(workingDir = directoryLocation)
  }
  cat("\n")
  catTheGlycos <- readline(prompt = "Would you like go ahead and process data (y/N)?")
  if(catTheGlycos == "y"){
    glycoCat(targets = pairing,
             graphIt = graphing, 
             colorBy = color_by)
  }
}

