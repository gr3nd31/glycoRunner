# glycoRunner
Script to directly compare whole protein gels to glycosylation gels. Ideally, the length of the lane extractions between the whole protein and the glycosylated protein should be the same, but this version of the script should be able to adjust lanes of varying sizes.
Additionally, extractions should be performed using a box ROI in ImageJ. While line extractions *can* work, they are extremely sensitive to random noise in the image.
The code has additional functionality to automatically call specific protein band peaks, but this is currently disabled until someone asks for it back. It *only* correlates data and makes pretty graphs now.
Relative glycosylation score is calculated by subtracting the whole protein values from the glycosylation protein values at each pixel. and performing lower quartile normalization.

### Example of the whole protein heatmap:
https://github.com/gr3nd31/glycoRunner/edit/main/Example_relGlyco.pdf

### Example of the intensity correlation plot:
https://github.com/gr3nd31/glycoRunner/edit/main/Example_signalCorr_relGlycoScore.pdf

## How to run:

1. Open R and source the glycorunn.R script
2. Setwd() to a directory containing the data files. Each lane or sample should have its own directory containing nothing but its specific files.
3. Call glycoRun() script, which takes the following parameters:
  - **ladderWP**: The csv file containing the whole protein ladder. Default is "ladder.wp.csv"
  - **LadderGP**: The csv file containing the glycoprotein ladder. This may be the same file as the ladderWP. Default is "ladder.gp.csv"
  - **ladderValuesWP**: An array of kilodalton sizes marked by the whole protein ladder. This list will be used to assign ladder peaks to sizes. Default is 'c(205, 150, 100, 75, 50, 37, 25, 20, 15, 10)'
  - **ladderValuesGP**: An array of kilodalton sizes marked by the whole protein ladder. This list will be used to assign ladder peaks to sizes. If the laddWP and ladderGP files are the same, this is ignored. Default is 'c(75, 25)'.
  - **wp**: The csv containing the whole protein lane extraction. This *MUST* be the same length as the ladderWP file. Default is "wp.csv"
  - **gp**: The csv containing the glycosylated protein lane extraction. This *MUST* be the same length as the ladderGP file. Default is "gp.csv"
  - **wp_background_file**: Depending on the quality of the gel or how it was imaged, it may be necessary to subtract non-specific background intensities from the whole protein lane. Therefore, you should include a control extraction csv from a nearby lane without any protein. If this file is not present, then background subtraction will not occur as this version no longer supports dynamic background adjustment.
  - **gp_background_file**: For the same reason as outlined for the wp_background_file, the gp data should also be corrected by using a control lane extraction from the gp image without any glycosylated signal.
  - **outputCSV**: The name of the output CSV.
  - **peakCSV**: The name of the CSV containing estimated peak information.
  - **givenThresholds**: When calling peaks, an intensity threshold can be given. This is an artifact of the original code that pulled protein peak data and has been left in the case that this functionality is requested again. Default is 'F'.
  - **ladderAttempts**: The script will try to find the number of peaks in the ladder data corresponding to the number of values given. If it can't, it will adjust the thresholds up or down to find the appropriate number of peaks. This parameter controls how many times this correction will be performd until it gives up. Default is '10'.
  - **upperOrLowerLadder**: If the script cannot find the appropriate number of peaks defined by the 'ladderAttempts' parameter, this parameter tells the script whether it should use the 'upper' or 'lower' part of the ladderValues array to assign peak sizes. Default is 'lower'.

4. Assuming you use the default file names, none of these paramters should need to be changed when running the script. It will generate two new directories ('csv' and 'figures') to save the generated data and figures.

## Considerations:
- While the script is fairly robust to small differences in how far the gels have run, the closer in length the whole protein and glycosylated protein extractions are, the better the alignment will be.
- The sample extraction and ladder extractions **MUST** be the same length.
- The line generated on the correlation graph is a slope of 1 with a y-intercept of the lower quartile value of the GP_values. Hypo- or hypergylcosylations should appear as distributions above or below the line.
