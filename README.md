# cyano_forecast
Updated INLA cyanobacteria forecast repository. 

This repository contains all code used in the INLA cyanobacteria study. It does not, however, include inpute files due to file size limitations. All input files are publicallly available for free download.  

Code and general workflow:  
1. Week assignments  
  a. generate_week_assignments_tibble.R  
    - Purpose: Generates a csv of dates in our study period with the week of the year it falls on. This was done to correspond with the CyAN data week numbers, wherein the first week of each year begins on the first Sunday of each year and counts until the first Sunday of the following year.  
    - Sbatch file: generate_week_assignments.sbatch  
2. Lake morphology  
  a. lake_morpho_code.R  
    - Purpose: Calculates lake morphology data for INLA model  
    - Sbatch file: lake_morpho.sbatch  
    - __NOTE__: While I configured this code so it _could_ run on atmos, it consistently fails on a different lake each time. Therefore, this code needs to be run on a local machine.  
3. CyAN data  
  a. parallel.step1plus2.R  
    - Purpose: Crops OLCI cyan images to conus boundary and masks out values that aren't resolvable lakes. Residual code exists in the file to mosaic individual cyan tiles, but this isn't necessary anymore. I've also recoded it so that the code is in the main working directory, OLCI_preprocessing-specific files remain in the OLCI_preprocessing subdirectory, and more generalized files are in the data directory. This helps prevent mix-ups with the resolvable lakes shapefile, etc.  
    - Sbatch file: step2plus2.sbatch  
b. cyanoCONUS_ice_step3_adjusted.R  
    - Purpose: Creates weekly ice masks as tifs and shapefiles and fills any holes within the mask. Also crops the mask to the CONUS boundary. As a side note, I don't think all the inputs are actually used in the code.  
    - Sbatch file: step3_adjusted.sbatch  
c. 
4. Ice data  

5. PRISM data  

6. Water temp data  

7. Data compilation  

8. INLA model  
