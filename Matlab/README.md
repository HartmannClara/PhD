# Pellet
MATLAB code for processing Results.csv files from experiment 1 (Pellet); 
# Pellet_date 
specific date files 

# Preprocessing
MATLAB code for preprocessing imaging files for the drive maze;

# Drivemaze.plotting
MATLAB code to plot Drive maze ethoplot, with imagig data (all cells or selective); includes frame cutting from start and errors;
averaging plots 

# aviTOtif
ImageJ macro to batch process several files within a folder structure
all it does is open the .avi file and save it as a .tiff under the same name

# Minian
Resorting of Raw Minian files 

# Missed_frames
Code to correct for missed frames in miniscope imaging data; also corrects event list from drivemaze accordingly

# DM_WIP_behavior2
script to analyze behavioral measures from the drivemaze, Work in progress

# A_DrivemazePipeline_ci_data
prouces ci_data file from raw traces imagej csv and drivemaze event list. corrects missed frames and saves in bsl-animals-day-session-data

# B__Shuffled_DrivemazePipelineAnalysis_events_averg
produces events_averg file from shuffled/shifted data

# C_DrivemazePipelineShufTresh
analyses shuffled data to produce histogram 5th and 95th percentile files for thresholding

# B_DrivemazePipelineAnalysis_events_averg
produces events_averg file after analysis of ci_data file

# C_DrivemazePipelinePlots
plotting from events_averg file, analyzing cells with shuffled thresholds