# Title: Stock Flow Consistent Model

# Purpose : This script serves to run all function that are written
#           and then stored in the folder R. These functions also 
#           include function for data cleaning and processing.
#           The script should use as little as possible of own code.
#           Everything should be reproducible every time the R is 
#           restarted. These main results regarding data and the model
#           are then saved. Outputs are created in a different script
#           which should be more flexible and versatile.

# Cite : citation() and used packages

# DataFile:

# Author: Alexander Toplitsch
# Contact details: alexander.toplitsch@s.wu.ac.at

# Date script created: We 11.10.2023 19:48 -------------
# Date script last modified:  ----

citation()
# Load Functions
source("R/packages.R")

# Run all necessary functions

packages(libs = c("tidyverse", "leaflet", "lubridate", "readxl", "sfcr", "scales"))





# For good practice
xfun::session_info()


