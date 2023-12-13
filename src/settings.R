#-------------------------- Packages and set up ------------------------------#
list_of_packages <- c("tidyverse", "FD", "janitor", "progress", "readxl",
                      "ggplot2", "ggthemes", "hrbrthemes", "writexl", "readr", "progress",
                      "RColorBrewer", "gcookbook", "janitor", "ggpmisc",
                      "vegan", "foreach", "doSNOW", "docstring", "renv", "ggcorrplot",
                      "reshape2", "conflicted", 
                      "remotes", "roxygen2", "abind", "prettyunits", "extrafont", "shiny",
                      "quantreg", "snow", "sp", "askpass", "gtools", 
                      "maptools", "minpack.lm" , "devtools", "geosphere",
                      "scales", "ecodist")

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(FD)
library(janitor) # Row to column names
library(progress) # Progress bar
library(readxl) # Import excel files
library(readr) # Import csv files
library(ggplot2)
library(ggthemes) # More themes for ggplot2
library(hrbrthemes) # More themes for ggplot2
library(writexl)
library(progress) # Progress bar
library(RColorBrewer) # More colors in ggplot2
library(gcookbook)
library(ggpmisc) # Correlation in ggplot2
library(vegan)
library(foreach) # Package for a parallel loop
library(doSNOW)
library(roxygen2) # Needed ffor docstring?
library(docstring) # Provide documentation for custom build functions
library(renv)
library(ggcorrplot) # Correlation plot
library(reshape2) # Needed for ggcorrplot
library(conflicted) # Function conflict_prefer
library(remotes)
#library(abind) # Script gave error otherwise
#library(prettyunits) # script gave error otherwise
#library(extrafont) # Script gave error otherwise
#library(shiny) # Script gave error otherwise
#library(quantreg) # Script gave error otherwise
#library(snow) # Script gave error otherwise
#library(maps) # Script gave error otherwise
#library(mapdata) # Script gave error otherwise
#library(sp) # Script gave error otherwise
#library(askpass) # Script gave error otherwise
#library(gtools) # Script gave error otherwise
#library(maptools) # Script gave error otherwise
#library(minpack.lm) # Script gave error otherwise
library(devtools)
library(taxonlookup)
library(geosphere) # Calculate distances between plots
library(scales) # Gradient colours in ggplot
library(ecodist) # Bray Curtis index

#remotes::install_github("traitecoevo/taxonlookup")

# Install and load austraits packages
#if("austraits" %notin% installed.packages()[,"Package"]) remotes::install_github("traitecoevo/austraits", dependencies = TRUE, upgrade = "never")
#if("ausplotsR" %notin% installed.packages()[,"Package"]) remotes::install_github("ternaustralia/ausplotsR", build_vignettes = TRUE, dependencies = TRUE, upgrade = "never")

#library(austraits)
#library(ausplotsR)

# Indicating preferences
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("recode", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("load", "base")
conflict_prefer("filter", "dplyr")
conflict_prefer("where", "dplyr")
conflict_prefer("check", "remotes")




# Display numbers without scientific notation
options(scipen = 15)   

# Set the seed to make reproducible
set.seed(19970606)



