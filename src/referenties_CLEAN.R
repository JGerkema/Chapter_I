install.packages("NCmisc")
library(NCmisc)

c("FD", "tidyverse", "FD", "janitor", "progress", "readxl", "readr", "ggplot2", 
  "ggthemes", "hrbrthemes", "writexl", "RColorBrewer", "vegan", "foreach", 
  "doSNOW", "roxygen2", "docstring", "renv", "conflicted", "AICcmodavg", "rstatix", 
  "ggpubr", "car", "devtools", "taxonlookup", "WorldFlora", "rtry", "lmerTest") %>%
  map(citation) %>%
  print(style = "bibtex") 



get_pkgs_info(pkgs = c("FD", "tidyverse", "janitor", "progress", "readxl", "readr", "ggplot2", 
                       "ggthemes", "hrbrthemes", "writexl", "RColorBrewer", "vegan", "foreach", 
                       "doSNOW", "roxygen2", "docstring", "renv", "conflicted", "AICcmodavg", "rstatix", 
                       "ggpubr", "car", "devtools", "taxonlookup", "WorldFlora", "rtry", "lmerTest", "ggh4x"), 
              out.dir = getwd(), out.format = "docx")


cite_packages(pkgs = c("ggplot2", "FD", "tidyverse", "janitor", "progress", "readxl", "readr",  
                       "ggthemes", "hrbrthemes", "writexl", "RColorBrewer", "vegan", "foreach", 
                       "doSNOW", "roxygen2", "docstring", "renv", "conflicted", "AICcmodavg", "rstatix", 
                       "ggpubr", "car", "devtools", "taxonlookup", "WorldFlora", "rtry", "lmerTest", "ggh4x"),
              out.format = "docx", out.dir = ".")



