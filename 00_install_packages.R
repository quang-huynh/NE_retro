
# Install DLMtool, then MSEtool, then mfSCA package (some miscellaneous functions for running the simulation)
install.packages("DLMtool_5.4.3.tar.gz", repos = NULL)
install.packages("MSEtool_2.0.0.tar.gz", repos = NULL)

devtools::install("mfSCA")

# Some other useful packages
install.packages("dplyr")
install.packages("ggplot2")