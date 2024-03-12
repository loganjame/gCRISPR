# init.R
#
# Example R code to install packages if not already installed
#
my_packages = c("shiny", "ggplot2","ggrepel", "reshape2", "dplyr" , "stringr" , "gridExtra", "DT")
install_if_missing = function(p) {
  if (p %in% rownames(installed.packages()) == FALSE) {
    install.packages(p)
  }
}
invisible(sapply(my_packages, install_if_missing))
