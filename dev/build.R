
library(roxygen2)
# setwd("~/git")
# devtools::create("recalc")

setwd("~/git/recalc")

devtools::document()

devtools::build_vignettes()

devtools::check()

#devtools::install()

# or from github, after push
library(devtools)
install_github("ianhussey/recalc")

library(recalc)

?recalc

detach("package:recalc", unload=TRUE)
