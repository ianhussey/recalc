
library(roxygen2)
#setwd("~/git/")
#devtools::create("recalc")
setwd("~/git/recalc")

devtools::document()

devtools::check(vignettes = FALSE)

#devtools::install()
# or from github, after push
devtools::install_github("ianhussey/recalc")

library(recalc)

?recalc

detach("package:recalc", unload=TRUE)

# once you have the package updated, you can use it to build the vignettes, check the whole thing, and reinstall again
devtools::build_vignettes()
devtools::check()

