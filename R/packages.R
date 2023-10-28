
###########################################################################
#-------------------------Call all necessary packages---------------------#
###########################################################################

packages <- function(libs){
  libs <- libs
  installed_libs <- libs %in% rownames(installed.packages())
  if (any(installed_libs == F)) {
    install.packages(libs[!installed_libs])
  }
  invisible(lapply(libs, library, character.only = T)) 
}


