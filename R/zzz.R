.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}
# .onUnload <- function (libpath) {
#   library.dynam.unload("mypackage", libpath)
# }

StartWelcomeMessage <- function(){

  paste("GDINA Package",
        " [Version ", utils::packageDescription("GDINA")$Version,
        "; ",utils::packageDescription("GDINA")$Date, "]\n",
        "More information: https://wenchao-ma.github.io/GDINA\n",
        sep="")
}
