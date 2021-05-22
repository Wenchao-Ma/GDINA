.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}
# .onUnload <- function (libpath) {
#   library.dynam.unload("mypackage", libpath)
# }

StartWelcomeMessage <- function(){

  paste("GDINA R Package ",
        "(version ", utils::packageDescription("GDINA")$Version,
        "; ",utils::packageDescription("GDINA")$Date, ")\n",
        "For tutorials, see https://wenchao-ma.github.io/GDINA\n",
        sep="")
}
