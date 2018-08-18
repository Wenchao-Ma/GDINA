.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}


StartWelcomeMessage <- function(){

  paste("GDINA Package",
        " [Version ", utils::packageDescription("GDINA")$Version,
        "; ",utils::packageDescription("GDINA")$Date, "]\n",
        "For examples and more information, \nsee https://wenchao-ma.github.io/GDINA\n",
        sep="")
}
