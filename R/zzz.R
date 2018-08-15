.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}


StartWelcomeMessage <- function(){

  paste("GDINA Package for Cognitive Diagnosis Modeling\n",
        "Version ", utils::packageDescription("GDINA")$Version,
        " (",utils::packageDescription("GDINA")$Date, ")\n",
        "For examples and more information, \nvisit https://wenchao-ma.github.io/GDINA\n",
        sep="")
}
