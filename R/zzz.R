.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}


StartWelcomeMessage <- function(){

  paste("==============================================\n",
        "GDINA Package for Cognitive Diagnosis Modeling\n",
        "        Version ", utils::packageDescription("GDINA")$Version,
        " (",utils::packageDescription("GDINA")$Date, ")\n",
        "==============================================\n",
        sep="")
}
