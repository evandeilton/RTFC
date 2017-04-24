CRANChoosen <- function()
{
  return(getOption("repos")["CRAN"] != "@CRAN@")
}

UsePackage <- function(package, defaultCRANmirror = "http://cran.at.r-project.org") 
{
  if(!InstalledPackage(package))
  {
    if(!CRANChoosen())
    {       
      chooseCRANmirror()
      if(!CRANChoosen())
      {
        options(repos = c(CRAN = defaultCRANmirror))
      }
    }
    
    suppressMessages(suppressWarnings(install.packages(package)))
    if(!InstalledPackage(package)) return(FALSE)
  }
  return(TRUE)
}


## Carrega os pacotes necessários, sem exibir logs dos mesmos.
rpacks <- function(..., install = TRUE){
    reqFun <- function(pack) {
        if(!suppressWarnings(suppressMessages(require(pack, character.only = TRUE)))) {
            message(paste0("Pacote nao carregado", pack, ": Tentando instalar e carregar!"))
            install.packages(pack)
            require(pack, character.only = TRUE)
        }
    }
    lapply(..., reqFun)
}

#
Rpacks <- function(...) {
	rpacks(c("plyr", "dplyr", "glmnet", "zoo","fBasics","forecast","tseries","TSA","lmtest","tsDyn","mondate", "reshape","relaimpo", "RColorBrewer"), install = TRUE)
}

# Chamada
#Rpacks()

packs <- c("plyr", "dplyr", "glmnet", "zoo","fBasics","forecast","tseries","TSA","lmtest","tsDyn","mondate","reshape","relaimpo", "RColorBrewer")

lapply(packs, FUN = function(X) {suppressMessages(do.call("require", list(X, quietly = TRUE)))})


