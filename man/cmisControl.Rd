\name{cmisControl}
\alias{cmisControl}
\title{
Controle de argumentos generico CMIS
}
\usage{
cmisControl(maxHorizon = 20, level = 95, cvMethod = "MAPE", tsfrequency = "day",
	OutType = "ts", outputFormat = "forecast", OutlierClean = TRUE, onlyfc = FALSE,
	userJDBC = "CMIS_OWNER", passJDBC = "CMIS_OWNER", driverJDBC = "ojdbc6.jar")
}

\arguments{
  \item{maxHorizon}{
Numero de pontos a projetar (horizonte de forecast)
}
  \item{level}{
Intervalor de confian�a fara os forecasts
}
  \item{cvMethod}{
Estat�stica de decis�o. Veja \code{\link{accuracy}}
}
  \item{tsfrequency}{
Frequencia dos dados, pode ser min, hour, day, month e year
}
  \item{OutType}{
Tipo de sa�da dos dados, ts ou xts
}
  \item{outputFormat}{
Tipo de saida da proje��o, forecast, df (data.frame) ou both (ambos)
}
  \item{OutlierClean}{
Tratar outliers
}
  \item{onlyfc}{
Para as fun��es wrapers, apenas forecast ou modelo completo
}
  \item{residlevel}{
Livel de signific�ncia para an�lise de res�duos. Padr�o � 0.05, mas pode ser definido pelo utilizador.
}
  \item{userJDBC}{
Usu�rio do banco de dados
}
  \item{passJDBC}{
Senha do usu�rio
}
  \item{driverJDBC}{
Driver RJDBC
}
}

\author{
LOPES, J. E.
}


\examples{
## Not run
# Monta controle 
#Control <- cmisControl(
#  maxHorizon = 45,
#  level = 95,
#  onlyfc = FALSE,
#  cvMethod = "MAPE",
#  tsfrequency = "hour",
#  outputFormat = "forecast",
#  driverJDBC = "ojdbc6.jar"
#)

#Coleta dados
#dd <- CMISForecastData(dados = NULL, 7, 1, 10000, connect = TRUE, Control = Control)
#str(dd)
}
