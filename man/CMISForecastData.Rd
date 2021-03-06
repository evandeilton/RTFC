\name{CMISForecastData}
\alias{CMISForecastData}
\title{
Coletor de dados do CMIS para forecasts
}

\usage{
CMISForecastData(dados = NULL, freq, idmini, idmfin, connect = FALSE, 
	Control = cmisControl(), ...)
}

\arguments{
  \item{dados}{
Se dados forem informados, n�o haver� conex�o com banco de dados.
}
  \item{freq}{
Frequencia dos dados: 1 - 24	Hora; 2 - 12	Mensal; 3 - 7	Di�rio
}
  \item{idmini}{
Faixa m�nima inicial de id_metrica_frequencia
}
  \item{idmfin}{
Faixa m�xima de id_metrica_frequencia
}
  \item{connect}{
Conecatar, sim ou n�o
}
  \item{Control}{
Controle gen�rica para passagem de argumentos.
}
  \item{\dots}{
Passagem de argumentos
}
}
\author{
LOPES, J. E.
}

\seealso{
\code{\link{cmis_JDBC_Connect}}
}
\examples{
# Not run
# Monta controle 
#Control <- cmisControl(
#  maxHorizon = 45,
# level = 95,
#  onlyfc = FALSE,
#  cvMethod = "MAPE",
#  tsfrequency = "hour",
#  outputFormat = "forecast",
#  driverJDBC = "ojdbc6.jar"
#)

#Coleta dados do CMIS
#dd <- CMISForecastData(dados = NULL, 7, 1, 10000, connect = FALSE,
#	Control = Control)
#str(dd)
}
