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
Se dados forem informados, não haverá conexão com banco de dados.
}
  \item{freq}{
Frequencia dos dados: 1 - 24	Hora; 2 - 12	Mensal; 3 - 7	Diário
}
  \item{idmini}{
Faixa mínima inicial de id_metrica_frequencia
}
  \item{idmfin}{
Faixa máxima de id_metrica_frequencia
}
  \item{connect}{
Conecatar, sim ou não
}
  \item{Control}{
Controle genérica para passagem de argumentos.
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
