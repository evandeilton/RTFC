\name{summary.multiforecast}
\alias{summary.multiforecast}
\title{
Funcao generica para objeto multiforecast
}
\description{
Exibe estatísticas e os graficos de forecasts multiplos
}
\usage{
summary.multiforecast(obj, digits = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
Objeto multiforecast
}
  \item{digits}{
Dígitos da saida
}
}

\author{
LOPES, J. E.
}

\examples{
## dados
data(diario)

## Controle
Control <- cmisControl(
  maxHorizon = 45,
  level = 95,
  onlyfc = FALSE,
  cvMethod = "MAPE",
  tsfrequency = "day",
  outputFormat = "forecast"
)

fit <- cmisMultiforecast(diario[,10], Control)
summary(fit)
}
