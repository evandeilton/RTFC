\name{plot.multiforecast}
\alias{plot.multiforecast}
\title{
Plota graficos de objetos da classe 'multiforecast'
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

fit <- cmisMultiforecast(diario[, 5], Control)
plot(fit)
}
