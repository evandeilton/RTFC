\name{switch.cvforecast}
\alias{switch.cvforecast}
\title{
Escolha de modelos de forecast
}
\description{
Com base em uma lista de opções faz ajuste de forecasts baseado em diversas metodologias
}
\usage{
switch.cvforecast(x, nmodelo, h, level = 95, onlyfc = FALSE)
}

\arguments{
  \item{x}{
Série temporal
}
  \item{nmodelo}{
Nome do método de forecast
}
  \item{h}{
Horizonte de forecast
}
  \item{level}{
Nível de significância para as bandas de confiança
}
  \item{onlyfc}{
Somente forecast ou modelos completo?
}
}

\value{
Objeto de forecast
}
\references{
Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast package for R", Journal of Statistical Software, 26(3).
}
\author{
LOPES, J. E.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{etsForecast}},\code{\link{auto.arimaForecast}}, \code{\link{HWsForecast}}, etc.
}
\examples{
data(diario)
fit <- switch.cvforecast(diario[,4], "etsForecast", h=20, level=95)
class(fit)
plot(fit)
}