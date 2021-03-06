\name{hwForecast}
\alias{hwForecast}
\title{
Forecast por HoltWinters
}
\usage{
hwForecast(x, h, level = 95, onlyfc = TRUE, ...)
}

\arguments{
  \item{x}{
S�rie temporal
}
  \item{h}{
Horizonte de forecast
}
  \item{level}{
Nivel de confian�a para os limites superiores e infeiores
}
  \item{onlyfc}{
Somente forecast ou modelo completo
}
  \item{\dots}{
Passagem de argumentos
}
}

\value{
Objeto de forecast
}

\seealso{
\code{\link{forecast}}
}
\examples{

## dados
data(diario)
y <- ConvertDataToTs(diario[,1:2], tsfrequency = "day", OutType = "ts")
hwForecast(y, h=20)
}
