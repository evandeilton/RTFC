\name{ForecastHorizon}
\alias{ForecastHorizon}

\title{
Define horizontes de forecast com dados com ou sem index de data
}

\usage{
ForecastHorizon(XtsData, tsfrequency, horizon)
}
\arguments{
  \item{XtsData}{
Dados formato xts ou data.frame formato data-valor
}
  \item{tsfrequency}{
Frequência min, hour, day, month year
}
  \item{horizon}{
Horizonte de forecast
}
}

\value{
Objeto de dados séries temporais
}

\author{
LOPES, J. E.
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Vetor
x <- ForecastHorizon(rnorm(100), 'day',20)
str(x)

## data.frame
data(diario)
y <- ConvertDataToTs(diario[,1:2], tsfrequency = "day", OutType = "xts")
y <- ForecastHorizon(y, 'day', 20)
}