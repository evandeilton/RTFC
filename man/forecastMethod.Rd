\name{forecastMethod}
\alias{forecastMethod}

\title{
Metodos de forecast implementados
}

\description{
Recebe como entrada um objeto série temporal e decide qual o melhor método de forecast entre vários. Os testes utilizados são de tendência e de linearidade.
}

\usage{
forecastMethod(x)
}

\arguments{
 \item{x}{Serie temporal para escolha do melhor modelo}
 \item{fcMethod}{Escolha do método de forecast definido pelo utilizador.}
}

\value{
Lista com cinco ou menos nomes dos melhores métodos de forecast recomendados para a série.
}
\references{
Em estudo...
}
\author{
LOPES, J. E.
}

\seealso{
\code{\link{nparTrend}}, \code{\link{linearityTest}}
}
\examples{
## vetor
x <- ForecastHorizon(rnorm(100), 'day',20)
str(x)
forecastMethod(x[[1]])

## data.frame
data(diario)
y <- ConvertDataToTs(diario[,1:2], tsfrequency = "day", OutType = "xts")
y <- ForecastHorizon(y, 'day', 20)
forecastMethod(y[[1]])
}
