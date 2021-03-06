\name{forecastMethod}
\alias{forecastMethod}

\title{
Metodos de forecast implementados
}

\description{
Recebe como entrada um objeto s�rie temporal e decide qual o melhor m�todo de forecast entre v�rios. Os testes utilizados s�o de tend�ncia e de linearidade.
}

\usage{
forecastMethod(x)
}

\arguments{
 \item{x}{Serie temporal para escolha do melhor modelo}
 \item{fcMethod}{Escolha do m�todo de forecast definido pelo utilizador.}
}

\value{
Lista com cinco ou menos nomes dos melhores m�todos de forecast recomendados para a s�rie.
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
