\name{linearityTest}
\alias{linearityTest}
\title{
Testes de linearidade
}

\usage{
linearityTest(x, Test)
}


\arguments{
  \item{x}{
Série temporal
}
  \item{Test}{
  Tipo de teste linearidade "terasvirta", "white", "keenan", "mcleodLi", "tsay", "tarTest"
  }
}

\details{
Em análise...
}

\value{
Resultado dos testes
}
\references{
Em estudo
}
\author{
LOPES, J. E.
}

\seealso{
\code{\link{terasvirta.test}}, \code{\link{white.test}}, \code{\link{Keenan.test}}, \code{\link{McLeod.Li.test}}, \code{\link{Tsay.test}}, \code{\link{tlrt}}
}
\examples{

## dados
data(diario)
y <- ConvertDataToTs(diario[,1:2], tsfrequency = "day", OutType = "ts")
linearityTest(y, "terasvirta")
}
