\name{nparTrend}
\alias{nparTrend}

\title{
Analise de tendencia por testes nao-parametricos
}

\usage{
nparTrend(x, ...)
}

\arguments{
  \item{x}{
S�rie temporal
}

  \item{npoints}{
Define o n�mero de pontos para fazer an�lise de tend�ncia (do mais recente para o mais antigo de acordo com a ordem dos dados).
}

  \item{\dots}{
Passagem de argumentos
}
}

\value{
Linear ou n�o linear
}
\references{
Em constru��o...
}
\author{
LOPES, J. E.
}
\seealso{
\code{\link{rkt}}
}
\examples{
## dados
data(diario)
y <- ConvertDataToTs(diario[,1:2], tsfrequency = "day", OutType = "ts")
nparTrend(y)
}