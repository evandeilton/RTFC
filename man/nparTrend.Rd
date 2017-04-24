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
Série temporal
}

  \item{npoints}{
Define o número de pontos para fazer análise de tendência (do mais recente para o mais antigo de acordo com a ordem dos dados).
}

  \item{\dots}{
Passagem de argumentos
}
}

\value{
Linear ou não linear
}
\references{
Em construção...
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