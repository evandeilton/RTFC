\name{ConvertDataToTs}

\title{
Converte dados em series temporais
}
\description{
Recebe como entrada data.frame, matrix ou vetor de dados tipo data-valor, converte em ts, mts ou xts.
}
\usage{
ConvertDataToTs(Data, tsfrequency = "month", OutType = "ts", OutlierClean = TRUE, ...)
}

\arguments{
  \item{Data}{
data.frame, matrix ou vetor de dados tipo data-valor
}
  \item{tsfrequency}{
Frequancia dos dados min, hour, day, month e year
}
  \item{OutType}{
Saida como ts ou xts
}
  \item{OutlierClean}{
Remover aoutliers, TRUE ou FALSE
}
  \item{\dots}{
Passagem de argumentos.
}
}

\value{
Objeto de dados para modelagem de series temporais.
}

\author{
LOPES, J. E.
}

\seealso{
\code{\link{xts}}
}

\examples{
x <- ConvertDataToTs(rnorm(100))
str(x)
}
