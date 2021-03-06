\name{trata.outliers.ts}
\alias{trata.outliers.ts}
\title{
Trata outliers em serie temporal
}
\description{
Trata outliers e missing em s�rie temporal
}
\usage{
trata.outliers.ts(dados.ts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados.ts}{
S�rie temporal j� tratada pela fun��o ts.dados
}
  \item{\dots}{
Passagem de par�metros
}
}

\value{
Lista com dados hist�ricos, tratados e outras informa��es
}

\author{
LOPES, J. E.

}

\seealso{
Veja tamb�m \code{\link{dataframe.outlier}} para lidar com data frames e \code{\link{tsclean}} para informa��es sobre metodologia utilizada.
}
\examples{
data(mensal)
dats <- ts.dados(mensal[,1:2], 10, 12)
td <- trata.outliers.ts(dats)
names(td)
}