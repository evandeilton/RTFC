\name{trata.outliers.ts}
\alias{trata.outliers.ts}
\title{
Trata outliers em serie temporal
}
\description{
Trata outliers e missing em série temporal
}
\usage{
trata.outliers.ts(dados.ts, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados.ts}{
Série temporal já tratada pela função ts.dados
}
  \item{\dots}{
Passagem de parâmetros
}
}

\value{
Lista com dados históricos, tratados e outras informações
}

\author{
LOPES, J. E.

}

\seealso{
Veja também \code{\link{dataframe.outlier}} para lidar com data frames e \code{\link{tsclean}} para informações sobre metodologia utilizada.
}
\examples{
data(mensal)
dats <- ts.dados(mensal[,1:2], 10, 12)
td <- trata.outliers.ts(dats)
names(td)
}