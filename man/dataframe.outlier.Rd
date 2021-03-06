\name{dataframe.outlier}
\alias{dataframe.outlier}
\title{
Trata outlier e missing em data.frame
}
\description{
Trata outliers em data.frames convertendo em s�ries temporais ou n�o
}
\usage{
dataframe.outlier(dados, n_passos_frente, freq, normalize = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados}{
Data.frame com vari�veis seguindo o modelo (data.frame (DATA, var1, var2, var3, ..., varn)
}
  \item{n_passos_frente}{
N�mero de passos utilizados na proje��o, necess�rio para a fun��o ts.dados trabalhar com datas
}
  \item{freq}{
Frequ�ncia dos dados. Mensal (12), di�rio(7), hor�rio(24)
}
  \item{normalize}{
Se TRUE faz normaliza��o dos dados
}
  \item{\dots}{
Passagem de par�metros
}
}

\value{
Lista com um objeto da classe mts, ts  ou matrix contendo dados tratados em forma de serie e um data.frame com os dados originais
\item{dados_ts}{Objeto mts, ts  ou matrix}
\item{dados}{data.frame dados originais tratados}
}

\author{
LOPES, J. L.
}

\seealso{
Esta fun��o utiliza \code{\link{ts.dados}}
}
\examples{
data(mensal)
Df <- dataframe.outlier(mensal[, 1:50], 10, 12)
head(Df)
}