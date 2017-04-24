\name{dataframe.outlier}
\alias{dataframe.outlier}
\title{
Trata outlier e missing em data.frame
}
\description{
Trata outliers em data.frames convertendo em séries temporais ou não
}
\usage{
dataframe.outlier(dados, n_passos_frente, freq, normalize = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados}{
Data.frame com variáveis seguindo o modelo (data.frame (DATA, var1, var2, var3, ..., varn)
}
  \item{n_passos_frente}{
Número de passos utilizados na projeção, necessário para a função ts.dados trabalhar com datas
}
  \item{freq}{
Frequência dos dados. Mensal (12), diário(7), horário(24)
}
  \item{normalize}{
Se TRUE faz normalização dos dados
}
  \item{\dots}{
Passagem de parâmetros
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
Esta função utiliza \code{\link{ts.dados}}
}
\examples{
data(mensal)
Df <- dataframe.outlier(mensal[, 1:50], 10, 12)
head(Df)
}