\name{ts.dados}
\alias{ts.dados}
\title{
Trata dados para a serie temporal
}
\description{
Trata dados para a série temporal com datas no formato adequado para ajustes
}
\usage{
ts.dados(dados_hist, n_passos_frente, freq, normalize = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados_hist}{
Vetor de dados históricos no formato data.frame(data, realizado)
}
  \item{n_passos_frente}{
Número de passos utilizados na projeção (quantos pontos serão projetados)
}
  \item{freq}{
Frequência dos dados: Mensal (12), diário(7), horário(24)
}
  \item{normalize}{
Se TRUE faz normalização dos dados. Vide \code{\link{scale}}
}
  \item{\dots}{
Passagem de parâmetros.
}
}
\details{
Esta função serve para transformar o vetor de dados brutos vindos das bases do CMIS de modo que a sazonalidade dos dados seja respeitada e o objeto final seja uma série temporal com dadtas compatíveis com o modo de operação interna do R. Modelos como HoltWinter sazonal exigem que exista pelo menos dois ciclos para executar ajustes, por exemplo.
}
\value{
Lista contendo várias estatísticas e dados dos procedimentos da função, como dados históricos com tratamento de outliers e missing dentre outras informações.
  \item{dados_historicos}{Dados históricos de entrada tratados}
  \item{historico_novo_ts }{Dados históricos em forma de objeto ts tratada}
  \item{data_previsao}{Datas de cada ponto projetado}
  \item{n_passos_frente}{Numero de passos à frente, passagem de argumento}
  \item{frequencia}{Frequencia dos dados, passagem de argumento}
}
\references{
Troyanskaya O., Cantor M., Sherlock G., Brown P., Hastie T., Tibshirani R., Botstein D., Altman R.B., (2001); Missing Value Estimation Methods for DNA microarrays Bioinformatics 17, 520-525.
}
\author{
LOPES, J. E.
}

\seealso{
Para mais informações sobre como lidar com outlier, missing e tratamento de datas e séries temporais, veja, \code{\link{interpNA}}, \code{\link{ISOdate}}, e \code{\link{tsoutliers}}
}
\examples{
data(mensal)
tsda <- ts.dados(mensal[,1:2], 10, 12)
tsda$dados_historicos
tsda$historico_novo_ts
}
