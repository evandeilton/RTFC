\name{ts.dados}
\alias{ts.dados}
\title{
Trata dados para a serie temporal
}
\description{
Trata dados para a s�rie temporal com datas no formato adequado para ajustes
}
\usage{
ts.dados(dados_hist, n_passos_frente, freq, normalize = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados_hist}{
Vetor de dados hist�ricos no formato data.frame(data, realizado)
}
  \item{n_passos_frente}{
N�mero de passos utilizados na proje��o (quantos pontos ser�o projetados)
}
  \item{freq}{
Frequ�ncia dos dados: Mensal (12), di�rio(7), hor�rio(24)
}
  \item{normalize}{
Se TRUE faz normaliza��o dos dados. Vide \code{\link{scale}}
}
  \item{\dots}{
Passagem de par�metros.
}
}
\details{
Esta fun��o serve para transformar o vetor de dados brutos vindos das bases do CMIS de modo que a sazonalidade dos dados seja respeitada e o objeto final seja uma s�rie temporal com dadtas compat�veis com o modo de opera��o interna do R. Modelos como HoltWinter sazonal exigem que exista pelo menos dois ciclos para executar ajustes, por exemplo.
}
\value{
Lista contendo v�rias estat�sticas e dados dos procedimentos da fun��o, como dados hist�ricos com tratamento de outliers e missing dentre outras informa��es.
  \item{dados_historicos}{Dados hist�ricos de entrada tratados}
  \item{historico_novo_ts }{Dados hist�ricos em forma de objeto ts tratada}
  \item{data_previsao}{Datas de cada ponto projetado}
  \item{n_passos_frente}{Numero de passos � frente, passagem de argumento}
  \item{frequencia}{Frequencia dos dados, passagem de argumento}
}
\references{
Troyanskaya O., Cantor M., Sherlock G., Brown P., Hastie T., Tibshirani R., Botstein D., Altman R.B., (2001); Missing Value Estimation Methods for DNA microarrays Bioinformatics 17, 520-525.
}
\author{
LOPES, J. E.
}

\seealso{
Para mais informa��es sobre como lidar com outlier, missing e tratamento de datas e s�ries temporais, veja, \code{\link{interpNA}}, \code{\link{ISOdate}}, e \code{\link{tsoutliers}}
}
\examples{
data(mensal)
tsda <- ts.dados(mensal[,1:2], 10, 12)
tsda$dados_historicos
tsda$historico_novo_ts
}
