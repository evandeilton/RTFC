\name{my_glm_vif}
\alias{my_glm_vif}
\title{
Analise de VIF (Variance Infraction Factor)
}
\description{
Realiza análise de VIF com intuito de estudar a multicolinearidade das variáveis explicativas em relação à explicada. Internamente são ajustados modelos lineares simples ou múltiplos para a análise.
}
\usage{
my_glm_vif(dados, vresposta, explicativas, nvif, regsimples = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados}{
data.frame no formato data.frame(data, var1, var2, var3, ..., varn)
}
  \item{vresposta}{
Variável reposta cujas explicativas terão análise de VIF
}
  \item{explicativas}{
Vetor de variáveis explicativas
}
  \item{nvif}{
Threshold de exclusão da variável explicativa pelo VIF da mesma. Ex.: VIF < 50 => exclusão
}
  \item{regsimples}{
Se TRUE, não realiza análise de multicolinearidade, apenas configura a saída para um formado compatível com outras funções internas do CMIS. Caso contrário análise de VIF é feita utilizando regressão linear multipla.
}
  \item{\dots}{
Passagem de argumentos
}
}
\details{
Nas métricas dos dados de modelagem, esta estatística é utilizada para eliminar dos modelos de regrssão linera multiplos a possibilidade de multicolinearidade que pode gerar estimativas falsas. Apenas as métricas que passarem neste teste segundo o threshold definido, seguirão em análise, as demais deixam de ser modeladas.
}
\value{
Lista contendo
  \item{dados}{Conjunto de dados do ajuste}
  \item{formula}{Formula final do modelo, reduzida ou não}
  \item{vif}{Estatística VIF}
  \item{vresposta}{Variáveis resposta}
  \item{vexplicativas}{Variáveis explicativas}
}
\references{
Davis CE, Hyde JE, Bangdiwala SI, Nelson JJ: An example of dependencies among variables in a conditional logistic regression. In Modern Statistical Methods in Chronic Disease Epidemiology, Eds SH Moolgavkar and RL Prentice, pp. 140-147. New York: Wiley; 1986.
}
\author{
LOPES, J. E.
}
\note{
Apenas modelos com mais de uma covariável deve passar por análise de VIF, pois combinações lineares similares podem ser encontradas entre variáveis explicativas causando perda de força do ajuste ou regressões espúrias.
}

\seealso{
Esta função depende da função \code{\link{vif}} deste pacote.
}

\examples{
data(mensal)
nomes <- names(mensal)
vresposta <- as.character(nomes[20])
explicativas <- as.character(nomes[3:10])
my_glm_vif(mensal, vresposta, explicativas, 100, regsimples=TRUE)
#my_glm_vif(mensal, vresposta, explicativas, 100, regsimples=FALSE)
}