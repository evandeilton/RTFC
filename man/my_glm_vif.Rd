\name{my_glm_vif}
\alias{my_glm_vif}
\title{
Analise de VIF (Variance Infraction Factor)
}
\description{
Realiza an�lise de VIF com intuito de estudar a multicolinearidade das vari�veis explicativas em rela��o � explicada. Internamente s�o ajustados modelos lineares simples ou m�ltiplos para a an�lise.
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
Vari�vel reposta cujas explicativas ter�o an�lise de VIF
}
  \item{explicativas}{
Vetor de vari�veis explicativas
}
  \item{nvif}{
Threshold de exclus�o da vari�vel explicativa pelo VIF da mesma. Ex.: VIF < 50 => exclus�o
}
  \item{regsimples}{
Se TRUE, n�o realiza an�lise de multicolinearidade, apenas configura a sa�da para um formado compat�vel com outras fun��es internas do CMIS. Caso contr�rio an�lise de VIF � feita utilizando regress�o linear multipla.
}
  \item{\dots}{
Passagem de argumentos
}
}
\details{
Nas m�tricas dos dados de modelagem, esta estat�stica � utilizada para eliminar dos modelos de regrss�o linera multiplos a possibilidade de multicolinearidade que pode gerar estimativas falsas. Apenas as m�tricas que passarem neste teste segundo o threshold definido, seguir�o em an�lise, as demais deixam de ser modeladas.
}
\value{
Lista contendo
  \item{dados}{Conjunto de dados do ajuste}
  \item{formula}{Formula final do modelo, reduzida ou n�o}
  \item{vif}{Estat�stica VIF}
  \item{vresposta}{Vari�veis resposta}
  \item{vexplicativas}{Vari�veis explicativas}
}
\references{
Davis CE, Hyde JE, Bangdiwala SI, Nelson JJ: An example of dependencies among variables in a conditional logistic regression. In Modern Statistical Methods in Chronic Disease Epidemiology, Eds SH Moolgavkar and RL Prentice, pp. 140-147. New York: Wiley; 1986.
}
\author{
LOPES, J. E.
}
\note{
Apenas modelos com mais de uma covari�vel deve passar por an�lise de VIF, pois combina��es lineares similares podem ser encontradas entre vari�veis explicativas causando perda de for�a do ajuste ou regress�es esp�rias.
}

\seealso{
Esta fun��o depende da fun��o \code{\link{vif}} deste pacote.
}

\examples{
data(mensal)
nomes <- names(mensal)
vresposta <- as.character(nomes[20])
explicativas <- as.character(nomes[3:10])
my_glm_vif(mensal, vresposta, explicativas, 100, regsimples=TRUE)
#my_glm_vif(mensal, vresposta, explicativas, 100, regsimples=FALSE)
}