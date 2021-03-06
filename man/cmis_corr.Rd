\name{cmis_corr}
\alias{cmis_corr}
\title{
Analise de relacao entre pares de variaveis
}
\description{
Faz analise de relacao entre todos os pares de combina��es de vari�veis de um data.frame 
}
\usage{
cmis_corr(tabela_y, tabela_x=NULL, method = "pearson", digits = 4, verbose = FALSE)
}
\arguments{
  \item{tabela_y}{
 data.frame ou matriz de dados numero 1}
  \item{tabela_x}{
 data.frame ou matriz de dados numero 2}
  \item{method}{
Os methodos de calculo de correla��o s�o: "pearson" (padr�o), "spearman" ou "kendall"}
  \item{digits}{N�mero de digitos da saida}
  \item{verbose}{Se TRUE, exibe no console o log das opera��es}  
  \item{...}{Passagem de argumentos.}
}



\details{
Esta fun��o calcula correla��es do tipo Pearson, Spearman e Kendal, para um conjunto de dados ou para duas matrizes pareadas de dados e exibe o resultado em formado de seis colunas, sendo duas para os pares de nomes das vari�veis de entrada e o restante para as estat�sticas entre elas. Existe um argumento chamado use na fun��o \code{cor} que trata a quest�o de valores missing, o padr�o definido aqui � correla��o apenas para "complete.cases", mas outras op��es podem ser escolhidas (vide help \code{cor}) Al�m de correla��o simples esta fun��o calcula tamb�m o R-Quadrado,R-Quadrado ajustado e estat�stica PRESS para a rela��o entre os pares de vari�veis. Para tanto, a fun��o \code{lm_press_stat} � invocada internamente. O data.frame de sa�da est� ordenado pelas maiores correla��es.

Internamente a fun��o efetura n(n-1)/2 regress�es simples, onde n � a quantidade de vari�veis no conjunto informado
}
\value{
Data frame com estat�sticas
}
\author{
LOPES, J. E.
}

\seealso{
\code{\link{lm_press_stat}}
}
\examples{
data(mensal)
cmis_corr(mensal[, 2:5])
cmis_corr(mensal[, 2:3], mensal[, 4:10], method="pearson")
#cmis_corr(mensal[, 2:3], mensal[, 4:10], method="kendall")
#cmis_corr(mensal[, 2:3], mensal[, 4:10], method="spearman")
}