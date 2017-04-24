\name{cmis_corr}
\alias{cmis_corr}
\title{
Analise de relacao entre pares de variaveis
}
\description{
Faz analise de relacao entre todos os pares de combinações de variáveis de um data.frame 
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
Os methodos de calculo de correlação são: "pearson" (padrão), "spearman" ou "kendall"}
  \item{digits}{Número de digitos da saida}
  \item{verbose}{Se TRUE, exibe no console o log das operações}  
  \item{...}{Passagem de argumentos.}
}



\details{
Esta função calcula correlações do tipo Pearson, Spearman e Kendal, para um conjunto de dados ou para duas matrizes pareadas de dados e exibe o resultado em formado de seis colunas, sendo duas para os pares de nomes das variáveis de entrada e o restante para as estatísticas entre elas. Existe um argumento chamado use na função \code{cor} que trata a questão de valores missing, o padrão definido aqui é correlação apenas para "complete.cases", mas outras opções podem ser escolhidas (vide help \code{cor}) Além de correlação simples esta função calcula também o R-Quadrado,R-Quadrado ajustado e estatística PRESS para a relação entre os pares de variáveis. Para tanto, a função \code{lm_press_stat} é invocada internamente. O data.frame de saída está ordenado pelas maiores correlações.

Internamente a função efetura n(n-1)/2 regressões simples, onde n é a quantidade de variáveis no conjunto informado
}
\value{
Data frame com estatísticas
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