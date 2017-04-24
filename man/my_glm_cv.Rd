\name{my_glm_cv}
\alias{my_glm_cv}
\title{
Analise de Cross Validacao
}
\description{
Faz análise de Cross Validação utilizando redução pelas metodologias LASSO e regressão ridge chamada de elasticnet regularization
}
\usage{
my_glm_cv(dados, vresposta, explicativas, nfolds = 5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dados}{
data.frame no formato data.frame(data, var1, var2, var3, ..., varn)
}
  \item{vresposta}{
Variável reposta
}
  \item{explicativas}{
Vetor de variáveis explicativas
}
  \item{nfolds}{
Número de repartições do conjunto de dados para a CrossValidation. Recomenda-se = 10 reamostragens para amostras grandes. PAdrão é 5.
}
  \item{\dots}{
Passagem de argumentos
}
}
\details{
Consultar as referências abaixo ou o help de algumas funções do pacote \code{glmnet} na seção veja também.
}
\value{
Lista contendo
  \item{dados}{Conjunto de dados do ajuste}
  \item{formula}{Formula final do modelo, reduzida ou não}
  \item{lambdacv}{Estatística de CrossValidation}
  \item{vresposta}{Variáveis resposta}
  \item{vexplicativas}{Variáveis explicativas}
}
\references{
Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, http://www.stanford.edu/~hastie/Papers/glmnet.pdf

Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
http://www.jstatsoft.org/v33/i01/
}
\author{
LOPES, J. E.
}

\seealso{
Sobre a metodologia elastic net regularization \link{cv.glmnet} e \link{glmnet} do pacote \code{glmnet}
}
\examples{
data(diario)
nomes <- names(diario)
vresposta <- as.character(nomes[20])
explicativas <- as.character(nomes[3:10])
Mcv <- my_glm_cv(diario, vresposta, explicativas, nfolds=5)
#my_glm_cv(diario, vresposta, explicativas, nfolds=10)
}
