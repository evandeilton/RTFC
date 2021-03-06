\name{elastic_net_reduction}
\alias{elastic_net_reduction}
\title{
Realiza reducao de variaveis via LASSO e Ridge usando CrossValidation
}
\description{
Framework para redu��o de vari�veis via LASSO e Ridge com uso da metodologia Elastic Net para o processo de regress�o din�mica.
}
\usage{
elastic_net_reduction(dados, vresposta, explicativas, nfolds = 10, ...)
}

\arguments{
  \item{dados}{
  Objeto da classe data.frame ou mts para ajustes iniciais. Precisa conter ao menos duas variaveis.
  }
  \item{vresposta}{
  Vetor com o nome da variavel resposta
  }
  \item{explicativas}{
  Vetor com no minimo um nome de covariavel
  }
  \item{nfolds}{
  Numero que indica a quantidade de resmostragem no processo Elastic Net Reduction. O padr�o � 10 reamostragens para uma boa CrossValidation.
  }
  \item{...}{
  Passagem de par�metros
  }  
}
\details{
Esta fun��o � capaz de eliminar variaveis com poss�veis problemas de multicolinearidade em modelos multiplos, bem como vari�veis fracas estat�sticamente na predi��o da vari�vel resposta dentro do processo de regress�o din�mica.
}
\value{
A sa�da � uma lista com informa��es que ser�o repassadas para o processo din�mico de regress�o.
}
\references{
Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast package for R", Journal of Statistical Software, 26(3).

Hyndman, R.J., Koehler, A.B., Snyder, R.D., and Grose, S. (2002) "A state space framework for automatic forecasting using exponential smoothing methods", International J. Forecasting, 18(3), 439-454.

Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, http://www.stanford.edu/~hastie/Papers/glmnet.pdf
Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
http://www.jstatsoft.org/v33/i01/

Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5) 1-13
http://www.jstatsoft.org/v39/i05/

}
\author{
LOPES, J. E.
}

\seealso{
Mais informa��es sobre a metodologia de forecast para os modelos ARIMA podem ser consultadas na fun��o \code{\link{auto.arima}} do pacote \code{forecast} e tamb�m \code{\link{cv.glmnet}}
}
\examples{

data(diario)
da <- dataframe.outlier(diario, n_passos_frente=10, freq=7)$dados_ts
nm <- dimnames(da)[[2]]
vresposta <- explicativas <- c()

for(i in 1:length(nm)){
	if (substr(nm[i], 1, 3) == "NEG") explicativas[i] <- nm[i]
	if (substr(nm[i], 1, 3) == "HDW") vresposta[i] <- nm[i]
}

vresposta <- as.character(na.omit(vresposta))
explicativas <- as.character(na.omit(explicativas))
elnet <- elastic_net_reduction(da, vresposta[1], explicativas)

}
