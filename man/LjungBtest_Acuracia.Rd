\name{LjungBtest_Acuracia}
\alias{LjungBtest_Acuracia}
\title{
Calcula as estatisticas de acuracia e bondade do forecast
}
\description{
Aplica o teste de autocorrelação de Box-Pierce e Ljung-Box e o teste de correlação serial de Durbin-Watson nos resíduos do modelo forecast. Calcula também as estatísticas de acurácia do modelo tais como MAPE, MAE, RMSE, etc.
}
\usage{
LjungBtest_Acuracia(model)
}

\arguments{
  \item{model}{
  Objeto do forecast
  }
  \item{...}{
  Passagem de parâmetros.
  }
}

\value{
Data.frame com as estatísticas calculadas

}
\references{
Box, G. E. P. and Pierce, D. A. (1970), Distribution of residual correlations in autoregressive-integrated moving average time series models. Journal of the American Statistical Association, 65, 1509-1526.

J. Durbin & G.S. Watson (1950), Testing for Serial Correlation in Least Squares Regression I. Biometrika 37, 409-428.

J. Durbin & G.S. Watson (1951), Testing for Serial Correlation in Least Squares Regression II. Biometrika 38, 159-178.

J. Durbin & G.S. Watson (1971), Testing for Serial Correlation in Least Squares Regression III. Biometrika 58, 1-19.

}
\author{
LOPES, J. E.
}

\seealso{
\code{\link{Box.test}}, \code{\link{dwtest} e \code{\link{accuracy}}}
}
\examples{
## Not run
#data(mensal)
#fit <- auto.arima(mensal[, 2])
#fore <- forecast(fit, h=12)
#LjungBtest_Acuracia(fore)
}
