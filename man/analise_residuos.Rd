\name{analise_residuos}
\alias{analise_residuos}
\title{
Analise de residuos em series temporais
}
\description{
Faz análise de resíduo de modelos ARIMA, ETS e HoltWinters
}
\usage{
analise_residuos(residuos, historico_novo, modelo)
}

\arguments{
  \item{residuos}{Rresiduos do modelo a ser testado
}
  \item{historico_novo}{
Dados históricos tratados
}
  \item{modelo}{
Objeto ARIMA, ETS ou Holt do ajuste do modelo testado.
}
}
\details{
São efetuados seis testes nos resíduos sendo eles: teste de independência de Box-Pierce e Ljung-Box, teste de nulidade da média de t-Student, teste de ruido branco de Portmanteau, teste de normalidade de Jarque Bera, teste de heterocedasticidade de Breusch-Pagan e por fim, teste de autocorrelação de Durbin-Watson. Um indicador de qualidade do resíduo é gerado ao final com pesos maiores para independência, heterocedasticidade e autocorrelação dos resíduos.
}
\value{
\item{analise_residual}{Data frame com estatísticas dos testes realizados sobre os resíduos}
\item{soma_peso}{Soma dos pesos atribuidos a cada teste}
\item{soma_peso_inv}{Indicador de bondade da análise de resíduo, equivale a \code{1/(colSums(analise_residual)[1])}}
}
\references{
Box, G. E. P. and Pierce, D. A. (1970), Distribution of residual correlations in autoregressive-integrated moving average time series models. Journal of the American Statistical Association, 65, 1509-1526.

Ljung, G. M. and Box, G. E. P. (1978), On a measure of lack of fit in time series models. Biometrika 65, 297-303.

Harvey, A. C. (1993) Time Series Models. 2nd Edition, Harvester Wheatsheaf, NY, pp. 44, 45.

J. B. Cromwell, W. C. Labys and M. Terraza (1994): Univariate Tests for Time Series Models, Sage, Thousand Oaks, CA, pages 20-22.

T.S. Breusch & A.R. Pagan (1979), A Simple Test for Heteroscedasticity and Random Coefficient Variation. Econometrica 47, 1287-1294

J. Durbin & G.S. Watson (1971), Testing for Serial Correlation in Least Squares Regression III. Biometrika 58, 1-19.

J. Racine & R. Hyndman (2002), Using R To Teach Econometrics. Journal of Applied Econometrics 17, 175-189.

}
\author{
LOPES, J. E.
}

\examples{
## Not run
require(forecast)
fit <- auto.arima(AirPassengers)
resids <- residuals(fit)
historico <- fit$x
analise_residuos(resids, historico, fit)
}
