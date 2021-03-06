\name{Mresid}
\alias{Mresid}
\title{
Diversos testes de residuos
}

\usage{
Mresid(forecast)
}

\arguments{
  \item{forecast}{
Objeto forecast
}
}
\details{
S�o efetuados seis testes nos res�duos sendo eles: teste de independ�ncia de Box-Pierce e Ljung-Box, teste de nulidade da m�dia de t-Student, teste de ruido branco de Portmanteau, teste de normalidade de Jarque Bera, teste de heterocedasticidade de Breusch-Pagan e por fim, teste de autocorrela��o de Durbin-Watson. Um indicador de qualidade do res�duo � gerado ao final com pesos maiores para independ�ncia, heterocedasticidade e autocorrela��o dos res�duos.
}
\value{
Resultados dos testes.
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

\seealso{
\code{\link{Box.test}}, \code{\link{t.test}}, \code{\link{LB.test}}
}

\examples{
## Not run
## dados
#data(diario)
#y <- ConvertDataToTs(diario[,1:2], tsfrequency = "day", OutType = "ts")
#fit <- auto.arima(y)
#Mresid(fit)
}
