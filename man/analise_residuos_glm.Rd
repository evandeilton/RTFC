\name{analise_residuos_glm}
\alias{analise_residuos_glm}
\title{
Faz analise de residuo em objetos lm ou glm
}
\description{
� uma vers�o modificada da fun��o \code{analise_residuos} para modelos lm e glm
}
\usage{
analise_residuos_glm(glm_obj)
}

\arguments{
  \item{glm_obj}{
Objeto de ajuste do tipo lm ou glm
}
}
\details{
S�o efetuados seis testes nos res�duos sendo eles: teste de independ�ncia de Box-Pierce e Ljung-Box, teste de nulidade da m�dia de t-Student, teste de ruido branco de Portmanteau, teste de normalidade de Jarque Bera, teste de heterocedasticidade de Breusch-Pagan e por fim, teste de autocorrela��o de Durbin-Watson. Um indicador de qualidade do res�duo � gerado ao final com pesos maiores para independ�ncia, heterocedasticidade e autocorrela��o dos res�duos.
}
\value{
\item{vetor}{Vetor com indicador de an�lise dos res�duos}
}
\references{
Box, G. E. P. and Pierce, D. A. (1970), Distribution of residual correlations in autoregressive-integrated moving average time series models. Journal of the American Statistical Association, 65, 1509-1526.

Ljung, G. M. and Box, G. E. P. (1978), On a measure of lack of fit in time series models. Biometrika 65, 297-303.

Harvey, A. C. (1993) Time Series Models. 2nd Edition, Harvester Wheatsheaf, NY, pp. 44, 45.

J. B. Cromwell, W. C. Labys and M. Terraza (1994): Univariate Tests for Time Series Models, Sage, Thousand Oaks, CA, pages 20-22.

T.S. Breusch & A.R. Pagan (1979), A Simple Test for Heteroscedasticity and Random Coefficient Variation. Econometrica 47, 1287-1294

J. Durbin & G.S. Watson (1971), Testing for Serial Correlation in Least Squares Regression III. Biometrika 58, 1-19.

}
\author{
LOPES, J. L.
}
\note{
Indicador de qualidade da an�lise de residuos: quanto mais pr�ximo de 1 melhor.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(forecast)
data(diario)
fit <- lm(HDW_REDE_602~NEG_NDA_257, data = diario)
analise_residuos_glm(fit)
}
