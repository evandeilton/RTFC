\name{ordem_diferencas}
\alias{ordem_diferencas}
\title{
Estima a ordem de diferencas
}
\description{
Estima o n�mero de diferen�as sazonais ou n�o que uma s�rie precisa para se tornar estacion�ria, utilizando ADF, KPSS e outros testes
}
\usage{
ordem_diferencas(x, alpha = 0.05, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vetor de dados
}
  \item{alpha}{
Signific�ncia dos testes
}
  \item{plot}{
Plota o gr�fico da an�lise
}
  \item{\dots}{
Passagem de par�metros
}
}
\details{
Esta fun��o utiliza testes de ra�z unit�ria para verificar o n�mero de diferen�as at� que a s�rie se torne estacion�ria. Entre os testes estat�sticos est�o do ADF (Augmented Dickey-Fuller), Phillips-Perron, KPSS, Canova-Hansen e Osborn-Chui-Smith-Birchenhall.
}
\value{
Lista contendo s�rie original, s�rie diferenciada, quando existir e o n�mero de diferen�as estimadas.
}
\references{
Canova F and Hansen BE (1995) "Are Seasonal Patterns Constant over Time? A Test for Seasonal Stability", Journal of Business and Economic Statistics 13(3):237-252.

Dickey DA and Fuller WA (1979), "Distribution of the Estimators for Autoregressive Time Series with a Unit Root", Journal of the American Statistical Association 74:427-431.

Kwiatkowski D, Phillips PCB, Schmidt P and Shin Y (1992) "Testing the Null Hypothesis of Stationarity against the Alternative of a Unit Root", Journal of Econometrics 54:159-178.

Osborn DR, Chui APL, Smith J, and Birchenhall CR (1988) "Seasonality and the order of integration for consumption", Oxford Bulletin of Economics and Statistics 50(4):361-377.

Osborn, D.R. (1990) "A survey of seasonality in UK macroeconomic variables", International Journal of Forecasting, 6:327-336.

Said E and Dickey DA (1984), "Testing for Unit Roots in Autoregressive Moving Average Models of Unknown Order", Biometrika 71:599-607.
}
\author{
Rob J Hyndman and Slava Razbash (fun��es ndiffs e nsdiffs)

LOPES, J. L. (ordem_diferencas)
}

\seealso{
Esta fun��o utliza a fun��o \code{\link{ndiffs}} e \code{\link{nsdiffs}} como suporte para os testes. Mais detalhes podem ser vistos na ajuda.
}
\examples{
data(mensal)
ordem_diferencas(ts(mensal[, 45], frequency=2))
}