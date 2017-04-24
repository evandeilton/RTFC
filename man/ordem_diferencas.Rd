\name{ordem_diferencas}
\alias{ordem_diferencas}
\title{
Estima a ordem de diferencas
}
\description{
Estima o número de diferenças sazonais ou não que uma série precisa para se tornar estacionária, utilizando ADF, KPSS e outros testes
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
Significância dos testes
}
  \item{plot}{
Plota o gráfico da análise
}
  \item{\dots}{
Passagem de parâmetros
}
}
\details{
Esta função utiliza testes de raíz unitária para verificar o número de diferenças até que a série se torne estacionária. Entre os testes estatísticos estão do ADF (Augmented Dickey-Fuller), Phillips-Perron, KPSS, Canova-Hansen e Osborn-Chui-Smith-Birchenhall.
}
\value{
Lista contendo série original, série diferenciada, quando existir e o número de diferenças estimadas.
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
Rob J Hyndman and Slava Razbash (funções ndiffs e nsdiffs)

LOPES, J. L. (ordem_diferencas)
}

\seealso{
Esta função utliza a função \code{\link{ndiffs}} e \code{\link{nsdiffs}} como suporte para os testes. Mais detalhes podem ser vistos na ajuda.
}
\examples{
data(mensal)
ordem_diferencas(ts(mensal[, 45], frequency=2))
}