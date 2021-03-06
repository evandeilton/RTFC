\name{cmisMultiforecast}
\alias{cmisMultiforecast}

\title{
Abordagem avancada para forecast por diveros metodos.
}
\description{
Al� de ARIMA, ETS e HoltWinters esta fun��o permite fazer foprecast por muitos outros m�todos de maneira inteligente utilizando testes de linearidade e de tend�ncia antes de aplicar a metodologia mais adequada de forecast. O resultado � um objeto com os melhores cinco ou menos modelos de onde se pode estrair estat�sticas de tamb�m gr�ficos.
}
\usage{
cmisMultiforecast(tsdata, Control = cmisControl(), ...)
}

\arguments{
  \item{tsdata}{
Data frame, ts, xts, matriz ou vetor de dados para forecast.
}
  \item{fcMethod}{
Escolha do m�todo de forecast pelo utilizador. Se NULL a escolha � autom�tica.
}
  \item{Control}{
Controle gen�rico para passagem de par�metros
}
  \item{\dots}{
Passagem de par�metros extras
}
}
\details{
Esta fun��o com apoio de \code{MultiForecast} utiliza v�rios m�todos de forecast implementados no pacote \code{forecast} para fazer modelagem de s�rie temporais. Foi criado para cada m�todo de forecast um wraper que permite ser chamado pela fun��o gen�rica \code{MultiForecast}. Estes wrapers s�o:

stsForecast = Modelos estruturais
hwForecast = Modelos HoltWinters
tbatsForecast = Modelos TBATS
auto.arimaForecast = Modelos ARIMA autom�ticos
sesForecast = Modelos de suaviza��o exponencial
meanForecast = Modelos de m�dias simples
holtForecast = Modelos Holt
batsForecast = Modelos BATS
etsForecast = Modelos de estados de espa�o (State Space)
arimaForecast  = Modelos ARIMA fixos
lmForecast = Modelos lineares
thetaForecast = Modelos theta
rwForecast = Modelos random walk (Passeio aleat�rio)
snaiveForecast = Modelos naive sazonais
naiveForecast = Modelos naive simples
nnetarForecast = Modelos de redes neurais (N�o possui limites de predi��o)
HWsForecast = Modelos HoltWinters sazonais
HWnsForecast = Modelos HoltWinters n�o sazonais
HWesForecast = Modelos HoltWinters com alizamento exponencial

No caso dos modelos ETS e ARIMA autom�ticos, a escolha dos melhores ajustes � feita pela estat�stica AIC, para os outros tipos de modelos o ajuste � feito sobre os dados sem crit�rio de bondade para a rodada inicial de ajustes. A escolha final dos melhores modelos � feita automaticamente utilizando a estat�stica definida no Controle. O padr�o � MAPE (Erro M�dio Absoluto Percentual). 

Em complemento o analista pode observar as estat�sticas de acur�cia e de qualidade do res�duos utilizando func��es especiais como \code{Acuracia} e \code{Mresid} no caso de uma an�lise fina individualmente, pois o objeto de sa�da cont�m os modelos escolhidos e permite acesso.
}
\value{
Lista simples ou aninhada com resultados de acordo com os argumentos de controle
}
\references{
R Core Team (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL http://www.R-project.org/.

Texts, Forecasting: principles and practice, < http://robjhyndman.com/talks/uwa/ >, acesso em 30/10/2014.

Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast package for R", Journal of Statistical Software, 26(3).

Hyndman, R.J., Koehler, A.B., Snyder, R.D., and Grose, S. (2002) "A state space framework for automatic forecasting using exponential smoothing methods",International J. Forecasting, 18(3), 439-454.

Hyndman, R.J., Akram, Md., and Archibald, B. (2008) "The admissible parameter space for exponential smoothing models". Annals of Statistical Mathematics, 60(2), 407-426.

Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008) Forecasting with exponential smoothing: the state space approach, Springer-Verlag.http://www.exponentialsmoothing.net.

Box, G. E. P., G. M. Jenkins and G. C. Reinsel (2008). Time series analysis: forecasting and control. 4th. Hoboken, NJ: John Wiley & Sons.

Brockwell, P. J. and R. A. Davis (2002). Introduction to time series and forecasting. 2nd ed. New York: Springer.

Chatfield, C. (2000). Time-series forecasting. Boca Raton: Chapman & Hall/CRC.

Pena, D., G.C. Tiao and R.S. Tsay, eds. (2001). A course in time series analysis. New York: John Wiley & Sons.

Shumway, R. H. and D. S. Stoffer (2011). Time series analysis and its applications: with R examples. 3rd ed. New York: Springer.
}
\author{
LOPES, J. E.
}

\seealso{
\code{\link{cmis_forecast}}, \code{\link{forecast}}, \code{\link{cmisControl}}
}
\examples{
# Not run
# Define controle para coleta e proje��es
#require(plyr)
#Control <- cmisControl(
#  maxHorizon = 45,
#  level = 95,
#  onlyfc = FALSE,
# cvMethod = "MAPE",
#  tsfrequency = "hour",
#  outputFormat = "forecast",
#  driverJDBC = #"C://Users//G0047743//Documents//sqldeveloper//jdbc//lib/ojdbc6.jar"
#)
#
#dd <- CMISForecastData(dados = NULL, 7, 1, 10000, connect = TRUE, #Control = Control)

## Fun�ao multiForecast (dados em forma de s�rie 'ts')
#f2 <- function(dados, ...){
#  llply(dados, function(X) {
#	X   <- X[,1:2]
#    fit <- Try_error(cmisMultiforecast(X, Control))
#    if (class(fit) != "try-error") fit else NULL
#  }, .progress = 'time')
#}

# Forecasts
#Forecasts <- f2(dd[1:5])
#plot(Forecasts[[1]])

# Estat�sticas
#ldply(Forecasts[[1]], Mresid)
#ldply(Forecasts[[1]], tsSummary)
}
