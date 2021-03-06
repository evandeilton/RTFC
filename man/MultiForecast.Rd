\name{MultiForecast}
\alias{MultiForecast}
\title{
Faz forecast com diversos metodos
}
\description{
Em analise...
}
\usage{
MultiForecast(x, Control = cmisControl(), fcMethod=NULL)
}

\arguments{
  \item{x}{
S�rie temporal
}
  \item{fcMethod}{
M�todo de forecast definida pelo utilizador. Se NULL, a escolha � feita automaticamente.
}

  \item{Control}{
Argumentos para a an�lise
}
}
\details{
Esta fun��o utiliza v�rios m�todos de forecast implementados no pacote \code{forecast} para fazer modelagem de s�rie temporais. Foi criado para cada m�todo de forecast um wraper que permite ser chamado pela fun��o gen�rica \code{MultiForecast}. Estes wrapers s�o:
	
	- stsForecast: Modelos estruturais
	
	- hwForecast: Modelos HoltWinters
	
	- tbatsForecast: Modelos TBATS
	
	- auto.arimaForecast: Modelos ARIMA autom�ticos
	
	- sesForecast: Modelos de suaviza��o exponencial
	
	- meanForecast: Modelos de m�dias simples
	
	- holtForecast: Modelos Holt
	
	- batsForecast: Modelos BATS
	
	- etsForecast: Modelos de estados de espa�o (State Space)
	
	- arimaForecast: Modelos ARIMA fixos
	
	- lmForecast: Modelos lineares
	
	- thetaForecast: Modelos theta
	
	- rwForecast: Modelos random walk (Passeio aleat�rio)
	
	- snaiveForecast: Modelos naive sazonais
	
	- naiveForecast: Modelos naive simples
	
	- nnetarForecast: Modelos de redes neurais (N�o possui limites de predi��o)
	
	- HWsForecast: Modelos HoltWinters sazonais
	
	- HWnsForecast: Modelos HoltWinters n�o sazonais
	
	- HWesForecast: Modelos HoltWinters com alizamento exponencial
	
No caso dos modelos ETS e ARIMA autom�ticos, a escolha dos melhores ajustes � feita pela estat�stica AIC, para os outros tipos de modelos o ajuste � feito sobre os dados sem crit�rio de bondade para a rodada inicial de ajustes. A escolha final dos melhores modelos � feita automaticamente utilizando a estat�stica definida no Controle. O padr�o � MAPE (Erro M�dio Absoluto Percentual). 

Em complemento o analista pode observar as estat�sticas de acur�cia e de qualidade do res�duos utilizando func��es especiais como \code{Acuracia} e \code{Mresid} no caso de uma an�lise fina individualmente, pois o objeto de sa�da cont�m os modelos escolhidos e permite acesso.
}
\value{
Lista com os melhores modelos de forecast escolhidos.
}
\references{
Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast package for R", Journal of Statistical Software, 26(3).
}
\author{
LOPES, J. E.
}

\seealso{
\code{\link{switch.cvforecast}}
}

\examples{
## Dados
data(diario)

## Controle
Control <- cmisControl(
  maxHorizon = 45,
  level = 95,
  onlyfc = FALSE,
  cvMethod = "MAPE",
  tsfrequency = "day",
  outputFormat = "forecast"
)

fit <- MultiForecast(diario[,2], Control)
# Analise de residuos (p.valor dos testes)
 sapply(fit, Mresid)
# Estat�stica dos modelos
 sapply(fit, tsSummary)
}
