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
Série temporal
}
  \item{fcMethod}{
Método de forecast definida pelo utilizador. Se NULL, a escolha é feita automaticamente.
}

  \item{Control}{
Argumentos para a análise
}
}
\details{
Esta função utiliza vários métodos de forecast implementados no pacote \code{forecast} para fazer modelagem de série temporais. Foi criado para cada método de forecast um wraper que permite ser chamado pela função genérica \code{MultiForecast}. Estes wrapers são:
	
	- stsForecast: Modelos estruturais
	
	- hwForecast: Modelos HoltWinters
	
	- tbatsForecast: Modelos TBATS
	
	- auto.arimaForecast: Modelos ARIMA automáticos
	
	- sesForecast: Modelos de suavização exponencial
	
	- meanForecast: Modelos de médias simples
	
	- holtForecast: Modelos Holt
	
	- batsForecast: Modelos BATS
	
	- etsForecast: Modelos de estados de espaço (State Space)
	
	- arimaForecast: Modelos ARIMA fixos
	
	- lmForecast: Modelos lineares
	
	- thetaForecast: Modelos theta
	
	- rwForecast: Modelos random walk (Passeio aleatório)
	
	- snaiveForecast: Modelos naive sazonais
	
	- naiveForecast: Modelos naive simples
	
	- nnetarForecast: Modelos de redes neurais (Não possui limites de predição)
	
	- HWsForecast: Modelos HoltWinters sazonais
	
	- HWnsForecast: Modelos HoltWinters não sazonais
	
	- HWesForecast: Modelos HoltWinters com alizamento exponencial
	
No caso dos modelos ETS e ARIMA automáticos, a escolha dos melhores ajustes é feita pela estatística AIC, para os outros tipos de modelos o ajuste é feito sobre os dados sem critério de bondade para a rodada inicial de ajustes. A escolha final dos melhores modelos é feita automaticamente utilizando a estatística definida no Controle. O padrão é MAPE (Erro Médio Absoluto Percentual). 

Em complemento o analista pode observar as estatísticas de acurácia e de qualidade do resíduos utilizando funcções especiais como \code{Acuracia} e \code{Mresid} no caso de uma análise fina individualmente, pois o objeto de saída contém os modelos escolhidos e permite acesso.
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
# Estatística dos modelos
 sapply(fit, tsSummary)
}
