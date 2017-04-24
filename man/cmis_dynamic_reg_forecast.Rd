\name{cmis_dynamic_reg_forecast}
\alias{cmis_dynamic_reg_forecast}
\title{
Ajusta regressao dinamica com CrossValidation via ARIMA automaticamente
}
\description{
Realiza ajuste de modelos de regressao dinamica utilizando ARIMA com reducaoo de variáveis via Elastic Net, Wald ou STEPWISE com escolha de modelos por CorssValidation de forma automatica.
}
\usage{
cmis_dynamic_reg_forecast(dados, vresposta, explicativas,  n_passos_frente, freq,
				nfolds = 10, nivel=95, p_aumento=NULL, trace = FALSE, ...)
}

\arguments{
  \item{dados}{
  Objeto mts ou data.frame com as variaveis de interesse
  }
  \item{vresposta}{
  Variavel resposta
  }
  \item{explicativas}{
  Vetor de variáveis explicativas
  }
  \item{n_passos_frente}{
  Horizonte de forecast, quantos pontos à frente deseja estimar
  }
  \item{nfolds}{
  Número máximo de reamostragens para redução elastic net.
  }  
  \item{sig}{
  Valor mínimo (p.valor) para eliminação de variáveis por Wald ou STEPWISE
  }
  
  \item{nivel}{
  Nível de significância para os limites de predição do forecast. Padrão é 0.95
  }
    \item{p_aumento}{
  Permite simular cenários de crescimento para as variáveis explicativas e seu efeito sovre a resposta. Aceita valores decimais. Exemplo: 50 por cento de aumento implica em 0.50
  }
  
  \item{trace}{Se TRUE, exibe informações do processo de estimação}

  \item{automatic}{Se TRUE, a escolha das covariáveis é feita automaticamente. Caso contrário, tenta ajustar um modelo saturado com todas as covariáveis de entrada. Padrão é TRUE.}
  
  \item{freq}{
  Frequencia da série. Diário - 7, mensal - 12 e horário - 24 pontos.
  }
  \item{nvif}{
  Fator de Inflação de Variância, utilizado para verificar multicolinearidade.
  }
  
  \item{...}{
  Passagem de argumentos.}
}

\details{
Esta função permite ajustes de regressão dinâmica utilizando ARIMA que por sua vez utiliza a ténica de Box-Jenkins para ajustes e forecast. Esta função é baseada na função \code{auto.arima} e antes do ajuste é feita uma eliminação automatizada de covariáveis fracas pela técnica de CrossValidation através da metodologia Elastic Net, que por sua vez utiliza algoritimos de reduão via LASSO e Ridge regression. O resultado é um forecast coerente e assertivo para os dados. É permitido ainda simulação de cenários de crescimento para as covariáveis adicionando percentual de crescimentos desejados a fim de medir este impacto na explicativa.

}
\value{
A saída é uma lista aninhada com várias entradas conforme abaixo.
\item{models}{Objeto dos dois modelos escolhidos automaticamente, sendo um simples e outro com o efeito das covariáveis e do aumento, caso seja definido}
\item{forecasts}{Forecast dos dois modelos escolhidos automaticamente, sendo um simples e outro com o efeito das covariáveis e do aumento, caso seja definido}
\item{acuracia}{Estatísticas dos modelos}
}
\references{
Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast package for R", Journal of Statistical Software, 26(3).


}
\author{
LOPES, J. E.
}

\seealso{
Mais informações sobre a metodologia de forecast para modelos ARIMA \code{\link{auto.arima}} e para redução por Elastic Net \code{\link{cv.glmnet}}.
\code{\link{model.select}} para redução Stepwise e \code{\link{fastbw}} para eliminação por Wald
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
  
  fit <- cmis_dynamic_reg_forecast(da, vresposta[1], explicativas, 20, 7, 
				nfolds = 10, p_aumento = 1, nivel = 95, trace = TRUE)
  plot(fit)
  summary(fit)
  
  ## Exemplo de looping varrendo uma base de dados.
  ## Pode demorar muito!!! ##
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
  
  out <- c()
  for (i in 1:length(vresposta[1:10])) {
    temp <- try(cmis_dynamic_reg_forecast(da, vresposta[i], explicativas, 20, 12,
				nfolds = 5, p_aumento = 0.05, nivel = 95, trace = TRUE, sig = 0.10))
    if(class(temp)!="try-error") out[[i]] <- temp
    try(plot(temp))
  }
}
