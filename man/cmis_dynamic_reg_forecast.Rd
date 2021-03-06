\name{cmis_dynamic_reg_forecast}
\alias{cmis_dynamic_reg_forecast}
\title{
Ajusta regressao dinamica com CrossValidation via ARIMA automaticamente
}
\description{
Realiza ajuste de modelos de regressao dinamica utilizando ARIMA com reducaoo de vari�veis via Elastic Net, Wald ou STEPWISE com escolha de modelos por CorssValidation de forma automatica.
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
  Vetor de vari�veis explicativas
  }
  \item{n_passos_frente}{
  Horizonte de forecast, quantos pontos � frente deseja estimar
  }
  \item{nfolds}{
  N�mero m�ximo de reamostragens para redu��o elastic net.
  }  
  \item{sig}{
  Valor m�nimo (p.valor) para elimina��o de vari�veis por Wald ou STEPWISE
  }
  
  \item{nivel}{
  N�vel de signific�ncia para os limites de predi��o do forecast. Padr�o � 0.95
  }
    \item{p_aumento}{
  Permite simular cen�rios de crescimento para as vari�veis explicativas e seu efeito sovre a resposta. Aceita valores decimais. Exemplo: 50 por cento de aumento implica em 0.50
  }
  
  \item{trace}{Se TRUE, exibe informa��es do processo de estima��o}

  \item{automatic}{Se TRUE, a escolha das covari�veis � feita automaticamente. Caso contr�rio, tenta ajustar um modelo saturado com todas as covari�veis de entrada. Padr�o � TRUE.}
  
  \item{freq}{
  Frequencia da s�rie. Di�rio - 7, mensal - 12 e hor�rio - 24 pontos.
  }
  \item{nvif}{
  Fator de Infla��o de Vari�ncia, utilizado para verificar multicolinearidade.
  }
  
  \item{...}{
  Passagem de argumentos.}
}

\details{
Esta fun��o permite ajustes de regress�o din�mica utilizando ARIMA que por sua vez utiliza a t�nica de Box-Jenkins para ajustes e forecast. Esta fun��o � baseada na fun��o \code{auto.arima} e antes do ajuste � feita uma elimina��o automatizada de covari�veis fracas pela t�cnica de CrossValidation atrav�s da metodologia Elastic Net, que por sua vez utiliza algoritimos de redu�o via LASSO e Ridge regression. O resultado � um forecast coerente e assertivo para os dados. � permitido ainda simula��o de cen�rios de crescimento para as covari�veis adicionando percentual de crescimentos desejados a fim de medir este impacto na explicativa.

}
\value{
A sa�da � uma lista aninhada com v�rias entradas conforme abaixo.
\item{models}{Objeto dos dois modelos escolhidos automaticamente, sendo um simples e outro com o efeito das covari�veis e do aumento, caso seja definido}
\item{forecasts}{Forecast dos dois modelos escolhidos automaticamente, sendo um simples e outro com o efeito das covari�veis e do aumento, caso seja definido}
\item{acuracia}{Estat�sticas dos modelos}
}
\references{
Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast package for R", Journal of Statistical Software, 26(3).


}
\author{
LOPES, J. E.
}

\seealso{
Mais informa��es sobre a metodologia de forecast para modelos ARIMA \code{\link{auto.arima}} e para redu��o por Elastic Net \code{\link{cv.glmnet}}.
\code{\link{model.select}} para redu��o Stepwise e \code{\link{fastbw}} para elimina��o por Wald
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
