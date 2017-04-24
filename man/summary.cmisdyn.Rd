\name{summary.cmisdyn}
\alias{summary.cmisdyn}
\title{
Summariza regressao dinamica
}
\description{
Recebe objeto cmisdyn e plota sumario das estatísticas do ajuste incluíndo gráfico.
}
\usage{
summary.cmisdyn(obj, plot=TRUE, digits = 10)
}

\arguments{
  \item{obj}{
  Objeto da classe cmisdyn
  }
  \item{plot}{
  Se TRUE exibe o gráfico do ajuste
  }
  \item{digits}{
  Total de digititos da saida
  }
}

\author{
LOPES, J. E.
}

\seealso{
\code{\link{plot.cmisdyn}}
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
summary(fit)
}
