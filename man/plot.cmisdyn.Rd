\name{plot.cmisdyn}
\alias{plot.cmisdyn}
\title{
Gera grafico do ajuste de regressao dinamica
}
\description{
Recebe um objeto do ajuste de regressão dinâmica e plota as estatísticas do ajuste e também da projeção antes e depois para a variável resposta.
}
\usage{
plot.cmisdyn(obj, ...)
}

\arguments{
  \item{obj}{
  Objeto da classe cmisdyn
  }
}

\author{
LOPES, J. E.
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

}
