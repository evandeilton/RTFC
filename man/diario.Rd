\name{diario}
\alias{Dados diarios}
\docType{data}
\title{
Conjunto de dados
}
\description{
Este é um conjunto de dados de métricas de IT incluíndo negócio e infraestrutura distribuidas entre, percentual de uso de CPU, de Memória física, redes, filesystem entre outras. Os nomes das métricas estão configurados do modo que as dunções \code{cmis_forecast} e \code{cmis_correlaco} conseguem trabalhar. Nota: apar trabalhar com estas duas funções citadas é mandatório entrar com os dados no mesmo estilo deste conjunto de dados aqui mostrado.
}
\usage{data("diario")}
\source{
Base de dados do CMIS aplicativo Siebel8, visão diaria
}

\examples{
data(diario)
names(diario)
projecao <- cmis_forecast(diario[,c(1,3)], 50, 7, foreplot = TRUE, saida_geral = FALSE)
projecao
regressao <- cmis_correlacao(diario[, 1:20], n_passos_frente = 1, freq = 7, nvif = 100,
	nfolds = 5, nmod = 5, r2 = 0.5, residteste = 0.5, sig = 0.1, prefixo.expl = "neg",
	prefixo.resp = "hdw", verbose = FALSE, saida_geral = FALSE, regsimples = FALSE) 
regressao
}
