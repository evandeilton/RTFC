\name{diario}
\alias{Dados diarios}
\docType{data}
\title{
Conjunto de dados
}
\description{
Este � um conjunto de dados de m�tricas de IT inclu�ndo neg�cio e infraestrutura distribuidas entre, percentual de uso de CPU, de Mem�ria f�sica, redes, filesystem entre outras. Os nomes das m�tricas est�o configurados do modo que as dun��es \code{cmis_forecast} e \code{cmis_correlaco} conseguem trabalhar. Nota: apar trabalhar com estas duas fun��es citadas � mandat�rio entrar com os dados no mesmo estilo deste conjunto de dados aqui mostrado.
}
\usage{data("diario")}
\source{
Base de dados do CMIS aplicativo Siebel8, vis�o diaria
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
