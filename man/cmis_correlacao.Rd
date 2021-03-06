\name{cmis_correlacao}
\alias{cmis_correlacao}
\title{
Faz analise Regressao/Correlacao automaticamente
}
\description{
Faz an�lise Regress�o/Correla��o entre m�tricas explicativas e explicadas automaticamente pelo m�todo do melhor conjunto de modelos, por�m utilizando an�lises de multicolinearidade, redu��o por CrossValidation e LASSO para selecionar as melhores vari�veis do conjunto de dados.
}
\usage{
cmis_correlacao(dados, n_passos_frente = 1, freq, nvif = 100, nfolds = 5, nmod = 10,
	r2 = 0.5, residteste = 0.8, sig = 0.15, prefixo.expl = "neg", prefixo.resp = "hdw",
	verbose = FALSE, saida_geral = FALSE, regsimples = TRUE, ...)
}
\arguments{
  \item{dados}{
Conjunto de dados com a primeira sendo a DATA no formato DD/MM/AAAA HH:MM:SS e outras vari�veis de neg�cio, infra, etc. m�nimo duas covari�veis
}
  \item{n_passos_frente}{
N�mero que indica o orizonte de proje��o.
M�trica mensal, a proje��o m�xima � de 18 meses (3 anos);
M�trica di�ria, a proje��o m�xima � de 180 dias (6 meses)
M�trica � por hora, a proje��o m�xima � de 336 dias (14 dias).

}
  \item{freq}{
Frequ�ncia da m�trica, di�ria � 7, mensal � 12 e hor�ria � 24.
}
  \item{nvif}{
Para modelos de regress�o m�ltipla, faz a escolhas das melhores m�tricas utilizando o crit�rio de multicolinearidade. Quanto menor for melhor. A literatura pede 10 como �timo, por�m nos dados de s�ries temporais certa flexibilidade foi imposta, deixando fixo em 300 para evitar perda de boas covari�veis j� que outras abordagem � realizada em conjunto na redu��o e escolha das m�tricas.
}
  \item{nfolds}{
Representa o n�mero de reamostragens a ser realizada para cada modelo a fim de decidir sobre a capacidade de predi��o do mesmo. Padr�o � 5 reamostragens.
}
  \item{nmod}{
Em modelos de regress�o m�ltipa indica o total de modelos a serem retornados por m�trica. Padr�o � 5.
}
  \item{r2}{
� o threshold padr�o para o coeficiente de deternina��o. 0.5 � o m�nimo padr�o. Quanto maior melhor.
}
  \item{residteste}{
Threshold para decis�o de escolha baseada na qualidade da an�lise de res�duos. Seis testes s�o aplicados e um indice de 1 a 10 � definido. Tomando-se o inverso deste n�mero, quanto mais pr�ximo de 0,10 melhor. 
}
  \item{sig}{
Alpha de signific�ncia para escolha dos betas estimados nas regress�es via teste estat�stico de Wald. O padr�o � 0,10. 
}
  \item{prefixo.expl}{
� o prefixo das vari�veis explicativas. Serve para qualquer combina��o de tr�s letras desde que as m�tricas do data frame de entrada sejam renomeadas, com exce��o da m�tricas de data, que deve ficar sempre com este nome. Padr�o � "neg", neg�cio.
}
  \item{prefixo.resp}{
� o prefixo das vari�veis resposta. Serve para qualquer combina��o de tr�s letras desde que as m�tricas do data frame de entrada sejam renomeadas, com exce��o da m�tricas de data, que deve ficar sempre com este nome. Padr�o � "hdw", infraestrutura.
}
  \item{verbose}{
Se TRUE exibe o progresso de escolha dos melhores modelos na tela. O tempo de processamento aumenta devido �s impress�es.
}
  \item{saida_geral}{
Se TRUE retorna uma lista aninhada, onde o primeiro elemento � o data.frame dos resultados no formato para o banco de dados SQL e o outro � uma lista dos modelos escolhidos no formato do R que podem ser manipulados com rotinas de loop do tipo sapply, ou lappy para an�lises diversas. Padr�o � FALSE.
}
  \item{regsimples}{
Se TRUE faz modelagem apenas para regress�o simples do tipo y = beta_0 + beta_1*X, caso contr�rio faz estima��o de todas as combina��es de modelos de regress�o linear sem itera��o, inclu�ndo simples e multipla. Padr�o � FALSE.
}
  \item{\dots}{
Argumentos opcionais para chamadas internas
}
}
\details{
O processo de modelagem pela fun��o \code{cmis_correlacao} envolve muita estat�stica e muitos detalhes te�ricos e computacionais. O tratamento inicial dos dados � feito pelas mesmas fun��es utilizadas para a parte do \code{cmis_forecast} que envolve tratamento de outliers, missing, m�tricas pobres de dados dentre outras. Em rela��o � modelagem, conceitos como redu��o de vari�veis via ANOVA tipo II, an�lise de CrossValidation e redu��o pela metodologia LASSO s�o empregadas em busca da melhor combina��o de t�cnicas para ajustes de modelos de forma mais inteligente evitando regress�es esp�rias e modelos estat�sticamente fracos. Em termos de decis�o, entram em cena as estat�sticas F, coeficiante de determina��o, estat�stica PRESS, ranking de res�duos e an�lise de AIC.

Os modelos abordados na an�lise de regress�o est�o divididos em duas classes: Modelos de Regress�o Linear Simples (MRLS) e Modelos de Regress�o Linear M�ltiplos (MRLM).

Quando o argumento de regress�o simples � verdadeiro, todas as possibilidades de modelos de regress�o simples para os dados s�o testadas e n�o h� redu��o de vari�veis, neste caso os modelos s�o checados por reamostragens internas via CrossValidation. Se os modelos s�o bons segundo a qualidade de seus res�duos e do coeficiente de determina��o, passa-se para a etapa de an�lise dos coeficientes estiamdos via estat�stcia Wald. No caso de modelos m�ltiplos uma verifica��o pr�via � feita nas covari�veis em busca de problemas de multicolinearidade que podem levar a modelos mal estimados e/ou regress�es esp�rias. Se alguma m�trica apresentar problemas de multicolinearidade ela � ent�o removida do modelo e o mesmo � reestimado e passa pela CrossValidation para um filtro mais sofisticado.

Em virtude de anomalias nos dados, os chamados pontos discrepantes e tamb�m valores faltantes (missing values), em muitas an�lises tais casos s�o tratados como exe��o e, em muitas situa��es, t�cnicas adequadas de tratamento devem ser utilizadas visando contornar sem deturpar o comportamento real dos dados. Diversas t�cnicas de detec��o de outliers s�o utilizadas na literatuta atual, contudo esta quest�o � delicada em fun��o do tipo de dado envolvido e de fatores inesperados como cat�strofes, por exemplo. A abordagem de tratamento de outliers no CMIS � feita em duas frentes que suporta dois tipos de situa��es suportadas pela abordagem de Hyndman.
 
I.	Dados n�o sazonais com ou sem missing: nestes casos interpola��o linear � empregada. Nesta abordagem, um modelo linear � ajustado com base no hist�rico e val�es estimados deste modelo s�o ent�o inseridos para completar a s�rie;

II.	Dados sazonais com ou sem missing: � feita uma decomposi��o da s�rie via stl(Seazonal Trend Loess ), que usa regress�o polinomial na s�rie dessazonalizada e em seguida resazonalisa a mesma com base na proje��o desta no per�do anterior. Embora complicada, esta t�cnica permite capturar comportamentos sazonais e evita a perda de informa��o da s�rie real.
}
\value{
\item{list[[1]]}{Um data.frame com as estat�sticas dos melhores modelos escolhidos pelo algoritimo no formato de tabela expandida para o servidor receber e levar para a tela do CMIS. As colunas do objeto s�o:}
\item{v.resposta}{Vari�vel resposta}
\item{explicativa}{Vari�vel (eis) explicativa (s)}
\item{coeficiente}{Coeficiente ou beta estimado}
\item{ic_2.5}{Intervalo de confian�a inferior de 0.05 de signific�ncia para o beta estimado}
\item{ic_97.5}{Intervalo de confian�a superior de 0.05 de signific�ncia para o beta estimado}
\item{p.valor}{valor p para cada beta estimado}
\item{aic}{Crit�rio de informa��o de Akaike do modelo}
\item{r2}{Coeficiente de determina��o do modelo}
\item{parcial.r2}{Coeficiente de determina��o parcial para modelos de regress�o linear multipla}
\item{rmse}{Ra�z do erro quadr�tico m�dio}
\item{resid.test}{Resultado do teste de res�duos}
\item{f.test}{Estat�stica F do modelo completo}
\item{press}{Estat�stica PRESS (PREdicted Residual Sum of Squares)}
\item{press.r2}{Coeficiente de determina��o do modelo segundo a estat�stica PRESS. Esta estat�stica serve para medir a capacidade preditiva do modelo para novas observa��es passadas no mesmo.}
\item{Id}{ID de ordena��o dos modelos, deve ser �nico por modelos, para evitar problemas de manipula��o dos modelos no servidor SQL.}
\item{list[[2]]}{Se saida_geral igual = TRUE, cont�m a lista completa dos melhores modelos escolhidos no formato de lista aninhada do R. Cada elemento da lista possui um ou mais objetos de modelos (lm), simples ou m�ltiplos e pode ser acessada para fins investigativos somente pelo R.}
}

\references{

R Core Team (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL  http://www.R-project.org/.

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

Alan Miller "Subset Selection in Regression" Chapman \& Hall

Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, http://www.stanford.edu/~hastie/Papers/glmnet.pdf

Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
http://www.jstatsoft.org/v33/i01/

Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5) 1-13
http://www.jstatsoft.org/v39/i05/

Buckland (1997) Model Selection: an Integral Part of Inference. Biometrics 10:41

Burnham & Anderson (2002) Model Selection and Multimodel Inference: an Information Theoretic Approach Calcagno \& de Mazancourt 2010 J. Stat. Soft. v34 i12. See http://www.jstatsoft.org/v34/i12

}
\author{
LOPES, J. L.
}

\examples{
data(mensal)
data(diario)

# Apenas regress�es simples com sa�da geral
fit0 <- cmis_correlacao(mensal, n_passos_frente = 1, freq = 12, nvif = 100, nfolds = 5, 
    nmod = 5, r2 = 0.5, residteste = 0.5, sig = 0.1, prefixo.expl = "neg", 
    prefixo.resp = "hdw", verbose = FALSE, saida_geral = TRUE, regsimples = TRUE) 
#fit0[[1]]
#fit0[[2]]

# Regress�es simples e multiplas com sa�da geral 
fit1 <- cmis_correlacao(diario, n_passos_frente = 1, freq = 7, nvif = 100, nfolds = 5, 
    nmod = 5, r2 = 0.5, residteste = 0.5, sig = 0.1, prefixo.expl = "neg", 
    prefixo.resp = "hdw", verbose = FALSE, saida_geral = TRUE, regsimples = FALSE)     
#fit1[[1]]
#fit1[[2]]    
}
