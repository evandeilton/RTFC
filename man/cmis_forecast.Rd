\name{cmis_forecast}
\alias{cmis_forecast}
\title{
Faz ajuste de modelos ARIMA, ETS e HoltWinters automaticamente
}
\description{
A fun��o recebe como argumentos de entrada um conjunto de dados com duas vari�veis: data e valor realizado conforme frequ�ncia, hor�ria, di�ria ou mensal e aplica tr�s tipos de metodologias de forecast em s�ries temporais, sendo eles ETS (Exponential smoothing state space model), HoltWinters (Holt-Winters Filtering) e ARIMA (Auto Regressive Integrated Moving Averege). Estas tr�s metodologias de forecast combinadas geram cerca de 34 tipos de modelos e o melhor deles � escolhido.
}
\usage{
cmis_forecast(dados, n_passos_frente, freq, foreplot = FALSE, id_metrica_frequencia = 0,
	primeiro_modelo = 0, parametro_1 = NA, parametro_2 = NA, parametro_3 = NA,
	parametro_4 = NA, parametro_5 = NA, parametro_6 = NA, parametro_7 = NA,
	modelo = 0, saida_geral = FALSE, ore_save = FALSE, ...)
}

\arguments{
  \item{dados}{
Conjunto de dados com duas vari�veis, a primeira sendo a DATA no formato DD/MM/AAAA HH:MM:SS e a outra contendo a m�trica realizada.
}
  \item{n_passos_frente}{
N�mero que indica o orizonte de proje��o.
M�trica mensal, a proje��o m�xima � de 18 meses (3 anos);
M�trica di�ria, a proje��o m�xima � de 180 dias (6 meses)
M�trica � por hora, a proje��o m�xima � de 336 dias (14 dias)

}
  \item{freq}{
Frequ�ncia da m�trica, di�ria � 7, mensal � 12 e hor�ria � 24
}
  \item{foreplot}{
Se TRUE, plota o gr�fico dos cinco melhores modelos ajustados
}
  \item{id_metrica_frequencia}{
Grava o nome correspondente nas tabelas do CMIS para identificar qual foi a proje��o calculada.
}
  \item{primeiro_modelo}{
Indica se o � o primeiro modelo ou j� existe um gravado na base do CMIS
}
  \item{parametro_1, ..., parametro_7}{
Par�metros de entrada de modelos j� presentes no CMIS. Caso um modelos j� exista no CMIS os par�tros deste passam pela fun��o CMIS_PROJECOES e testes s�o realizados com os mesmos. Eles tamb�m s�o par�metros de sa�da para o novo modelo ajustado para os dados.
}
  \item{modelo}{
Esta flag indica se j� existe ou n�o um modelo ativo no CMIS
}
  \item{saida_geral}{
Se TRUE Exibe v�rias estat�sticas do melhor modelo de s�rie temporal ajustado, como o objeto do modelo, a tabela de bonda do ajuste e escolha, os par�metros, os dados hist�ricos e a proje��o dos dados. Se FALSE, exibe apenas a proje��o com limites superior e inferior.
}
  \item{ore_save}{
Se TRUE salva a os outliers tratados, a proje��o e os par�metros para o melhor modelo escolhido no servidor CMIS.
}
  \item{\dots}{
Argumentos adicionais para fun��es internas
}
}
\details{
O c�lculo das proje��es passa basicamente em seis etapas: receber os dados, tratar as datas, tratar outliers, estimar modelos iniciais, escolher melhores modelos, salvar resultados.
Para cada etapa fun��es espec�ficas foram desenvolvidas e o conjunto integrado de todas resulta na fun��o geral \code{cmis_forecast}. Tratamento de outliers e informa��es detalhadas sobre os tipos de modelos de previs�o ou forecast ser�o discutidos na se��o modelos abordados.

TRATAMENTO DE OUTLIERS

Em virtude de anomalias nos dados, os chamados pontos discrepantes e tamb�m valores faltantes (missing values), em muitas an�lises tais casos s�o tratados como exe��o e, em muitas situa��es, t�cnicas adequadas de tratamento devem ser utilizadas visando contornar sem deturpar o comportamento real dos dados. Diversas t�cnicas de detec��o de outliers s�o utilizadas na literatuta atual, contudo esta quest�o � delicada em fun��o do tipo de dados envolvidos e de fatores inesperados como cat�strofes, por exemplo. No CMIS a vis�o estat�stica � abordada em conjunto com t�cnicas de visualiza��o em gr�ficos o que nos d� margem para maior robustez nas an�lises. 
A abordagem de tratamento de outliers no CMIS � feita em duas frentes que suporta dois tipos de situa��es suportadas pela abordagem de Hyndman .
 
I.	Dados n�o sazonais com ou sem missing: nestes casos interpola��o linear � empregada. Nesta abordagem, um modelo linear � ajustado com base no hist�rico e val�es estimados deste modelo s�o ent�o inseridos para completar a s�rie;

II.	Dados sazonais com ou sem missing: � feita uma decomposi��o da s�rie via stl(Seazonal Trend Loess ), que usa regress�o polinomial na s�rie dessazonalizada e em seguida resazonalisa a mesma com base na proje��o desta no per�do anterior. Embora complicada, esta t�cnica permite capturar comportamentos sazonais e evita a perda de informa��o da s�rie real.

TIPOS DE MODELOS ABORDADOS

Os modelos abordados no c�lculo das proje��es est�o divididos em tr�s grandes classes: ETS, HoltWinters e ARIMA, contudo estes modelos possuem casos especiais, abaixo segue um resumo dos principais pontos de cada tipo:

Modelos ETS

A termina��o ETS foi introduzida por Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008) em seu livro e no pacote forecast onde o autor consolidou os trinta tipos de modelos divididos entre os com alisamento ou suaviza��o exponenciais e os de espa�os de estados com inova��es em uma �nica equa��o e criou uma rotina em R capaz fazer proje��o a partir do melhor dos modelos ajustados atrav�s pela equa��o geral. Esta classe de modelos possui horizonte de proje��o e vida curtos. Em termos de tend�ncia, podem ser: N (nenhuma),A(aditiva), Ad (aditiva suavizada), M (multiplicativa), Ma (multiplicativa suavizada) e em termos de sazonalidade: N (nenhuma),	A (aditiva)	e M (multiplicativa) segundo a taxonomia empregada por Rob J Hyndman.

Modelos HoltWinters

O m�todo de HoltWinters � uma extens�o do m�todo de Holt que por sua vez utiliza m�dias m�veis exponenciais para fazer proje��es. Tal m�todo consiste em extender o m�todo de Holt de modo a lidar com varia��es tanto na tend�ncia quanto na sazonalidade dos dados. Estes modelos tamb�m possuem vers�es aditivas e multiplicativas, semelhantemente aos modelos ETS vistos anteriormente. Estes modelos possuem horizonte de proje��o e de vida de curto a m�dio.

Modelos ARIMA

Os Modelos de Alisamento Exponencial e os modelos ARIMA s�o os mais utilizados para fazer proje��es devido a robustez dos mesmos na modelagem de s�ries temporais a principal diferen�a entres as duas aboradagens � que a primeira se baseia na descri��o da tend�ncia e da sazonalidade dos dados enquanto o segundo nos ajuda a descrever melhor o fen�meno de autocorrela��o dos dados. Para mais informa��es sobre estes modelos veja as refr�ncias. Modelos ARIMA s�o compostos por dois componentes fundamentais: um processo estoc�stico Auto Regresivo (AR) e um processo M�dias M�veis (MA), do ingl�s Moving Average mais uma ordem de diferen�a. O horizonte de proje��o e a vida destes modelos � longa pras duas abordagens, sozonais e n�o sazonais.

ESTRAT�GIA DE SELE��O DE MODELOS

O Pacote forecast cont�m uma s�rie de m�todos de proje��o/forecast automatizados onde in�meros modelos s�o testados e os melhores s�o escolhidos para cada categoria. No CMIS o desafio � escolher n�o entre modelos de mesmo tipo, ex. ARIMA, mas entre modelos dos tipos ARIMA, ETS e HoltWinters. Embora para cada um destes modelos o algoritmo de escolha baseado em AICc seja plaus�vel, no nosso caso precisamos de m�todos que avaliem a qualidade de cada modelo de proje��o.
No mercado, existem tr�s principais abordagens de decis�o sobre a bondade do ajuste que s�o elas: Mean Absolute Percent Error (MAPE), Mean Absolute Error (MAE) e Root Mean Squared Error (RMSE). Uma t�cnica bastante utilizada � a de rankeamento onde cada estat�st�stica � ordenada do menor para o maior valor, e geralmente o modelo com menor MAPE (< 40 por cento) � escolhido como o mais adequado.
Para o CMIS, uma abordagem complementar desta t�cnica foi desenvolvida criando-se a soma de rankings que al�m das estat�sticas MAPE, MAE e RMSE adiciona mais duas estat�sticas que s�o o inverso do ranking de seis testes estat�sticos sobre os res�duos do modelo testado e o Mean Absolute Scaled Error (MASE). Assim, o modelo escolhido ser� aquele com a menor das somas de cinco rankings. A an�lise de res�duo do modelo ajustado � parte fundamental na qualidade do ajuste serve para verificar se o mesmo atende aos pressupostos assumidos para os res�duos e � feita com base em seis testes, vide \code{analise_residuos}.

S�o atribu�dos pesos para cada teste realizado nos res�duos. Quanto maior o peso mais importante � o teste, nesta vis�o os testes de independ�ncia, homocedasticidade e autocorrela��o s�o os mais importantes.

A regra para compor o ranking dos testes para os res�duos � composta da seguinte maneira (para cada modelo ajustado):

I.	� fixado o n�vel de confian�a � de 0.95 e ent�o realiza-se o teste e confere o p-valor do mesmo;

II.	Se a hip�tese de que o pressuposto foi atendido for verdadeira, recebe o peso de acordo com a tabela de pesos, sen�o 0. Faze-se a soma dos pesos e cria-se um sinalizador com o inverso da soma dos pesos que aqui damos o nome de 1/RESIDUOS (o m�ximo � 10, passou nos 6 testes);

A regra para compor o ranking geral � a seguinte:

I.	Ajusta-se cinco tipos de modelos iniciais, ARIMA, ETS, HoltWinters Sazonal, HoltWinters n�o sazonal e HoltWinters com Alisamento Exponencial e armazena as estat�sticas MAPE, MAE, RMSE, MASE e RESIDUOS;

II.	Ranqueia-se estas estat�sticas do menor para o maior valor para cada modelo e faz-se a soma destes valores e comp�e-se uma tabela.

III.	A menor soma corresponde ao modelo escolhido para a proje��o.
A Tabela 8 exemplifica a escolha de um modelo ARIMA dentre cinco possibilidades para uma m�trica qualquer atrav�s da t�cnica de kankeamento das estat�sticas de bondade e de an�lise de res�duos. Neste caso o modelo ETS ficou se segundo lugar e assim por diante.
Nesta tabela a estat�stica res�duos ficou com o valor 1/0,125 = 8 para o modelo escolhido, sabe-se que a soma m�xima � 10, contudo em todos os testes a bondade do ajuste foi efetiva e uma viola��o de 2 pontos na soma dos resultados dos testes de res�duos n�o desclassifica o modelo. Outra raz�o � que nosso interesse � mais forte na proje��o e n�o nos valores preditos para a m�trica. Pode haver situa��es em que a soma de rankings total ficar� acima de 5, nestas situa��es a menor soma ainda ser� escolhida, pois nem todas as m�tricas s�o plaus�veis de possu�rem os melhores modelos em virtude de situa��es diversas como outliers ou grandes varia��es e at� mesmo situa��es onde os dados s�o constantes ao longo do tempo.


}
\value{
Quando o argumento saida_geral = TRUE um objeto tipo lista � exibido com as seguintes sa�das:

\item{modelo}{Melhor modelo de s�rie temporal ajustado dentre de 34 tipos: ETS (Exponential smoothing state space model), HoltWinters (Holt-Winters Filtering) e ARIMA (Auto Regressive Integrated Moving Averege)}

\item{par_modelo}{Par�metros do modelo. Se for ARIMA, retornam 7 par�metros, caso HoltWinters ou ETS, menos par�metros s�o retornados.}

\item{bondade}{Tabela fundamental para a decis�o entre os modelos testados. Nela constam os rankings das estat�sticas de bondade utilizadas na tomada de decis�o. Entre elas est�o as utilizadas pelo mercado de forecast mundial: MAPE(Mean Absolute Percent Error), MAE (Mean Absolute Error) e RMSE(Root Mean Squared Error) al�m da estat�stica baseada na an�lise dos res�duos do modelo criada por Capacidade que � o inverso da soma de rankings das estat�sticas de teste dos res�duos. Mais detalhes adiante.}

\item{previsao}{Tabela da proje��o com limite superio e inferior a 0.95 de certeza. Esta tabela ser� salva no CMIS para exibi��o no dashboard.}

\item{projecao}{Tabela completa que cont�m al�m da proje��o os limites superiores e inferiores para 0.85 e 0.95 de certeza.}

\item{dados_historicos}{Os dados originais de entrada.}

\item{outliers}{Tabela contendo aquelas observa��es tratadas como autliers no processo de estima��o estat�stica.}

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
}


\author{
Desenvolvimento inicial por

GUILHERME, O.T.(jose.torres2@gvt.com.br)

Aprimoramentos feito por

LOPES, J. L. (jose.evandeilton@gvt.com.br)
}

\examples{

data(mensal)
fit <- cmis_forecast(mensal[,c(1,3)], 50, 12, foreplot = TRUE, saida_geral = TRUE)
}
