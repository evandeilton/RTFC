\name{cmis_forecast}
\alias{cmis_forecast}
\title{
Faz ajuste de modelos ARIMA, ETS e HoltWinters automaticamente
}
\description{
A função recebe como argumentos de entrada um conjunto de dados com duas variáveis: data e valor realizado conforme frequência, horária, diária ou mensal e aplica três tipos de metodologias de forecast em séries temporais, sendo eles ETS (Exponential smoothing state space model), HoltWinters (Holt-Winters Filtering) e ARIMA (Auto Regressive Integrated Moving Averege). Estas três metodologias de forecast combinadas geram cerca de 34 tipos de modelos e o melhor deles é escolhido.
}
\usage{
cmis_forecast(dados, n_passos_frente, freq, foreplot = FALSE, id_metrica_frequencia = 0,
	primeiro_modelo = 0, parametro_1 = NA, parametro_2 = NA, parametro_3 = NA,
	parametro_4 = NA, parametro_5 = NA, parametro_6 = NA, parametro_7 = NA,
	modelo = 0, saida_geral = FALSE, ore_save = FALSE, ...)
}

\arguments{
  \item{dados}{
Conjunto de dados com duas variáveis, a primeira sendo a DATA no formato DD/MM/AAAA HH:MM:SS e a outra contendo a métrica realizada.
}
  \item{n_passos_frente}{
Número que indica o orizonte de projeção.
Métrica mensal, a projeção máxima é de 18 meses (3 anos);
Métrica diária, a projeção máxima é de 180 dias (6 meses)
Métrica é por hora, a projeção máxima é de 336 dias (14 dias)

}
  \item{freq}{
Frequência da métrica, diária é 7, mensal é 12 e horária é 24
}
  \item{foreplot}{
Se TRUE, plota o gráfico dos cinco melhores modelos ajustados
}
  \item{id_metrica_frequencia}{
Grava o nome correspondente nas tabelas do CMIS para identificar qual foi a projeção calculada.
}
  \item{primeiro_modelo}{
Indica se o é o primeiro modelo ou já existe um gravado na base do CMIS
}
  \item{parametro_1, ..., parametro_7}{
Parâmetros de entrada de modelos já presentes no CMIS. Caso um modelos já exista no CMIS os parâtros deste passam pela função CMIS_PROJECOES e testes são realizados com os mesmos. Eles também são parâmetros de saída para o novo modelo ajustado para os dados.
}
  \item{modelo}{
Esta flag indica se já existe ou não um modelo ativo no CMIS
}
  \item{saida_geral}{
Se TRUE Exibe várias estatísticas do melhor modelo de série temporal ajustado, como o objeto do modelo, a tabela de bonda do ajuste e escolha, os parâmetros, os dados históricos e a projeção dos dados. Se FALSE, exibe apenas a projeção com limites superior e inferior.
}
  \item{ore_save}{
Se TRUE salva a os outliers tratados, a projeção e os parâmetros para o melhor modelo escolhido no servidor CMIS.
}
  \item{\dots}{
Argumentos adicionais para funções internas
}
}
\details{
O cálculo das projeções passa basicamente em seis etapas: receber os dados, tratar as datas, tratar outliers, estimar modelos iniciais, escolher melhores modelos, salvar resultados.
Para cada etapa funções específicas foram desenvolvidas e o conjunto integrado de todas resulta na função geral \code{cmis_forecast}. Tratamento de outliers e informações detalhadas sobre os tipos de modelos de previsão ou forecast serão discutidos na seção modelos abordados.

TRATAMENTO DE OUTLIERS

Em virtude de anomalias nos dados, os chamados pontos discrepantes e também valores faltantes (missing values), em muitas análises tais casos são tratados como exeção e, em muitas situações, técnicas adequadas de tratamento devem ser utilizadas visando contornar sem deturpar o comportamento real dos dados. Diversas técnicas de detecção de outliers são utilizadas na literatuta atual, contudo esta questão é delicada em função do tipo de dados envolvidos e de fatores inesperados como catástrofes, por exemplo. No CMIS a visão estatística é abordada em conjunto com técnicas de visualização em gráficos o que nos dá margem para maior robustez nas análises. 
A abordagem de tratamento de outliers no CMIS é feita em duas frentes que suporta dois tipos de situações suportadas pela abordagem de Hyndman .
 
I.	Dados não sazonais com ou sem missing: nestes casos interpolação linear é empregada. Nesta abordagem, um modelo linear é ajustado com base no histórico e valões estimados deste modelo são então inseridos para completar a série;

II.	Dados sazonais com ou sem missing: é feita uma decomposição da série via stl(Seazonal Trend Loess ), que usa regressão polinomial na série dessazonalizada e em seguida resazonalisa a mesma com base na projeção desta no perído anterior. Embora complicada, esta técnica permite capturar comportamentos sazonais e evita a perda de informação da série real.

TIPOS DE MODELOS ABORDADOS

Os modelos abordados no cálculo das projeções estão divididos em três grandes classes: ETS, HoltWinters e ARIMA, contudo estes modelos possuem casos especiais, abaixo segue um resumo dos principais pontos de cada tipo:

Modelos ETS

A terminação ETS foi introduzida por Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008) em seu livro e no pacote forecast onde o autor consolidou os trinta tipos de modelos divididos entre os com alisamento ou suavização exponenciais e os de espaços de estados com inovações em uma única equação e criou uma rotina em R capaz fazer projeção a partir do melhor dos modelos ajustados através pela equação geral. Esta classe de modelos possui horizonte de projeção e vida curtos. Em termos de tendência, podem ser: N (nenhuma),A(aditiva), Ad (aditiva suavizada), M (multiplicativa), Ma (multiplicativa suavizada) e em termos de sazonalidade: N (nenhuma),	A (aditiva)	e M (multiplicativa) segundo a taxonomia empregada por Rob J Hyndman.

Modelos HoltWinters

O método de HoltWinters é uma extensão do método de Holt que por sua vez utiliza médias móveis exponenciais para fazer projeções. Tal método consiste em extender o método de Holt de modo a lidar com variações tanto na tendência quanto na sazonalidade dos dados. Estes modelos também possuem versões aditivas e multiplicativas, semelhantemente aos modelos ETS vistos anteriormente. Estes modelos possuem horizonte de projeção e de vida de curto a médio.

Modelos ARIMA

Os Modelos de Alisamento Exponencial e os modelos ARIMA são os mais utilizados para fazer projeções devido a robustez dos mesmos na modelagem de séries temporais a principal diferença entres as duas aboradagens é que a primeira se baseia na descrição da tendência e da sazonalidade dos dados enquanto o segundo nos ajuda a descrever melhor o fenômeno de autocorrelação dos dados. Para mais informações sobre estes modelos veja as refrências. Modelos ARIMA são compostos por dois componentes fundamentais: um processo estocástico Auto Regresivo (AR) e um processo Médias Móveis (MA), do inglês Moving Average mais uma ordem de diferença. O horizonte de projeção e a vida destes modelos é longa pras duas abordagens, sozonais e não sazonais.

ESTRATÉGIA DE SELEÇÃO DE MODELOS

O Pacote forecast contém uma série de métodos de projeção/forecast automatizados onde inúmeros modelos são testados e os melhores são escolhidos para cada categoria. No CMIS o desafio é escolher não entre modelos de mesmo tipo, ex. ARIMA, mas entre modelos dos tipos ARIMA, ETS e HoltWinters. Embora para cada um destes modelos o algoritmo de escolha baseado em AICc seja plausível, no nosso caso precisamos de métodos que avaliem a qualidade de cada modelo de projeção.
No mercado, existem três principais abordagens de decisão sobre a bondade do ajuste que são elas: Mean Absolute Percent Error (MAPE), Mean Absolute Error (MAE) e Root Mean Squared Error (RMSE). Uma técnica bastante utilizada é a de rankeamento onde cada estatístística é ordenada do menor para o maior valor, e geralmente o modelo com menor MAPE (< 40 por cento) é escolhido como o mais adequado.
Para o CMIS, uma abordagem complementar desta técnica foi desenvolvida criando-se a soma de rankings que além das estatísticas MAPE, MAE e RMSE adiciona mais duas estatísticas que são o inverso do ranking de seis testes estatísticos sobre os resíduos do modelo testado e o Mean Absolute Scaled Error (MASE). Assim, o modelo escolhido será aquele com a menor das somas de cinco rankings. A análise de resíduo do modelo ajustado é parte fundamental na qualidade do ajuste serve para verificar se o mesmo atende aos pressupostos assumidos para os resíduos e é feita com base em seis testes, vide \code{analise_residuos}.

São atribuídos pesos para cada teste realizado nos resíduos. Quanto maior o peso mais importante é o teste, nesta visão os testes de independência, homocedasticidade e autocorrelação são os mais importantes.

A regra para compor o ranking dos testes para os resíduos é composta da seguinte maneira (para cada modelo ajustado):

I.	É fixado o nível de confiança é de 0.95 e então realiza-se o teste e confere o p-valor do mesmo;

II.	Se a hipótese de que o pressuposto foi atendido for verdadeira, recebe o peso de acordo com a tabela de pesos, senão 0. Faze-se a soma dos pesos e cria-se um sinalizador com o inverso da soma dos pesos que aqui damos o nome de 1/RESIDUOS (o máximo é 10, passou nos 6 testes);

A regra para compor o ranking geral é a seguinte:

I.	Ajusta-se cinco tipos de modelos iniciais, ARIMA, ETS, HoltWinters Sazonal, HoltWinters não sazonal e HoltWinters com Alisamento Exponencial e armazena as estatísticas MAPE, MAE, RMSE, MASE e RESIDUOS;

II.	Ranqueia-se estas estatísticas do menor para o maior valor para cada modelo e faz-se a soma destes valores e compõe-se uma tabela.

III.	A menor soma corresponde ao modelo escolhido para a projeção.
A Tabela 8 exemplifica a escolha de um modelo ARIMA dentre cinco possibilidades para uma métrica qualquer através da técnica de kankeamento das estatísticas de bondade e de análise de resíduos. Neste caso o modelo ETS ficou se segundo lugar e assim por diante.
Nesta tabela a estatística resíduos ficou com o valor 1/0,125 = 8 para o modelo escolhido, sabe-se que a soma máxima é 10, contudo em todos os testes a bondade do ajuste foi efetiva e uma violação de 2 pontos na soma dos resultados dos testes de resíduos não desclassifica o modelo. Outra razão é que nosso interesse é mais forte na projeção e não nos valores preditos para a métrica. Pode haver situações em que a soma de rankings total ficará acima de 5, nestas situações a menor soma ainda será escolhida, pois nem todas as métricas são plausíveis de possuírem os melhores modelos em virtude de situações diversas como outliers ou grandes variações e até mesmo situações onde os dados são constantes ao longo do tempo.


}
\value{
Quando o argumento saida_geral = TRUE um objeto tipo lista é exibido com as seguintes saídas:

\item{modelo}{Melhor modelo de série temporal ajustado dentre de 34 tipos: ETS (Exponential smoothing state space model), HoltWinters (Holt-Winters Filtering) e ARIMA (Auto Regressive Integrated Moving Averege)}

\item{par_modelo}{Parâmetros do modelo. Se for ARIMA, retornam 7 parâmetros, caso HoltWinters ou ETS, menos parâmetros são retornados.}

\item{bondade}{Tabela fundamental para a decisão entre os modelos testados. Nela constam os rankings das estatísticas de bondade utilizadas na tomada de decisão. Entre elas estão as utilizadas pelo mercado de forecast mundial: MAPE(Mean Absolute Percent Error), MAE (Mean Absolute Error) e RMSE(Root Mean Squared Error) além da estatística baseada na análise dos resíduos do modelo criada por Capacidade que é o inverso da soma de rankings das estatísticas de teste dos resíduos. Mais detalhes adiante.}

\item{previsao}{Tabela da projeção com limite superio e inferior a 0.95 de certeza. Esta tabela será salva no CMIS para exibição no dashboard.}

\item{projecao}{Tabela completa que contém além da projeção os limites superiores e inferiores para 0.85 e 0.95 de certeza.}

\item{dados_historicos}{Os dados originais de entrada.}

\item{outliers}{Tabela contendo aquelas observações tratadas como autliers no processo de estimação estatística.}

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
