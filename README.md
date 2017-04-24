# O pacote RTFC
RTFC é a sigla para R Trend, Forecast and Correlation. A ideia base deste pacote é fazer uma abordagem
estatística baseada em forecast de séries tempoarais para um CMIS (Common management Information Service).
Em outras palabras busca-se criar um módulo de forecast automatizado que pode ser integrado a um CMIS
genérico para trabalhar forncendo previsões para deterninadas métricas ou séries temporais

# Descrição
Faz análise Regressão/Correlação entre métricas explicativas e explicadas automaticamente pelo método 
do melhor conjunto de modelos, porém utilizando análises de multicolinearidade, redução por CrossValidation 
e LASSO para selecionar as melhores variáveis do conjunto de dados.

# Metodologia
O processo de modelagem pela função _cmis_correlacao_ envolve muita estatística e muitos detalhes 
teóricos e computacionais. O tratamento inicial dos dados é feito pelas mesmas funções utilizadas para 
a parte do _cmis_forecast_ que envolve tratamento de outliers, missing, métricas pobres de dados 
dentre outras. Em relação à modelagem, conceitos como redução de variáveis via ANOVA tipo II, análise 
de CrossValidation e redução pela metodologia LASSO são empregadas em busca da melhor combinação de 
técnicas para ajustes de modelos de forma mais inteligente evitando regressões espúrias e modelos 
estatísticamente fracos. Em termos de decisão, entram em cena as estatísticas F, coeficiante de 
determinação, estatística PRESS, ranking de resíduos e análise de AIC.

Os modelos abordados na análise de regressão estão divididos em duas classes: Modelos de Regressão
Linear Simples (MRLS) e Modelos de Regressão Linear Múltiplos (MRLM).

Quando o argumento de regressão simples é verdadeiro, todas as possibilidades de modelos de regressão
simples para os dados são testadas e não há redução de variáveis, neste caso os modelos são checados 
por reamostragens internas via CrossValidation. Se os modelos são bons segundo a qualidade de seus 
resíduos e do coeficiente de determinação, passa-se para a etapa de análise dos coeficientes estiamdos
via estatístcia Wald. No caso de modelos múltiplos uma verificação prévia é feita nas covariáveis em
busca de problemas de multicolinearidade que podem levar a modelos mal estimados e/ou regressões espúrias.
Se alguma métrica apresentar problemas de multicolinearidade ela é então removida do modelo e o mesmo
é reestimado e passa pela CrossValidation para um filtro mais sofisticado.

Em virtude de anomalias nos dados, os chamados pontos discrepantes e também valores faltantes (missing values),
em muitas análises tais casos são tratados como exeção e, em muitas situações, técnicas adequadas de 
tratamento devem ser utilizadas visando contornar sem deturpar o comportamento real dos dados. 
Diversas técnicas de detecção de outliers são utilizadas na literatuta atual, contudo esta questão é delicada
em função do tipo de dado envolvido e de fatores inesperados como catástrofes, por exemplo. A abordagem 
de tratamento de outliers no CMIS é feita em duas frentes que suporta dois tipos de situações 
suportadas pela abordagem de Hyndman.
 
I.	Dados não sazonais com ou sem missing: nestes casos interpolação linear é empregada. Nesta abordagem,
um modelo linear é ajustado com base no histórico e valões estimados deste modelo são então inseridos 
para completar a série;

II.	Dados sazonais com ou sem missing: é feita uma decomposição da série via stl(Seazonal Trend Loess),
que usa regressão polinomial na série dessazonalizada e em seguida resazonalisa a mesma com base na projeção 
desta no perído anterior. Embora complicada, esta técnica permite capturar comportamentos sazonais e evita 
a perda de informação da série real.

# A função core
```{R}
cmis_correlacao(dados, n_passos_frente = 1, freq, nvif = 100, nfolds = 5, nmod = 10,
	r2 = 0.5, residteste = 0.8, sig = 0.15, prefixo.expl = "neg", prefixo.resp = "hdw",
	verbose = FALSE, saida_geral = FALSE, regsimples = TRUE, ...)
}
```
## Argumentos
* ***dados :***
Conjunto de dados com a primeira sendo a DATA no formato DD/MM/AAAA HH:MM:SS e outras variáveis de negócio, infra, etc.
mínimo duas covariáveis.

* ***n_passos_frente :***
Número que indica o orizonte de projeção.
Métrica mensal, a projeção máxima é de 18 meses (3 anos);
Métrica diária, a projeção máxima é de 180 dias (6 meses)
Métrica é por hora, a projeção máxima é de 336 dias (14 dias).

* ***freq :***
Frequência da métrica, diária é 7, mensal é 12 e horária é 24.

* ***vif :***
Para modelos de regressão múltipla, faz a escolhas das melhores métricas utilizando o critério de multicolinearidade.
Quanto menor for melhor. A literatura pede 10 como ótimo, porém nos dados de séries temporais certa flexibilidade foi
imposta, deixando fixo em 300 para evitar perda de boas covariáveis já que outras abordagem é realizada em conjunto 
na redução e escolha das métricas.

* ***nfolds :***
Representa o número de reamostragens a ser realizada para cada modelo a fim de decidir sobre a capacidade de predição
do mesmo. Padrão é 5 reamostragens.

* ***nmod :***
Em modelos de regressão múltipa indica o total de modelos a serem retornados por métrica. Padrão é 5.

* ***r2 :***
É o threshold padrão para o coeficiente de deterninação. 0.5 é o mínimo padrão. Quanto maior melhor.

* ***residteste :***
Threshold para decisão de escolha baseada na qualidade da análise de resíduos. Seis testes são aplicados
e um indice de 1 a 10 é definido. Tomando-se o inverso deste número, quanto mais próximo de 0,10 melhor. 

* ***sig :***
Alpha de significância para escolha dos betas estimados nas regressões via teste estatístico de Wald. O padrão é 0,10. 

* ***prefixo.expl :***
É o prefixo das variáveis explicativas. Serve para qualquer combinação de três letras desde que as métricas
do data frame de entrada sejam renomeadas, com exceção da métricas de data, que deve ficar sempre com este nome.
Padrão é "neg", negócio.

* ***prefixo.resp :***
É o prefixo das variáveis resposta. Serve para qualquer combinação de três letras desde que as métricas do data
frame de entrada sejam renomeadas, com exceção da métricas de data, que deve ficar sempre com este nome.
Padrão é "hdw", infraestrutura.

* ***verbose :***
Se TRUE exibe o progresso de escolha dos melhores modelos na tela. O tempo de processamento aumenta devido às impressões.

* ***saida_geral :***
Se TRUE retorna uma lista aninhada, onde o primeiro elemento é o data.frame dos resultados no formato para o
banco de dados SQL e o outro é uma lista dos modelos escolhidos no formato do R que podem ser manipulados
com rotinas de loop do tipo sapply, ou lappy para análises diversas. Padrão é FALSE.

## Exemplos
```{R}
data(mensal)
data(diario)
# Apenas regressões simples com saída geral
fit0 <- cmis_correlacao(mensal, n_passos_frente = 1, freq = 12, nvif = 100, nfolds = 5, 
    nmod = 5, r2 = 0.5, residteste = 0.5, sig = 0.1, prefixo.expl = "neg", 
    prefixo.resp = "hdw", verbose = FALSE, saida_geral = TRUE, regsimples = TRUE) 
fit0[[1]]
fit0[[2]]

# Regressões simples e multiplas com saída geral 
fit1 <- cmis_correlacao(diario, n_passos_frente = 1, freq = 7, nvif = 100, nfolds = 5, 
    nmod = 5, r2 = 0.5, residteste = 0.5, sig = 0.1, prefixo.expl = "neg", 
    prefixo.resp = "hdw", verbose = FALSE, saida_geral = TRUE, regsimples = FALSE)     
#fit1[[1]]
#fit1[[2]]    
}
```

# Referências
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

