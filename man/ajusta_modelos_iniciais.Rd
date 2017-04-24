\name{ajusta_modelos_iniciais}
\alias{ajusta_modelos_iniciais}
\title{
Ajusta ARIMA, ETS e HoltWinters
}
\description{
Recebe data.frame no formato data.frame(data, metrica), ajusta 35 tipos de modelos de séries temporais e escolhe um ARIMA, um ETS e três HoltWinters, sazonal, não sazonal e outro com alisamento exponencial.
}
\usage{
ajusta_modelos_iniciais(obj.dados, ...)
}

\arguments{
  \item{obj.dados}{
data.frame para ajustes iniciais
}
}
\details{
Atualmente esta função faz ajustes de 35 tipos de modelos ARIMA, ETS e HolWinters e escolhe 5 finais. A função faz ainda a análise de resíduos e a projeção de acordo com as informações presentes no objeto de dados (obj.dados).

}
\value{
A saída é uma lista aninhada com várias entradas conforme abaixo.
\item{modelo_escolhido}{Objeto do modelo escolhido automaticamente (melhor dos 5)}
\item{gof}{Estatística da bondade do modelo escolhido}
\item{modelos_testados}{Sublista com informações dos cinco melhores modelos ajustados para a série}
\item{dados_historicos}{São os dados históricos}
\item{data_previsao}{Objeto com as datas projetadas}
\item{outliers}{Objeto contendo os outliers encontrados, caso haja algum}
\item{escolha}{Nome do modelo escolhido}
}
\references{
Hyndman, R.J. and Khandakar, Y. (2008) "Automatic time series forecasting: The forecast package for R", Journal of Statistical Software, 26(3).

Hyndman, R.J., Koehler, A.B., Snyder, R.D., and Grose, S. (2002) "A state space framework for automatic forecasting using exponential smoothing methods", International J. Forecasting, 18(3), 439-454.

Hyndman, R.J., Akram, Md., and Archibald, B. (2008) "The admissible parameter space for exponential smoothing models". Annals of Statistical Mathematics, 60(2), 407-426.

Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008) Forecasting with exponential smoothing: the state space approach, Springer-Verlag. http://www.exponentialsmoothing.net.

C. C. Holt (1957) Forecasting trends and seasonals by exponentially weighted moving averages, ONR Research Memorandum, Carnegie Institute of Technology 52.

P. R. Winters (1960) Forecasting sales by exponentially weighted moving averages, Management Science 6, 324-342.
}
\author{
LOPES, J. E.
}
\note{
 Modelos de redes neurais podem ser ajustados, mas atualmente o pacote forecast não possui implementação de bandas de confiança para este tipo de modelo. Sendo assim, ele foi desativado nesta versão.
}

\seealso{
Mais informações sobre a metodologia de forecast para os modelos ARIMA, ETS e HoltWinters podem ser consultadas nas funções \code{\link{auto.arima}}, \code{\link{ets} e \code{\link{HoltWinters}}}
}
\examples{
data(mensal)
dats <- ts.dados(mensal[, 1:2], 10, 12)
modelos <- ajusta_modelos_iniciais(dats)
}
