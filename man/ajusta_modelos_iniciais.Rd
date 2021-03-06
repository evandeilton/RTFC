\name{ajusta_modelos_iniciais}
\alias{ajusta_modelos_iniciais}
\title{
Ajusta ARIMA, ETS e HoltWinters
}
\description{
Recebe data.frame no formato data.frame(data, metrica), ajusta 35 tipos de modelos de s�ries temporais e escolhe um ARIMA, um ETS e tr�s HoltWinters, sazonal, n�o sazonal e outro com alisamento exponencial.
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
Atualmente esta fun��o faz ajustes de 35 tipos de modelos ARIMA, ETS e HolWinters e escolhe 5 finais. A fun��o faz ainda a an�lise de res�duos e a proje��o de acordo com as informa��es presentes no objeto de dados (obj.dados).

}
\value{
A sa�da � uma lista aninhada com v�rias entradas conforme abaixo.
\item{modelo_escolhido}{Objeto do modelo escolhido automaticamente (melhor dos 5)}
\item{gof}{Estat�stica da bondade do modelo escolhido}
\item{modelos_testados}{Sublista com informa��es dos cinco melhores modelos ajustados para a s�rie}
\item{dados_historicos}{S�o os dados hist�ricos}
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
 Modelos de redes neurais podem ser ajustados, mas atualmente o pacote forecast n�o possui implementa��o de bandas de confian�a para este tipo de modelo. Sendo assim, ele foi desativado nesta vers�o.
}

\seealso{
Mais informa��es sobre a metodologia de forecast para os modelos ARIMA, ETS e HoltWinters podem ser consultadas nas fun��es \code{\link{auto.arima}}, \code{\link{ets} e \code{\link{HoltWinters}}}
}
\examples{
data(mensal)
dats <- ts.dados(mensal[, 1:2], 10, 12)
modelos <- ajusta_modelos_iniciais(dats)
}
