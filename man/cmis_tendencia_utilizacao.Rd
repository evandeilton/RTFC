\name{cmis_tendencia_utilizacao}
\alias{cmis_tendencia_utilizacao}
\title{
Indicador de tendencia de utilizacao e crescimento
}
\description{
Analisa tendência de utilização e crescimento de séries
}
\usage{
cmis_tendencia_utilizacao(metrica, thx = 80, forma = "tcm", plotm = FALSE)
}
\arguments{
  \item{metrica}{
Vetor de dados ou ts simples para calculo
}
  \item{thx}{
Threshold que indica o limitante de utilização de recurso sugerido
}
  \item{forma}{
Método de cálculo da tendência. "tcm" = Utiliza Taxa de Cresciento Médio. "tcme" = utiliza Taxa de Crescimento Exponencial
}
  \item{plotm}{
Plota o gráfico dos resultados
}
}
\details{
Para o cálculo do indicador de utilização: Toma-se os 50 por cento de pontos mais recentes da média móvel exponencial ponderada e tira-se destes as taxas de utilização em relação ao treshold. Isso é feito para dar mais peso aos dados mais recentes da variável analisada.
}
\value{
\item{data.frame}{data.frame com os resultados da análise. Este objeto contém as informações de tedência, taxas de utilização e média ponderada de utilização.}
}

\author{
LOPES, J. L.
}

\seealso{
Mais detalhes em \code{\link{WMA}}
}
\examples{
# Not run
#data(mensal)
#cmis_tendencia_utilizacao(mensal[, 5])
#cmis_tendencia_utilizacao(mensal[, 56])
#cmis_tendencia_utilizacao(mensal[, 62])
#cmis_tendencia_utilizacao(mensal[, 113])
#cmis_tendencia_utilizacao(AirPassengers)
}