\name{cmis_tendencia_utilizacao}
\alias{cmis_tendencia_utilizacao}
\title{
Indicador de tendencia de utilizacao e crescimento
}
\description{
Analisa tend�ncia de utiliza��o e crescimento de s�ries
}
\usage{
cmis_tendencia_utilizacao(metrica, thx = 80, forma = "tcm", plotm = FALSE)
}
\arguments{
  \item{metrica}{
Vetor de dados ou ts simples para calculo
}
  \item{thx}{
Threshold que indica o limitante de utiliza��o de recurso sugerido
}
  \item{forma}{
M�todo de c�lculo da tend�ncia. "tcm" = Utiliza Taxa de Cresciento M�dio. "tcme" = utiliza Taxa de Crescimento Exponencial
}
  \item{plotm}{
Plota o gr�fico dos resultados
}
}
\details{
Para o c�lculo do indicador de utiliza��o: Toma-se os 50 por cento de pontos mais recentes da m�dia m�vel exponencial ponderada e tira-se destes as taxas de utiliza��o em rela��o ao treshold. Isso � feito para dar mais peso aos dados mais recentes da vari�vel analisada.
}
\value{
\item{data.frame}{data.frame com os resultados da an�lise. Este objeto cont�m as informa��es de ted�ncia, taxas de utiliza��o e m�dia ponderada de utiliza��o.}
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