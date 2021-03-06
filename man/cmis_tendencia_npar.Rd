\name{cmis_tendencia_npar}
\alias{cmis_tendencia_npar}
\title{
Faz analise nao parametrica de tendencia
}
\description{
Faz an�lise de tend�ncia de s�ries utilizando testes n�o param�tricos
}
\usage{
cmis_tendencia_npar(x, metodo = "stl", plot = TRUE, ...)
}

\arguments{
  \item{x}{
Vetor de dados ou ts simples
}
  \item{metodo}{
M�todo de c�lculo, utiliza decomposi��o da s�rie em STL com "stl" e decomposi��o por m�dias m�veis para "mm"
}
  \item{plot}{
Plota o gr�fico da an�lise
}
  \item{\dots}{
Argumentos extras para depend�ncias internas ou passagem de par�metros.
}
}
\details{
Faz extra��o da tend�ncia utilizando decomposi��o em STL ou M�dias M�veis e depois Calcula a dire��o desta tend�ncia utilizando O Slope de Sen e o teste de Mann-Kendal que s�o testes n�o param�tricos. No aobjeto de sa�da h� um elemento chamado sinal que representa a tend�ncia da s�rie analisada. Este n�mero � um indicador criado com base nas estat�sticas de Sen e tau de Mann-Kendal.
}
\value{
  \item{Vetor}{Vetor de resultados contendo o \code{slope}, \code{tau} de Mann-Kendal, \code{score}, \code{p.valor} de tau e \code{sinal},  um n�mero que indica a tend�ncia. A regra do indicador � aseguinte:
  
    2 = Tend�nca de crescimento cr�tico
    
    1 = Tend�nca de crescimento simples
    
    0 = Sem tend�nca de crescimento
    
   -1 = Tend�nca de queda simples
    
   -2 = Tend�nca de queda cr�tica
   }
}

\references{
Marchetto A., Rogora M., Arisci S. 2013 Trend analysis of atmospheric deposition data: a comparison of statistical approaches. Atmospheric Environment 64, 95-102

Helsel D.R., Frans L.M. 2006 The regional Kendall test for trend: Environmental Science and Technology 40, 4066-4073

Helsel D.R., Mueller D.K., Slack J.R. 2006 Computer program for the Kendall family of trend tests U.S. Geological Survey Scientific Investigations Report 2005-5275, 4 pp.

Hirsch R.M., Slack J.R., Smith R.A. Techniques of trend ananlyis for monthly water quality data. Water Resources Research 18, 107-121

Hirsch R.M., Slack J.R. 1984 A nonparametric test for seasonal data with serial dependance. Water Resources Research 20, 727-732

Libiseller C., Grimvall A. 2002 Perfomance of partial Mann-Kendall tests for trend detection in the presence of covariates. Environmetrics 13, 71-84

Mann H.B. 1945. Nonparametric tests against trend. Econometrica 13, 245-249
}
\author{
LOPES, J. L.
}

\seealso{
Para mais informa��es sobre os testes e sobre a teoria utilizada, vide \code{\link{rkt}} ou as refer�ncias.
}
\examples{
## Not run
#data(mensal)
#cmis_tendencia_npar(mensal[, 5])
#cmis_tendencia_npar(mensal[, 56])
#cmis_tendencia_npar(mensal[, 62])
#cmis_tendencia_npar(mensal[, 113])
#cmis_tendencia_npar(AirPassengers)
}