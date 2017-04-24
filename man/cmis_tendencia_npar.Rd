\name{cmis_tendencia_npar}
\alias{cmis_tendencia_npar}
\title{
Faz analise nao parametrica de tendencia
}
\description{
Faz análise de tendência de séries utilizando testes não paramétricos
}
\usage{
cmis_tendencia_npar(x, metodo = "stl", plot = TRUE, ...)
}

\arguments{
  \item{x}{
Vetor de dados ou ts simples
}
  \item{metodo}{
Método de cálculo, utiliza decomposição da série em STL com "stl" e decomposição por médias móveis para "mm"
}
  \item{plot}{
Plota o gráfico da análise
}
  \item{\dots}{
Argumentos extras para dependências internas ou passagem de parâmetros.
}
}
\details{
Faz extração da tendência utilizando decomposição em STL ou Médias Móveis e depois Calcula a direção desta tendência utilizando O Slope de Sen e o teste de Mann-Kendal que são testes não paramétricos. No aobjeto de saída há um elemento chamado sinal que representa a tendência da série analisada. Este número é um indicador criado com base nas estatísticas de Sen e tau de Mann-Kendal.
}
\value{
  \item{Vetor}{Vetor de resultados contendo o \code{slope}, \code{tau} de Mann-Kendal, \code{score}, \code{p.valor} de tau e \code{sinal},  um número que indica a tendência. A regra do indicador é aseguinte:
  
    2 = Tendênca de crescimento crítico
    
    1 = Tendênca de crescimento simples
    
    0 = Sem tendênca de crescimento
    
   -1 = Tendênca de queda simples
    
   -2 = Tendênca de queda crítica
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
Para mais informações sobre os testes e sobre a teoria utilizada, vide \code{\link{rkt}} ou as referências.
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