\name{descritiva}
\alias{descritiva}
\title{
Analise descritiva de dados
}
\description{
Faz análise descritiva em objetos com variáveis numéricas
}
\usage{
descritiva(x, basic = TRUE, desc = TRUE, norm = FALSE, p = 0.95, dig = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Vetor, data.frame ou matriz de dados
}
  \item{basic}{
Se TRUE, exibe estatísticas descritivas simples
}
  \item{desc}{
Se TRUE, exibe outras estatísticas e intervalos de confiança 
}
  \item{norm}{
Se TRUE, faz normalização dos dados antes de fazer a descritiva
}
  \item{p}{
Nível de significância dos testes estatístos internos
}
  \item{dig}{
Número de digitos no data.frame de saida
}
}

\value{
Data frame com estatísticas da análise
}
\author{
Frédéric Ibanez (ibanez@obs-vlfr.fr), Philippe Grosjean (phgrosjean@sciviews.org)
}
\note{
Esta função foi adaptada e renomeada do pacote pastecs. O nome original é \code{stat.desc}
}

\examples{
data(mensal)
descritiva(mensal[,50:60], basic = TRUE)
descritiva(mensal[,50:60], basic = FALSE, dig = 3)
}