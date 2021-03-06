\name{descritiva}
\alias{descritiva}
\title{
Analise descritiva de dados
}
\description{
Faz an�lise descritiva em objetos com vari�veis num�ricas
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
Se TRUE, exibe estat�sticas descritivas simples
}
  \item{desc}{
Se TRUE, exibe outras estat�sticas e intervalos de confian�a 
}
  \item{norm}{
Se TRUE, faz normaliza��o dos dados antes de fazer a descritiva
}
  \item{p}{
N�vel de signific�ncia dos testes estat�stos internos
}
  \item{dig}{
N�mero de digitos no data.frame de saida
}
}

\value{
Data frame com estat�sticas da an�lise
}
\author{
Fr�d�ric Ibanez (ibanez@obs-vlfr.fr), Philippe Grosjean (phgrosjean@sciviews.org)
}
\note{
Esta fun��o foi adaptada e renomeada do pacote pastecs. O nome original � \code{stat.desc}
}

\examples{
data(mensal)
descritiva(mensal[,50:60], basic = TRUE)
descritiva(mensal[,50:60], basic = FALSE, dig = 3)
}