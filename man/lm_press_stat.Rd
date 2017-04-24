\name{lm_press_stat}
\alias{lm_press_stat}
\title{
Estatistica PRESS
}
\description{
Retorna as estatísticas R-quadrado, R-quadrado parcial e PRESS de modelos lineares
}
\usage{
lm_press_stat(obj)
}
\arguments{
  \item{package}{
Objeto das classe \code{lm} model
}
}
\details{
Utiliza resíduos de uma ANOVA feita a partir do objeto lm e retorna as estatísticas PRESS, R-Quadrado PRESS, R-Quadrado e R-Quadrado Ajustado
}

\author{
LOPES, J. L.
}

\examples{
data(diario)
fit <- lm(HDW_REDE_602~NEG_NDA_257, data = diario)
Lp <- lm_press_stat(fit)
}