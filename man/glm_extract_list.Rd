\name{glm_extract_list}
\alias{glm_extract_list}
\title{
Extrai estatistica de objetos LM ou GLM
}
\description{
Extrai estat�stica de objetos LM ou GLM sem intervalo de confian�a
}
\usage{
glm_extract_list(fit, r2, residteste, ...)
}
\arguments{
  \item{fit}{
Objeto lm ou glm
}
  \item{r2}{
Threshold do coeficiente de determina��o ou r-quadrado
}
  \item{residteste}{
Threshold do indicador de qualidade do res�duo gerado por \code{analise_residuos_glm}
}
  \item{\dots}{
Passagem de par�metros
}
}
\details{
Esta fun��o retorna diversas estatisticas do modelo testado de acordo com os argumentos de threshold. Entre elas est�o R^2, VIF, RMSE, PRESS, etc.

Nota: as estatistica retornadas dependem dos argumentos de threshold, logo, caso o modelo n�o apresente bons R^2 e/ou res�duos nada ser� retornado. O R^2 parcial � calculado com aporio da fun��o \code{calc.relimp}, do pacote relaimpo.
}
\value{
Data frame com coeficientes e estat�sticas do modelo
}
\author{
LOPES, J. L.
}

\seealso{
Esta fun��o depende das fun��es \code{analise_residuos_glm}, \code{my_lm_press} deste pacote mais \code{calc.relimp} do pacote \code{relaimp}
}
\examples{
# Objeto com dados para modelagem
data(mensal)
glmdata <- my_glm_data(mensal, 1, 12, "neg","hdw")
damod <- glmdata$ds_modelagem
vresposta <- glmdata$respostas
explicativas <- glmdata$explicativas
mreg <- paste(vresposta, " ~ 1")

fun  <- as.formula(paste(c(mreg[20], explicativas[1:5]), collapse=" + "))
fit  <- step(lm(fun, data = glmdata$ds_modelagem))
Glme <- glm_extract_list(fit, r2 = 0.5, residteste = 0.5)
}