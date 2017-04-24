\name{glm_extract_list}
\alias{glm_extract_list}
\title{
Extrai estatistica de objetos LM ou GLM
}
\description{
Extrai estatística de objetos LM ou GLM sem intervalo de confiança
}
\usage{
glm_extract_list(fit, r2, residteste, ...)
}
\arguments{
  \item{fit}{
Objeto lm ou glm
}
  \item{r2}{
Threshold do coeficiente de determinação ou r-quadrado
}
  \item{residteste}{
Threshold do indicador de qualidade do resíduo gerado por \code{analise_residuos_glm}
}
  \item{\dots}{
Passagem de parâmetros
}
}
\details{
Esta função retorna diversas estatisticas do modelo testado de acordo com os argumentos de threshold. Entre elas estão R^2, VIF, RMSE, PRESS, etc.

Nota: as estatistica retornadas dependem dos argumentos de threshold, logo, caso o modelo não apresente bons R^2 e/ou resíduos nada será retornado. O R^2 parcial é calculado com aporio da função \code{calc.relimp}, do pacote relaimpo.
}
\value{
Data frame com coeficientes e estatísticas do modelo
}
\author{
LOPES, J. L.
}

\seealso{
Esta função depende das funções \code{analise_residuos_glm}, \code{my_lm_press} deste pacote mais \code{calc.relimp} do pacote \code{relaimp}
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