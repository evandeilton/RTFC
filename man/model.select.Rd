\name{model.select}
\alias{model.select}
\title{
Selecao automatica de modelos lm
}
\description{
Sele��o autom�tica de modelos lm simples e multiplos por redu��o stepwise utilizando p-valor dos betas estimados
}
\usage{
model.select(model, keep, sig, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
Objeto de ajuste lm. Este objeto deve conter o modelo completo (com maior n�mero de preditoras poss�vel, para a redu��o de vari�veis em cada itera��o)
}
  \item{keep}{
lista de termos do modelo que devem permanecer em cada itera��o
}
  \item{sig}{
Threshold de decis�o sobre o coeficiente, se ele entra ou n�o no modelo. O Padr�o � 0.05, mas atualmente usa-se 0.10.
}
  \item{verbose}{
Se TRUE, exibe log de todas as itera��es
}
  \item{\dots}{
Passagem de par�metros
}
}
\details{
O modelo de entrada deve ser aquele com o maior n�mero de preditoras, ou seja, modelo completo. Desta forma inicia a elimina��o. Sempre que o p-valor de algum coeficiente for superior a sig, este sair� do modelo e o modelo passar� a ser reajustado. Quando n�o restarem mais covari�veis a sair ou entrar no modelo o script reajusta o modelo final e exibe no formato lm. O processo de elimina��o segue utiliza os p-valores gerados atrav�s de uma ANOVA tipo 2 para cada beta estimado.
}
\value{
Objeto lm com ajuste reduzido via stepwise
}
\references{
Brown, P. J. (1994) Measurement, Regression and Calibration Oxford.

http://stackoverflow.com/a/3701896

http://rstudio-pubs-static.s3.amazonaws.com/17812_139f1aec12d94206b25577ffd3730de1.html
}
\author{
Joris Meys
}
\note{
Esta fun��o foi adapdata da vers�o feita por Joris Meys publicada no stackoverflow
}

\examples{
# Objeto com dados para modelagem
data(mensal)
glmdata <- my_glm_data(mensal, 1, 12, "neg","hdw")
damod <- glmdata$ds_modelagem
vresposta <- glmdata$respostas[1]
explicativas <- glmdata$explicativas
mreg <- paste(vresposta, " ~ 1")
fun  <- as.formula(paste(c(mreg, explicativas), collapse=" + "))
lmfit <- lm(fun, data = damod)
summary(lmfit)
ms <- model.select(lmfit, sig = 0.05, verbose = FALSE)
summary(ms)
}