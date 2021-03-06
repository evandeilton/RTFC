\name{my_lm_press}
\alias{my_lm_press}
\title{
Calcula estatistica PRESS
}
\description{
Estat�stica de for�a da predi��o do modelo (PRESS - PREdicted Residual Sum of Squares) para modelos lm e glm
}
\usage{
my_lm_press(obj, wts = c())
}
\arguments{
  \item{obj}{
Objeto do modelo ajustado (lm ou glm)
}
  \item{wts}{
Pesos no caso de regress�o ponderada (pouco usado), se nulo regress�o n�o ponderada que � o padr�o.
}
}
\details{
A estat�stica PRESS vem sendo utilizada para medir o qu�o forte um modelo � na predi��o de valores extras. Ela � importante para forecasts lineares que s� fazem sentido quando as extrapola��es s�o feitas via cen�rios conhecidos ou no caso do CMIS, projetados por s�ries temporais robustas como ARIMA.
}
\value{
data.frame com as estat�sticas da an�lise
}
\references{
https://github.com/WinVector/Examples/blob/master/PRESS/
http://www.win-vector.com/blog/2014/09/estimating-generalization-error-with-the-press-statistic/
http://www.r-bloggers.com/estimating-generalization-error-with-the-press-statistic/
}
\author{
Nina Zumel

Adapta��o:

LOPES, J. L.
}

\examples{
data(mensal)
glmdata <- my_glm_data(mensal, 1, 12, "neg","hdw")
damod <- glmdata$ds_modelagem
vresposta <- glmdata$respostas[5]
explicativas <- glmdata$explicativas[1:3]
mreg <- paste(vresposta, " ~ 1")
fun  <- as.formula(paste(c(mreg, explicativas), collapse=" + "))
lmfit <- lm(fun, data = damod)
lmfit <- model.select(lmfit, sig = 0.05, verbose = FALSE)
Lmfit <- my_lm_press(lmfit)
}