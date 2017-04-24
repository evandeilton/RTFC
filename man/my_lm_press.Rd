\name{my_lm_press}
\alias{my_lm_press}
\title{
Calcula estatistica PRESS
}
\description{
Estatística de força da predição do modelo (PRESS - PREdicted Residual Sum of Squares) para modelos lm e glm
}
\usage{
my_lm_press(obj, wts = c())
}
\arguments{
  \item{obj}{
Objeto do modelo ajustado (lm ou glm)
}
  \item{wts}{
Pesos no caso de regressão ponderada (pouco usado), se nulo regressão não ponderada que é o padrão.
}
}
\details{
A estatística PRESS vem sendo utilizada para medir o quão forte um modelo é na predição de valores extras. Ela é importante para forecasts lineares que só fazem sentido quando as extrapolações são feitas via cenários conhecidos ou no caso do CMIS, projetados por séries temporais robustas como ARIMA.
}
\value{
data.frame com as estatísticas da análise
}
\references{
https://github.com/WinVector/Examples/blob/master/PRESS/
http://www.win-vector.com/blog/2014/09/estimating-generalization-error-with-the-press-statistic/
http://www.r-bloggers.com/estimating-generalization-error-with-the-press-statistic/
}
\author{
Nina Zumel

Adaptação:

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