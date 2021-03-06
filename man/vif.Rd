\name{vif}
\alias{vif}
\title{
Calcula o VIF
}
\description{
Recebe objeto lm ou glm e calcula o VIF (Variance Infraction Factor) das vari�veis envolvidas no modelo.
}
\usage{
vif(obj, digits = 5)
}

\arguments{
  \item{obj}{
Objeto lm ou glm
}
  \item{digits}{
N�mero de casas decimais da sa�da
}
}

\value{
Vetor com o valor do VIF das vari�veis envolvidas
}

\references{
Fox, J. and Monette, G. (1992) Generalized collinearity diagnostics. JASA, 87, 178-183.

Fox, J. (2008) Applied Regression Analysis and Generalized Linear Models, Second Edition. Sage.

Fox, J. and Weisberg, S. (2011) An R Companion to Applied Regression, Second Edition, Sage.
}

\author{
GIOLO, S. R.
LOPES, J. E.
}
\note{
Se houver algum error ou n�o for poss�vel calcular VIF, o n�mero 9999 retornar�
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
#lmfit <- model.select(lmfit, sig = 0.05, verbose = TRUE)
vif(lmfit)
}