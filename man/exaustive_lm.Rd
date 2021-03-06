\name{exaustive_lm}
\alias{exaustive_lm}
\title{
Ajuste e escolha de modelos lineares simples e multiplos automaticamente
}
\description{
Utiliza diversos crit�rios estat�sticos para ajustar e escolher os melhore modelos lm de atrav�s de um conjunto de s�ries temporais de forma automatizada
}
\usage{
exaustive_lm(vresposta, explicativas, dados, r2, nmod, residteste, sig, verbose,
			regsimples = FALSE, ...)
}
\arguments{
    \item{vresposta}{
  Vari�vel reposta
  }
    \item{explicativas}{
  Vetor de vari�veis explicativas
  }
    \item{dados}{
  Data.frame no formato data.frame(data, var1, var2, var3, ..., varn)
  }
    \item{r2}{
  Threshold do coeficiente de determina��o ou r-quadrado
  }
    \item{nmod}{
  N�mero de melhores modelos a retornar para a vari�vel explicada em fun��o das explicativas
  }
    \item{residteste}{
  Threshold do indicador de qualidade do res�duo gerado por \code{analise_residuos_glm}
  }
    \item{sig}{
  Threshold de decis�o sobre o coeficiente, se ele entra ou n�o no modelo. O Padr�o � 0.05 mas atualmente usa-se 0.10.
  }
    \item{verbose}{
  Se TRUE, exibe log de todas as itera��es
  }
    \item{regsimples}{
  Se TRUE, ajusta apenas modelos de regress�o linear simples, caso contr�rio, simples e m�ltiplos
  }
    \item{\dots}{
  Passagem de par�metros internos
  }
}
\details{
Esta fun��o utiliza a alimina��o stepwise por p-valor para chegar aos melhores modelos para cada vari�vel explicativa e retorna alguns deles segundo exist�ncia ou defini��o do threshold nmod (n�mero de modelos). Esta fun��o integra v�rias outras do processo de ajuste, an�lise de res�duos e decis�o final. A sa�da cont�m os melhores modelos em formato de objeto lm, assim como data.frame com as principais estat�sticas no formato para o banco de dados. O n�mero de modelos cresce na escala de 2^n-1 para cada variavel explicativa adicional. Devido a isto, este procedimento se torna muito lento e consome muita mem�ria para conjuntos de vari�veis muito grande.
}

\value{
Lista contendo os seguintes elementos
  \item{objects}{Lista dos objetos lm resultados da an�lise, apenas os melhores escolhidos}
  \item{formulas}{As formulas matem�ticas e cada modelo}
  \item{vresposta}{Vari�vel resposta}
  \item{modelos}{data.frame com os modelos em formado de tabela para o banco de dados}
}

\references{
Alan Miller "Subset Selection in Regression" Chapman \& Hall

Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, http://www.stanford.edu/~hastie/Papers/glmnet.pdf

Journal of Statistical Software, Vol. 33(1), 1-22 Feb 2010
http://www.jstatsoft.org/v33/i01/

Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5) 1-13
http://www.jstatsoft.org/v39/i05/

Buckland (1997) Model Selection: an Integral Part of Inference. Biometrics 10:41

Burnham & Anderson (2002) Model Selection and Multimodel Inference: an Information Theoretic Approach Calcagno \& de Mazancourt 2010 J. Stat. Soft. v34 i12. See http://www.jstatsoft.org/v34/i12
}
\author{
LOPES, J. L.
}

\seealso{
Esta fun��o depende de \code{\link{model.select}}, \code{\link{analise_residuos_glm}} e  \code{\link{glm_extract_list_ci}}
}
\examples{
# Objeto com dados para modelagem (demorado!)

#data(mensal)
#glmdata <- my_glm_data(mensal, 1, 12, "neg","hdw")
#damod <- glmdata$ds_modelagem
#vresposta <- glmdata$respostas[1]
#explicativas <- glmdata$explicativas

## Infraestrutura em fun��o de neg�cio
#modelos <- exaustive_lm(vresposta, explicativas, damod, r2=0.50, nmod=10,
#			residteste = 0.3, sig = 0.05, verbose=FALSE)
}