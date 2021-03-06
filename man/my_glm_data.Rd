\name{my_glm_data}
\alias{my_glm_data}
\title{
Trata dados para regressoes
}
\description{
Trata dados para entrada no algoritimo das regress�es
}
\usage{
my_glm_data(dados, n_passos_frente, freq, prefixo.expl, prefixo.resp, ...)
}

\arguments{
  \item{dados}{
data.frame com vari�veis seguindo o modelo (data.frame (DATA, var1, var2, var3, ..., varn)
}
  \item{n_passos_frente}{
N�mero de passos utilizados na proje��o, necess�rio para a fun�ao ts.dados trabalhar com datas
}
  \item{freq}{
Frequ�ncia dos dados. Mensal (12), di�rio(7), hor�rio(24)
}
  \item{prefixo.expl}{
String tamanho 3 com o prefixo das vari�veis explicativas na forma (aaa). Ex: negocio = (neg), aplica��o = (apl). Estes nomes precisam ser definidos no data.frame antes da entrada na fun��o
}
  \item{prefixo.resp}{
String tamanho 4 com o prefixo das vari�veis respostas na forma (bbb). Ex: infraestrutura = (inf), hardware = (hdw). Estes nomes precisam ser definidos no data.frame antes da entrada na fun��o.
}
  \item{\dots}{
Passagem de par�metros
}
}
\details{
Esta fun��o faz uso da fun��o \code{dataframe.outlier} para tratar os dados e retornar um data.frame limpo pronto para modelagem.
}
\value{
Lista com variaveis resposta, explicativas e data.frame com os dados limpos para modelagem.
\item{explicativas}{vetor de variaveis explicativas}
\item{respostas}{vetor de vari�veis resposta}
\item{ds_modelagem}{data.frame limpo e pronto para modelagem}
}

\author{
LOPES, J. L.
}

\examples{
data(mensal)
glmdata <- my_glm_data(mensal, 1, 12, "neg","hdw")
names(glmdata)
}
