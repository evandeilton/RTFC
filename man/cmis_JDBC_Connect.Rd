\name{Extracao de dados via RJDBC}
\alias{cmis_JDBC_Connect}
\title{
Conecta a base do CMIS via RJDB e coleta dados
}
\description{
Conecta a base Oracle do CMIS via RJDB e coleta dados utilizando o driver RJDBC
}
\usage{
cmis_JDBC_Connect(local_jdbcdriver, user, pass, JavaRAM = "-XmX2048m", 
	setenv = FALSE, JAVA_HOME = NULL, autocommit = FALSE, ...)
}

\arguments{
  \item{local_jdbcdriver}{
Local do driver jdbc
}
  \item{user}{ 
Nome do usuario a logar no banco de dados
}
  \item{pass}{
Senha de acesso ao banco de dados
}
  \item{JavaRAM}{
Memoria RAM alocada para a sessão
}
  \item{setenv}{
Configura o ambiente de trabalho
}
  \item{JAVA_HOME}{
Variável de ambiente do JAVA_HOME
}
  \item{autocommit}{
Comit automatico (cuidado!)
}
  \item{\dots}{
Passagem de argumentos.
}
}

\author{
LOPES, J. L.
}

\seealso{
\code{\link{dbConnect}}, \code{\link{JDBC}}
}
\examples{
## Not run
# O driver 'ojdbc6.jar' é parte do Java, encontre ele caso o R não ache automaticamente.
#local_jdbcdriver <- "C://Users//G0047743//Documents//sqldeveloper//jdbc//lib/ojdbc6.jar"
#user <- "CMIS_OWNER"
#pass <- "CMIS_OWNER"
#connect <- cmis_JDBC_Connect(local_jdbcdriver, user, pass)
}
