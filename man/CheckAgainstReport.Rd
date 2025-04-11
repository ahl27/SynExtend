\name{CheckAgainstReport}
\alias{CheckAgainstReport}
\title{
Pull an assembly from the NCBI FTP site.
}
\description{
This function is designed to work internally to functions within \code{SynExtend} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
CheckAgainstReport(FTP_ADDRESS,
                   CHECK_ADDRESS,
                   RETRY = 5L)
}
\arguments{
  \item{FTP_ADDRESS}{
  Character; the ftp address of an ncbi assembly.
}
  \item{CHECK_ADDRESS}{
  Character; the ftp address of an ncbi assembly report.
}
  \item{RETRY}{
  Integer; the number of times to retry an assembly download should it not pull correctly.
}
}
\details{
  On occasion, \code{\link{readDNAStringSet}} fails to completely pull assemblies from the ncbi ftp site. It is not clear why, though it is infrequent but replicable at large scale. \code{CheckAgainstReport} checks the captured \code{DNAStringSet} against the reported assembly size and string widths.
}
\value{
  A \code{DNAStringSet}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{readDNAStringSet}}
}
\examples{
#
}
