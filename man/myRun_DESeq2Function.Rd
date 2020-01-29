\name{myRun_DESeq2Function}
\alias{myRun_DESeq2Function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{myRun_DESeq2Function
%%  ~~function to do ... ~~
}
\description{Run_DESeq2 function, which also return the dds object
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
myRun_DESeq2Function(countDF, targets, cmp, independent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countDF}{ =   date.frame, containing raw  read counts}
  \item{targets}{ =   targets, data.frame}
  \item{cmp}{ =       character matrix, contrasts defined in two columns}
%%     ~~Describe \code{x} here~~
}

\details{if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list containing data.frame with DESeq2 results & dds DESeq2 object
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
