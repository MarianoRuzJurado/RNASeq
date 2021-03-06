\name{myrun_edgeR}
\alias{myrun_edgeR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{myrun_edgeR
%%  ~~function to do ... ~~
}
\description{run_edgeR function, which also return normalized counts
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
myrun_edgeR(countDF, targets, cmp, independent = TRUE, paired = NULL,
                     mdsplot = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{countDF}{ =   date.frame, containing raw  read counts}
  \item{targets}{ =   targets, data.frame}
  \item{cmp}{ =       character matrix, contrasts defined in two columns}
  \item{paired}{ =  defines pairs for paired analysis}
%%     ~~Describe \code{x} here~~
}

\details{if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list containing data.frame with edgeR results and normalized counts

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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
