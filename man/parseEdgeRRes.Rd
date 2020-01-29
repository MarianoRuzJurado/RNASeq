\name{parseEdgeRRes}
\alias{parseEdgeRRes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{parseEdgeRRes
%%  ~~function to do ... ~~
}
\description{function, which parse myrun_edgeR results into a collectResults list and creates folders for edgeR results
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
parseEdgeRRes(cmp,myRes,collectResults=collectResults,padj =padjDEgenes ,lfc=lFCDEgenes, baseDir=baseDir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cmp}{ =       character matrix, contrasts defined in two columns}
  \item{collectResults}{ =      list, with results of myRun_edgeR}
  \item{baseDir}{ =       path, work directory}
  \item{myRes}{ =       list, containing ojects out of myRun_edgeR}
  \item{padj}{ =        global enviromnent threshhold for padj-value}
  \item{lfc}{ =         global enviromnent threshhold for log2foldchange-value}
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
list
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
