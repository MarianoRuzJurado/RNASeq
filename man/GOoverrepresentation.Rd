\name{GOoverrepresentation}
\alias{GOoverrepresentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{GOoverrepresentation
%%  ~~function to do ... ~~
}
\description{runs a GO over-representation test
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GOoverrepresentation(collectResults,
                      OrgDb,
                      pvalueCutoff,
                      pAdjustMethod,
                      qvalueCutoff,
                      minGSSize,
                      maxGSSize)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collectResults}{ =    a list, with sub lists}
  \item{OrgDb}{ =   organism db (e.g org.Hs.eg.db)}
  \item{pvalueCutoff}{ =    Cutoff value of pvalue}
  \item{pAdjustMethod}{ =   one method of the following: "holm", "hochberg", "hommel","bonferroni","BH","BY","fdr","none"}
  \item{qvalueCutoff}{ =    Cutoff value of qvalue}
  \item{minGSSize}{ =   minimal size of genes annotated by ontology term for testing}
  \item{maxGSSize}{ =   maximal size of genes annotated for testing}
%%     ~~Describe \code{x} here~~

}
\details{if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
%%  ~~ If necessary, more details than the description above ~~
}
\value{.xls files with results for molecular function, biological process and celular component and results are added to your collectResults list

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
