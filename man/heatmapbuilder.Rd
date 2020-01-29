\name{heatmapbuilder}
\alias{heatmapbuilder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{heatmapbuilder
%%  ~~function to do ... ~~
}
\description{Creates heatmaps based on a list to visualize differential expression in sample sets. Differential expressed genes are outfiled in a .xls-file and heatmaps as .pdf-file.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatmapbuilder(collectResults)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{collectResults}{ =    a list, with sub lists
%%     ~~Describe \code{x} here~~
}
}
\details{collectResults should contain: used differential expression tools (like DESeq2), contrasts of your samples, normalizedcounts, paths to your results folder and the differential expressed genes
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
Your list should look something like this:

collectResults$DETool$CommonResults$counts
                     $contrasts$folder
                               $diffGenes
                               $DF

DETool= Tool used (e.g. DESeq2), CommonResults = sub list with rawcounts and normalized counts, DF=Data Frame (e.g. output of run_DESeq2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
