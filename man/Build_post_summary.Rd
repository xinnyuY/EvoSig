\name{Build_post_summary}
\alias{Build_post_summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Build_post_summary(input, output = NA, mergeType = F, typefile = "NA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
%%     ~~Describe \code{input} here~~
}
  \item{output}{
%%     ~~Describe \code{output} here~~
}
  \item{mergeType}{
%%     ~~Describe \code{mergeType} here~~
}
  \item{typefile}{
%%     ~~Describe \code{typefile} here~~
}
}
\details{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (input, output = NA, mergeType = F, typefile = "NA") 
{
    library(magrittr)
    library(dplyr)
    sample_list <- dir(input)
    colnames <- c("samplename", "Tumor_Sample_Barcode", "n_mutations", 
        "ccf_0-1_percentage", "ccf_0-2_percentage", "Ncluster", 
        "purity", "ccube_purity", "ccf_mean_cluster1", "ccf_mean_cluster2", 
        "ccf_mean_cluster3", "ccf_mean_cluster4", "ccf_mean_cluster5")
    post_summary_analyse <- function(samplename) {
        library(dplyr)
        suppressWarnings(rm(ssm))
        ssm <- load_ccf(samplename, input = input)
        ccf <- unique(ssm$ccube_ccf_mean)
        ccf_mean_order <- sort(ccf, decreasing = T)
        Ncluster <- length(ccf)
        post_summary <- data.frame(samplename <- samplename) \%>\% 
            mutate(Tumor_Sample_Barcode = unique(ssm$Tumor_Sample_Barcode), 
                n_mutations = nrow(ssm), ccf_01_percentage = mean(ssm$ccube_ccf <= 
                  1), ccf_02_percentage = mean(ssm$ccube_ccf <= 
                  2), Ncluster = Ncluster, purity = unique(ssm$purity), 
                ccube_purity = ifelse(exists("ssm$ccube_purity"), 
                  unique(ssm$ccube_purity), NA), ccf_mean_cluster1 = ifelse(Ncluster >= 
                  1, ccf_mean_order[1], 0), ccf_mean_cluster2 = ifelse(Ncluster >= 
                  2, ccf_mean_order[2], 0), ccf_mean_cluster3 = ifelse(Ncluster >= 
                  3, ccf_mean_order[3], 0), ccf_mean_cluster4 = ifelse(Ncluster >= 
                  4, ccf_mean_order[4], 0), ccf_mean_cluster5 = ifelse(Ncluster >= 
                  5, ccf_mean_order[5], 0)) \%>\% set_colnames(colnames)
        return(post_summary)
    }
    post_summary <- do.call(rbind, lapply(sample_list, post_summary_analyse))
    if (mergeType == T & !is.na(typefile)) {
        cancertype <- read.csv(file = typefile)[, -1]
        colnames(cancertype)[1] <- "samplename"
        post_summary <- left_join(post_summary, cancertype, by = "samplename")
        if ("Types" \%in\% colnames(post_summary)) 
            colnames(post_summary)[which(colnames(post_summary) == 
                "Types")] <- "cancertype"
    }
    if (!is.na(output)) {
        write.csv(post_summary, file = paste0(output, "post_summary_", 
            length(sample_list), "_", Sys.Date(), ".csv"))
        save(post_summary, file = paste0(output, "post_summary_", 
            length(sample_list), "_", Sys.Date(), ".RData"))
    }
    return(post_summary)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
