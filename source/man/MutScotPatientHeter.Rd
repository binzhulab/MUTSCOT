
\name{MutScotPatientHeter}
\alias{MutScotPatientHeter}
\title{Test for heterogeneity}
\description{
  Test genes for subtype heterogeneity                  
}
\usage{
MutScotPatientHeter(obj, patientInfo, testGenes, options=NULL)
}
\arguments{
  \item{obj}{ Object returned from \code{\link{MutScot}}. No default. }
  \item{patientInfo}{ Data frame containing sample ids and numeric columns 
                 of adjusted continuous covariates. The first column  must be the ids
                and there must be at least one covariate (see details). 
                 No default. }
  \item{testGenes}{ Character vector of gene names to test. No default. }
  \item{options}{ List of options (see details). The default is NULL.}
}

\details{
  The first column of patient ids in the \code{patientInfo} data frame are matched
  with \code{obj$patientIds}. 

  \bold{Options} \cr
  \tabular{rrr}{
    Option   \tab Description                       \tab Default \cr
    MCT      \tab number of Monte-Carlo simulations \tab 1e4 
  } % END: tabular
}

\value{
  A list with the following names and descriptions: 

  \itemize{
    \item \code{pvalues} Data frame of gene-specific p-values for missense, nonsense, splicing, indel, 
                          and the test of subtype heterogeneity. A row of \code{NAs} will appear
                          for any gene not found in \code{obj$genes}.
    \item \code{glm_miss} List of the glm fitting of missense mutation for each gene.
    \item \code{glm_nsen} List of the glm fitting of nonsense mutation for each gene.
    \item \code{glm_splc} List of the glm fitting of splicing site mutation for each gene.
    \item \code{glm_indel} List of the glm fitting of indel mutation for each gene.
  } % END: itemize
}
%\references{
%  Choo-Wosoba H., Albert P.S., Zhu B.
%  A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.
%}

\author{Ho-Hsiang Wu and Bin Zhu <bin.zhu@nih.gov>}

\examples{
  data(data, package="MUTSCOT")
  data(MutSigCV.out, package="MUTSCOT")
  data(PatInfo, package="MUTSCOT")
  
  # Not run due to the time it takes
  #retM    <- MutScot(data, MutSigCV.out=MutSigCV.out)
  #x       <- retM$significantGenes
  #tmp     <- (x[, "Benchmark"]==TRUE) & (x[, "MutScot"]==TRUE)
  #TestSet <- as.vector(x[tmp, "Gene"])
  #retH    <- MutScotPatientHeter(retM, PatInfo, TestSet)
}


\keyword{ mutation } 
\seealso{\code{\link{MutScot}}}
