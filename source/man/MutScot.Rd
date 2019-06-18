
\name{MutScot}
\alias{MutScot}
\title{MutScot}
\description{
  MutScot                        
}
\usage{
MutScot(data, MutSigCV.out=NULL, options=NULL)
}
\arguments{
  \item{data}{ Data frame contining five columns for sample id,
    chromosome, position (hg19), reference allele, and alternate allele.
    The columns must be in the order given above. No default }
  \item{MutSigCV.out}{ Information provided by MutSigCV algorithm (see details).
              The default is NULL.}
  \item{options}{ List of options (see details). The default is NULL.}
}

\details{
  See \code{data(MutSigCV)} for an example format of \code{MutSigCV.out}. Currently, the only columns 
  of this object that are used are: \code{gene}, \code{X}, \code{x}, and \code{q}.

  \bold{Options} \cr
  \tabular{lll}{
    Option         \tab Description                                                        \tab Default \cr
    maxMutPerPat   \tab maximum number of mutations per patient                            \tab 300 \cr 
    multitesting   \tab method used for multiple tests (see \code{\link[stats]{p.adjust}}) \tab "fdr"  \cr
    alpha.q        \tab significance level                                                 \tab 0.1 \cr
    benchmarkGenes \tab vector of genes for summary                                        \tab genes returned from \cr
                   \tab                                                                    \tab \code{\link{getBenchmarkGenes}}
  } % END: tabular
}
\value{
  A list of class "MutScot" with the following names and descriptions: 

  \itemize{
    \item \code{driverGenes} Vector of identified driver genes.
    \item \code{significantGenes} Data frame of significant genes among the methods determined by \code{options$alpha.q}.
    \item \code{MutScotGenes} Data frame of p-values and adjusted p-values for the MutScot method.
    \item \code{countMat} Matrix of TP, FP, F1 counts by method.
    \item \code{dndsout} Output from the dndscv algorithm.
    \item \code{patientIds} Vector of patient ids.
    \item \code{genes} Vector of gene names.
    \item \code{mutationRates} List containing relative SNV mutation rates of patients (\code{PatRat}),
               weight parameter of patients w.r.t SNV mutations (\code{wt}),
               relative Indel mutation rates of patients (\code{PatRatI}),
               weight parameter of patients w.r.t Indel mutations (\code{wtI}),
               annotated Indel mutations (\code{Indels_record}) ,
               annotated observed/expected Indel mutations (\code{geneindels}),
               IVW estimated gene-specific background mutation rates (\code{Lambda_g}).
    \item \code{obsExp} List containing gene-specific observed and expected missense mutation counts (\code{miss}),
                        gene-specific observed and expected nonsense mutation counts (\code{nsen}),
                        gene-specific observed and expected splicing site mutation counts (\code{splc}),
                        gene-specific observed and expected Indel mutation counts (\code{indel}).
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
  # Not run due to the time it takes
  #MutScot(data, MutSigCV.out=MutSigCV.out)
}
\keyword{ mutation } 
\seealso{\code{\link{MutScotPatientHeter}}}
