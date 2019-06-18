
\name{MutScot-package}
\alias{MUTSCOT}
\docType{package}
\title{
 MutScot
}
\description{
 The MutScot package for identifying driver genes
}
\details{
Identifying cancer driver genes is essential for understanding mechanisms of carcinogenesis and designing therapeutic strategies. 
Consequently, a set of driver genes has been identified for each cancer types, assumed to be identical across subtypes. 
This assumption may not hold, and the sets of driver genes are possibly distinct across cancer subtypes. 
The MutScot package includes a statistical framework that identifies driver genes and utilizes patient information
 to investigate subtype heterogeneity of driver genes. The function \code{\link{MutScot}} identifies driver genes, and 
 the function \code{\link{MutScotPatientHeter}} is used for conducting a subtype heterogeneity test.
}
%\references{
%  Choo-Wosoba H., Albert P.S., Zhu B.
%  A Hidden Markov Modeling Approach for Identifying Tumor Subclones in Next-Generation Sequencing Studies.
%}
\author{Ho-Hsiang Wu and Bin Zhu <bin.zhu@nih.gov>}
\keyword{package}