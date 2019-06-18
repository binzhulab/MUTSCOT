Sub192_Pat_Rate = function(np, 
                                   ng, 
                                   Pat,
                                   dndsout, 
                                   RefCDS, 
                                   substmodel,
                                   MutSigCV.out,
                                   gene_L, known_cancergenes, covs
        ){ 
          
          Lambda_g = matrix(0, ncol=1,   nrow=ng)
          
          n_non_sil_pm = matrix(0, ncol=192, nrow=np)
          
          n_sil_pm = matrix(0, ncol=192, nrow=np)
          
          n_ind_pm = matrix(0, ncol=1,   nrow=np)
          
          
          E_non_sil_pm = matrix(0, ncol=192, nrow=1)
          
          E_sil_pm = matrix(0, ncol=192, nrow=1)
          
          E_ind_pm = matrix(0, ncol=1,   nrow=1)
          
          
          ######################### silent model  
          
          n_silent = subset( dndsout$annotmuts, dndsout$annotmuts$impact== "Synonymous" )
          
          n_silent_pat = n_silent[,1]
          
          trinucsubsind=trinuc()
          
          n_silent_192 = trinucsubsind[paste(n_silent$ref3_cod, n_silent$mut3_cod, sep = ">")]
          
          n_silent = cbind.data.frame(n_silent_pat, n_silent_192)
          
          n_silent$n_silent_pat = factor(n_silent$n_silent_pat, levels=as.factor(Pat))
          
          n_silent$n_silent_192 = factor(n_silent$n_silent_192, levels=c(1:192))
          
          n_sil_pm = (as.matrix(table(n_silent$n_silent_pat, n_silent$n_silent_192), ncol=192, nrow=np))
          
          n_sil_pm = n_sil_pm
          
          rm(list = c("n_silent", "n_silent_pat", "n_silent_192") ) 
          gc()
          
          ######################### Non-silent model  
          
          n_non_silent = subset( dndsout$annotmuts, dndsout$annotmuts$impact%in%c("Missense", "Nonsense", "Essential_Splice") )
          
          n_non_silent_pat = n_non_silent[,1]
          
          n_non_silent_192 = trinucsubsind[paste(n_non_silent$ref3_cod, n_non_silent$mut3_cod, sep = ">")]
          
          n_non_silent = cbind.data.frame(n_non_silent_pat, n_non_silent_192)
          
          n_non_silent$n_non_silent_pat = factor(n_non_silent$n_non_silent_pat, levels=as.factor(Pat))
          
          n_non_silent$n_non_silent_192 = factor(n_non_silent$n_non_silent_192, levels=c(1:192))
          
          n_non_sil_pm = (as.matrix(table(n_non_silent$n_non_silent_pat, n_non_silent$n_non_silent_192), ncol=192, nrow=np))
          
          n_non_sil_pm = n_non_sil_pm
          
          rm(list = c("n_non_silent", "n_non_silent_pat", "n_non_silent_192") ) 
          gc()
          
          ######################### Indel model observation, n_ind_pm
          
          n_Ind = subset( dndsout$annotmuts, dndsout$annotmuts$impact== "no-SNV" )
          
          IndL_ptLoc = paste0(n_Ind$sampleID, n_Ind$pos)
          
          rmIndl = which(IndL_ptLoc%in%IndL_ptLoc[which(duplicated(IndL_ptLoc))])
          
          if(length(rmIndl)>0){n_Ind=n_Ind[-rmIndl,]}
          
          Indels_record <- n_Ind
          
          n_Ind$sampleID =factor(n_Ind$sampleID, levels=as.factor(Pat))
          
          n_ind_pm = (as.matrix(table(n_Ind$sampleID), ncol=1, nrow=np))
          
          rm(list = c("rmIndl", "IndL_ptLoc", "n_Ind"))
          gc()
          
          
          #################  Aggregate expectation of silent mutations of p by m
          
          ###### Silent mutation 
          
          par = dndsout$mle_submodel
          
          poissmodel = dndsout$poissmodel$model
          
          parmle = setNames(par[, 2], par[, 1])
          
          mutrates = sapply(substmodel[, 1], function(x) prod(parmle[base::strsplit(x, split = "\\*")[[1]]]))
          
          ####### Indel mutation
          
          indels = subset(dndsout$annotmuts, dndsout$annotmuts$impact=="no-SNV")
          
          geneindels = as.data.frame(array(0, dim = c(length(RefCDS),  8)))
          
          colnames(geneindels) = c( "gene_name",                # gene name
                                    "n_ind",                    # number of indel mutation
                                    "n_induniq",                # number of unique indel mutation, different patient but same indel will not be counted
                                    "n_indused",                # == n_ind if no want to remove unique, ==n_induniq instead
                                    "cds_length",               # length of coding DNA sequence of gene
                                    "excl",                     # the list of 558 cancer genes in the Cancer Gene Census version 73 are excluded when fitting the negative binomial model
                                    "exp_unif",                 # prediction with no cv
                                    "exp_indcv")                # prediction with cv
          
          geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
          
          geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[,1]])
          
          geneindels[is.na(geneindels[, 2]), 2] = 0
          
          geneindels$n_induniq = as.numeric(table(unique(indels[,-1])$gene)[geneindels[, 1]])
          
          geneindels[is.na(geneindels[, 3]), 3] = 0
          
          geneindels$n_indused = geneindels[, 3]
          
          geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
          
          geneindels$excl = (geneindels[, 1] %in% known_cancergenes)
          
          if(sum(geneindels[!geneindels$excl, "n_indused"]) ==0){
            
            geneindels$excl = F
            
          }
          
          geneindels$exp_unif = sum(geneindels[!geneindels$excl,"n_indused"])/sum(geneindels[!geneindels$excl,"cds_length"])*geneindels$cds_length
          
          nbrdf = cbind(geneindels[, c("n_indused", "exp_unif")], covs)[!geneindels[, 6], ]     # combine the cv, remove the known cancer driver gene
          
          model = suppressWarnings(glm.nb(n_indused ~ offset(log(exp_unif))+., data = nbrdf))   # Fit the negative binomial nb model
          
          if(!is.null(model$th.warn) | nrow(dndsout$genemuts) < 500){   # if less mutation then fit no model
            
            model = glm.nb(n_indused ~ offset(log(exp_unif)) - 1, data = nbrdf) # Fit no model
            
          }
        
          nbrdf_all = cbind(geneindels[, c("n_indused", "exp_unif")], covs) # Put back known driver gene for prediction
          
          theta_indels = model$theta
          
          nbregind = model
          
          geneindels$exp_indcv = exp(predict(model, nbrdf_all)) 
          
          rm(nbrdf_all, model, nbregind, theta_indels)
          gc()
          
          xL=sapply(RefCDS, function(x) x$L)
          
          genemuts <- as.matrix(dndsout$genemuts[, -1, drop=FALSE])
          shape    <- dndsout$nbreg$theta 
          if (length(MutSigCV.out)) {
            kmatch   <- match(gene_L, MutSigCV.out$gene)
            tmp      <- is.na(kmatch)
            if (any(tmp)) kmatch[tmp] <- 0 
          } else {
            kmatch <- rep(0, ng)
          }
          ivec1    <- 1:192
          ivec2    <- 193:384
          ivec3    <- 385:576
          ivec4    <- 577:768
          sumXL1   <- colSums(xL[ivec1, ])*np 
          lambda   <- rep(NA, ng)
       
          for(j in 1:ng){ ###  Gene  Loop  
        
            y          <- genemuts[j, ] 
            exp_rel    <- y[5:8]/y[5]                                                    
            scale      <- y[9]/shape 
            indneut_PS <- 1
            sumy       <- sum(y[indneut_PS])
            sume       <- sum(exp_rel[indneut_PS])
            opt_t_PS   <- mle_tcv(n_neutral=sumy, exp_rel_neutral=sume, shape=shape, scale=scale)
            #mrfold_PS  <- max(1e-10, opt_t_PS/y[5])
            
            k <- kmatch[j]
            if (k) {
              mutsig <- MutSigCV.out[k,] 
              NBG    <- mutsig$X 
              BG     <- mutsig$x/NBG
            } else {
              BG     <- 0
              NBG    <- 1
            }

            sumXLj    <- sumXL1[j]
            r1        <-  y[1]/sumXLj    
            r1Flag    <- r1 != 0         
            vr1       <- ((r1*(1-r1))/sumXLj)^r1Flag
            r2        <- BG
            r2Flag    <- r2 != 0
            vr2       <- ((r2*(1-r2))/NBG)^r2Flag
            r3        <- opt_t_PS/sumXLj 
            vr3       <- var_tcv(n_neutral=sumy, exp_rel_neutral=sume, shape=shape, scale=scale)*(1/sumXLj)^2 
            lambda[j] <- sum(r1/vr1, r2/vr2, r3/vr3)/sum(r1Flag/vr1, r2Flag/vr2, (r3!=0)/vr3)              
          }

          Lambda_g[, 1]     <- lambda
          tmp               <- matrix(lambda, byrow=TRUE, nrow=length(ivec2), ncol=ncol(xL))
          E_non_sil_pm[1, ] <- rowSums(tmp*(xL[ivec2, ] + xL[ivec3, ] + xL[ivec4, ])*np)
          E_sil_pm[1, ]     <- rowSums(tmp*xL[ivec1, ]*np)
          E_ind_pm[1,1]     <- sum(as.numeric(geneindels$exp_indcv))  

          rm(lambda, xL, sumXL1, ivec1, ivec2, ivec3, ivec4)
          gc()           

          ##### Estimating Tau_m and Delta_p
          
          Xm = colSums(n_sil_pm)
          Ym = colSums(n_non_sil_pm)
          names(Xm)=NULL
          names(Ym)=NULL
          
          Sm = E_sil_pm/np
          Nm = E_non_sil_pm/np
          
          MutRat <- (Xm+Ym+Xm[192:1]+Ym[192:1])/( (Sm+Nm+Sm[192:1]+Nm[192:1])*np )
          
          
          Xp = rowSums(n_sil_pm)
          Yp = rowSums(n_non_sil_pm)
          
          summ = sum(MutRat*(Sm+Nm))
          
          PatRat <-  (Xp+Yp)/summ
          
          
          ##### Indel Model
          XI = n_ind_pm[,1]
          NI = E_ind_pm/np      
          PatRatI <- XI/as.numeric(NI) 
          
          
          ######################## assign weight to substitution model
          wt = rep(1, np)
          
          percentile = ecdf(PatRat)
          
          quantl = seq(ceiling(percentile(1)*10),10,1)/10
   
          for(qq in 1:(length(quantl)-1) ){
            
            eval(parse(text=paste0("Hgp",qq,"=which(PatRat>quantile(PatRat,",quantl[qq],") & PatRat<=quantile(PatRat,",quantl[qq+1],"))")))
            
            eval(parse(text=paste0("Hw",qq,"=1/mean(PatRat[Hgp",qq,"])")))
            
            eval(parse(text=paste0("wt[Hgp",qq,"]=Hw",qq)))
          }
          
          wt <- wt
          
          
          ######################## assign weight to Indel model
          wtI = rep(1, np)
          
          percentileI = ecdf(PatRatI)
          
          quantlI = seq(ceiling(percentileI(1)*10),10,1)/10
          
          for(qq in 1:(length(quantlI)-1) ){
            
            eval(parse(text=paste0("HgpI",qq,"=which(PatRatI>quantile(PatRatI,",quantlI[qq],") & PatRatI<=quantile(PatRatI,",quantlI[qq+1],"))")))
            
            eval(parse(text=paste0("HwI",qq,"=1/mean(PatRatI[HgpI",qq,"])")))
            
            eval(parse(text=paste0("wtI[HgpI",qq,"]=HwI",qq)))
          }
          
          wtI <- wtI
          
          ret <- list(Indels_record=Indels_record, geneindels=geneindels, Lambda=Lambda_g, 
                      MutRat=MutRat, PatRat=PatRat, PatRatI=PatRatI, wt=wt, wtI=wtI)

          ret
          
        } # END: Sub192_Pat_Rate


dNdS_model <- function(np, gene_L, Pat, substmodel, Indels_record, geneindels, 
                 RefCDS, dndsout, Lambda_g, MutRat, PatRat, PatRatI, wt, wtI){
          
          ng <- length(gene_L)
          Fun_imp= c("silent", "miss", "nsen", "splc", "Indel")
          
          Fun_imp_complete = c("Synonymous", "Missense", "Nonsense", "Essential_Splice")
      
          for(v in 2: 5 ){
            
            eval(parse(text=paste0("count_ALL_", Fun_imp[v] ,"=lapply(1:ng, matrix, data= 0, nrow=1, ncol=2)")))
            
          }
        
          par        = dndsout$mle_submodel
          
          poissmodel = dndsout$poissmodel$model
          
          parmle     = setNames(par[, 2], par[, 1])
          
          mutrates   = sapply(substmodel[, 1], function(x) prod(parmle[base::strsplit(x, split = "\\*")[[1]]]))
          
          xL=sapply(RefCDS, function(x) x$L)
      
          fn2Vec     <- ((2-1)*192+1):(2*192)
          fn3Vec     <- ((3-1)*192+1):(3*192)
          fn4Vec     <- ((4-1)*192+1):(4*192)
          sumXL2     <- colSums(xL[fn2Vec, ])
          sumXL3     <- colSums(xL[fn3Vec, ])
          sumXL4     <- colSums(xL[fn4Vec, ])
          levels     <- as.factor(Pat)  
          geneind    <- dndsout$annotmuts$gene
          sampleID0  <- factor(dndsout$annotmuts$sampleID, levels=levels)
          impact0    <- dndsout$annotmuts$impact
          Igeneind   <- Indels_record$geneind 
          IsampleID0 <- factor(Indels_record$sampleID, levels=levels)
          wt3        <- wt^3
          sumPwt3    <- sum(PatRat*wt3) 
          wtI3       <- wtI^3

          tmp    <- as.matrix(xL[fn2Vec, , drop=FALSE])*matrix(data=MutRat, nrow=length(fn2Vec), ncol=ng)
          cm2Vec <- Lambda_g[,1]*colSums(tmp)*sumPwt3
          tmp    <- as.matrix(xL[fn3Vec, , drop=FALSE])*matrix(data=MutRat, nrow=length(fn3Vec), ncol=ng)
          cm3Vec <- Lambda_g[,1]*colSums(tmp)*sumPwt3
          tmp    <- as.matrix(xL[fn4Vec, , drop=FALSE])*matrix(data=MutRat, nrow=length(fn4Vec), ncol=ng)
          cm4Vec <- Lambda_g[,1]*colSums(tmp)*sumPwt3
          rm(fn2Vec, fn3Vec, fn4Vec)
          gc()

          imVec   <- as.numeric(geneindels$exp_indcv)/np*sum(PatRatI*wtI3)
          kmatch  <- match(gene_L, geneind)
          tmp     <- is.na(kmatch)
          if (any(tmp)) kmatch[tmp] <- 0
          Ikmatch <- match(gene_L, Igeneind)
          tmp     <- is.na(Ikmatch)
          if (any(tmp)) Ikmatch[tmp] <- 0

          for(j in 1:ng){ ###  Gene  Loop              
            if (Lambda_g[j, 1] > 0) {
              if (kmatch[j]) {

                #tmp        <- geneind == j
                tmp        <- geneind == gene_L[j]
                sampleID   <- sampleID0[tmp]   
                impact     <- impact0[tmp]     

                if (sumXL2[j] > 0) {
                  tmp                 <- impact == Fun_imp_complete[2]  
                  non_silent_ID       <- sampleID[tmp]
                  Cm11                <- sum(table(non_silent_ID)*wt3)
                  Cm12                <- cm2Vec[j]
                  count_ALL_miss[[j]] <- count_ALL_miss[[j]] + c(Cm11, Cm12)
                } 
                if (sumXL3[j] > 0) {
                  tmp                 <- impact == Fun_imp_complete[3]  
                  non_silent_ID       <- sampleID[tmp]
                  Cm11                <- sum(table(non_silent_ID)*wt3)
                  Cm12                <- cm3Vec[j]
                  count_ALL_nsen[[j]] <- count_ALL_nsen[[j]] + c(Cm11, Cm12)
                } 
                if (sumXL4[j] > 0) {
                  tmp                 <- impact == Fun_imp_complete[4]  
                  non_silent_ID       <- sampleID[tmp]
                  Cm11                <- sum(table(non_silent_ID)*wt3)
                  Cm12                <- cm4Vec[j]
                  count_ALL_splc[[j]] <- count_ALL_splc[[j]] + c(Cm11, Cm12)
                } 
              }
              ######################### Indel model 
           
              #n_Ind                <- Indels_record[Igeneind %in% j, , drop=FALSE] 
              #n_Ind_ID             <- factor(n_Ind$sampleID, levels=levels)
              if (Ikmatch[j]) {
                #n_Ind_ID             <- IsampleID0[Igeneind == j]
                n_Ind_ID             <- IsampleID0[Igeneind == gene_L[j]] 
                Ind11                <- sum(as.numeric(table(n_Ind_ID)*wtI3))
                Ind12                <- imVec[j] 
                count_ALL_Indel[[j]] <- count_ALL_Indel[[j]] + c(Ind11, Ind12)
              }
            }           
          } 
          

          ret <- list(miss=count_ALL_miss,  nsen=count_ALL_nsen,
                      splc=count_ALL_splc,  indel=count_ALL_Indel)

          ret

} # END: dNdS_model

getBenchmarkGenes <- function() {

  ret <- c("AKT1", "ARID1A", "BRCA1", "CASP8", "CBFB", "CDH1", "CDKN1B", "CHD4", 
           "CTCF", "ERBB2", "FBXW7", "FOXA1", "GATA3", "GPS2", "KMT2C", "KRAS", 
           "MAP2K4", "MAP3K1", "NCOR1", "NF1", "PIK3CA", "PIK3R1", "PTEN", 
           "PTPRD", "RB1", "RUNX1", "SF3B1", "TBX3", "TP53")
  ret

} # END: getBenchmarkGenes

MutScot <- function(data, MutSigCV.out=NULL, options=NULL) {
            
          options <- default.list(options, 
                     c("multitesting", "maxMutPerPat", "alpha.q", "benchmarkGenes"), 
                     list("fdr", 300, 0.1, NULL))              
          if (!is.data.frame(data)) stop("ERROR: data must be a data frame")
          if (ncol(data) != 5) stop("ERROR: data must have five columns")
          
          if (is.null(MutSigCV.out))  message("MutSigCV.out not input")
          
          #### Filter out hyper mutators (mutations per patient >300)
          Total_mut_count_pat = table(data$sampleID)
          
          Hyper = names(Total_mut_count_pat)[which(Total_mut_count_pat>= options$maxMutPerPat)]
          
          data <- subset(data, !(data$sampleID%in%Hyper) )
          
          rm(list=c("Hyper", "Total_mut_count_pat"))
          gc()

          message(" Running dndscv algorithm..") 
          dndsout <- dndscv(data)

          # gene ID
          gene_L = sapply(RefCDS, function(x) x$gene_name)
          # number of genes
          ng <- length(RefCDS)

          # Patient ID after removing hyper mutators
          Pat = unique(dndsout$annotmuts$sampleID)
          # sample size
          np = length(Pat)
          
          # calculating patient/192 mutation rate  
          message(" calculating patient/192 mutation rate ..")
          ret192.1 <- Sub192_Pat_Rate(np, 
                          ng, 
                          Pat,
                          dndsout, 
                          RefCDS, 
                          substmodel,
                          MutSigCV.out,
                          gene_L, known_cancergenes, covs) 

          # processing observed-to-expected matrix
          message(" processing observed-to-expected matrix ..")

          ret1 <- dNdS_model(np, 
                     gene_L, 
                     Pat, 
                     substmodel, 
                     ret192.1$Indels_record, 
                     ret192.1$geneindels, 
                     RefCDS, 
                     dndsout, 
                     ret192.1$Lambda, 
                     ret192.1$MutRat,
                     ret192.1$PatRat, 
                     ret192.1$PatRatI, 
                     ret192.1$wt, 
                     ret192.1$wtI)
          
          # testsing driver genes 
          message(" testing driver genes ..")
          pval <- Test(np, 
               ng, 
               ret1$miss, 
               ret1$nsen, 
               ret1$splc, 
               ret1$indel)       

          # summary results
          message(" summary results ..")  
          MutScot.q <- p.adjust(pval[, 1], method=options$multitesting)
          MutScot.q <- data.frame(Gene=gene_L, Pvalue=pval[, 1], Pvalue.adj=MutScot.q, stringsAsFactors=FALSE)                     
          MutScot.q <- MutScot.q[order(pval[, 1]), ]
          rownames(MutScot.q) <- NULL
 
          g.mutscot <- as.character(MutScot.q[MutScot.q$Pvalue.adj < options$alpha.q,1])
          g.dndscv  <- as.character(subset(dndsout$sel_cv, qglobal_cv < options$alpha.q)[,1])
          if (!is.null(MutSigCV.out)) { 
            g.mutsig <- as.character(subset(MutSigCV.out, q < options$alpha.q)[,1]) 
          } else {
            g.mutsig <- NULL
          }
          g.benchmark <- options[["benchmarkGenes", exact=TRUE]]
          if (is.null(g.benchmark)) g.benchmark <- getBenchmarkGenes()
          genes <- sort(unique(c(g.mutscot, g.dndscv, g.mutsig, g.benchmark)))
          genes <- data.frame(Gene=genes, stringsAsFactors=FALSE)
          genes[, "Benchmark"]  <- genes[, "Gene"] %in% g.benchmark
          genes[, "dNdScv"]     <- genes[, "Gene"] %in% g.dndscv
          genes[, "MutSigCV"]   <- genes[, "Gene"] %in% g.mutsig
          genes[, "MutScot"]    <- genes[, "Gene"] %in% g.mutscot
          rownames(genes)       <- NULL
                   
                CntMatrix = matrix(0, ncol=4, nrow = 4)
                colnames(CntMatrix)=c("Benchmark", "dNdScv", "MutSigCV", "MutScot")
                row.names(CntMatrix)=c("Total", "TP", "FP", "F1")
                DL <- genes
                CntMatrix[1,]= colSums(DL[,2:5], na.rm=T)
                CntMatrix[2,]= c(colSums(DL[,2,drop=F], na.rm=T), length(which(DL[,2]==1&DL[,3]==1)), 
                                 length(which(DL[,2]==1&DL[,4]==1)), length(which(DL[,2]==1&DL[,5]==1)))
                CntMatrix[3,]= CntMatrix[1,]-CntMatrix[2,]
                CntMatrix[4,]= c(100, round(F1(CntMatrix[2,2], CntMatrix[3,2], CntMatrix[1,1]),2)*100, 
                                 round(F1(CntMatrix[2,3], CntMatrix[3,3], CntMatrix[1,1]),2)*100, 
                                 round(F1(CntMatrix[2,4], CntMatrix[3,4], CntMatrix[1,1]),2)*100 )
              
        
  ret <- list(driverGenes=g.mutscot, significantGenes=genes, MutScotGenes=MutScot.q, countMat=CntMatrix, 
              dndsout=dndsout, patientIds=Pat, genes=gene_L, mutationRates=ret192.1, obsExp=ret1)
  class(ret) <- "MutScot"

  ret
        
} # END: MutScot 

checkPatientInfo <- function(x, ids) {

  if (is.null(x)) return(NULL)
  if (!is.data.frame(x)) stop("ERROR: patientInfo must be a data frame")
  nc <- ncol(x)
  if (nc < 2) stop("ERROR: patientInfo must contain at least two columns")
  for (i in 2:nc) {
    if (!is.numeric(x[, i])) stop(paste("ERROR: column ", i, " of patientInfo must be numeric", sep=""))
  }
  tmp <- ids %in% x[, 1]
  if (!any(tmp)) stop("ERROR: no matching ids in patientInfo")
  if (!all(tmp)) warning("Not all ids found in patientInfo")
  ids <- ids[tmp]
  tmp <- match(ids, x[, 1])
  x   <- x[tmp, , drop=FALSE]
 
  x

} # END: checkPatientInfo

MutScotPatientHeter <- function(obj, patientInfo, testGenes, options=NULL) {

            if (class(obj) != "MutScot") stop("ERROR: obj must be of class MutScot")
            testGenes       <- unique(testGenes)
            if (!length(testGenes)) stop("ERROR: testGenes must be specified")
            
            options         <- default.list(options, c("MCT"), list(1e4))
            MCT             <- options$MCT
            if (MCT < 100) stop("ERROR: options$MCT must be at least 100")
            Pat             <- obj$patientIds
            
            # Returned object from cgeckPatientInfo will have the same order as Pat
            patientInfo     <- checkPatientInfo(patientInfo, Pat)
            PatRat          <- obj$mutationRates$PatRat
            PatRatI         <- obj$mutationRates$PatRatI  
            wt              <- obj$mutationRates$wt
            wtI             <- obj$mutationRates$wtI 
            
            # Pat contains all ids from call to MutScot, subset objects if needed
            tmp             <- Pat %in% patientInfo[, 1] 
            if (!all(tmp)) {
              Pat           <- Pat[tmp]
              PatRat        <- PatRat[tmp]
              PatRatI       <- PatRatI[tmp]  
              wt            <- wt[tmp]
              wtI           <- wtI[tmp]
            }
            
            count_ALL_miss  <- obj$obsExp$miss
            count_ALL_nsen  <- obj$obsExp$nsen
            count_ALL_splc  <- obj$obsExp$splc
            count_ALL_Indel <- obj$obsExp$indel        
            np              <- length(Pat)
            gene_L          <- obj$genes
                                      
            Fun_imp      <- c("silent", "miss", "nsen", "splc", "Indel")
            pvFN_ALL     <- matrix(1, nrow=length(testGenes), ncol=4)
            list_of_glm2 <- list() 
            list_of_glm3 <- list()
            list_of_glm4 <- list()
            list_of_glm5 <- list()
            
            ######################### Non-silent model  
            tmp       <- obj$dndsout$annotmuts$gene %in% testGenes
            n_pat_mis <- subset( obj$dndsout$annotmuts, (obj$dndsout$annotmuts$impact%in%c("Missense"))&(tmp) )[,c("sampleID", "gene", "impact")]
            n_pat_non <- subset( obj$dndsout$annotmuts, (obj$dndsout$annotmuts$impact%in%c("Nonsense"))&(tmp) )[,c("sampleID", "gene", "impact")]
            n_pat_spl <- subset( obj$dndsout$annotmuts, (obj$dndsout$annotmuts$impact%in%c("Essential_Splice"))&(tmp) )[,c("sampleID", "gene", "impact")]
            n_pat_Ind <- subset( obj$dndsout$annotmuts, (obj$dndsout$annotmuts$impact== "no-SNV")&(tmp) )[,c("sampleID", "gene", "impact")]
            
            tmp              <- as.factor(Pat)
            misSampleID      <- factor(n_pat_mis$sampleID, levels=tmp)
            nonSampleID      <- factor(n_pat_non$sampleID, levels=tmp)
            splSampleID      <- factor(n_pat_spl$sampleID, levels=tmp)
            IndSampleID      <- factor(n_pat_Ind$sampleID, levels=tmp)
            wtI3             <- wtI^3
            wt3              <- wt^3 
            FunDat           <- as.data.frame(matrix(0, ncol=3, nrow=np))
            colnames(FunDat) <- paste0("V", c(1,2,3) )
            ivec             <- 1:np
            sumWP            <- sum(wt3*PatRat)
            sumWPI           <- sum(wtI3*PatRatI)
            Z                <- as.matrix(patientInfo[,-1])
            WZ               <- cbind(1, Z) 
            NR               <- nrow(Z)
            NC               <- ncol(Z)
                 
            for(L in 1: length(testGenes) ){
               
              cpL1 <- testGenes[L] 
              j    <- which(gene_L==cpL1)
              if (length(j)) {
                n_pat_mis_g <- misSampleID[n_pat_mis$gene==cpL1]
                n_pat_non_g <- nonSampleID[n_pat_non$gene==cpL1]
                n_pat_spl_g <- splSampleID[n_pat_spl$gene==cpL1]
                n_pat_Ind_g <- IndSampleID[n_pat_Ind$gene==cpL1]
              
                for(fn2 in 2:5){ 
                
                  if (fn2 < 5) {
                    PatRatt <- PatRat
                    wtt     <- wt3    
                    denom   <- sumWP
                    if (fn2==2) {
                      fnCt  <- table(n_pat_mis_g)
                      count <- count_ALL_miss[[j]]
                    } else if (fn2==3) {
                      fnCt  <- table(n_pat_non_g)                
                      count <- count_ALL_nsen[[j]]  
                    } else if (fn2==4) {
                      fnCt  <- table(n_pat_spl_g)
                      count <- count_ALL_splc[[j]]
                    }
                  } else {
                    PatRatt <- PatRatI
                    wtt     <- wtI3
                    count   <- count_ALL_Indel[[j]]  
                    fnCt    <- table(n_pat_Ind_g)
                    denom   <- sumWPI
                  }

                  FunDat[ivec ,1] <- fnCt
                  FunDat[ivec ,2] <- count[,2]*PatRatt/denom
                  FunDat[ivec ,3] <- wtt
                              
                  if ( sum(FunDat$V1) >= sum(FunDat$V2) ){
                           
                    glm0   = try(glm(V1~ -1 + V3 + offset(log(V2)), data=FunDat, family=poisson(link="log")  ), silent = T)
                  
                    if (!("try-error" %in% class(glm0))) {
                   
                      EP      <- fitted.values(glm0)                    
                      tmp     <- FunDat$V1-EP
                      Score1  <- t(wtt)%*%tmp
                      Score2  <- t(Z)%*%tmp
                      Score3  <- c(Score1, Score2)
                      WZ[, 1] <- wtt 
               
                      rc       <- as.integer(-1)
                      rv       <- as.numeric(rep(0, (NC+1)*(NC+1)))
                      tmp      <- .C("C_getScoreMat", as.numeric(wtt), as.numeric(Z), as.numeric(EP),
                                      as.integer(NR), as.integer(NC), ret_code=rc, ret_vec=rv)
                      if (tmp$ret_code) next                          
                      tmp      <- matrix(tmp$ret_vec, nrow=NC+1, ncol=NC+1, byrow=FALSE)
                      inVScore <-  try(solve(tmp), silent=TRUE)
                      if ("try-error" %in% class(inVScore)) next  
                      TestStat  = t(Score3)%*%inVScore%*%Score3 
                               
                      eval(parse(text=paste0("list_of_glm", fn2 ,"[[L]] = (try(glm(V1~ -1+offset(wtt)+offset(log(V2))+Z, data=FunDat, family=poisson(link=\"log\")  ), silent = T))" )))
                    
                      rc  <- as.integer(-1)
                      rp  <- as.numeric(-1)     
                      tmp <- .C("C_getP", as.numeric(TestStat), as.integer(MCT), as.numeric(wtt), as.numeric(EP), as.numeric(Z), 
                                 as.integer(NR), as.integer(NC), as.numeric(inVScore), ret_code=rc, ret_p=rp)
                      if (tmp$ret_code) next                 
                      pvFN_ALL[L,(fn2-1)] <- tmp$ret_p  
                     
                    }
                  
                  } else {
                    eval(parse(text=paste0("list_of_glm", fn2 ,"[[L]] <- NULL" )))
             
                  }
                }# fn2 loop end  
                  
              } else {
                # test gene not found
                list_of_glm2[[L]] <- NULL
                list_of_glm3[[L]] <- NULL
                list_of_glm4[[L]] <- NULL
                list_of_glm5[[L]] <- NULL
                pvFN_ALL[L, ]     <- NA 
              }     
            }# gene_L loop end
            
            Sol = apply( pvFN_ALL, 1, function (x) -2*sum(log(x)) )
 
            Sol = matrix(Sol, ncol=1)
            
            Solp = matrix(apply(Sol, 1, function(x) (1-pchisq(x, df=8)) ), ncol=1)
             
            pmat           <- data.frame(Gene=testGenes, x=pvFN_ALL, x2=Solp)
            colnames(pmat) <- c("Gene", "miss.p", "nsen.p", "splc.p", "indel.p", "subtype.het.p")
            rownames(pmat) <- NULL 
    
          ret <- list(pvalues=pmat,
                      glm_miss=list_of_glm2,
                      glm_nsen=list_of_glm3,
                      glm_splc=list_of_glm4,
                      glm_indel=list_of_glm5)
          
} #END: MutScotPatientHeter

F1 <- function(TP,FP,TTP){
          
          Recall=TP/TTP
          
          Precis=TP/(TP+FP)
          
          F1S=2*(Precis*Recall)/(Precis+Recall)
          
          return(F1S)
          
} # END: F1

trinuc <- function(tri=192){
          
          nt = c("A", "C", "G", "T")
          
          trinucs = paste(rep(nt, each = 16, times = 1), rep(nt, each = 4, times = 4), rep(nt, each = 1, times = 16), sep = "")
          
          trinucinds = setNames(1:64, trinucs)
          
          trinucsubs = NULL
          
          for (j in 1:length(trinucs)) {
            
            trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j], 1, 1), setdiff(nt, substr(trinucs[j], 2, 2)), substr(trinucs[j], 3, 3), sep = ""), sep = ">"))
            
          }
          
          trinucsubsind = setNames(1:192, trinucsubs)
          
          return(trinucsubsind) 
          
} # END: trinuc

mle_tcv <- function(n_neutral, exp_rel_neutral, shape, scale){
          
          tml = (n_neutral + shape - 1)/(exp_rel_neutral + (1/scale))
          
          if (shape <= 1) {
            
            tml = max(shape * scale, tml)
            
          }
          
          return(tml)
          
} # END: mle_tcv 
        
var_tcv <- function(n_neutral, exp_rel_neutral, shape, scale){
          
          tml = (n_neutral  + shape) / (exp_rel_neutral + (1/scale))^2
          
          if (shape <= 1) {
            
            tml = max(shape * scale^2, tml)
            
          }
          
          return(tml)
          
} # END: var_tcv

Test <- function(np, ng, count_ALL_miss, count_ALL_nsen, count_ALL_splc, count_ALL_Indel){
          
          count <- NULL
          Fun_imp= c("silent", "miss", "nsen", "splc", "Indel")
          
          pval_PS   = matrix(1, nrow=ng, ncol=1)
          
          pvFN_ALL  = matrix(0, nrow=ng, ncol=4)
          
          for(fn2 in 2:5){  # Functional Lopp
            
            eval(parse(text=paste0("count = count_ALL_", Fun_imp[fn2] )))
            
            for(j in 1:ng){ # Gene Loop  
              
              obC = round(sum(count[[j]][,1]),0)
              
              exC = sum(count[[j]][,2]) 
        
              if(obC > exC){
                
                pvFN_ALL[j,(fn2-1)] = 1 - ppois(obC, lambda=exC) + dpois(obC, lambda=exC)  
                
              }else{
                
                pvFN_ALL[j,(fn2-1)] = 1  
                
              }
              
            }
            
          }
          
          PPP = apply( pvFN_ALL, 1, function (x) -2*sum(log(x)) )
          
          PPP = matrix(PPP, ncol=1)
          
          pval_PS[, 1]= matrix(apply(PPP, 1, function(x) (1-pchisq(x, df=8)) ), ncol=1)
          
          pval <- pval_PS   

          pval
          
} # END: Test  

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list

