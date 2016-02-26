fdrGSEA <- function(gsStatsAll,gsStatsAllPerm,nGenes,signMethod) {
   
   # Pre-allocate:
   res <- list()
   pValuesAllUpAdj <- vector()
   pValuesAllDnAdj <- vector()
   
   # For each contrast/comparison:
   for(iContrast in 1:ncol(gsStatsAll)) {
   
      # Get random background matrix (for all gene-sets 'sizes' and permutations):
      randBgMat <- gsStatsAllPerm[[iContrast]]
      
      # For each gene-set (or gene-set size):
      randBgMatNorm <- randBgMat
      for(iGeneSetSize in 1:nrow(randBgMat)) {
         
         # Normalize background:
         geneSetBg <- randBgMat[iGeneSetSize,]
         posMean <- max(c(mean(geneSetBg[geneSetBg >= 0]),0),na.rm=TRUE)
         negMean <- max(c(abs(mean(geneSetBg[geneSetBg <= 0])),0),na.rm=TRUE)
         geneSetBg[geneSetBg > 0] <- geneSetBg[geneSetBg > 0]/posMean
         geneSetBg[geneSetBg < 0] <- geneSetBg[geneSetBg < 0]/negMean
         randBgMatNorm[iGeneSetSize,] <- geneSetBg
      }
      
      # For each gene set:
      NES <- rep(NA,nrow(gsStatsAll))
      randBgMatNormFull <- as.data.frame(matrix(nrow=nrow(gsStatsAll),ncol=ncol(randBgMat)))
      
      message("before for")
      message(date())
      
      
      for(iGeneSet in 1:nrow(gsStatsAll)) {
         
         # Create whole background matrix (for all 'gene-sets' and permutations)
         if(signMethod == "geneperm") {
            sizeGS <- nGenes[iGeneSet,iContrast]
            randBgMatNormFull[iGeneSet,] <- randBgMatNorm[as.character(sizeGS),]
         } else {
            randBgMatNormFull[iGeneSet,] <- randBgMatNorm[iGeneSet,]
         }
         
         # Normalize ES:
         ES <- gsStatsAll[iGeneSet,iContrast]
         if(signMethod == "geneperm") {
            tmp <- randBgMat[as.character(sizeGS),]
         } else {
            tmp <- randBgMat[iGeneSet,]
         }
         
         #if(ES > 0 & max(tmp) <= 0) stop("normalization of enrichment scores failed, increase the number of permutations")
         #else if(ES <= 0 & min(tmp) >= 0) stop("normalization of enrichment scores failed, increase the number of permutations")
         
         if(ES > 0) {
            NES[iGeneSet] <- ES/max(c(mean(tmp[tmp >= 0]),0),na.rm=TRUE) # if no backgound, set mean to 0
         } else if(ES < 0) {
            NES[iGeneSet] <- ES/max(c(abs(mean(tmp[tmp <= 0])),0),na.rm=TRUE) # if no backgound, set mean to 0
         } else {
            NES[iGeneSet] <- 0  
         }
      }
       
      message("afer for")
      message(date())
      
      # Negated for edgeCases of findInterval below
      randBgNormUpNeg <- sort(-randBgMatNormFull[randBgMatNormFull >= 0])
      randBgNormDn <- sort(randBgMatNormFull[randBgMatNormFull <= 0])
      NES.up <- sort(-NES[NES >= 0])
      NES.dn <- sort(NES[NES <= 0])
      
      
      # For each gene set:
      FDRup <- rep(NA,nrow(gsStatsAll))
      FDRdn <- rep(NA,nrow(gsStatsAll))
      
      {
          iGeneSet <- which(NES > 0)
          iGeneSet <- iGeneSet[order(-NES[iGeneSet])]
          tmp1 <- findInterval(-NES[iGeneSet], randBgNormUpNeg) / length(randBgNormUpNeg)
          tmp2 <- findInterval(-NES[iGeneSet], NES.up) / length(NES.up)
          FDRup[iGeneSet] <- tmp1/tmp2
      }
      
      {
          iGeneSet <- which(NES < 0)
          iGeneSet <- iGeneSet[order(NES[iGeneSet])]
          tmp1 <- findInterval(NES[iGeneSet], randBgNormDn) / length(randBgNormDn)
          tmp2 <- findInterval(NES[iGeneSet], NES.dn) / length(NES.dn)
          FDRdn[iGeneSet] <- tmp1/tmp2
      }
      
      {
          iGeneSet <- which(NES == 0)
          FDRup[iGeneSet] <- 1  
          FDRdn[iGeneSet] <- 1
      }
      
      FDRup[which(FDRup > 1)] <- 1
      FDRdn[which(FDRdn > 1)] <- 1
      
      message("afer for2")
      message(date())
      # Save FDR in matrix:
      pValuesAllUpAdj <- cbind(pValuesAllUpAdj,FDRup)
      pValuesAllDnAdj <- cbind(pValuesAllDnAdj,FDRdn)
   }
   
   # Return results:
   res$pValuesAllUpAdj <- pValuesAllUpAdj
   res$pValuesAllDnAdj <- pValuesAllDnAdj
   return(res)
}
