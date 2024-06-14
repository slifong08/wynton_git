#BiocManager::install("DiffBind")
library(DiffBind)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

setwd("~/Library/CloudStorage/GoogleDrive-dcheng1@vt.edu/My Drive/UCSF/Ultrasound/4-Ultrasound/ATAC-seq/Diffbind")

cells=c('hob',"k562", 'bj', 'hepg2') 

for (cell in cells){
  print(cell)
  
  atac <- dba(minOverlap=2,
              sampleSheet=paste(cell,".csv",sep=''),
              config = data.frame(AnalysisMethod=DBA_DESEQ2,  #SF - this is default
                                  bUsePval=TRUE,  # SF - NOT default, should be False if using FDR
                                  doBlacklist=TRUE, #SF - this is default
                                  doGreylist=TRUE),#Greylist FALSE for ChIP data #SF - this is default
              bRemoveM = TRUE, #SF - this is default
              bRemoveRandom = TRUE)  #SF - this is default
  
  atac
  
  atac_count <- dba.count(atac,
                          score=DBA_SCORE_RPKM_FOLD, # not default
                          bScaleControl = TRUE, # dfault
                          summits = 200) #ATAC summit = 200bp (deafult); K27ac summit 400bp
  
  atac_count
  
  atac_contrast <- dba.contrast(atac_count,
                                categories = DBA_TREATMENT,
                                minMembers = 3) # use 2 for BJ as it only has 2 replicates
  
  atac_analyze <- dba.analyze(atac_contrast,
                              method=DBA_ALL_METHODS)
  
  atac.deseq2 <- dba.report(atac_analyze,
                            contrast = 1,
                            method = DBA_DESEQ2)
  #plot(atac.deseq2)
  
  sum(atac.deseq2$Fold>0)
  sum(atac.deseq2$Fold<0)
  
  
  peakAnno <- annotatePeak(atac.deseq2, 
                           tssRegion=c(-2000, 2000),
                           TxDb=txdb, 
                           annoDb="org.Hs.eg.db")
  
  plotAnnoBar(peakAnno)
  
  write.csv(data.frame(peakAnno), file=paste(cell,"deseq2.csv",sep='_'), row.names=FALSE)
}

