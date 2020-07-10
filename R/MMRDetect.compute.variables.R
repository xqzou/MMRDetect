# Calculate variables for MMRDetect


MMRDetect.compute.variables <- function(subs, indels, tissue_type,MMR_sig_indel, tissue_subsig96, conversion_mtx){
  
  
  sub_summary <- data.frame(table(subs$Sample))
  names(sub_summary) <- c("Sample","sub_num")
  indel_classied <- indel_classifier(indels)
  
  indel_classied_rep <- indel_classied[!indel_classied$indeltype_short%in%c("[-]Mh","[->1]Other","[+]Mh","[+>1]Other","Complex","[-C]Rep=1","[-T]Rep=1","[+C]Rep=0","[+T]Rep=0", "[+C]Rep=1","[+T]Rep=1"),]
  indel_classied_rep_summary <- data.frame(table(indel_classied_rep$Sample))
  names(indel_classied_rep_summary) <- c("Sample","RepIndel_num")
  
  muts_summary <- merge(sub_summary,indel_classied_rep_summary,by="Sample")
  
  write.table(muts_summary,paste0("muts_summary_",tissue_type,".txt"),sep = "\t", col.names = T, row.names = F, quote = F)
  #sample_highburden <- muts_summary[muts_summary$sub_num>=190 & muts_summary$RepIndel_num>=80,"Sample"]
  
  
  # Generate catalouge for subs_highburden
  sub_catalouge <- GenCatalogue(subs,"Sample")
  sub_catalouge <- sub_catalouge[,-2]
  
  # Step 2:  Signature fitting for subs using Andrea's tissue-specific signatures 
  selected_tissueSig <- tissue_subsig96[,c(1,which(sub("_[^_]+$","",names(tissue_subsig96))==tissue_type))]
  
  mut_sig <- merge(selected_tissueSig,sub_catalouge,by="MutationType")
  mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
  mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
  sig_cat <- mut_sig[,2:dim(selected_tissueSig)[2]]
  mut_cat <- mut_sig[,(dim(selected_tissueSig)[2]+1):dim(mut_sig)[2]]
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  
  tissue_exposure <- t(a$E_median_filtered)
  write.table(tissue_exposure,paste0("TissueSigexposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  tissue_exposure <- read.table("TissueSigexposure_Cervix.txt", sep = "\t", header = T, as.is = T)
  # conversion matrix
  conversion_mtx_tissue <- conversion_mtx[rownames(conversion_mtx)%in% names(selected_tissueSig),] 
  conversion_mtx_tissue <- conversion_mtx_tissue[colnames(tissue_exposure),]
  
  # convert tissue-specific signatures to reference sigs
  refsig_exposure <- as.matrix(tissue_exposure) %*% as.matrix(conversion_mtx_tissue)
  
  write.table(refsig_exposure,paste0("exposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  # calculate cosine similarity of indel profile and  MLH1/MSH2/MSH6 PMS2 signatures 
  Sample_MMR <- refsig_exposure[,c("RefSig.MMR1","RefSig.MMR2")]
  
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  for(i in 2:dim(MMR_sig_indel)[2]){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==rownames(Sample_MMR)[j],]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_RepIndel2(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,rownames(Sample_MMR)[j])
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  
  MMRsig_2 <- merge(Del_rep_mean[,-2],Ins_rep_mean[,-2],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,indel_classied_rep_summary,by="Sample")
  
  Sample_MMR <- as.data.frame(Sample_MMR)
  Sample_MMR$Sample <- rownames(Sample_MMR)
  Sample_MMR$MMR_sum <- Sample_MMR$RefSig.MMR1+Sample_MMR$RefSig.MMR2
  MMRsig_2 <- merge(MMRsig_2,Sample_MMR,by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  return(MMRsig_2) 
}
