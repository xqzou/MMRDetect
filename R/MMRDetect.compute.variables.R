


#' Calculate cosine similarities between indel profile of a sample and KO, input is a list of classified indels
#'
#' @param SingleSample_classified_indels A list of classified indels
#' @param Sig indel signature of MMRD KOs
#' @return cosine similarities between indel profile of a sample and KO
#' @export
Generate_CossimVector_SingleSample_RepIndel <- function(SingleSample_classified_indels, Sig=MMRKO.indelsig){
  mut_catalogue <-  gen_indelmuttype_MMRD(SingleSample_classified_indels,"Sample","indeltype_short")
  return(Calculae_Cossim_catalogue_RepIndel(mut_catalogue, Sig))
}

Generate_CossimVector_SingleSample_CTRepIndel <- function(SingleSample_classified_indels, Sig=MMRKO.indelsigCT){
  mut_catalogue <-  gen_indelmuttype_MMRD(SingleSample_classified_indels,"Sample","indeltype_short")
  return(Calculae_Cossim_catalogue_CTRepIndel(mut_catalogue, Sig))
}


#' Calculate cosine similarities between indel profile of a sample and KO, input is indel profile of 45 channels
#'
#' @param SingleSample_classified_indels A list of classified indels
#' @param Sig indel signature of MMRD KOs
#' @return cosine similarities between indel profile of a sample and KO
#' @export
Calculae_Cossim_catalogue_RepIndel <- function(singlesample_indel_catalogue, Sig){
  total_subs_sig <- merge(singlesample_indel_catalogue,Sig, by="indelsubtype")
  
  cossim_del <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Del",3],total_subs_sig[total_subs_sig$type=="Del",4]))
  cossim_ins <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Ins",3],total_subs_sig[total_subs_sig$type=="Ins",4]))
  cossim <- abs(cos_similarity(total_subs_sig[,3],total_subs_sig[,4]))
  
  return(c(cossim_del, cossim_ins,cossim))
}

#' Calculate cosine similarities between indel profile of a sample and KO, input is indel profile of 45 channels
#'
#' @param SingleSample_classified_indels A list of classified indels
#' @param Sig indel signature of MMRD KOs
#' @return cosine similarities between indel profile of a sample and KO

#' @export
Calculae_Cossim_catalogue_CTRepIndel <- function(singlesample_indel_catalogue, Sig){
  total_subs_sig <- merge(singlesample_indel_catalogue,Sig, by="indelsubtype")
  total_subs_sig$type <- substr(total_subs_sig$indelsubtype, 2,3)
  cossim_Cdel <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="-C",3],total_subs_sig[total_subs_sig$type=="-C",4]))
  cossim_Cins <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="+C",3],total_subs_sig[total_subs_sig$type=="+C",4]))
  cossim_Tdel <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="-T",3],total_subs_sig[total_subs_sig$type=="-T",4]))
  cossim_Tins <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="+T",3],total_subs_sig[total_subs_sig$type=="+T",4]))
  cossim <- abs(cos_similarity(total_subs_sig[,3],total_subs_sig[,4]))
  
  return(c(cossim_Cdel, cossim_Cins, cossim_Tdel, cossim_Tins, cossim))
}



#' Calculate cosine similarities between indel profile of a sample and KO, input is indel profile of 45 channels
#'
#' @param SingleSample_classified_indels A list of classified indels
#' @param Sig indel signature of MMRD KOs
#' @return cosine similarities between C/T insertion and deletion profiles of a sample and KO

#' @export
Calculae_Cossim_catalogue_RepIndel <- function(singlesample_indel_catalogue, Sig){
  total_subs_sig <- merge(singlesample_indel_catalogue,Sig, by="indelsubtype")
  
  cossim_del <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Del",3],total_subs_sig[total_subs_sig$type=="Del",4]))
  cossim_ins <- abs(cos_similarity(total_subs_sig[total_subs_sig$type=="Ins",3],total_subs_sig[total_subs_sig$type=="Ins",4]))
  cossim <- abs(cos_similarity(total_subs_sig[,3],total_subs_sig[,4]))
  
  return(c(cossim_del, cossim_ins,cossim))
}


#' Calculate variables for MMRDetect
#'
#' @param subs A list of substitutions
#' @param indels A list of indels
#' @return Signature exposures
#' @examples
#' @export
MMRDetect.compute.variables <- function(subs, indels, tissue_type,MMR_subsig96,MMR_sig_indel, tissue_subsig96){
  
  
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
  
  
  # Compute similarity of tissue-specific signatures with MMR KO sigs
  selected_tissueSig <- tissue_subsig96[,c(1,which(sub("_[^_]+$","",names(tissue_subsig96))==tissue_type))]
  
  MMR1Sig <- c("Breast_A", "Colorectal_F", "Liver_E","Stomach_H", "Uterus_C", "Uterus_J")
  MMR2Sig <- c("Biliary_E", "Breast_D", "Colorectal_E", "Ovary_F", "Pancreas_H", "Stomach_B", "Uterus_D", "Kidney_D")
  
  tissue_MMR1Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR1Sig)]
  tissue_MMR2Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR2Sig)]
  
  
  if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)>0){
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)==0){
    # add PMS2
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)>0){
    # add MLH1
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR1")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)==0){
    
    # add MLH1 and PMS2 
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    names(selected_tissueSig)[dim(selected_tissueSig)[2]-1] <- c("MMR1")
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  write.table(a$E_median_filtered,paste0("exposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  
  # Using tissue-specific signatures to fit, but remove the ones that have high similarity with MMR KO signatures.
  #Sig_MMR <- selected_tissueSig[,which(names(selected_tissueSig) %in% ko_pc_cossim[ko_pc_cossim$similarity>=0.85,"CosmicSig"])]
  
  # Choose the samples that have exposure in MMR signatures
  #Exposure <- read.table("exposure_breast741_250breastsamples.txt", sep="\t",header = T, as.is=T)
  Exposure <-a$E_median_filtered
  sample_exposure <- as.data.frame(Exposure)
  MMRsig_sample <- sample_exposure[which(rownames(sample_exposure) %in% c(tissue_MMR1Sig, tissue_MMR2Sig,"MMR2","MMR1")),]
  #MMRsig_sample <- sample_exposure[which(rownames(sample_exposure) %in% ko_pc_cossim_MMR$CosmicSig),]
  # MMRsig_sample <- MMRsig_sample[,which(colSums(MMRsig_sample)>0)]
  
  MMRsig_sample$AndreaSig <- rownames(MMRsig_sample)
  #MMRsig_sample <- merge(MMRsig_sample,ko_pc_cossim_MMR_dcast[,c("CosmicSig","KO_sig")],by="CosmicSig")
  MMRsig_sample_melt <- melt(MMRsig_sample,c("AndreaSig"))
  names(MMRsig_sample_melt) <- c("AndreaSig","Sample","exposure")
  MMRsig_sample_melt_dcast <- dcast(MMRsig_sample_melt,Sample~AndreaSig,value.var="exposure")
  
  if(dim(MMRsig_sample_melt_dcast)[2]==2){
    MMRsig_sample_melt_dcast$MMR_sum <- MMRsig_sample_melt_dcast[,-1]
  }else{
    MMRsig_sample_melt_dcast$MMR_sum <- rowSums(MMRsig_sample_melt_dcast[,-1])
    
  }
  # MMRsig_sample_melt_dcast_realMMR <- MMRsig_sample_melt_dcast[MMRsig_sample_melt_dcast$MMR_sum>=190,]
  #write.table(MMRsig_sample_melt_dcast,"MMRexposure_breast741_175breastsamples.txt",sep = "\t",row.names = F, col.names = T, quote = F)
  
  
  # Step 3. Measure the cosine similarity between RepDel, RepIns and RepIndel
  # Using indel profile to discriminate MLH1/MSH2/MSH6 and PMS2 signatures 
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  Sample_MMR <- MMRsig_sample_melt_dcast
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==as.character(Sample_MMR[j,"Sample"]),]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  Indel_rep_MMR1 <- ddply(cossim_allsample[cossim_allsample$MMRgene!="PMS2",],c("Sample"),summarise,N=length(Sample),Indel_rep_MMR1_mean=mean(Indel_rep),Indel_rep_MMR1_sd=sd(Indel_rep))
  Indel_rep_MMR2 <- cossim_allsample[cossim_allsample$MMRgene=="PMS2",c("Sample","Indel_rep")]
  
  MMRsig_2 <- merge(Del_rep_mean[,-2],Ins_rep_mean[,-2],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,Indel_rep_MMR1[,-2],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,Indel_rep_MMR2,by="Sample")
  MMRsig_2 <- merge(MMRsig_2,indel_classied_rep_summary,by="Sample")
  names(MMRsig_2)[8]="Indel_rep_MMR2"
  
  MMRsig_2 <- merge(MMRsig_2,MMRsig_sample_melt_dcast,by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  return(MMRsig_2) 
}


#' Calculate variables for MMRDetect for indel catalogues
#'
#' @param indels A list of indels
#' @param tissue_type tissue types of the samples
#' @param MMR_sig_indel experimental MMRD signatures
#' @return Signature exposures
#' @exampl
#' @export
MMRDetect.compute.indelcatalogue.similarity <- function(indel_cat, tissue_type,MMR_sig_indel=MMRKO_indelsig){
  
  # Measure the cosine similarity between RepDel, RepIns and RepIndel
  # Using indel profile to discriminate MLH1/MSH2/MSH6 and PMS2 signatures 
  
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  Sample_MMR <- data.frame("Sample"=names(indel_cat)[3:dim(indel_cat)[2]], "indel_num"=colSums(indel_cat[,3:dim(indel_cat)[2]]))
  names(Sample_MMR) <- c("Sample","indel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_cat[,c("indelsubtype", "type",Sample_MMR[j,"Sample"])]
     
        cossim <- Calculae_Cossim_catalogue_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))

      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
 
  MMRsig_2 <- merge(Del_rep_mean[,c("Sample","Del_rep_mean")],Ins_rep_mean[,c("Sample","Ins_rep_mean")],by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  
  return(MMRsig_2) 
}

MMRDetect.compute.Repindelcatalogue.similarity <- function(indel_cat, tissue_type,MMR_sig_indel=MMRKO_indelsig){
  
  # Measure the cosine similarity between RepDel, RepIns and RepIndel
  # Using indel profile to discriminate MLH1/MSH2/MSH6 and PMS2 signatures 
  
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL

  Sample_MMR <- data.frame("Sample"=names(indel_cat)[3:dim(indel_cat)[2]], "RepIndel_num"=colSums(indel_cat[!indel_cat$indelsubtype %in%c("[-]Mh","[->1]Other","[+]Mh","[+>1]Other","Complex","[-C]Rep=1","[-T]Rep=1","[+C]Rep=0","[+T]Rep=0", "[+C]Rep=1","[+T]Rep=1"),3:dim(indel_cat)[2]]))
  names(Sample_MMR) <- c("Sample","RepIndel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_cat[,c("indelsubtype", "type",Sample_MMR[j,"Sample"])]
      
      cossim <- Calculae_Cossim_catalogue_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
      cossim_allsample <- rbind(cossim_allsample,cossim)
      cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
      cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
      
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  MMRsig_2 <- merge(Del_rep_mean[,c("Sample","Del_rep_mean")],Ins_rep_mean[,c("Sample","Ins_rep_mean")],by="Sample")
  MMRsig_2 <- merge(Sample_MMR,MMRsig_2, by="Sample")
  
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  
  return(MMRsig_2) 
}

#' Calculate variables for MMRDetect for indels: Sample, Chrom, Pos, Ref, Alt 
#'
#' @param indels A list of indels
#' @param tissue_type tissue types of the samples
#' @param MMR_sig_indel experimental MMRD signatures
#' @return Signature exposures
#' @examples
#' @export
MMRDetect.compute.indel.similarity <- function(indels, tissue_type,MMR_sig_indel=MMRKO_indelsig){
  
  indel_classied <- indel_classifier(indels)
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  Sample_MMR <- data.frame(table(indels$Sample))
  names(Sample_MMR) <- c("Sample","indel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==as.character(Sample_MMR[j,"Sample"]),]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  MMRsig_2 <- merge(Del_rep_mean[,c("Sample","Del_rep_mean")],Ins_rep_mean[,c("Sample","Ins_rep_mean")],by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)

  
  return(MMRsig_2) 
}


MMRDetect.compute.Repindel.similarity <- function(indels, tissue_type,MMR_sig_indel=MMRKO_indelsig){
  
  indel_classied <- indel_classifier(indels)
  
  indel_classied_rep <- indel2610_classied[!indel2610_classied$indeltype_short%in%c("[-]Mh","[->1]Other","[+]Mh","[+>1]Other","Complex","[-C]Rep=1","[-T]Rep=1","[+C]Rep=0","[+T]Rep=0", "[+C]Rep=1","[+T]Rep=1"),]
  Sample_MMR <- data.frame(table(indel_classied_rep$Sample))
  names(Sample_MMR) <- c("Sample","RepIndel_num")
  
  
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
 # Sample_MMR <- data.frame(table(indels$Sample))
#  names(Sample_MMR) <- c("Sample","indel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==as.character(Sample_MMR[j,"Sample"]),]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_RepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("Del_rep","Ins_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  Del_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Del_rep_mean=mean(Del_rep),Del_rep_sd=sd(Del_rep))
  Ins_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),Ins_rep_mean=mean(Ins_rep),Ins_rep_sd=sd(Ins_rep))
  
  MMRsig_2 <- merge(Del_rep_mean[,c("Sample","Del_rep_mean")],Ins_rep_mean[,c("Sample","Ins_rep_mean")],by="Sample")
  MMRsig_2 <- merge(Sample_MMR,MMRsig_2, by="Sample")
  write.table(MMRsig_2,paste0("sample_feature_summary_",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  return(MMRsig_2) 
}


MMRDetect.compute.indel.exposure_andrea_notworking <- function(indels, MMR_sig_indel=MMRKO_indelsig, outputname="exposure_indelMMRD.txt"){
 
  indel_classied <- indel_classifier(indels)
  indel_cat <- gen_indelmuttype_MMRD(indel_classied,"Sample","indeltype_short")
  Sig_tofit <- MMR_sig_indel
  mut_sig <- merge(Sig_tofit,indel_cat[,-2],by="indelsubtype")
  mut_sig <- mut_sig[-dim(mut_sig)[2]]
  sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
  mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  write.table(a$E_median_filtered,paste0("indel_sigexposure",".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  
  
  # Choose the samples that have exposure in MMR signatures
  Exposure <-a$E_median_filtered
  sample_exposure <- as.data.frame(Exposure)
  MMRsig_sample <- sample_exposure
  MMRsig_sample$MMRKOsig <- rownames(MMRsig_sample)
  MMRsig_sample_melt <- melt(MMRsig_sample,c("MMRKOsig"))
  names(MMRsig_sample_melt) <- c("MMRKOsig","Sample","exposure")
  MMRsig_sample_melt_dcast <- dcast(MMRsig_sample_melt,Sample~MMRKOsig,value.var="exposure")
  MMRsig_sample_melt_dcast$indelMMR_sum <- rowSums(MMRsig_sample_melt_dcast[,-1])
  write.table(MMRsig_sample_melt_dcast,outputname,sep = "\t",row.names = F, col.names = T, quote = F)
  
  
  return(MMRsig_sample_melt_dcast) 
}


# Using signature.tools.lib to calculate MMRD signature exposure, not working.
MMRDetect.compute.indelcatalogue.exposure_andrea_notworking <- function(indel_cat, MMR_sig_indel=MMRKO_indelsig, outputname="exposure_indelMMRD.txt"){
  Sig_tofit <- MMR_sig_indel
  mut_sig <- merge(Sig_tofit,indel_cat[,-2],by="indelsubtype")
  rownames(mut_sig) <- mut_sig$indelsubtype
  #mut_sig <- mut_sig[-dim(mut_sig)[2]]
  sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
  mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
  
#  a <- signature.tools.lib::SignatureFit_withBootstrap_Analysis("./",mut_cat,sig_cat,type_of_mutations="generic")
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat[,1:10],sig_cat, threshold_p.value = 0.01, method = "SA")

  write.table(a$E_median_filtered,paste0("indel_sigexposure",".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  Exposure <-a$E_median_filtered
  sample_exposure <- as.data.frame(Exposure)
  MMRsig_sample <- sample_exposure
  MMRsig_sample$MMRKOsig <- rownames(MMRsig_sample)
  MMRsig_sample_melt <- melt(MMRsig_sample,c("MMRKOsig"))
  names(MMRsig_sample_melt) <- c("MMRKOsig","Sample","exposure")
  MMRsig_sample_melt_dcast <- dcast(MMRsig_sample_melt,Sample~MMRKOsig,value.var="exposure")
  MMRsig_sample_melt_dcast$indelMMR_sum <- rowSums(MMRsig_sample_melt_dcast[,-1])
  write.table(MMRsig_sample_melt_dcast,outputname,sep = "\t",row.names = F, col.names = T, quote = F)
  
  
  return(MMRsig_sample_melt_dcast) 
}

# Detect the maximum exposure of signature in genome
MaxSigExposure <- function(mcat, sigx, nboot,p=0.9){

    m_quantile <- apply(mcat, 2, function(x) Bootstrap_boundary(x,nboot,p))
    m_count <- m_quantile/sigx
    m_mincount <- apply(m_count, 2, min)
    return(m_mincount)
  
}

# using bootstraping to identify the quantiles at give probabilities
# input: one genome vector
# return a vector of
Bootstrap_boundary <- function(x,nboot, p){
  
  Repx <- matrix(rep(x,nboot),ncol = nboot)
  bootstrapx <- apply(Repx, 2, function(x) rmultinom(1, sum(x), x))
  return(apply(bootstrapx, 1, function(x) quantile(x, probs = p)))
  
}



MMRDetect.compute.indelcatalogue.exposure <- function(indel_cat, MMR_sig_indel=MMRKO_indelsig, outputname="exposure_indelMMRD.txt"){
  Sig_tofit <- MMR_sig_indel
  mut_sig <- merge(Sig_tofit,indel_cat[,-2],by="indelsubtype")
  rownames(mut_sig) <- mut_sig$indelsubtype
  indelexposure <- data.frame("Sample"=names(indel_cat[,-c(1,2)]))
  # Use +-T at 7,8,9
  mut_sig <- mut_sig[c("[-T]Rep=7","[-T]Rep=8","[-T]Rep=9","[+T]Rep=7","[+T]Rep=8","[+T]Rep=9"),]
  sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
  mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
  
  mut_cat_nonzero <- mut_cat[,which(colSums(mut_cat)>0)]
  
  a <- apply(sig_cat, 2, function(x) MaxSigExposure(mut_cat_nonzero,x,nboot=100,p=0.9))
  
  MMRsig_sample <- as.data.frame(a)
  MMRsig_sample$indelMMR <- apply(MMRsig_sample,1,max)
  MMRsig_sample$Sample <- rownames(MMRsig_sample)
  
  indelexposure <- merge(indelexposure,MMRsig_sample,by="Sample",all.x=T)
  
  write.table(MMRsig_sample,outputname,sep = "\t",row.names = T, col.names = T, quote = F)
  
  
  return(MMRsig_sample) 
}

#' Calculate variables for MMRDetect
#'
#' @param subs A list of substitutions
#' @param indels A list of indels
#' @return Signature exposures
#' @examples
#' @export
MMRDetect.compute.sub.exposure <- function(subs, tissue_type,MMR_sig_sub,tissue_subsig96){
  
 # Generate catalouge for subs_highburden
  sub_catalouge <- GenCatalogue(subs,"Sample")
  sub_catalouge <- sub_catalouge[,-2]
  
  # Step 2:  Signature fitting for subs using Andrea's tissue-specific signatures 
  
  selected_tissueSig <- tissue_subsig96[,c(1,which(sub("_[^_]+$","",names(tissue_subsig96))==tissue_type))]
  
  MMR1Sig <- c("Breast_A", "Colorectal_F", "Liver_E","Stomach_H", "Uterus_C", "Uterus_J")
  MMR2Sig <- c("Biliary_E", "Breast_D", "Colorectal_E", "Ovary_F", "Pancreas_H", "Stomach_B", "Uterus_D", "Kidney_D")
  
  tissue_MMR1Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR1Sig)]
  tissue_MMR2Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR2Sig)]
  
  
  if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)>0){
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)==0){
    # add PMS2
    selected_tissueSig <- merge(selected_tissueSig,MMR_sig_sub[,c("MutationType","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)>0){
    # add MLH1
    selected_tissueSig <- merge(selected_tissueSig,MMR_sig_sub[,c("MutationType","MLH1")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR1")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)==0){
    
    # add MLH1 and PMS2 
    selected_tissueSig <- merge(selected_tissueSig,MMR_sig_sub[,c("MutationType","MLH1","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    names(selected_tissueSig)[dim(selected_tissueSig)[2]-1] <- c("MMR1")
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  write.table(a$E_median_filtered,paste0("exposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  

  # Choose the samples that have exposure in MMR signatures
  Exposure <-a$E_median_filtered
  sample_exposure <- as.data.frame(Exposure)
  MMRsig_sample <- sample_exposure[which(rownames(sample_exposure) %in% c(tissue_MMR1Sig, tissue_MMR2Sig,"MMR2","MMR1")),]

  MMRsig_sample$AndreaSig <- rownames(MMRsig_sample)
  MMRsig_sample_melt <- melt(MMRsig_sample,c("AndreaSig"))
  names(MMRsig_sample_melt) <- c("AndreaSig","Sample","exposure")
  MMRsig_sample_melt_dcast <- dcast(MMRsig_sample_melt,Sample~AndreaSig,value.var="exposure")
  
  if(dim(MMRsig_sample_melt_dcast)[2]==2){
    MMRsig_sample_melt_dcast$MMR_sum <- MMRsig_sample_melt_dcast[,-1]
  }else{
    MMRsig_sample_melt_dcast$MMR_sum <- rowSums(MMRsig_sample_melt_dcast[,-1])
    
  }
  write.table(MMRsig_sample_melt_dcast,paste0("subsig_exposure_",tissue_type,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
   
  return(MMRsig_sample_melt_dcast) 
}


#' Calculate variables for MMRDetect
#'
#' @param subs A dataframe of substitution 96-channel catalogue
#' @param indels A list of indels
#' @return Signature exposures
#' @examples
#' @export
MMRDetect.compute.subcatalogue.exposure <- function(sub_cat, tissue_type,MMR_subsig96,tissue_subsig96){
  
  
  # Generate catalouge for subs_highburden
  sub_catalouge <- sub_cat

  
  # Step 2:  Signature fitting for subs using Andrea's tissue-specific signatures 
  
  selected_tissueSig <- tissue_subsig96[,c(1,which(sub("_[^_]+$","",names(tissue_subsig96))==tissue_type))]
  
  MMR1Sig <- c("Breast_A", "Colorectal_F", "Liver_E","Stomach_H", "Uterus_C", "Uterus_J")
  MMR2Sig <- c("Biliary_E", "Breast_D", "Colorectal_E", "Ovary_F", "Pancreas_H", "Stomach_B", "Uterus_D", "Kidney_D")
  
  tissue_MMR1Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR1Sig)]
  tissue_MMR2Sig <- names(selected_tissueSig)[which(names(selected_tissueSig) %in%MMR2Sig)]
  
  
  if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)>0){
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)>0 & length(tissue_MMR2Sig)==0){
    # add PMS2
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)>0){
    # add MLH1
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR1")
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
    
  }else if(length(tissue_MMR1Sig)==0 & length(tissue_MMR2Sig)==0){
    
    # add MLH1 and PMS2 
    selected_tissueSig <- merge(selected_tissueSig,MMR_subsig96[,c("MutationType","MLH1","PMS2")])
    names(selected_tissueSig)[dim(selected_tissueSig)[2]] <- c("MMR2")
    names(selected_tissueSig)[dim(selected_tissueSig)[2]-1] <- c("MMR1")
    
    Sig_tofit <- selected_tissueSig
    mut_sig <- merge(Sig_tofit,sub_catalouge,by="MutationType")
    mut_sig$Mutation <- substr(mut_sig$MutationType,3,5)
    mut_sig <- mut_sig[order(mut_sig$Mutation),-dim(mut_sig)[2]]
    sig_cat <- mut_sig[,2:dim(Sig_tofit)[2]]
    mut_cat <- mut_sig[,(dim(Sig_tofit)[2]+1):dim(mut_sig)[2]]
    
  }
  
  a <- signature.tools.lib::SignatureFit_withBootstrap(mut_cat,sig_cat)
  write.table(a$E_median_filtered,paste0("exposure_",tissue_type,".txt"),sep = "\t",row.names = T, col.names = T, quote = F)
  
  
  
  # Choose the samples that have exposure in MMR signatures
  Exposure <-a$E_median_filtered
  sample_exposure <- as.data.frame(Exposure)
  MMRsig_sample <- sample_exposure[which(rownames(sample_exposure) %in% c(tissue_MMR1Sig, tissue_MMR2Sig,"MMR2","MMR1")),]
  
  MMRsig_sample$AndreaSig <- rownames(MMRsig_sample)
  MMRsig_sample_melt <- melt(MMRsig_sample,c("AndreaSig"))
  names(MMRsig_sample_melt) <- c("AndreaSig","Sample","exposure")
  MMRsig_sample_melt_dcast <- dcast(MMRsig_sample_melt,Sample~AndreaSig,value.var="exposure")
  
  if(dim(MMRsig_sample_melt_dcast)[2]==2){
    MMRsig_sample_melt_dcast$MMR_sum <- MMRsig_sample_melt_dcast[,-1]
  }else{
    MMRsig_sample_melt_dcast$MMR_sum <- rowSums(MMRsig_sample_melt_dcast[,-1])
    
  }
  write.table(MMRsig_sample_melt_dcast,paste0("subsig_exposure_",tissue_type,".txt"),sep = "\t",row.names = F, col.names = T, quote = F)
  
  
  return(MMRsig_sample_melt_dcast) 
}


cos_similarity <- function(v1,v2){
  v1v2 <- sum(v1*v2)
  v1_length <- sqrt(sum(v1*v1))
  v2_length <- sqrt(sum(v2*v2))
  return(v1v2/v1_length/v2_length)
}

MMRDetect.compute.subcatalogue.similarity <- function(mcat,scat){
  mut_sig <- merge(scat,mcat,by="MutationType")
  sig_cat <- mut_sig[,2:dim(scat)[2]]
  mut_cat <- mut_sig[,(dim(scat)[2]+1):dim(mut_sig)[2]]
  mut_cat <- mut_cat/colSums(mut_cat)[col(mut_cat)]
  a <- apply(sig_cat, 2, function(x) SigCossim(mut_cat,x))
  subsim <- as.data.frame(a)
  subsim$Sample <- rownames(subsim)
  subsim$Sample <- sub(".DNA","-DNA",subsim$Sample)
  subsim$maxcossim <- apply(subsim[,-which(names(subsim)=="Sample")],1,max)
  
  return(subsim)
  
}

SigCossim <- function(mcat, sigx){
  
  return(apply(mcat, 2, function(x) cos_similarity(x,sigx)))
}




#' Calculate variables for MMRDetect for indel catalogues
#'
#' @param indels A list of indels
#' @param tissue_type tissue types of the samples
#' @param MMR_sig_indel experimental MMRD signatures
#' @return Signature exposures
#' @exampl
#' @export
MMRDetect.compute.indelcatalogue.CTsimilarity <- function(indel_cat, tissue_type,MMR_sig_indel=MMRKO_indelsigCT){
  
  # Measure the cosine similarity between RepDel, RepIns and RepIndel
  # Using indel profile to discriminate MLH1/MSH2/MSH6 and PMS2 signatures 
  
  indel_classied <- indel_cat
  #write.table(indel_classied,"indel_classied.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  # Measure the cosine similarity between RepDel, RepIns
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  Sample_MMR <- data.frame(table(indels$Sample))
  names(Sample_MMR) <- c("Sample","indel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==as.character(Sample_MMR[j,"Sample"]),]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_CTRepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("CDel_rep","CIns_rep","TDel_rep","TIns_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  CDel_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),CDel_rep_mean=mean(CDel_rep),Del_rep_sd=sd(CDel_rep))
  CIns_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),CIns_rep_mean=mean(CIns_rep),Ins_rep_sd=sd(CIns_rep))
  
  TDel_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),TDel_rep_mean=mean(TDel_rep),Del_rep_sd=sd(TDel_rep))
  TIns_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),TIns_rep_mean=mean(TIns_rep),Ins_rep_sd=sd(TIns_rep))
  
  
  MMRsig_2 <- merge(CDel_rep_mean[,c("Sample","CDel_rep_mean")],CIns_rep_mean[,c("Sample","CIns_rep_mean")],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,TDel_rep_mean[,c("Sample","TDel_rep_mean")],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,TIns_rep_mean[,c("Sample","TIns_rep_mean")],by="Sample")
  
  
  write.table(MMRsig_2,paste0("summary_indel_cossim",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  
  return(MMRsig_2) 
}


#' Calculate variables for MMRDetect for indels: Sample, Chrom, Pos, Ref, Alt 
#'
#' @param indels A list of indels
#' @param tissue_type tissue types of the samples
#' @param MMR_sig_indel experimental MMRD signatures
#' @return Signature exposures
#' @examples
#' @export
MMRDetect.compute.indel.CTsimilarity <- function(indels, tissue_type,MMR_sig_indel=MMRKO_indelsigCT){
  
  indel_classied <- indel_classifier(indels)
  #write.table(indel_classied,"indel_classied.txt", sep = "\t", col.names = T, row.names = F, quote = F)
  # Measure the cosine similarity between RepDel, RepIns
  cossim_allsample <- NULL
  cossim_MMR <- NULL
  cossim_sample <- NULL
  Sample_MMR <- data.frame(table(indels$Sample))
  names(Sample_MMR) <- c("Sample","indel_num")
  for(i in 2:5){
    for(j in 1:dim(Sample_MMR)[1]){
      current_sample_indel <- indel_classied[indel_classied$Sample==as.character(Sample_MMR[j,"Sample"]),]
      if(dim(current_sample_indel)[1]>1){
        cossim <- Generate_CossimVector_SingleSample_CTRepIndel(current_sample_indel,MMR_sig_indel[,c(1,i)])    
        cossim_allsample <- rbind(cossim_allsample,cossim)
        cossim_MMR <- c(cossim_MMR,names(MMR_sig_indel)[i])
        cossim_sample <- c(cossim_sample,as.character(Sample_MMR[j,"Sample"]))
        
      }
      
    }
    
  }
  cossim_allsample <- data.frame(cossim_allsample)
  names(cossim_allsample) <- c("CDel_rep","CIns_rep","TDel_rep","TIns_rep","Indel_rep")
  cossim_allsample$Sample <- cossim_sample
  cossim_allsample$MMRgene <- cossim_MMR
  
  CDel_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),CDel_rep_mean=mean(CDel_rep),Del_rep_sd=sd(CDel_rep))
  CIns_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),CIns_rep_mean=mean(CIns_rep),Ins_rep_sd=sd(CIns_rep))
  
  TDel_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),TDel_rep_mean=mean(TDel_rep),Del_rep_sd=sd(TDel_rep))
  TIns_rep_mean <- ddply(cossim_allsample,c("Sample"),summarise,N=length(Sample),TIns_rep_mean=mean(TIns_rep),Ins_rep_sd=sd(TIns_rep))
  
  
  MMRsig_2 <- merge(CDel_rep_mean[,c("Sample","CDel_rep_mean")],CIns_rep_mean[,c("Sample","CIns_rep_mean")],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,TDel_rep_mean[,c("Sample","TDel_rep_mean")],by="Sample")
  MMRsig_2 <- merge(MMRsig_2,TIns_rep_mean[,c("Sample","TIns_rep_mean")],by="Sample")
  
  
  write.table(MMRsig_2,paste0("summary_indel_cossim",tissue_type,".txt"),sep="\t", col.names = T, row.names = F, quote = F)
  
  
  return(MMRsig_2) 
}
