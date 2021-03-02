# Train MMRDetect

MMRDetectTC.train <- function(mutationVariable, classification, cancerType = NULL) {
  
  
  ## match the data with classification
  trainset = mutationVariable
  trainset = merge(trainset, classification[,c("Sample","MSI_status")], by="Sample")
  
  if(nrow(trainset)<50){
    warning('training set size < 50')  
  }
  
  ## build model with trainset
  trainset$MSI_status<-as.factor(trainset$MSI_status)
  tree.model = J48(MSI_status ~ . , data=trainset[,c("CDel_rep_mean","CIns_rep_mean","TDel_rep_mean","TIns_rep_mean","MMR_sum","MSI_status")])
  
  tree.eval = evaluate_Weka_classifier(tree.model, numFolds = 5, seed = 1)
  cat('5 fold cross validation result: ',unname(tree.eval$detail[1]), '\n', sep='')
  tree.model
  
}

MMRDetect.train <- function(mutationVariable, classification, cancerType = NULL) {
  
  
  ## match the data with classification
  trainset = mutationVariable
  trainset = merge(trainset, classification[,c("Sample","MSI_status")], by="Sample")
  
  if(nrow(trainset)<50){
    warning('training set size < 50')  
  }
  
  ## build model with trainset
  trainset$MSI_status<-as.factor(trainset$MSI_status)
  tree.model = J48(MSI_status ~ . , data=trainset[,c("Del_rep_mean","Ins_rep_mean","MMR_sum","RepIndel_num","MSI_status")])
  
  tree.eval = evaluate_Weka_classifier(tree.model, numFolds = 5, seed = 1)
  cat('5 fold cross validation result: ',unname(tree.eval$detail[1]), '\n', sep='')
  tree.model
  
}


MMRDetect.train <- function(mutationVariable, classification, cancerType = NULL) {
  
  
  ## match the data with classification
  trainset = mutationVariable
  trainset = merge(trainset, classification[,c("Sample","MSI_status")], by="Sample")
  
  if(nrow(trainset)<50){
    warning('training set size < 50')  
  }
  
  ## build model with trainset
  trainset$MSI_status<-as.factor(trainset$MSI_status)
  tree.model = J48(MSI_status ~ . , data=trainset[,c("Del_rep_mean","Ins_rep_mean","MMR_sum","MSI_status")])
  
  tree.eval = evaluate_Weka_classifier(tree.model, numFolds = 5, seed = 1)
  cat('5 fold cross validation result: ',unname(tree.eval$detail[1]), '\n', sep='')
  tree.model
  
}
