MMRDetectTC.classify <- function(mutationVariable, classifier = MMRDclassifier) {
  classifyset = mutationVariable[,c("CDel_rep_mean","CIns_rep_mean","TDel_rep_mean","TIns_rep_mean","MMR_sum")]
  
  
  classify.result = predict(classifier, newdata=classifyset, type="class")
  result = as.data.frame(cbind(rownames(classifyset), 
                               as.character(classify.result)))
  colnames(result) = c("Sample", "MSI_status")
  result
} 

MMRDetect.classify <- function(mutationVariable, classifier = MMRDclassifier) {
  classifyset = mutationVariable[,c("Del_rep_mean","Ins_rep_mean","MMR_sum","RepIndel_num")]
  
  
  classify.result = predict(classifier, newdata=classifyset, type="class")
  result = as.data.frame(cbind(rownames(classifyset), 
                               as.character(classify.result)))
  colnames(result) = c("Sample", "MSI_status")
  result
} 



MMRDetect.classify <- function(mutationVariable, classifier = MMRDclassifier) {
  classifyset = mutationVariable[,c("Del_rep_mean","Ins_rep_mean","MMR_sum")]
  
  
  classify.result = predict(classifier, newdata=classifyset, type="class")
  result = as.data.frame(cbind(rownames(classifyset), 
                               as.character(classify.result)))
  colnames(result) = c("Sample", "MSI_status")
  result
} 
