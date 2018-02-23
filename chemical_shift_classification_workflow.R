####CODE TO REPRODUCE STUDY

######
#Helper functions to source
calculation=function(initial_dataset,ind,ind2,metadata) {
  ctrl <- trainControl(method = "boot632",number=500,sampling="up",classProbs=T)
  ctrl2 <- trainControl(method = "boot632",number=500,sampling="up",classProbs=T,summaryFunction =twoClassSummary)
  
  variable_importance_subset=rep(0,ncol(initial_dataset))
  names(variable_importance_subset)=colnames(initial_dataset)
  classification_data=matrix(NA,3,3)
  ctrlrfe <- rfeControl(functions = rfFuncs,
                        method = "boot",
                        verbose = FALSE)
  iii=1
  dataset=initial_dataset
  colnames(dataset)=make.names(colnames(dataset))
  #Both kinds of information
  dataset$sample_type=factor(make.names(metadata))
  set.seed(iii);rfe_fr <- rfe(sample_type ~ .,data = dataset,
                              rfeControl = ctrlrfe)
  set.seed(iii);plsFit <- train(sample_type ~ .,data = dataset[,colnames(dataset) %in%c('sample_type',rfe_fr$optVariables)],method = "rf",trControl = ctrl)
  set.seed(iii);plsFit2 <- train(sample_type ~ .,data = dataset[,colnames(dataset) %in%c('sample_type',rfe_fr$optVariables)],method = "rf",trControl = ctrl2)
  classification_data[,1]=c(apply(plsFit$results[,2:3],2,max),max(plsFit2$results[,2]))
  #Variable importance information
  variable_importance_subset[colnames(dataset)[-ncol(dataset)] %in% rfe_fr$optVariables]=varImp(plsFit)$importance[,1]
  
  #Only quantification information
  dataset=initial_dataset[,ind]
  colnames(dataset)=make.names(colnames(dataset))
  dataset$sample_type=factor(make.names(metadata))
  set.seed(iii);rfe_fr <- rfe(sample_type ~ .,data = dataset,
                              rfeControl = ctrlrfe)
  set.seed(iii);plsFit <- train(sample_type ~ .,data = dataset[,colnames(dataset) %in%c('sample_type',rfe_fr$optVariables)],method = "rf",trControl = ctrl)
  set.seed(iii);plsFit2 <- train(sample_type ~ .,data = dataset[,colnames(dataset) %in%c('sample_type',rfe_fr$optVariables)],method = "rf",trControl = ctrl2)
  classification_data[,2]=c(apply(plsFit$results[,2:3],2,max),max(plsFit2$results[,2]))
  #Only chemical shift information
  dataset=initial_dataset[,ind2]
  colnames(dataset)=make.names(colnames(dataset))
  
  dataset$sample_type=factor(make.names(metadata))
  set.seed(iii);rfe_fr <- rfe(sample_type ~ .,data = dataset,
                              rfeControl = ctrlrfe)
  set.seed(iii);plsFit <- train(sample_type ~ .,data = dataset[,colnames(dataset) %in%c('sample_type',rfe_fr$optVariables)],method = "rf",trControl = ctrl)
  set.seed(iii);plsFit2 <- train(sample_type ~ .,data = dataset[,colnames(dataset) %in%c('sample_type',rfe_fr$optVariables)],method = "rf",trControl = ctrl2)
  classification_data[,3]=c(apply(plsFit$results[,2:3],2,max),max(plsFit2$results[,2]))
  # loc=rbind(loc,as.vector(classification_data))
  
  indicators=classification_data
  colnames(indicators)=c("Both kinds of information","Quantification information","Chemical shift information")
  rownames(indicators)=c("Accuracy", "kappa", "AUROC")
  var_imp=round(variable_importance_subset,3)
  var2=var_imp[order(var_imp,decreasing = T)]
  gr=list(indicators=indicators,var_imp=var_imp)
  return(gr)
}

preparation_dataset= function(fitting_error,signal_area_ratio,quantification,
                              chemical_shift,metadata) {
  #Removing all NA columns and converting quantifications with suboptimal quality indicatiors into NA values.
  ind=which(is.na(colSums(fitting_error)))
  quantification[signal_area_ratio<3]=NA
  quantification[,-ind][fitting_error[,-ind]>0.15]=NA
  chemical_shift[signal_area_ratio<3]=NA
  chemical_shift[,-ind][fitting_error[,-ind]>0.15]=NA
  chemical_shift[chemical_shift==Inf]=NA

  #Converting outliers for each signal for every group (smoker vs non-smoker) into NA
  for (j in unique(metadata)) {
    for (i in seq(ncol(quantification))) {
      quantification[(quantification[which(metadata==j),i] %in% boxplot.stats(quantification[which(metadata==j),i])$out),i]=NA
    }
    for (i in seq(ncol(chemical_shift))) {
      chemical_shift[(chemical_shift[which(metadata==j),i] %in% boxplot.stats(chemical_shift[which(metadata==j),i])$out),i]=NA
    }
  }



  #Remvoving signals with too many NA values
  ind=apply(quantification,2,function(x)length(which(is.na(x))))
  quantification_clean=quantification[,ind<(0.4*nrow(quantification))]
  ind=apply(chemical_shift,2,function(x)length(which(is.na(x))))
  chemical_shift_clean=chemical_shift[,ind<(0.4*nrow(quantification))]


  #Preparing specific colnames for signal chemical shifts
  colnames_chemical_shift_part1=sapply(colnames(chemical_shift_clean),function(x)substr(x,1,nchar(x)-2))
  colnames_chemical_shift_part2=round(colMeans(chemical_shift_clean,na.rm=T),3)
  colnames_chemical_shift=paste(colnames_chemical_shift_part1,' - chemical shift (',colnames_chemical_shift_part2,' ppm)',sep='')
  colnames(chemical_shift_clean)=colnames_chemical_shift

  #Choosing for every metabolite the signal with lowest number of NAs, and removing the TSP quntification
  ind=c()
  def=sapply(colnames(quantification_clean),function(x)substr(x,1,nchar(x)-2))
  def2=which(duplicated(def)==F)
  def4=apply(quantification_clean,2,function(x)length(which(is.na(x))))
  for (i in def2) {
    def3=which(def==def[i])
    ind=c(ind,def3[which.min(def4[def3])])
  }
  quantification_clean=quantification_clean[,ind]
  quantification_clean=quantification_clean[,-1]

  #Preparing specific colnames for metabolite quantifications
  colnames_quantification=sapply(colnames(quantification_clean),function(x)substr(x,1,nchar(x)-2))
  colnames_quantification=paste(colnames_quantification,'quantification',sep=' - ')
  colnames(quantification_clean)=colnames_quantification

  #Imputing NAs
  quantification_clean=missForest(data.matrix(quantification_clean))$ximp
  chemical_shift_clean=missForest(data.matrix(chemical_shift_clean))$ximp
  #Cleaning non-informative chemical shifts and performing dimensionality reduction through PCA
  ab=psych::alpha(scale(chemical_shift_clean),check.keys=T)
  chemical_shift_clean=chemical_shift_clean[,which(ab$alpha.drop$`S/N`<ab$total$`S/N`)]
  chemical_shift_pcs=prcomp(scale(chemical_shift_clean))$x[,1:5]
  colnames(chemical_shift_pcs)=paste(colnames(chemical_shift_pcs)," - Chemical shift",sep='')
  #Merging both quantification and  chemical shift datasets
  combined_dataset=as.data.frame(cbind(scale(quantification_clean),chemical_shift_pcs))
  li=list(combined_dataset=combined_dataset,ind=1:ncol(quantification_clean),ind2=(ncol(quantification_clean)+1):ncol(combined_dataset))
  return(li)
}

######
# Previous requirements
library(readr)
library(data.table)
library(caret)
library(readxl)
library(missForest)

#Change the working directory by the one where the script is located.
setwd("C:/Users/Usuario/Documents/GitHub/chemical_shift_classification")
# Download this zip: https://www.dropbox.com/s/9q9mygus57g0kfn/MTBLS_data_profiling.rar?dl=0 and copy its unzipped .RData files into the working directory.


###WORKFLOW FOR MTBLS374 DATASET
load("MTBLS374_console_profiling.RData")
#Reading metadata extracted from MTBLS374 website and converting into factor associated with MTBLS profiling dataset
metadata=fread("s_BoEfRTP2 Serum NMR.txt")
metadata=metadata[which(!is.na(metadata$`Factor Value[smoking status]`)),]
metadata=factor(metadata$`Factor Value[smoking status]`)

combined_mtbls374=preparation_dataset(profiling_data$fina$fitting_error,
                                      profiling_data$fina$signal_area_ratio,
                                      profiling_data$fina$quantification,
                                      profiling_data$fina$chemical_shift,metadata)

mtbls374_data_classification=calculation(combined_mtbls374$combined_dataset,
                                         combined_mtbls374$ind,
                                         combined_mtbls374$ind2,
                                         metadata)
write.csv(round(mtbls374_data_classification$var_imp[order(mtbls374_data_classification$var_imp,decreasing=T)],3),file="var_imp_mtbls374.csv")
write.csv(round(mtbls374_data_classification$indicators,3),file="classification_data_mtbls374.csv")

###WORKFLOW FOR MTBLS237 DATASET
set.seed(1)

#Reading metadata
metadata <- factor(as.vector(as.data.frame(read_csv("Metadata.csv"))[,3]))
metadata=droplevels(metadata[metadata!=4])

#Reading outputted profiling information (quantification quantification, chemical shift and quality indicators) of MTBLS374 dataset
load("MTBLS237_console_profiling.RData")


combined_mtbls237=preparation_dataset(profiling_data$fina$fitting_error,
                                      profiling_data$fina$signal_area_ratio,
                                      profiling_data$fina$quantification,
                                      profiling_data$fina$chemical_shift,metadata)

mtbls237_data_1=calculation(combined_mtbls237$combined_dataset[metadata!=3,],
                            combined_mtbls237$ind,
                            combined_mtbls237$ind2,
                            droplevels(metadata[metadata!=3]))

mtbls237_data_2=calculation(combined_mtbls237$combined_dataset[metadata!=2,],
                            combined_mtbls237$ind,
                            combined_mtbls237$ind2,
                            droplevels(metadata[metadata!=2]))


mtbls237_data_3=calculation(combined_mtbls237$combined_dataset[metadata!=1,],
                            combined_mtbls237$ind,
                            combined_mtbls237$ind2,
                            droplevels(metadata[metadata!=1]))

dummy=round(rbind(mtbls237_data_1$indicators,
                  mtbls237_data_2$indicators,
                  mtbls237_data_3$indicators),3)
write.csv(dummy,file="classification_data_mtbls237.csv")
dummy=cbind(mtbls237_data_1$var_imp,mtbls237_data_2$var_imp,mtbls237_data_3$var_imp)
dummy=round(dummy[order(rowMeans(dummy),decreasing=T),],3)
write.csv(dummy,file="var_imp_mtbls237.csv")

######

###WORKFLOW FOR MTBLS1 DATASET
load("MTBLS1_console_profiling.RData")
#Reading metadata extracted from MTBLS1 website and converting into factor associated with MTBLS profiling dataset
metadata=factor(rep(c(1,2),c(48,84)))

combined_mtbls1=preparation_dataset(profiling_data$fina$fitting_error,
                                      profiling_data$fina$signal_area_ratio,
                                      profiling_data$fina$quantification,
                                      profiling_data$fina$chemical_shift,metadata)

mtbls1_data_classification=calculation(combined_mtbls1$combined_dataset,
                                         combined_mtbls1$ind,
                                         combined_mtbls1$ind2,
                                         metadata)
write.csv(round(mtbls1_data_classification$var_imp[order(mtbls1_data_classification$var_imp,decreasing=T)],3),file="var_imp_mtbls1.csv")
write.csv(round(mtbls1_data_classification$indicators,3),file="classification_data_mtbls1.csv")
