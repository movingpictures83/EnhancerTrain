#!/usr/bin/env Rscript
######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
######################################################################################

##  Training enhancer prediction model from annotated regions (known enhancers
##  and known negative regions)
##  Please email Yupeng He (yupeng.he.bioinfo at gmail) for feedback, question and bug. 


## option parsing
suppressPackageStartupMessages(library("REPTILE",verbose=FALSE))


dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {
	#option_parser = get_option_parser_training()
#options <- optparse::parse_args(option_parser)
###No input
#if(length(commandArgs(TRUE)) == 0){
#    optparse::print_help(option_parser)
#    quit("no")
#}
##End of option parsing
	######################################################################################
## Input from command line options
##Data information file

data_info_file = paste(pfix, parameters["data_info_file", 2], sep="/")# options$data_info_file
##Epimark file of DMRs
DMR_epimark_file = paste(pfix, parameters["DMR_epimark_file", 2], sep="/")#options$DMR_epimark_file
##Epimark file of annotated regions
annotated_region_epimark_file = paste(pfix, parameters["annotated_region_epimark_file", 2], sep="/")#options$annotated_region_epimark_file
##Label file
label_file = paste(pfix, parameters["label_file", 2], sep="/")#options$label_file
##Samples in which the activities of annotated regions are used for training
samples_for_training = unlist(strsplit(parameters["samples_for_training", 2],","))
##Reference sample(s)
#ref_sample = options$ref_sample
#if(!is.null(ref_sample)){
#    ref_sample = unlist(strsplit(options$ref_sample,","))
#}
ref_sample = NULL
##Prefix of output files
output_prefix = outputfile
##Classifier family to use
classifier_family = ""
num_trees = 1
incl_dev = 0
#classifier_family = options$classifier_family
##Numbers of tree in random forest classifier
#num_trees = options$num_trees
##Whether to compute deviation as additional features
#incl_dev = options$incl_dev

## Read data informaiton
data_info = read.table(data_info_file,header=T,stringsAsFactors=F)[,1:2]

epimark_DMR = NULL
label_DMR = NULL
epimark_region = NULL
label_region = NULL
for(query_sample in samples_for_training){
    if(!is.null(DMR_epimark_file)){
        ## Read the epigenomic signature of DMRs
        epimark_DMR_sample <- read_epigenomic_data(data_info,
                                                   DMR_epimark_file,
                                                   query_sample=query_sample,
                                                   ref_sample=ref_sample,
                                                   incl_dev=incl_dev
                                                   )
        DMR_region_id = sapply(strsplit(rownames(epimark_DMR_sample),":"),function(x) return(x[2]))
        rownames(epimark_DMR_sample) <- paste0(rownames(epimark_DMR_sample),
                                               ":",
                                               query_sample,"_",DMR_region_id)
        
        label_DMR_sample = read_label(label_file,query_sample)
        label_DMR_sample = label_DMR_sample[DMR_region_id]
        
        ##Ignore rows with NAs
        label_DMR_sample = label_DMR_sample[rowSums(is.na(epimark_DMR_sample))==0]
        epimark_DMR_sample = epimark_DMR_sample[rowSums(is.na(epimark_DMR_sample))==0,]
        epimark_DMR_sample = epimark_DMR_sample[!is.na(label_DMR_sample),]
        label_DMR_sample = label_DMR_sample[!is.na(label_DMR_sample)]
        
        epimark_DMR = rbind(epimark_DMR,epimark_DMR_sample)
        names(label_DMR_sample) = paste0(query_sample,"_",names(label_DMR_sample))
        label_DMR = c(label_DMR,label_DMR_sample)
    }
    ## Read the epigenomic signature of query regions
    epimark_region_sample <- read_epigenomic_data(data_info,
                                                  annotated_region_epimark_file,
                                                  query_sample=query_sample,
                                                  ref_sample=ref_sample,
                                                  incl_dev=incl_dev
                                                  )
    region_id = rownames(epimark_region_sample)
    rownames(epimark_region_sample) <- paste0(rownames(epimark_region_sample),
                                              ":",
                                              query_sample,"_",rownames(epimark_region_sample))
    
    label_region_sample = read_label(label_file,query_sample)
    label_region_sample = label_region_sample[region_id]

    ##Ignore rows with NAs
    label_region_sample = label_region_sample[rowSums(is.na(epimark_region_sample))==0]
    epimark_region_sample = epimark_region_sample[rowSums(is.na(epimark_region_sample))==0,]
    epimark_region_sample = epimark_region_sample[!is.na(label_region_sample),]
    label_region_sample = label_region_sample[!is.na(label_region_sample)]
    
    epimark_region = rbind(epimark_region,epimark_region_sample)
    names(label_region_sample) = paste0(query_sample,"_",names(label_region_sample))
    label_region = c(label_region,label_region_sample)
}
## Model training
if(!is.null(DMR_epimark_file)){
    label_DMR = factor(label_DMR)
    label_region = factor(label_region)
    reptile <- reptile_train(epimark_region,label_region,
                             epimark_DMR,label_DMR,
                             family=classifier_family,
                             ntree=num_trees,nodesize=1)    
}else{
    label_region = factor(label_region)
    reptile <- reptile_train(epimark_region,label_region,
                             NULL,NULL,
                             family=classifier_family,
                             ntree=num_trees,nodesize=1)
}

## Output
save(
    reptile,
    file = paste0(output_prefix,".reptile")
    )
}
