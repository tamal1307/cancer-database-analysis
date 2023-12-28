#data downloading steps
#refer to this website: https://www.costalab.org/wp-content/uploads/2022/11/handout_day5.html

library("TCGAbiolinks")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library('DESeq2')
library(ComplexHeatmap)
library(tidyverse)

#or you could use this functions to load the packages all at once

#pkgs <- c('TCGAbiolinks', 'limma', 'edgeR', 'factoextra', 'caret', 
#            'RColorBrewer', 'gProfileR', 'genefilter')
#invisible(lapply(pkgs, function(x) library(x, character.only=TRUE)))


###skip this part and directly load the pre-captured data###
  

#Since there are many lets look at the projects from tcga
  myproject <- 'TCGA-project' #use cancer abbreviation in place of 'project'
  GDCprojects = getGDCprojects()
  TCGAbiolinks:::getProjectSummary(myproject) 

#we have to keep only primary tumor samples. As sample siye in others in really low.
#editing the query to have only primary tissue.
  query_TCGA = GDCquery(
    project = myproject,
    data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
    experimental.strategy = "RNA-Seq",
    workflow.type = "STAR - Counts", 
    access = 'open', sample.type = 'Primary Tumor')
  GDCdownload(query = query_TCGA)

#To visualize the query results in a more readable way, we can use the command getResults.
  project_res = getResults(query_TCGA) # make results as table
  unique(project_res$sample_type) # first 6 types of tissue.
  summary(factor(project_res$sample_type)) # summary of distinct tissues types present in this study

tcga_data = GDCprepare(query_TCGA)
  # Save the data as a file, if you need it later, you can just load this file
  # instead of having to run the whole pipeline again
saveRDS(object = tcga_data,
          file = "tcga_data.RDS",
          compress = FALSE)
  
###start fro here if tcga_data had already been downloaded and prepared earlier###
tcga_data = readRDS(file = "tcga_data.RDS")

#writing a function to convert ENSG ids to gene symbols

convertids2<- function(mygenes, tcga_data){
  #mygenes <-  mygenes[(mygenes %in% tmp$gene_name)]
  tmp <- tcga_data@rowRanges
  return(tmp[match(mygenes,names(tmp)),]$gene_name)
}

#getting an overview of the data
dim(tcga_data)
colnames(colData(tcga_data)) #this contains clinical data
table(tcga_data@colData$vital_status)
table(tcga_data@colData$tissue_or_organ_of_origin)#

#convert the ENSG id to gene symbols
rownames(tcga_data) <- convertids2(rownames(tcga_data), tcga_data)

#quality control of clinical data
clinical.info <- colData(tcga_data)
clinical.info <- clinical.info[!is.na(clinical.info $ days_to_death), ]#removing missing death info

#to include preffered stages (here stage I+II+III only)
#hash (#) the line if want to include all
#you could include specific clinical criterias in place of 'ajcc_pathologic_stage'

clinical.info <- clinical.info[clinical.info$ajcc_pathologic_stage %in% c("Stage I",
                                                                      "Stage IA",
                                                                      "Stage II",
                                                                      'Stage IIA',
                                                                      'Stage IIB',
                                                                      'Stage IIC',
                                                                      'Stage III',
                                                                      'Stage IIIA',
                                                                      'Stage IIIB',
                                                                      'Stage IIIC'),]

dim(clinical.info)

#modifying specific columns to create binary operators with if else.
#here is an example for TCGA-COAD data to specify left and right colon samples

#clinic.coad$tumor_site <- ifelse(clinic.coad$tissue_or_organ_of_origin ==  'Descending colon', 'left_colon',
                             ifelse(clinic.coad$tissue_or_organ_of_origin ==  'Sigmoid colon', 'left_colon', 
                                    ifelse(clinic.coad$tissue_or_organ_of_origin ==  'Splenic flexure of colon', 'left_colon',
                                           ifelse(clinic.coad$tissue_or_organ_of_origin ==  'Cecum', 'right_colon',
                                                  ifelse(clinic.coad$tissue_or_organ_of_origin ==  'Ascending colon', 'right_colon',
                                                         ifelse(clinic.coad$tissue_or_organ_of_origin ==  'Hepatic flexure of colon', 'right_colon', 'NA'))))))

#this line has to be hashed when including stage IV patients
clinical.info<- clinical.info[clinical.info$ajcc_pathologic_m %in% c ('M0'),] #only non-metastasis samples included

#creating binary operators for vital status(dead=1 or alive=0)
clinical.info$vital_status <- ifelse(clinical.info$vital_status == "Dead", 1, 0)
table(clinical.info$vital_status)
dim(clinical.info)

#obtaining feature barcode matrix
matx <- assay(tcga_data)
rownames(matx) <- convertids2(rownames(matx), tcga_data)
dim(matx)

#matching pateint ids bethwwen clinical and expression data
matx <- matx[,(colnames(matx)%in%rownames(clinical.info))]
dim(matx)
match(colnames(matx),rownames(clinical.info)

#reading the clinical data into object abc for further analysis
clinical.info -> abc

#rename the 'gene symbol' with target gene to stratify patients according to expression level
abc$gene_exp <- matx['gene symbol',]
dim(clinical.info
abc$gene_exp <- ifelse(abc$gene_exp >= median(matx['gene symbol',]), 'high', 'low' )

dim(abc)
abc$survival_time <- (abc$days_to_death)/30
fit.a <- surv_fit(Surv(survival_time, vital_status)~ gene_exp, data = abc)
fit.a

#you have to change the title according to the selection criteria applied above
#in the clinical.info object

plot.a <- ggsurvplot(fit.a, break.x.by = 12, pval = TRUE, 
                     risk.table = TRUE, risk.table.title = "",
                     surv.scale = "percent", palette = c("red", "blue"), 
                     title = "
                     TCGA patients, 
                     primary site: all,
                     stage: all
                     M: all
                     ", xlab = "Time (months)")

plot.a





























