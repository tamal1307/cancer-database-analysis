setwd("working/directory")

library(UCSCXenaTools)
library('DESeq2')
library(ComplexHeatmap)

# new function to look match ensembl ids with gene symbols
convertids2<- function(mygenes, tcga_data){
  #mygenes <-  mygenes[(mygenes %in% tmp$gene_name)]
  #tmp is a temporary file
  tmp <- tcga_data@rowRanges
  return(tmp[match(mygenes,names(tmp)),]$gene_name)
}

data(XenaData)

head(XenaData)
# The options in XenaFilter function support Regular Expression
#grep('icgc', XenaData$XenaCohorts, ignore.case = T )
#grep('donor', XenaData$XenaDatasets, ignore.case = T)

#XenaGenerate(subset = XenaHostNames=="icgcHub") %>% 
#  XenaFilter(filterDatasets = "PACA-AU") %>% 
#  XenaFilter(filterDatasets = "normalised") -> df_todo

#XenaGenerate(subset = XenaHostNames=="icgcHub") %>% 
#  XenaFilter(filterDatasets = "exp_seq.all_projects.donor.USonly.xena.tsv") -> df_todo2

XenaGenerate(subset = XenaHostNames=="icgcHub") %>% 
  XenaFilter(filterCohorts = "ICGC") %>%
  XenaFilter(filterDatasets = "donor") %>% 
  XenaFilter(filterDatasets = "donor/donor.all_projects.overallSurvival")-> df_OS

df_OS@datasets

#download
#XenaQuery(df_todo) %>%
#  XenaDownload() -> xe_download # this downloads normaliyed counts

#XenaQuery(df_todo2) %>%
#  XenaDownload() -> xe_download2 #this downloads data for all icgc cohorts

XenaQuery(df_OS) %>%
  XenaDownload() ->xe_OS_download #this downloads OS data for total icgcHub

#prepare data for R analysis
#cli = XenaPrepare(xe_download)
#class(cli)

#downloading OS data
cli_OS <- XenaPrepare(xe_OS_download)
class(cli_OS)
names(cli_OS)

#loading expression data with raw counts
cancer_exp <- read.table('exp-raw-counts.txt', header=T, row.names=1)
print('this file has donor id as column names')

#have to match this two. using the donor.tsv file in icgc folder
#donor.tsv has to be downloaded from ICGC website
abc= read.table('donor.tsv', sep = '\t', header=T)

#now i have donor ids for the whole icgc database and 
# submitter id & donor id match file for the cancer cohort. 
cli_OS<- cli_OS[(cli_OS$icgc_donor_id %in% abc$icgc_donor_id),]
print('PACA-AU ids sorted from icgcHub')

#cli_OS contains donor id and OS info of ICGC project patient extracted from the icgcHub
#now have to extract the icgc submitter id for the donor ids
#colnames(abc)
#cli_OS$submitted_donor_id <- abc$submitted_donor_id 
#print('addition of submitter id info to cli_OS data')

#compare between colnames(cancer_exp) and submitter_donor_id of cli_OS
cli_OS <- cli_OS[match(colnames(cancer_exp), cli_OS$icgc_donor_id),]
cli_OS<- cli_OS[!is.na(cli_OS$OS.time), ]
cli_OS<- cli_OS[!is.na(cli_OS$OS.time), ]
print('omitted patient info with NA')
cancer_exp <- cancer_exp[,match(cli_OS$icgc_donor_id, colnames(cancer_exp))]
dim(cli_OS)
dim(cancer_exp)


#add new paramiter 'dead' to the cdata, with conditions early and late 
hist(cli_OS$OS.time)
quantile(cli_OS$OS.time) #obtain the median survival time from here
cli_OS$died<- c(ifelse((cli_OS$OS.time<427), #median surviavl time is 427 days for the example cohort
                      'early', 'late')) 



#the following steps are commons for any DESeq2 and hitmap analysis
#define new dds for early vs late death analysis
dds <- DESeqDataSetFromMatrix(countData = cancer_exp,
                              colData = cli_OS,
                              design= ~ died)
dds <- DESeq(dds)

#saving the RDS file for future analysis
#saveRDS(dds, file = "my_data.rds")

res<- results(dds, contrast = c('died', 'early', 'late'))
head(res[order(res$pvalue),])

####create heatmap
##variance stabilizing transformation for heatmap
vsd<- vst(dds)

#getting genes that are upregulated
res<- as.data.frame(res)
dim(res)
mygenes<- rownames(res [!is.na(res$padj) & res$padj < 0.01 & res$log2FoldChange > 1, ])
mygenes

#genes downregulated in early dying patients
downgenes <- rownames(res [!is.na(res$padj) & res$padj < 0.01 & res$log2FoldChange < -1, ])
downgenes
#create a new matrix with gene names, load convertids2 function from sig-tcga-paad.R
mt<- assay(vsd) [mygenes, ]
#rownames(mt) <- convertids2(rownames(mt), tcga_data)
head(mt)
rownames(mt)
dim(mt)

mt_down <- assay(vsd) [downgenes, ]
#rownames(mt_down) <- convertids2(rownames(mt_down), tcga_data)

#write.table(rownames(mt), 'test-icgc.txt', row.names = F, col.names = F, quote = F)

#scale-wise normalisation (z standardisation)
#each sample will be tested for its distance from normalised mean (how many sd away?)
mt<- t(scale(t(mt)))
head(mt)
apply(mt, 1, mean)
dim(mt)
#annotation for heatmap
haC <- HeatmapAnnotation(died= colData(dds)$died,
                         days= cli_OS$OS.time)

#plot the heatmap
Heatmap(mt,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 1),
        top_annotation = haC)

