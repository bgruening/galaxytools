rm(list=ls())
options(stringAsfactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

############################################################################################################
############################################################################################################
########################Prepare all the parameters here#####################################################
############################################################################################################
############################################################################################################

#######################options in IMA.methy450R#############################################################
######################Load data#############################################################################
libPaths = "/home/galaxy/R/x86_64-pc-linux-gnu-library/2.15"                              ### Specify the location of your R library
MethyFileName = args[1] #"/home/galaxy/install/IMA/inst/extdata/SampleMethFinalReport.txt" ### Specfiy the original methylation data produced by the GenomeStudio
PhenoFileName = args[2] #"/home/galaxy/install/IMA/inst/extdata/SamplePhenotype.txt"                 ### Specify the phenotype for each sample
############################################################################################################

###########output file######################################################################################
siteleveltest = args[3] #"./sitelevle.test.xls" ### Specify the path and name for the site-level testing result
list11excel = args[4] #"./list11excel.xls"     ### Specify the path and name for region-level analysis results
list11Rdata = args[5] #"./list11.Rdata"        ### Specify the path of Rdata file which stores the region-level analysis results
qcPlot = args[6] ### Path to the QC file
############################################################################################################

#################Preprocessing:IMA.methy450PP ##############################################################
samplefilterdetectP = as.double(args[14])   ### The cutoff for sample-level detection Pvalue
samplefilterperc = as.double(args[15])      ### The percent of loci with detection Pvalue less than "samplefilterdetectP" in each sample
sitefilterdetectP = as.double(args[16])     ### The cutoff for site-level detection Pvalue
sitefilterperc = as.double(args[17])         ### The percent of samples with detection Pvalue less than "sitefilterdetectP" for each site
na.omit = args[8]               ### Remove the sites containing missing beta value
XYchrom = args[9]                ### Remove the sites on chromosome X
peakcorrection = args[10]       ### If TRUE, peak correction is performed
cat(args[10])
normalization = args[11]        ### If TRUE, quantile normalization performed
transfm = args[12]              ### If FALSE, no transform is performed; if "arcsinsqr", arcsin square root transformation is performed; if "logit", logit transformation is performed;
locidiff = args[13]             ### If FALSE, don't filter sites by the difference of group beta value. Otherwise, remove the sites with beta value difference smaller than the specified value
locidiffgroup = c("g1","g2") ### Specify which two groups are considered to check the loci difference (if "locidiff" is not true)
snpfilter = FALSE            ### If FALSE, keep the loci whose methylation level are measured by probes containing SNP(s) at/near the targeted CpG site; otherwise, filter out the list of SNP containing loci by specifying the snp file name and location
                             ### A list of SNP-containing probes (based on dbSNP v132) could be accessed by the command: snpfilter = system.file("extdata/snpsites.txt",package ="IMA")
##############################################################################################################

############sitetest/regionwrapper############################################################################
testmethod = args[7]       ### Other options of differential testing methods: "wilcox"/"pooled"/"satterthwaite" for the comparison between two group
concov = "OFF"             ### If "ON", covariates is continuous variable
gcase = "g2"               ### Specify the case group index in the sample.txt file (if "concov" is "ON")
gcontrol = "g1"            ### Specify the control group index in the sample.txt file (if "concov" is "ON")
Padj = "BH"                ### Options for multiple testing correction. The user can choose the methods provided by p.adjust function of R stat package 
indexmethod ="mean"        ### Options for deriving an index of overall methylation value of each region. mean/median/tbrm: "tbrm" is Tukey's Biweight robust average 
paired = FALSE             ### If ture, the differential test methods would change to the corresponding paired-test methods
##############################################################################################################

####################################output the differential sites#############################################
rawpcut = NULL             ### cut off for raw pvalue 
adjustpcut = NULL          ### cut off for adjusted pvalue
betadiffcut = NULL         ### cut off for beta value difference
##############################################################################################################

##############################################################################################################
###############################End of the parameter specification#############################################
##############################################################################################################


################################ Analysis Routes #############################################################
##############################################################################################################

.libPaths(libPaths) ## Specify your R library
library(IMA)        ## load the IMA package
data =IMA.methy450R(fileName = MethyFileName,columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName,qcfile = qcPlot) ## load the data
dataf = IMA.methy450PP(data,na.omit = na.omit,normalization=normalization,peakcorrection = peakcorrection,transfm = transfm,samplefilterdetectP = samplefilterdetectP,samplefilterperc = samplefilterperc,sitefilterdetectP = sitefilterdetectP,locidiff = locidiff, locidiffgroup = locidiffgroup,XYchrom = XYchrom,snpfilter = snpfilter) ## QC filtering

sitetest = sitetest(dataf,gcase=gcase,gcontrol=gcontrol,concov=concov,testmethod = testmethod,Padj=Padj,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut,paired = paired) ## site-level testing with the "BH" adjustment
write.table(sitetest,file=siteleveltest,row.names=TRUE) ## saving the reults (note that writeXLS won't work on the data exceeds 65535 rows or 256 columns)

regionswrapper(dataf,indexmethod =indexmethod,gcase = gcase,gcontrol=gcontrol,testmethod = testmethod,Padj=Padj,concov = concov,list11excel=list11excel,list11Rdata=list11Rdata,rawpcut = rawpcut,adjustpcut = adjustpcut,betadiffcut = betadiffcut,paired = paired)   ## region-level testing for all 11 categories of annotated regions
############################### End of Analysis Routes #########################################################
################################################################################################################

