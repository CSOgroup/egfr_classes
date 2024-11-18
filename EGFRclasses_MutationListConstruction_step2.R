
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_utils.R")
CommonDataDir = "/mnt/ed2/daniele/Common_Data/"
DataDir = "/mnt/ndata/daniele/alfredo_egfr/Data/"
MainDir = "/mnt/ndata/daniele/alfredo_egfr/Processed/EGFRclasses_pipeline/Repro/"
OutDir = paste0(MainDir,"EGFR_classes/")
dir.create(OutDir)
ordered_classes = c( "common","uncommon","compound","T790M","ex20ins" )
colorz_classes = c("steelblue4","tomato3","mediumpurple4","darksalmon","orange2")

################## FUNCTIONS ##################
strReverse = function(x){
	sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}
###############################################

#################### Main #####################

### Tagging mutations with CBioPortal
load(file = paste0( OutDir,"ec_annotated.RData" ) )
cba = dan.read(paste0(OutDir,"LollipopPlot_ec_exons18_21_uncommon_table_CBioPortal_annotation.tsv"))
for (rn in rownames(cba)){
	ann = unlist(strsplit(cba[rn,"Annotation"],split=";"))
	cba[rn,"OncoKB"] = gsub("OncoKB: ","",gsub(",.*","",ann[1]) )
	cba[rn,"CIViC"] = gsub("CIViC: ","",ann[2])
	cba[rn,"CIViC"] = sum(as.numeric(gsub(".*:","",unlist(strsplit(cba[rn,"CIViC"],split=",")))))
	cba[rn,"MyCancerGenome"] = gsub("MyCancerGenome: ","",ann[3])
	cba[rn,"CancerHotspot"] = gsub("CancerHotspot: ","",ann[4])
	cba[rn,"3DHotspot"] = gsub("3DHotspot: ","",ann[5])
}
cba[is.na(cba$CIViC),"CIViC"] = 0
cba = cba[,c("Protein.Change","Mutation.Type","OncoKB","CIViC" )]
keep_oncokb = cba$OncoKB %in% c( "Likely Oncogenic","Oncogenic","Resistance" )
keep_civic = cba$CIViC>0
cba = cba[keep_oncokb | keep_civic,]
rownames(cba) = paste0("p.",cba$Protein.Change)
ec$Oncogenic = NA
ec[intersect(rownames(ec),rownames(cba)),"Oncogenic"] = "yes"
save(ec,file = paste0( OutDir,"ec_annotated_oncogenic.RData" ) )

### Intersecting with alphamissense
# am = dan.read(paste0(DataDir,"AlphaMissense_hg19.tsv"),header=F)
# am_egfr = am[grepl("ENST00000275493",am$V7),]
# save(am_egfr,file = paste0(DataDir,"AlphaMissense_hg19_EGFR.RData"))

load(file = paste0(DataDir,"AlphaMissense_hg19_EGFR.RData"))
ec$AlphaMissense = NA
am_egfr = am_egfr[am_egfr$V8 %in% gsub("p.","",ec$Mutation ), ]
am_egfr = am_egfr[!duplicated(am_egfr$V8),]
rownames(am_egfr) = paste0( "p.",am_egfr$V8 )
ec[intersect(rownames(ec),rownames(am_egfr)),"AlphaMissense"] = am_egfr[intersect(rownames(ec),rownames(am_egfr)),"V10"]
boo = dan.read(paste0(DataDir,"boostDM-EGFR.LUAD-prediction/app/tables/EGFR.LUAD.prediction.tsv"))
boo = boo[boo$boostDM_score>0.5,]
boo = boo[!duplicated(boo$aachange),]
rownames(boo) = paste0( "p.",boo$aachange )
intersect(rownames(boo),rownames(ok))
commonz = intersect(rownames(boo),ec$Mutation)
ec$boostdm_driver = NA
ec[commonz,"boostdm_driver"] = "yes"
save(ec,file = paste0( OutDir,"ec_annotated_oncogenic_alphamissense_boostdm.RData" ) )
load(file = paste0( OutDir,"ec_annotated_oncogenic_alphamissense_boostdm.RData" ) )

ec = ec[(ec$Class!="uncommon") | ((!is.na(ec$Oncogenic)) | ((ec$AlphaMissense=="pathogenic") %in% c(T) )), ]
save(ec,file = paste0( OutDir,"ec_exons18_21_oncogenic_alphamissense_filtered.RData" ) )




