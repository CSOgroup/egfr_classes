
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_utils.R")
CommonDataDir = "/mnt/ed2/daniele/Common_Data/"
DataDir = "/mnt/ndata/daniele/alfredo_egfr/Data/"
dir.create(DataDir)
OutDir = "/mnt/ndata/daniele/alfredo_egfr/Data/using_AlphaMissense/"
dir.create(OutDir)
ordered_classes = c( "common","uncommon","compound","T790M","ex20ins" )
colorz_classes = c("steelblue4","tomato3","mediumpurple4","darksalmon","orange2")

################## FUNCTIONS ##################

###############################################

#################### Main #####################

# load(file = paste0( DataDir,"ec_exons18_21_oncogenic_alphamissense_filtered.RData" ) )
load(file = "/mnt/ndata/daniele/alfredo_egfr/Processed/EGFR_classes/ec_exons18_21_oncogenic_alphamissense_filtered.RData" ) # new version, including china

mutdf = dan.df(0,c("EGFR_status","Dataset","Count","Percentage" ))
scdf = dan.df(0,c("EGFR_class","Dataset","Count","Percentage" ))
pcdf = dan.df(0,c("EGFR_class","Dataset","Count","Percentage" ))
egfr_classes = c(unique(ec$Class),"compound")
egfr_classes = ordered_classes

## Maybe clinical tables patient-level and sample-level
######### TCGA
Clinical_extended = as.data.frame(t(read.table(paste0(CommonDataDir,"Clinical_extended/LUAD.Clinical.txt"), quote = '', sep = "\t", header = FALSE, row.names = 1)), stringsAsFactors = FALSE)
Clinical_extended[,"Patient"] <- toupper(as.character(Clinical_extended[,"bcr_patient_barcode"]))
rownames(Clinical_extended) = Clinical_extended[,"Patient"]
load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
# maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation","Nonsense_Mutation","Splice_Site"),]
maf_ss$Patient = substr(maf_ss$Tumor_Sample_Barcode,1,12)
maf_ss$Sample = substr(maf_ss$Tumor_Sample_Barcode,1,16)

nonp = length(unique(intersect( maf_ss[(maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Nonsense_Mutation","Splice_Site")),"Patient"],Clinical_extended$Patient)))
op1 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & (!(maf_ss$HGVSp_Short %in% rownames(ec))) ,"Patient"],Clinical_extended$Patient)))
op2 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & (maf_ss$HGVSp_Short %in% rownames(ec)) ,"Patient"],Clinical_extended$Patient)))
wtp = length(unique(intersect( maf_ss$Patient,Clinical_extended$Patient)))-nonp-op1-op2
tdf = data.frame(EGFR_status=c( "WT or synonymous","Nonsense or splice site","Outside exons 18-21","Inside exons 18-21" ), Dataset="TCGA",Count=c(wtp,nonp,op1,op2), Percentage=round(c(wtp,nonp,op1,op2)*100/sum(c(wtp,nonp,op1,op2)),1),stringsAsFactors=F )
mutdf = rbind(mutdf,tdf)

maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"),]
maf_ss = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
maf_ss = maf_ss[maf_ss$HGVSp_Short %in% rownames(ec),]
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
dtable(mafe$Variant_Classification)
# Subsetting Clin to EGFR mutant only
Clin = Clinical_extended[Clinical_extended$Patient %in% mafe$Patient,]
for (rn in rownames(Clin)){
	mafep = mafe[mafe$Patient==rn,]
	for (rn2 in 1:length(rownames(mafep))){
		Clin[rn,paste0("egfr_mut_",rn2)] = mafep[rn2,"HGVSp_Short"]
		Clin[rn,paste0("egfr_class_",rn2)] = ec[mafep[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin$egfr_class_1)
dtable(Clin$egfr_class_2)
dtable(Clin$egfr_class_3)
## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin$egfr_class_consensus = NA
for (rn in rownames(Clin)){
	all_muts = c( Clin[rn,"egfr_class_1"],Clin[rn,"egfr_class_2"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "TCGA", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
pcdf = rbind(pcdf, tdf)
save(Clin, file = paste0( OutDir,"Clin_TCGA_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_TCGA.RData" ))


# Sample-level
mafe2 = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Sample, maf_ss$HGVSp_Short)),]
Clin2 = mafe2[!duplicated(mafe2$Sample),c( "Patient","Sample","HGVSp_Short","Variant_Classification")]
rownames(Clin2) = Clin2$Sample
for (rn in rownames(Clin2)){
	mafep = mafe2[mafe2$Sample==rn,]
	for (rn2 in 1:length(rownames(mafep))){
		Clin2[rn,paste0("egfr_mut_",rn2)] = mafep[rn2,"HGVSp_Short"]
		Clin2[rn,paste0("egfr_class_",rn2)] = ec[mafep[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin2$egfr_class_1)
dtable(Clin2$egfr_class_2)
dtable(Clin2$egfr_class_3)
## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin2$egfr_class_consensus = NA
for (rn in rownames(Clin2)){
	all_muts = c( Clin2[rn,"egfr_class_1"],Clin2[rn,"egfr_class_2"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin2[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin2[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin2[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin2[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin2[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin2$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "TCGA", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin2$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin2$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
scdf = rbind(scdf, tdf)
## adding clinical variables
Clinical_extended$Patient = NULL
for (rn in rownames(Clin2)){
	Clin2[rn,colnames(Clinical_extended)] = Clinical_extended[Clin2[rn,"Patient"],colnames(Clinical_extended)]
}
save(Clin2, file = paste0( OutDir,"Clin2_TCGA_SampleLevel_complete.RData" ))
Clin2 = Clin2[Clin2$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin2, file = paste0( OutDir,"Clin2_TCGA_SampleLevel.RData" ))
## Relationship between sample-level and patient-level
load(file = paste0( OutDir,"Clin_TCGA.RData" ))
load(file = paste0( OutDir,"Clin2_TCGA_SampleLevel.RData" ))
### NO DIFFERENCE, great










######### ChenEAS
chen_indir = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoSsignature_onChen/"
load(file = paste0(chen_indir,"ClinChen.RData"))
maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
maf$Patient = maf$Tumor_Sample_Barcode

maf_ss = maf
nonp = length(unique(intersect( maf_ss[(maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("nonsense_mutation","Splice_Site")),"Patient"],ClinChen$Patient)))
op1 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins", "missense"))) & (!(maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],ClinChen$Patient)))
op2 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins", "missense"))) & ((maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],ClinChen$Patient)))
wtp = length(unique(intersect( maf_ss$Patient,ClinChen$Patient)))-nonp-op1-op2
tdf = data.frame(EGFR_status=c( "WT or synonymous","Nonsense or splice site","Outside exons 18-21","Inside exons 18-21" ), Dataset="ChenEAS",Count=c(wtp,nonp,op1,op2), Percentage=round(c(wtp,nonp,op1,op2)*100/sum(c(wtp,nonp,op1,op2)),1),stringsAsFactors=F )
mutdf = rbind(mutdf,tdf)

table(maf$Variant_Classification)
maf_ss = maf[(maf$Variant_Classification %in% c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins", "missense")),]
maf_ss = maf_ss[maf_ss$HGVSp_Short %in% rownames(ec),]
mafe = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
mafe = mafe[!duplicated(paste0(mafe$Hugo_Symbol, mafe$Patient, mafe$HGVSp_Short)),]
Clin = ClinChen
rownames(Clin) = Clin$Patient.ID
Clin$Patient = Clin$Patient.ID
Clin = Clin[Clin$Patient %in% mafe$Patient,]
for (rn in rownames(Clin)){
	mafep = mafe[mafe$Patient==rn,]
	for (rn2 in 1:length(rownames(mafep))){
		Clin[rn,paste0("egfr_mut_",rn2)] = mafep[rn2,"HGVSp_Short"]
		Clin[rn,paste0("egfr_class_",rn2)] = ec[mafep[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin$egfr_class_1)
dtable(Clin$egfr_class_2)
dtable(Clin$egfr_class_3)
dtable(Clin$egfr_class_4)
## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin$egfr_class_consensus = NA
for (rn in rownames(Clin)){
	all_muts = c( Clin[rn,"egfr_class_1"],Clin[rn,"egfr_class_2"],Clin[rn,"egfr_class_3"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "ChenEAS", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
pcdf = rbind(pcdf, tdf)
scdf = rbind(scdf, tdf)
dtable(Clin$egfr_class_1,Clin$egfr_class_2,Clin$egfr_class_consensus)
dtable(Clin$egfr_class_1,Clin$egfr_class_2,Clin$egfr_class_3,Clin$egfr_class_consensus)
Clin$Sample = Clin$Patient
save(Clin, file = paste0( OutDir,"Clin_Chen_complete.RData" ))
Clin2 = Clin
save(Clin2, file = paste0( OutDir,"Clin2_Chen_SampleLevel_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_Chen.RData" ))
Clin2 = Clin
save(Clin2, file = paste0( OutDir,"Clin2_Chen_SampleLevel.RData" ))
# no difference between samples and patients





#### Genie, new release (14) downloaded by me.
Clin2 = dan.read(paste0(DataDir,"genie_14.0/data_clinical_sample.txt"))
Clin2$Patient = Clin2$PATIENT_ID
Clin2$Sample = Clin2$SAMPLE_ID
Clin2 = Clin2[Clin2$CANCER_TYPE %in% c( "Non-Small Cell Lung Cancer" ),]
Clin2 = Clin2[Clin2$CANCER_TYPE_DETAILED=="Lung Adenocarcinoma",]
rownames(Clin2) = Clin2$Sample

Clin = dan.read(paste0(DataDir,"genie_14.0/data_clinical_patient.txt"))
Clin$Patient = Clin$PATIENT_ID
rownames(Clin) = Clin$Patient
Clin = Clin[Clin$Patient %in% Clin2$Patient,]

load(file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]

nonp = length(unique(intersect( maf_ss[(maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Nonsense_Mutation","Splice_Site")),"Patient"],Clin$Patient)))
op1 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & (!(maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],Clin$Patient)))
op2 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & ((maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],Clin$Patient)))
wtp = length(unique(intersect( maf_ss$Patient,Clin$Patient)))-nonp-op1-op2
tdf = data.frame(EGFR_status=c( "WT or synonymous","Nonsense or splice site","Outside exons 18-21","Inside exons 18-21" ), Dataset="Genie",Count=c(wtp,nonp,op1,op2), Percentage=round(c(wtp,nonp,op1,op2)*100/sum(c(wtp,nonp,op1,op2)),1),stringsAsFactors=F )
mutdf = rbind(mutdf,tdf)

maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")) & (maf_ss$Hugo_Symbol=="EGFR"),]
maf_ss = maf_ss[maf_ss$HGVSp_Short %in% rownames(ec),]

# patient-level
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
# Subsetting Clin to EGFR mutant only
Clin = Clin[Clin$Patient %in% mafe$Patient,]
for (rn in rownames(Clin)){
	mafep = mafe[mafe$Patient==rn,]
	for (rn2 in 1:length(rownames(mafep))){
		Clin[rn,paste0("egfr_mut_",rn2)] = mafep[rn2,"HGVSp_Short"]
		Clin[rn,paste0("egfr_class_",rn2)] = ec[mafep[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin$egfr_class_1)
dtable(Clin$egfr_class_2)
dtable(Clin$egfr_class_3)
dtable(Clin$egfr_class_4)
dtable(Clin$egfr_class_5)
dtable(Clin$egfr_class_6)

## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin$egfr_class_consensus = NA
for (rn in rownames(Clin)){
	all_muts = c( Clin[rn,"egfr_class_1"],Clin[rn,"egfr_class_2"],Clin[rn,"egfr_class_3"],Clin[rn,"egfr_class_4"],Clin[rn,"egfr_class_5"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "Genie", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
pcdf = rbind(pcdf, tdf)
save(Clin, file = paste0( OutDir,"Clin_Genie_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_Genie.RData" ))
# Sample-level
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Sample, maf_ss$HGVSp_Short)),]
Clin2 = Clin2[Clin2$Sample %in% mafe$Sample,]
for (rn in rownames(Clin2)){
	mafe2p = mafe[mafe$Sample==rn,]
	for (rn2 in 1:length(rownames(mafe2p))){
		Clin2[rn,paste0("egfr_mut_",rn2)] = mafe2p[rn2,"HGVSp_Short"]
		Clin2[rn,paste0("egfr_class_",rn2)] = ec[mafe2p[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin2$egfr_class_1)
dtable(Clin2$egfr_class_2)
dtable(Clin2$egfr_class_3)
dtable(Clin2$egfr_class_4)
dtable(Clin2$egfr_class_5)
dtable(Clin2$egfr_class_6)

## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin2$egfr_class_consensus = NA
for (rn in rownames(Clin2)){
	all_muts = c( Clin2[rn,"egfr_class_1"],Clin2[rn,"egfr_class_2"],Clin2[rn,"egfr_class_3"],Clin2[rn,"egfr_class_4"],Clin2[rn,"egfr_class_5"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin2[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin2[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin2[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin2[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin2[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin2$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "Genie", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin2$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin2$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
scdf = rbind(scdf, tdf)
save(Clin2, file = paste0( OutDir,"Clin2_Genie_SampleLevel_complete.RData" ))
Clin2 = Clin2[Clin2$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin2, file = paste0( OutDir,"Clin2_Genie_SampleLevel.RData" ))
## Relationship between sample-level and patient-level
load(file = paste0( OutDir,"Clin_Genie.RData" ))
load(file = paste0( OutDir,"Clin2_Genie_SampleLevel.RData" ))
load(file = paste0( "/mnt/ndata/daniele/alfredo_egfr/Data/using_AlphaMissense/","Clin_Genie_complete.RData" ))
load(file = paste0( "/mnt/ndata/daniele/alfredo_egfr/Data/using_AlphaMissense/","Clin2_Genie_SampleLevel_complete.RData" ))
all(Clin$Patient %in% Clin2$Patient)
all(Clin2$Patient %in% Clin$Patient)
for (rn in rownames(Clin2)){
	Clin2[rn,"PatientLevel_egfr_class_consensus"] = Clin[Clin2[rn,"Patient"],"egfr_class_consensus"]
}
dtable(Clin2$egfr_class_consensus, Clin2$PatientLevel_egfr_class_consensus)




##### Adding clinical information to Genie from various sources.
# Order them from earlier to latest paper! 

load(file = paste0( OutDir,"Clin_Genie_complete.RData" ))
Clin$Sex = tolower(Clin$SEX)
Clin[Clin$Sex=="unknown","Sex"] = NA
Clin$Race = tolower(Clin$PRIMARY_RACE)
Clin[ tolower(Clin$PRIMARY_RACE) %in% c("unknown","not collected","other"), "Race"] = NA
Clin$Ethnicity = tolower(Clin$ETHNICITY)
Clin[ Clin$Ethnicity=="non-spanish/non-hispanic", "Ethnicity"] = "non-hispanic"
Clin[ Clin$Ethnicity=="spanish/hispanic", "Ethnicity"] = "hispanic"
Clin[ Clin$Ethnicity %in% c("unknown","not collected"), "Ethnicity"] = NA

################### MSK 2017 ###################
msk = dan.read( paste0(DataDir,"lung_msk_2017_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
	if (!is.na(mean(tmsk$Tumor.Purity, na.rm=T))) { Clin[rn,"Purity"] = mean(tmsk$Tumor.Purity, na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Age = msk$Diagnosis.Age
msk$Smoking = ifelse(msk$Smoking.History=="Never","never smoker","ever smoker")
msk$Stage = NA
msk[msk$Stage.At.Diagnosis %in% c( "IA","IB" ),"Stage"] = "I"
msk[msk$Stage.At.Diagnosis %in% c( "IIA","IIB" ),"Stage"] = "II"
msk[msk$Stage.At.Diagnosis %in% c( "IIIA","IIIB","III" ),"Stage"] = "III"
msk[msk$Stage.At.Diagnosis %in% c( "IV" ),"Stage"] = "IV"
varz = c( "Age", "Smoking", "Stage" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### MSK-IMPACT June 2017 ###################
msk = dan.read( paste0(DataDir,"msk_impact_2017_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
	if (!is.na(mean(tmsk$Tumor.Purity, na.rm=T))) { Clin[rn,"Purity"] = mean(tmsk$Tumor.Purity, na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Smoking = ifelse( msk$Smoking.History=="Unknown",NA,ifelse( msk$Smoking.History=="Never","never smoker","ever smoker" ) )
msk$OS_months = msk$Overall.Survival..Months.
msk$OS_status = as.numeric(substr(msk$Overall.Survival.Status,1,1))
varz = c( "Smoking", "OS_months", "OS_status" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### MSK NSCLC PD1 2018 ###################
msk = dan.read( paste0(DataDir,"nsclc_pd1_msk_2018_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Age = msk$Diagnosis.Age
msk$Smoking = ifelse(msk$Smoker=="Never","never smoker","ever smoker")
varz = c( "Age", "Smoking" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### JTO, December 2020 ###################
jto = dan.read(paste0(DataDir,"luad_mskcc_2020_JTO_patterns_clinical_data.tsv"))
rownames(jto) = paste0("GENIE-MSK-", jto$Patient.ID ) # JTO patients are unique (1 sample per patient), great
commonz = (intersect(rownames(jto),rownames(Clin)))
length(commonz)
jto = jto[rownames(jto) %in% commonz,]
jto$Age = jto$Age.At.Surgery
jto$FGA = jto$Fraction.Genome.Altered
jto$nMut = jto$Mutation.Count
jto$Stage = jto$Pathologic.Stage
jto$Stage = ifelse((jto$Stage==1) %in% c(T),"I",ifelse((jto$Stage==2) %in% c(T),"II",ifelse((jto$Stage==3) %in% c(T),"III",ifelse((jto$Stage==4) %in% c(T),"IV",NA) )))
jto$TMB = tolower(jto$Tumor.Mutation.Burden)
jto$TMB_nonsynonymous = tolower(jto$TMB..nonsynonymous.)
jto$OS_months = tolower(jto$Overall.Survival..Months.)
jto$OS_status = tolower(jto$Overall.Survival.Status)
jto$Pattern = tolower(jto$Predominant.Histologic.Subtype)
jto$Smoking = tolower(jto$Smoking.History)
msk$Smoking = ifelse(msk$Smoker=="Never","never smoker","ever smoker")
varz = c( "Age","FGA","Smoking","nMut","Stage","Pattern","TMB","TMB_nonsynonymous","OS_months","OS_status" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(jto[rn,var])) { Clin[rn,var] = jto[rn,var] }
	}
}

################### MSK NPJPO Jul 2021 ###################
msk = dan.read( paste0(DataDir,"luad_msk_npjpo_2021_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Age = msk$Age.at.Resection
msk$Smoking = tolower(msk$Smoking.Status)
msk$Stage = ifelse(msk$Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage=="Stage I","I","II")
varz = c( "Age", "Smoking", "Stage" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### MSK-IMPACT Feb 2022 ###################
msk = dan.read( paste0(DataDir,"luad_mskimpact_2021_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
	if (!is.na(mean(tmsk$Tumor.Purity, na.rm=T))) { Clin[rn,"Purity"] = mean(tmsk$Tumor.Purity, na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Ethnicity = NA
msk[!is.na(msk$Ethnicity.Category),"Ethnicity"] = ifelse(msk[!is.na(msk$Ethnicity.Category),"Ethnicity.Category"]=="Non-Spanish; Non-Hispanic","non-hispanic","hispanic")
msk$Race = NA
msk[ (msk$Race.Category=="ASIAN-FAR EAST/INDIAN SUBCONT") %in% c(T),"Race" ] = "asian"
msk[ (msk$Race.Category=="BLACK OR AFRICAN AMERICAN") %in% c(T),"Race" ] = "black"
msk[ (msk$Race.Category=="WHITE") %in% c(T),"Race" ] = "white"
msk$Age = msk$Age.at.Which.Sequencing.was.Reported..Years.
msk$OS_months = msk$Overall.Survival..Months.
msk$OS_status = as.numeric(substr(msk$Overall.Survival.Status,1,1))
varz = c( "Age", "Race","Ethnicity","OS_months", "OS_status" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### MSK MET Feb 2022 ###################
msk = dan.read( paste0(DataDir,"msk_met_2021_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
	if (!is.na(mean(tmsk$Tumor.Purity, na.rm=T))) { Clin[rn,"Purity"] = mean(tmsk$Tumor.Purity, na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Race = NA
msk[ (msk$Race.Category=="Asian-far east/indian subcont") %in% c(T),"Race" ] = "asian"
msk[ (msk$Race.Category=="Black or african american") %in% c(T),"Race" ] = "black"
msk[ (msk$Race.Category=="White") %in% c(T),"Race" ] = "white"
msk[ (msk$Race.Category=="Native american-am ind/alaska") %in% c(T),"Race" ] = "native american"
msk[ (msk$Race.Category=="Native hawaiian or pacific isl") %in% c(T),"Race" ] = "pacific islander"
msk$Age = msk$Age.at.Sequencing
msk$OS_months = msk$Overall.Survival..Months.
msk$OS_status = as.numeric(substr(msk$Overall.Survival.Status,1,1))
varz = c( "Age", "Race","OS_months", "OS_status" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### MSK Mind Oct 2022 ###################
msk = dan.read( paste0(DataDir,"lung_msk_mind_2020_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
	if (!is.na(mean(tmsk$Tumor.Purity, na.rm=T))) { Clin[rn,"Purity"] = mean(tmsk$Tumor.Purity, na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Age = msk$Age.at.Which.Sequencing.was.Reported..Years.
msk$Smoking = ifelse( msk$What.is.the.patient.s.smoking.status.=="Never smoker","never smoker","ever smoker")
varz = c( "Age","Smoking" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### MSK NSCLC CTDX Nov 2022 ###################
msk = dan.read( paste0(DataDir,"nsclc_ctdx_msk_2022_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
	if (!is.na(mean(tmsk$Tumor.Purity, na.rm=T))) { Clin[rn,"Purity"] = mean(tmsk$Tumor.Purity, na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Age = msk$Patient.Current.Age
msk$Smoking = NA
msk[(msk$Smoking.Status=="False") %in% c(T),"Smoking"] = "never smoker"
msk[(msk$Smoking.Status=="True") %in% c(T),"Smoking"] = "ever smoker"
msk$Stage = "IV"
msk$OS_months = msk$Overall.Survival..Months.
msk$OS_status = as.numeric(substr(msk$Overall.Survival.Status,1,1))
varz = c( "Age", "Smoking", "Stage", "OS_months", "OS_status" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### MSK NSCLC BM Aug 2023 ###################
msk = dan.read( paste0(DataDir,"bm_nsclc_mskcc_2023_clinical_data.tsv") )
msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
commonz = intersect(msk$Patient,Clin$Patient)
length(commonz)
msk = msk[msk$Patient %in% commonz,]
## For sample-level measurements (FGA, TMB, ...) we take the mean
for (rn in commonz){
	tmsk = msk[msk$Patient==rn,]
	if (!is.na(mean(tmsk$Fraction.Genome.Altered, na.rm=T))) { Clin[rn,"FGA"] = mean(tmsk$Fraction.Genome.Altered, na.rm=T) }
	if (!is.na(mean(tmsk$Mutation.Count, na.rm=T))) { Clin[rn,"nMut"] = mean(tmsk$Mutation.Count, na.rm=T) }
	if (!is.na(mean(tmsk$TMB..nonsynonymous., na.rm=T))) { Clin[rn,"TMB_nonsynonymous"] = mean(tmsk$TMB..nonsynonymous., na.rm=T) }
}
msk = msk[!duplicated(msk$Patient),]
rownames(msk) = msk$Patient
msk$Smoking = NA
msk[(msk$Smoking.Status=="Never") %in% c(T),"Smoking"] = "never smoker"
msk[(msk$Smoking.Status %in% c("Current","Former") %in% c(T)),"Smoking"] = "ever smoker"
msk$Stage = "IV"
msk$OS_months = msk$Overall.Survival..Months.
msk$OS_status = as.numeric(substr(msk$Overall.Survival.Status,1,1))
varz = c( "Smoking", "Stage", "OS_months", "OS_status" )
for (rn in commonz){
	for (var in varz){
		if (!is.na(msk[rn,var])) { Clin[rn,var] = msk[rn,var] }
	}
}

################### BPC, 1 Sept 2023 ###################
bpc = read.csv(paste0(DataDir,"genie_BPC_NSCLC/cancer_level_dataset_index.csv"),stringsAsFactors=F)
bpc = bpc[!duplicated(bpc$record_id),]
rownames(bpc) = bpc$record_id
length(intersect(rownames(bpc),rownames(Clin)))
dtable(bpc$ca_lung_cigarette)
bpc$Smoking = ifelse( bpc$ca_lung_cigarette=="Unknown",NA,ifelse( bpc$ca_lung_cigarette=="Never used","never smoker","ever smoker" ) )
bpc$Age = bpc$age_dx
bpc$Stage = NA
bpc[bpc$stage_dx=="Stage I","Stage"] = "I"
bpc[bpc$stage_dx=="Stage II","Stage"] = "II"
bpc[bpc$stage_dx=="Stage III","Stage"] = "III"
bpc[bpc$stage_dx=="Stage IV","Stage"] = "IV"
bpc$OS_days = bpc$tt_os_dx_days
bpc$OS_months = bpc$tt_os_dx_mos
bpc$OS_status = bpc$os_dx_status
varz = c( "Age","Smoking","Stage","OS_days", "OS_months", "OS_status" )
for (rn in intersect(rownames(bpc),rownames(Clin))){
	for (var in varz){
		if (!is.na(bpc[rn,var])) { Clin[rn,var] = bpc[rn,var] }
	}
}

save(Clin,file = paste0( OutDir,"Clin_Genie_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_Genie.RData" ))












##### Adding pattern to TCGA and Chen
load(file = paste0( OutDir,"Clin_TCGA_complete.RData" ))
Clint = Clin
load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/TCGA_histology_formatted_permissive.RData"))
rownames(Clin) = Clin$Patient
Clint$Pattern = NA
Clint[intersect(rownames(Clint),rownames(Clin)),"Pattern"] = Clin[intersect(rownames(Clint),rownames(Clin)),"Pattern"]
Clin = Clint
save(Clin,file = paste0( OutDir,"Clin_TCGA_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_TCGA.RData" ))

load(file = paste0( OutDir,"Clin2_TCGA_SampleLevel_complete.RData" ))
Clint = Clin2
load(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/TCGA_histology_formatted_permissive.RData"))
rownames(Clin) = Clin$Patient
Clint$Pattern = NA
rownames(Clint) = Clint$Patient
Clint[intersect(rownames(Clint),rownames(Clin)),"Pattern"] = Clin[intersect(rownames(Clint),rownames(Clin)),"Pattern"]
rownames(Clint) = Clint$Sample
Clin2 = Clint
save(Clin2,file = paste0( OutDir,"Clin2_TCGA_SampleLevel_complete.RData" ))
Clin2 = Clin2[Clin2$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin2, file = paste0( OutDir,"Clin2_TCGA_SampleLevel.RData" ))







##### Zhang (lung cancer in never-smokers)
Clin = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_clinical_patient.txt"))
Clin$Patient = Clin$PATIENT_ID
Clin = Clin[Clin$HISTOLOGY=="Adenocarcinomas",]

Clin2 = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_clinical_sample.txt"))
Clin2$Patient = Clin2$PATIENT_ID
Clin2$Sample = Clin2$SAMPLE_ID
Clin2 = Clin2[Clin2$CANCER_TYPE %in% c( "Non-Small Cell Lung Cancer" ),]
Clin2 = Clin2[Clin2$CANCER_TYPE_DETAILED=="Lung Adenocarcinoma",]
rownames(Clin2) = Clin2$Sample
Clin2 = Clin2[Clin2$Patient %in% Clin$Patient,]

Clin = Clin[Clin$Patient %in% Clin2$Patient,]
rownames(Clin) = Clin$Patient

maf = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_mutations.txt"))
maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]

nonp = length(unique(intersect( maf_ss[(maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Nonsense_Mutation","Splice_Site")),"Patient"],Clin$Patient)))
op1 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & (!(maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],Clin$Patient)))
op2 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & ((maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],Clin$Patient)))
wtp = length(unique(intersect( maf_ss$Patient,Clin$Patient)))-nonp-op1-op2
tdf = data.frame(EGFR_status=c( "WT or synonymous","Nonsense or splice site","Outside exons 18-21","Inside exons 18-21" ), Dataset="Zhang",Count=c(wtp,nonp,op1,op2), Percentage=round(c(wtp,nonp,op1,op2)*100/sum(c(wtp,nonp,op1,op2)),1),stringsAsFactors=F )
mutdf = rbind(mutdf,tdf)

maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")) & (maf_ss$Hugo_Symbol=="EGFR"),]
intersect(maf_ss$HGVSp_Short, rownames(ec))

maf_ss = maf_ss[maf_ss$HGVSp_Short %in% rownames(ec),]

# patient-level
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
# Subsetting Clin to EGFR mutant only
Clin = Clin[Clin$Patient %in% mafe$Patient,]
for (rn in rownames(Clin)){
	mafep = mafe[mafe$Patient==rn,]
	for (rn2 in 1:length(rownames(mafep))){
		Clin[rn,paste0("egfr_mut_",rn2)] = mafep[rn2,"HGVSp_Short"]
		Clin[rn,paste0("egfr_class_",rn2)] = ec[mafep[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin$egfr_class_1)
dtable(Clin$egfr_class_2)
dtable(Clin$egfr_class_3)


## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin$egfr_class_consensus = NA
for (rn in rownames(Clin)){
	all_muts = c( Clin[rn,"egfr_class_1"],Clin[rn,"egfr_class_2"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "Zhang", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
pcdf = rbind(pcdf, tdf)
save(Clin, file = paste0( OutDir,"Clin_Zhang_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_Zhang.RData" ))
# Sample-level
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Sample, maf_ss$HGVSp_Short)),]
Clin2 = Clin2[Clin2$Sample %in% mafe$Sample,]
for (rn in rownames(Clin2)){
	mafe2p = mafe[mafe$Sample==rn,]
	for (rn2 in 1:length(rownames(mafe2p))){
		Clin2[rn,paste0("egfr_mut_",rn2)] = mafe2p[rn2,"HGVSp_Short"]
		Clin2[rn,paste0("egfr_class_",rn2)] = ec[mafe2p[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin2$egfr_class_1)
dtable(Clin2$egfr_class_2)
dtable(Clin2$egfr_class_3)

## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin2$egfr_class_consensus = NA
for (rn in rownames(Clin2)){
	all_muts = c( Clin2[rn,"egfr_class_1"],Clin2[rn,"egfr_class_2"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin2[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin2[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin2[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin2[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin2[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin2$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "Zhang", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin2$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin2$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
scdf = rbind(scdf, tdf)
save(Clin2, file = paste0( OutDir,"Clin2_Zhang_SampleLevel_complete.RData" ))
Clin2 = Clin2[Clin2$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin2, file = paste0( OutDir,"Clin2_Zhang_SampleLevel.RData" ))
## Relationship between sample-level and patient-level
load(file = paste0( OutDir,"Clin_Zhang.RData" ))
load(file = paste0( OutDir,"Clin2_Zhang_SampleLevel.RData" ))
all(Clin$Patient %in% Clin2$Patient)
for (rn in rownames(Clin2)){
	Clin2[rn,"PatientLevel_egfr_class_consensus"] = Clin[Clin2[rn,"Patient"],"egfr_class_consensus"]
}
dtable(Clin2$egfr_class_consensus, Clin2$PatientLevel_egfr_class_consensus)
## Exactly the same

load(file = paste0( OutDir,"Clin_Zhang_complete.RData" ))
load(file = paste0( OutDir,"Clin2_Zhang_SampleLevel_complete.RData" ))
# c( "Stage","Sex","Race","Ethnicity","Smoking","Pattern","Age" )
Clin$Stage = Clin$TUMOR_STAGE
Clin$Ethnicity = NA
Clin$Race = NA
Clin$Smoking = "never smoker"
Clin$Pattern = NA
Clin$Stage[Clin$Stage %in% c( "IA","IB" )] = "I"
Clin$Sex = tolower(Clin$SEX)
Clin$Age = as.numeric(Clin$AGE)
rownames(Clin2) = Clin2$Patient
Clin$TMB_nonsynonymous = Clin2[rownames(Clin),"TMB_NONSYNONYMOUS"]
Clin$Purity = Clin2[rownames(Clin),"TUMOR_PURITY"]
save(Clin,file = paste0( OutDir,"Clin_Zhang_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_Zhang.RData" ))






##### TRACERx 421
Clin = readRDS(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_all_patient_df.rds"))
Clin$Patient = Clin$cruk_id
rownames(Clin) = Clin$cruk_id
Clin = Clin[grepl("LUAD",Clin$histology_multi_full_genomically.confirmed),]
Clin2 = readRDS(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_all_tumour_df.rds"))
Clin2 = as.data.frame(Clin2,stringsAsFactors=F)
Clin2$Patient = Clin2$cruk_id
Clin2$Sample = Clin2$tumour_id_muttable_cruk
rownames(Clin2) = Clin2$Sample
Clin2 = Clin2[Clin2$Patient %in% Clin$Patient,]
Clin = Clin[Clin$Patient %in% Clin2$Patient,]
rownames(Clin) = Clin$Patient
library(fst)
maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
mafr = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221123_TRACERx421_mutation_table_region.fst"))
# mafr = mafr[mafr$tumour_id %in% Clin2$Sample,]
# expr = read_fst(paste0(DataDir,"TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_counts_mat.fst"))
maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
maf_ss$Sample = maf_ss$tumour_id
maf_ss$Patient = maf_ss$patient_id
mafe = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
mafe = mafe[((mafe$exonic.func %in% c("frameshift substitution","frameshift insertion", "nonframeshift insertion", "nonframeshift substitution", "nonsynonymous","nonsynonymous SNV")) %in% c(T)) & (mafe$Hugo_Symbol=="EGFR"),]

nonp = length(unique(intersect( maf_ss[(maf_ss$Hugo_Symbol=="EGFR") & ((substr(maf_ss$AAChange,nchar(maf_ss$AAChange),nchar(maf_ss$AAChange))=="*") %in% c(T)  ),"Patient"],Clin$Patient)))
op1 = length(unique(intersect( mafe[!((mafe$NucleotideChange %in% rownames(ec)) | ( (mafe$AAChange %in% rownames(ec)) %in% c(T) )),"Patient"],Clin$Patient)))
op2 = length(unique(intersect( mafe[(mafe$NucleotideChange %in% rownames(ec)) | ( (mafe$AAChange %in% rownames(ec)) %in% c(T) ),"Patient"],Clin$Patient)))

op2 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") ) & ((maf_ss$NucleotideChange %in% rownames(ec)) | ( (maf_ss$AAChange %in% rownames(ec)) %in% c(T) )),"Patient"],Clin$Patient)))

wtp = length(unique(intersect( maf_ss$Patient,Clin$Patient)))-nonp-op1-op2
tdf = data.frame(EGFR_status=c( "WT or synonymous","Nonsense or splice site","Outside exons 18-21","Inside exons 18-21" ), Dataset="TRACERx421",Count=c(wtp,nonp,op1,op2), Percentage=round(c(wtp,nonp,op1,op2)*100/sum(c(wtp,nonp,op1,op2)),1),stringsAsFactors=F )
mutdf = rbind(mutdf,tdf)

maf_ss = mafe
maf_ss = maf_ss[(maf_ss$NucleotideChange %in% rownames(ec)) | ( (maf_ss$AAChange %in% rownames(ec)) %in% c(T) ),]

# patient-level
maf_ss$mut = maf_ss$AAChange
maf_ss[is.na(maf_ss$mut),"mut"] = maf_ss[is.na(maf_ss$mut),"NucleotideChange"]
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$mut)),]

# Subsetting Clin to EGFR mutant only
Clin = Clin[Clin$Patient %in% mafe$Patient,]
for (rn in rownames(Clin)){
	mafep = mafe[mafe$Patient==rn,]
	for (rn2 in 1:length(rownames(mafep))){
		Clin[rn,paste0("egfr_mut_",rn2)] = ifelse(is.na(mafep[rn2,"AAChange"]),mafep[rn2,"NucleotideChange"],mafep[rn2,"AAChange"])
		Clin[rn,paste0("egfr_class_",rn2)] = ifelse(mafep[rn2,"AAChange"] %in% rownames(ec),ec[mafep[rn2,"AAChange"],"Class"],ec[mafep[rn2,"NucleotideChange"],"Class"])
	}
}
dtable(Clin$egfr_class_1)
dtable(Clin$egfr_class_2)
dtable(Clin$egfr_class_3)
## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin$egfr_class_consensus = NA
for (rn in rownames(Clin)){
	all_muts = c( Clin[rn,"egfr_class_1"],Clin[rn,"egfr_class_2"])
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "TRACERx421", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
pcdf = rbind(pcdf, tdf)
save(Clin, file = paste0( OutDir,"Clin_TRACERx421_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_TRACERx421.RData" ))
# Sample-level
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Sample, maf_ss$mut)),]
Clin2 = Clin2[Clin2$Sample %in% mafe$Sample,]
Clin2save=Clin2
for (rn in rownames(Clin2)){
	mafe2p = mafe[mafe$Sample==rn,]
	for (rn2 in 1:length(rownames(mafe2p))){
		Clin2[rn,paste0("egfr_mut_",rn2)] = ifelse(is.na(mafe2p[rn2,"AAChange"]),mafe2p[rn2,"NucleotideChange"],mafe2p[rn2,"AAChange"])
		Clin2[rn,paste0("egfr_class_",rn2)] = ifelse(mafe2p[rn2,"AAChange"] %in% rownames(ec),ec[mafe2p[rn2,"AAChange"],"Class"],ec[mafe2p[rn2,"NucleotideChange"],"Class"])
	}
}
dtable(Clin2$egfr_class_1)
dtable(Clin2$egfr_class_2)
dtable(Clin2$egfr_class_3)

## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin2$egfr_class_consensus = NA
for (rn in rownames(Clin2)){
	all_muts = c( Clin2[rn,"egfr_class_1"],Clin2[rn,"egfr_class_2"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin2[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin2[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin2[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin2[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin2[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin2$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "TRACERx421", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin2$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin2$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
scdf = rbind(scdf, tdf)
save(Clin2, file = paste0( OutDir,"Clin2_TRACERx421_SampleLevel_complete.RData" ))
Clin2 = Clin2[Clin2$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin2, file = paste0( OutDir,"Clin2_TRACERx421_SampleLevel.RData" ))
## Relationship between sample-level and patient-level
load(file = paste0( OutDir,"Clin_TRACERx421.RData" ))
load(file = paste0( OutDir,"Clin2_TRACERx421_SampleLevel.RData" ))
all(Clin$Patient %in% Clin2$Patient)
for (rn in rownames(Clin2)){
	Clin2[rn,"PatientLevel_egfr_class_consensus"] = Clin[Clin2[rn,"Patient"],"egfr_class_consensus"]
}
dtable(Clin2$egfr_class_consensus, Clin2$PatientLevel_egfr_class_consensus)

load(file = paste0( OutDir,"Clin_TRACERx421_complete.RData" ))
load(file = paste0( OutDir,"Clin2_TRACERx421_SampleLevel_complete.RData" ))
c( "Stage","Sex","Race","Ethnicity","Smoking","Pattern","Age" )
Clin$Race = "white"
Clin[Clin$ethnicity %in% c( "Indian","Middle eastern" ),"Race"] = "asian"
Clin[Clin$ethnicity %in% c( "Black" ),"Race"] = "black"
Clin$Ethnicity = "non-hispanic"
Clin[Clin$ethnicity=="South American","Ethnicity"] = "hispanic"
Clin$Stage = "I"
Clin$Stage[Clin$pathologyTNM %in% c("IIA","IIB")] = "II"
Clin$Stage[Clin$pathologyTNM=="IIIA"] = "III"
Clin$Sex = tolower(Clin$sex)
Clin$smoking_status_merged = as.character(Clin$smoking_status_merged)
Clin$Smoking = "never smoker"
Clin[ Clin$smoking_status_merged %in% c( "Ex-Smoker","Smoker" ),"Smoking"] = "ever smoker"
Clin$Pattern = Clin$LUAD_pred_subtype
Clin$Age = as.numeric(Clin$age)
Clin$TMB_nonsynonymous = NA
save(Clin, file=paste0( OutDir,"Clin_TRACERx421_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_TRACERx421.RData" ))








##### China Origimed
Clin = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_clinical_patient.txt"))
Clin2 = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_clinical_sample.txt"))
Clin2 = Clin2[Clin2$CANCER_TYPE_DETAILED %in% c( "Lung Adenocarcinoma" ),]
Clin = Clin[Clin$PATIENT_ID %in% Clin2$PATIENT_ID,]
rownames(Clin) = Clin$PATIENT_ID
rownames(Clin2) = Clin2$PATIENT_ID
Clin$Patient = Clin$PATIENT_ID
Clin2$Sample = Clin2$SAMPLE_ID
Clin = cbind(Clin,Clin2[rownames(Clin),c( "Sample","AJCC_PATHOLOGIC_TUMOR_STAGE","SAMPLE_TYPE","TUMOR_PURTITY","GDNA.NG.","TMB_NONSYNONYMOUS" )])
Clin2 = Clin
rownames(Clin2) = Clin2$Sample
maf = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_mutations_extended.txt"))
maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]

nonp = length(unique(intersect( maf_ss[(maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Nonsense_Mutation","Splice_Site")),"Patient"],Clin$Patient)))
op1 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & (!(maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],Clin$Patient)))
op2 = length(unique(intersect( maf_ss[((maf_ss$Hugo_Symbol=="EGFR") & (maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"))) & ((maf_ss$HGVSp_Short %in% rownames(ec))),"Patient"],Clin$Patient)))
wtp = length(unique(intersect( maf_ss$Patient,Clin$Patient)))-nonp-op1-op2
tdf = data.frame(EGFR_status=c( "WT or synonymous","Nonsense or splice site","Outside exons 18-21","Inside exons 18-21" ), Dataset="Origimed",Count=c(wtp,nonp,op1,op2), Percentage=round(c(wtp,nonp,op1,op2)*100/sum(c(wtp,nonp,op1,op2)),1),stringsAsFactors=F )
mutdf = rbind(mutdf,tdf)

maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")) & (maf_ss$Hugo_Symbol=="EGFR"),]
intersect(maf_ss$HGVSp_Short, rownames(ec))

maf_ss = maf_ss[maf_ss$HGVSp_Short %in% rownames(ec),]

# patient-level
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
# Subsetting Clin to EGFR mutant only
Clin = Clin[Clin$Patient %in% mafe$Patient,]
for (rn in rownames(Clin)){
	mafep = mafe[mafe$Patient==rn,]
	for (rn2 in 1:length(rownames(mafep))){
		Clin[rn,paste0("egfr_mut_",rn2)] = mafep[rn2,"HGVSp_Short"]
		Clin[rn,paste0("egfr_class_",rn2)] = ec[mafep[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin$egfr_class_1)
dtable(Clin$egfr_class_2)
dtable(Clin$egfr_class_3)
dtable(Clin$egfr_class_4)

## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin$egfr_class_consensus = NA
for (rn in rownames(Clin)){
	all_muts = c( Clin[rn,"egfr_class_1"],Clin[rn,"egfr_class_2"],Clin[rn,"egfr_class_3"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "Origimed", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
pcdf = rbind(pcdf, tdf)
save(Clin, file = paste0( OutDir,"Clin_Origimed_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_Origimed.RData" ))
# Sample-level
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Sample, maf_ss$HGVSp_Short)),]
Clin2 = Clin2[Clin2$Sample %in% mafe$Sample,]
for (rn in rownames(Clin2)){
	mafe2p = mafe[mafe$Sample==rn,]
	for (rn2 in 1:length(rownames(mafe2p))){
		Clin2[rn,paste0("egfr_mut_",rn2)] = mafe2p[rn2,"HGVSp_Short"]
		Clin2[rn,paste0("egfr_class_",rn2)] = ec[mafe2p[rn2,"HGVSp_Short"],"Class"]
	}
}
dtable(Clin2$egfr_class_1)
dtable(Clin2$egfr_class_2)
dtable(Clin2$egfr_class_3)
dtable(Clin2$egfr_class_4)

## Assigning compounds. They are compounds if: common+uncommon, uncommon+uncommon
Clin2$egfr_class_consensus = NA
for (rn in rownames(Clin2)){
	all_muts = c( Clin2[rn,"egfr_class_1"],Clin2[rn,"egfr_class_2"],Clin2[rn,"egfr_class_3"] )
	all_muts = all_muts[!is.na(all_muts)]
	all_muts_unf = all_muts
	all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
	if (length(all_muts)>0){
		if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
		Clin2[rn,"egfr_class_consensus"] = "common" 
		} else {
			Clin2[rn,"egfr_class_consensus"] = "compound" 
		}
		if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ Clin2[rn,"egfr_class_consensus"] = "uncommon" }
	}
	if (any(tolower(all_muts_unf)=="t790m")) { Clin2[rn,"egfr_class_consensus"] = "T790M" }
	if (any(all_muts_unf=="ex20ins")) { Clin2[rn,"egfr_class_consensus"] = "ex20ins" }
}
dtable(Clin2$egfr_class_consensus)
tdf = data.frame(row.names = egfr_classes, EGFR_class = egfr_classes, Dataset = "Origimed", Count=0, stringsAsFactors=F)
tdf[names(dtable(Clin2$egfr_class_consensus)),"Count"] = as.numeric(dtable(Clin2$egfr_class_consensus))
tdf$Percentage = round(tdf$Count*100/sum(tdf$Count),1)
scdf = rbind(scdf, tdf)
save(Clin2, file = paste0( OutDir,"Clin2_Origimed_SampleLevel_complete.RData" ))
Clin2 = Clin2[Clin2$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin2, file = paste0( OutDir,"Clin2_Origimed_SampleLevel.RData" ))
## Relationship between sample-level and patient-level
load(file = paste0( OutDir,"Clin_Origimed.RData" ))
load(file = paste0( OutDir,"Clin2_Origimed_SampleLevel.RData" ))
all(Clin$Patient %in% Clin2$Patient)
for (rn in rownames(Clin2)){
	Clin2[rn,"PatientLevel_egfr_class_consensus"] = Clin[Clin2[rn,"Patient"],"egfr_class_consensus"]
}
dtable(Clin2$egfr_class_consensus, Clin2$PatientLevel_egfr_class_consensus)
## Exactly the same

load(file = paste0( OutDir,"Clin_Origimed_complete.RData" ))
load(file = paste0( OutDir,"Clin2_Origimed_SampleLevel_complete.RData" ))
# c( "Stage","Sex","Race","Ethnicity","Smoking","Pattern","Age" )
Clin = Clin[Clin$AJCC_PATHOLOGIC_TUMOR_STAGE!=0,]
Clin$Stage = Clin$AJCC_PATHOLOGIC_TUMOR_STAGE
Clin$Ethnicity = NA
Clin$Race = "asian"
dtable(Clin$SMOKE.STATUS)
Clin$Smoking = NA
Clin[Clin$SMOKE.STATUS=="Smoker","Smoking"] = "ever smoker"
Clin[Clin$SMOKE.STATUS=="Nonsmoker","Smoking"] = "never smoker"
Clin$Pattern = NA
Clin$Sex = tolower(Clin$SEX)
Clin$Age = as.numeric(Clin$DIAGNOSIS.AGE)
Clin$TMB_nonsynonymous = Clin[rownames(Clin),"TMB_NONSYNONYMOUS"]
Clin$Purity = Clin[rownames(Clin),"TUMOR_PURITY"]
save(Clin,file = paste0( OutDir,"Clin_Origimed_complete.RData" ))
Clin = Clin[Clin$egfr_class_consensus %in% c( "common","compound","uncommon" ),]
save(Clin, file = paste0( OutDir,"Clin_Origimed.RData" ))
load(file = paste0( OutDir,"Clin_Origimed_complete.RData" ))
Clin2 = Clin
save(Clin2,file = paste0( OutDir,"Clin2_Origimed_SampleLevel_complete.RData" ))
load(file = paste0( OutDir,"Clin_Origimed.RData" ))
Clin2 = Clin
save(Clin2,file = paste0( OutDir,"Clin2_Origimed_SampleLevel.RData" ))









save(scdf, file = paste0( OutDir,"EGFRclasses_SampleLevel_barplot.RData" ))
save(pcdf, file = paste0( OutDir,"EGFRclasses_PatientLevel_barplot.RData" ))
save(mutdf, file = paste0( OutDir,"EGFRstatus_barplot.RData" ))

load( file = paste0( OutDir,"EGFRclasses_SampleLevel_barplot.RData" ))
load( file = paste0( OutDir,"EGFRclasses_PatientLevel_barplot.RData" ))
load( file = paste0( OutDir,"EGFRstatus_barplot.RData" ))

mutdf$EGFR_status = factor(mutdf$EGFR_status, levels = c( "WT or synonymous","Nonsense or splice site","Outside exons 18-21","Inside exons 18-21" ))
mutdf$TotalPatients = NA
mutdf[(mutdf$EGFR_status=="WT or synonymous") & (mutdf$Dataset=="TCGA"),"TotalPatients"] = paste0("N = ",sum(mutdf[(mutdf$Dataset=="TCGA"),"Count"]))
mutdf[(mutdf$EGFR_status=="WT or synonymous") & (mutdf$Dataset=="Genie"),"TotalPatients"] = paste0("N = ",sum(mutdf[(mutdf$Dataset=="Genie"),"Count"]))
mutdf[(mutdf$EGFR_status=="WT or synonymous") & (mutdf$Dataset=="ChenEAS"),"TotalPatients"] = paste0("N = ",sum(mutdf[(mutdf$Dataset=="ChenEAS"),"Count"]))
mutdf[(mutdf$EGFR_status=="WT or synonymous") & (mutdf$Dataset=="Zhang"),"TotalPatients"] = paste0("N = ",sum(mutdf[(mutdf$Dataset=="Zhang"),"Count"]))
mutdf[(mutdf$EGFR_status=="WT or synonymous") & (mutdf$Dataset=="TRACERx421"),"TotalPatients"] = paste0("N = ",sum(mutdf[(mutdf$Dataset=="TRACERx421"),"Count"]))
mutdf[(mutdf$EGFR_status=="WT or synonymous") & (mutdf$Dataset=="Origimed"),"TotalPatients"] = paste0("N = ",sum(mutdf[(mutdf$Dataset=="Origimed"),"Count"]))
colorz = c("gray55","gold","darkorange1","firebrick")
pdf(paste0(OutDir,"EGFRstatus_barplot.pdf"),4.5,2.5)
p = ggplot(data=mutdf, aes(x=Dataset, y=Percentage, fill=EGFR_status)) +
  geom_bar(stat="identity", colour="black", size = 0.1) + ylab("Percentage of patients") + xlab("Dataset") + scale_fill_manual(name = "EGFR mutation",values=colorz) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# scale_x_discrete(labels=rownames(ldf)) +
  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = 6/.pt)
  p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
print(p)
dev.off()

pcdf$EGFR_class = factor(pcdf$EGFR_class, levels = ordered_classes)
colorz = colorz_classes
pcdf$TotalPatients = NA
for (dataset in unique(pcdf$Dataset)){ pcdf[(pcdf$EGFR_class=="common") & (pcdf$Dataset==dataset),"TotalPatients"] = paste0("N = ",sum(pcdf[(pcdf$Dataset==dataset),"Count"])) }
pdf(paste0(OutDir,"EGFRclasses_PatientLevel_barplot.pdf"),7.5,4.2)
p = ggplot(data=pcdf, aes(x=Dataset, y=Percentage, fill=EGFR_class)) +
  geom_bar(stat="identity", colour="black", size = 0.1) + ylab("Percentage of patients") + xlab("Dataset") + scale_fill_manual(name = "EGFR mutation class",values=colorz) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# scale_x_discrete(labels=rownames(ldf)) +
  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)
print(p)
dev.off()

## Piechart
load( file = paste0( OutDir,"EGFRclasses_PatientLevel_barplot.RData" ))
pie = aggregate(Count~EGFR_class,pcdf,FUN='sum')
library(ggplot2)
library(dplyr)
pie$percent = 100*pie$Count/sum(pie$Count)
pie$percent_rounded = round(pie$percent,1)
pie$label = paste0( pie$Count,"\n(",pie$percent_rounded,"%)" )
# Basic piechart
pie$EGFR_class = factor(pie$EGFR_class, levels=ordered_classes)
plot = ggplot(pie, aes(x="", y=Count, fill=EGFR_class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  geom_text(aes(x=1.25,label = label),size = 6/.pt,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual( values=colorz_classes )
pdf( paste0(OutDir,"Pie_PooledDatasets.pdf"),3,3 )
plot = plot + theme(text = element_text(size=6), legend.key.size = unit(10, 'pt'))
print(plot)
dev.off()

#### Now the same but only on the classes of interest (common, uncommon, compound)
ordered_classes = c( "common","uncommon","compound" )
colorz_classes = c("steelblue4","tomato3","mediumpurple4")
load( file = paste0( OutDir,"EGFRclasses_PatientLevel_barplot.RData" ))
pcdf = pcdf[pcdf$EGFR_class %in% ordered_classes,]
for (d in unique(pcdf$Dataset)){
	total = sum(pcdf[pcdf$Dataset==d,"Count"])
	for (cl in ordered_classes){
		pcdf[(pcdf$Dataset==d) & (pcdf$EGFR_class==cl),"Percentage"] = pcdf[(pcdf$Dataset==d) & (pcdf$EGFR_class==cl),"Count"]*100/total
	}
}
pcdf$Percentage = round(pcdf$Percentage,1)

pcdf$EGFR_class = factor(pcdf$EGFR_class, levels = ordered_classes)
colorz = colorz_classes
pcdf$TotalPatients = NA
for (dataset in unique(pcdf$Dataset)){ pcdf[(pcdf$EGFR_class=="common") & (pcdf$Dataset==dataset),"TotalPatients"] = paste0("N = ",sum(pcdf[(pcdf$Dataset==dataset),"Count"])) }
pdf(paste0(OutDir,"EGFRclasses_PatientLevel_barplot_onlyCUC.pdf"),4,2.5)
p = ggplot(data=pcdf, aes(x=Dataset, y=Percentage, fill=EGFR_class)) +
  geom_bar(stat="identity", colour="black", size = 0.1) + ylab("Percentage of patients") + xlab("Dataset") + scale_fill_manual(name = "EGFR mutation class",values=colorz) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# scale_x_discrete(labels=rownames(ldf)) +
  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = 6/.pt)
p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
print(p)
dev.off()

## Piechart
pie = aggregate(Count~EGFR_class,pcdf,FUN='sum')
library(ggplot2)
library(dplyr)
pie$percent = 100*pie$Count/sum(pie$Count)
pie$percent_rounded = round(pie$percent,1)
pie$label = paste0( pie$Count,"\n(",pie$percent_rounded,"%)" )
# Basic piechart
pie$EGFR_class = factor(pie$EGFR_class, levels=ordered_classes)
plot = ggplot(pie, aes(x="", y=Count, fill=EGFR_class)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  geom_text(aes(x=1.25,label = label),size = 6/.pt,
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual( values=colorz_classes )
pdf( paste0(OutDir,"Pie_PooledDatasets.pdf"),3,3)
plot = plot + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), legend.key.size = unit(10, 'pt'))
print(plot)
dev.off()










