
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

ec = dan.df(0,c("Mutation", "Class", "First_Codon", "Last_Codon", "Domain", "Variant_Classification") )

### Genie
load( file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")) & (maf_ss$Hugo_Symbol=="EGFR"),]
maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
mafe = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Sample, maf_ss$HGVSp_Short)),]
ecnew = data.frame(Mutation = mafe$HGVSp_Short, Class = NA, Variant_Classification = mafe$Variant_Classification, stringsAsFactors=F)
ecnew = ecnew[!duplicated(ecnew$Mutation),]
stringz = gsub("\\_.*","",ecnew$Mutation )
stringz = gsub("\\*.*","",stringz)
stringz = as.numeric(gsub("\\D", "", stringz))
ecnew$First_Codon = stringz
aa = sapply(1:length(ecnew$Mutation),function(x) gsub(stringz[x],"",ecnew$Mutation[x]))
stringz = as.numeric(gsub("\\D", "", aa))
ecnew$Last_Codon = stringz
ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"Last_Codon" ] = ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"First_Codon" ]
ecnew = ecnew[order(ecnew$First_Codon),]
ecnew$Domain = "others"
ecnew[(ecnew$First_Codon>=688) & (ecnew$First_Codon<=875),"Domain" ] = "tirosine_kinase_18_21"
rownames(ecnew) = ecnew$Mutation
ec = rbind(ec,ecnew[,colnames(ec)])
ec = ec[order(ec$First_Codon),]

### TCGA
load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")) & (maf_ss$Hugo_Symbol=="EGFR"),]
mafe = maf_ss[(!duplicated(maf_ss$HGVSp_Short)) & (!(maf_ss$HGVSp_Short %in% ec$Mutation)),]
ecnew = data.frame(Mutation = mafe$HGVSp_Short, Class = NA, Variant_Classification = mafe$Variant_Classification, stringsAsFactors=F)
ecnew = ecnew[!duplicated(ecnew$Mutation),]
stringz = gsub("\\_.*","",ecnew$Mutation )
stringz = gsub("\\*.*","",stringz)
stringz = as.numeric(gsub("\\D", "", stringz))
ecnew$First_Codon = stringz
aa = sapply(1:length(ecnew$Mutation),function(x) gsub(stringz[x],"",ecnew$Mutation[x]))
stringz = as.numeric(gsub("\\D", "", aa))
ecnew$Last_Codon = stringz
ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"Last_Codon" ] = ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"First_Codon" ]
ecnew = ecnew[order(ecnew$First_Codon),]
ecnew$Domain = "others"
ecnew[(ecnew$First_Codon>=688) & (ecnew$First_Codon<=875),"Domain" ] = "tirosine_kinase_18_21"
rownames(ecnew) = ecnew$Mutation
ec = rbind(ec,ecnew[,colnames(ec)])
ec = ec[order(ec$First_Codon),]

### ChenEAS
maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
maf_ss = maf[(maf$Variant_Classification %in% c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins", "missense")) & (maf$Hugo_Symbol=="EGFR"),]
mafe = maf_ss[(!duplicated(maf_ss$HGVSp_Short)) & (!(maf_ss$HGVSp_Short %in% ec$Mutation)),]
ecnew = data.frame(Mutation = mafe$HGVSp_Short, Class = NA, Variant_Classification = mafe$Variant_Classification, stringsAsFactors=F)
ecnew = ecnew[!duplicated(ecnew$Mutation),]
stringz = gsub("\\_.*","",ecnew$Mutation )
stringz = gsub("\\*.*","",stringz)
stringz = as.numeric(gsub("\\D", "", stringz))
ecnew$First_Codon = stringz
aa = sapply(1:length(ecnew$Mutation),function(x) gsub(stringz[x],"",ecnew$Mutation[x]))
stringz = as.numeric(gsub("\\D", "", aa))
ecnew$Last_Codon = stringz
ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"Last_Codon" ] = ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"First_Codon" ]
ecnew = ecnew[order(ecnew$First_Codon),]
ecnew$Domain = "others"
ecnew[(ecnew$First_Codon>=688) & (ecnew$First_Codon<=875),"Domain" ] = "tirosine_kinase_18_21"
rownames(ecnew) = ecnew$Mutation
ecnew[ecnew$Variant_Classification=="frame_shift_del","Variant_Classification"] = "Frame_Shift_Del"
ecnew[ecnew$Variant_Classification=="frame_shift_ins","Variant_Classification"] = "Frame_Shift_Ins"
ecnew[ecnew$Variant_Classification=="in_frame_del","Variant_Classification"] = "In_Frame_Del"
ecnew[ecnew$Variant_Classification=="in_frame_ins","Variant_Classification"] = "In_Frame_Ins"
ecnew[ecnew$Variant_Classification=="missense","Variant_Classification"] = "Missense_Mutation"
ec = rbind(ec,ecnew[,colnames(ec)])
ec = ec[order(ec$First_Codon),]

### Zhang
Clin = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_clinical_patient.txt"))
Clin = Clin[Clin$HISTOLOGY=="Adenocarcinomas",]
maf = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_mutations.txt"))
Clin2 = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_clinical_sample.txt"))
Clin2 = Clin2[Clin2$PATIENT_ID %in% Clin$PATIENT_ID,]
maf_ss = maf[((maf$Tumor_Sample_Barcode %in% Clin2$SAMPLE_ID) & ( maf$Hugo_Symbol=="EGFR" )) & ( maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation") ),]
mafe = maf_ss[(!duplicated(maf_ss$HGVSp_Short)) & (!(maf_ss$HGVSp_Short %in% ec$Mutation)),]
ecnew = data.frame(Mutation = mafe$HGVSp_Short, Class = NA, Variant_Classification = mafe$Variant_Classification, stringsAsFactors=F)
ecnew = ecnew[!duplicated(ecnew$Mutation),]
stringz = gsub("\\_.*","",ecnew$Mutation )
stringz = gsub("\\*.*","",stringz)
stringz = as.numeric(gsub("\\D", "", stringz))
ecnew$First_Codon = stringz
aa = sapply(1:length(ecnew$Mutation),function(x) gsub(stringz[x],"",ecnew$Mutation[x]))
stringz = as.numeric(gsub("\\D", "", aa))
ecnew$Last_Codon = stringz
ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"Last_Codon" ] = ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"First_Codon" ]
ecnew = ecnew[order(ecnew$First_Codon),]
ecnew$Domain = "others"
ecnew[(ecnew$First_Codon>=688) & (ecnew$First_Codon<=875),"Domain" ] = "tirosine_kinase_18_21"
rownames(ecnew) = ecnew$Mutation
ec = rbind(ec,ecnew[,colnames(ec)])
ec = ec[order(ec$First_Codon),]

### TRACERx421
Clin = readRDS(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_all_patient_df.rds"))
Clin$Patient = Clin$cruk_id
rownames(Clin) = Clin$cruk_id
Clin = Clin[grepl("LUAD",Clin$histology_multi_full_genomically.confirmed),]
Clin2 = readRDS(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_all_tumour_df.rds"))
Clin2$Patient = Clin2$cruk_id
Clin2$Sample = Clin2$tumour_id_muttable_cruk
rownames(Clin2) = Clin2$Sample
Clin2 = Clin2[Clin2$Patient %in% Clin$Patient,]
Clin = Clin[Clin$Patient %in% Clin2$Patient,]
rownames(Clin) = Clin$Patient
library(fst)
maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
maf_ss$Sample = maf_ss$tumour_id
maf_ss$Patient = maf_ss$patient_id
mafe = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
maf_ss = mafe[((mafe$exonic.func %in% c("frameshift substitution","frameshift insertion", "nonframeshift insertion", "nonframeshift substitution", "nonsynonymous","nonsynonymous SNV")) %in% c(T)) & (mafe$Hugo_Symbol=="EGFR"),]
maf_ss$nucleotide_start = as.numeric(gsub(".*?([0-9]+).*", "\\1", sub("_.*", "", maf_ss$NucleotideChange)))
maf_ss$AAChange_numeric = as.numeric(gsub(".*?([0-9]+).*", "\\1", maf_ss$AAChange))
# removing variants already present in ec
maf_ss = maf_ss[is.na(maf_ss$AAChange) | (!(maf_ss$AAChange %in% ec$Mutation)),]
# keeping only those between 688 and 875. But some indels don't have AAchange, so I'll check them one by one when constructing ecnew
maf_ss = maf_ss[is.na(maf_ss$AAChange_numeric) | ( (maf_ss$AAChange_numeric>=688) & (maf_ss$AAChange_numeric<=875) ),]
maf_ss = maf_ss[!duplicated(paste0(maf_ss$start,"_",maf_ss$stop)),]
ecnew = data.frame(row.names=c( "c.2299_2300insCCAGCGTGT","c.2300_2301insCAGCGTGGA","c.2296_2297insTGGCCAGCG","c.2236_2236delins-AATTAAGAGAAGCAACAT","c.2234_2234delins-GGAATTAAGAGAAGC","c.2237_2237delins-ATTAAGAGAAGC","c.2239_2239delins-TAAGAGAAGCAACATCTC","c.2235_2235delins-GAATTAAGAGAAGCA" ), 
	Mutation=c( "c.2299_2300insCCAGCGTGT","c.2300_2301insCAGCGTGGA","c.2296_2297insTGGCCAGCG","c.2236_2236delins-AATTAAGAGAAGCAACAT","c.2234_2234delins-GGAATTAAGAGAAGC","c.2237_2237delins-ATTAAGAGAAGC","c.2239_2239delins-TAAGAGAAGCAACATCTC","c.2235_2235delins-GAATTAAGAGAAGCA" ), 
	Class=c("ex20ins","ex20ins","ex20ins","common","common","common","common","common"), 
	First_Codon=c(767,767,766,746,745,746,747,745), Last_Codon=c(767,767,766,746,745,746,747,745), Domain="tirosine_kinase_18_21", 
	Variant_Classification=c("In_Frame_Ins","In_Frame_Ins","In_Frame_Ins","In_Frame_Del","In_Frame_Del","In_Frame_Del","In_Frame_Del","In_Frame_Del"),stringsAsFactors=F)
ec = rbind(ec,ecnew)

###### Adding variants from China Pan-can study
Clin = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_clinical_patient.txt"))
Clin2 = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_clinical_sample.txt"))
Clin2 = Clin2[Clin2$CANCER_TYPE_DETAILED %in% c( "Lung Adenocarcinoma" ),]
Clin = Clin[Clin$PATIENT_ID %in% Clin2$PATIENT_ID,]
rownames(Clin) = Clin$PATIENT_ID
rownames(Clin2) = Clin2$PATIENT_ID
Clin = cbind(Clin,Clin2[rownames(Clin),c( "SAMPLE_ID","AJCC_PATHOLOGIC_TUMOR_STAGE","SAMPLE_TYPE","TUMOR_PURTITY","GDNA.NG.","TMB_NONSYNONYMOUS" )])
Clin2 = Clin
maf = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_mutations_extended.txt"))
maf_ss = maf[((maf$Tumor_Sample_Barcode %in% Clin2$SAMPLE_ID) & ( maf$Hugo_Symbol=="EGFR" )) & ( maf$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation") ),]
mafe = maf_ss[(!duplicated(maf_ss$HGVSp_Short)) & (!(maf_ss$HGVSp_Short %in% ec$Mutation)),]
ecnew = data.frame(Mutation = mafe$HGVSp_Short, Class = NA, Variant_Classification = mafe$Variant_Classification, stringsAsFactors=F)
ecnew = ecnew[!duplicated(ecnew$Mutation),]
stringz = gsub("\\_.*","",ecnew$Mutation )
stringz = gsub("\\*.*","",stringz)
stringz = as.numeric(gsub("\\D", "", stringz))
ecnew$First_Codon = stringz
aa = sapply(1:length(ecnew$Mutation),function(x) gsub(stringz[x],"",ecnew$Mutation[x]))
stringz = as.numeric(gsub("\\D", "", aa))
ecnew$Last_Codon = stringz
ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"Last_Codon" ] = ecnew[is.na(ecnew$Last_Codon) | ( (ecnew$Last_Codon<=ecnew$First_Codon) %in% c(T) ),"First_Codon" ]
ecnew = ecnew[order(ecnew$First_Codon),]
ecnew$Domain = "others"
ecnew[(ecnew$First_Codon>=688) & (ecnew$First_Codon<=875),"Domain" ] = "tirosine_kinase_18_21"
rownames(ecnew) = ecnew$Mutation
ec = rbind(ec,ecnew[,colnames(ec)])
ec = ec[order(ec$First_Codon),]

save(ec, file = paste0(OutDir,"ec_all_unannotated.RData"))
load( file = paste0(OutDir,"ec_all_unannotated.RData"))
ec = ec[ec$Domain=="tirosine_kinase_18_21",]
### annotation rules:  
# common: p.L858R or In_Frame_Del overlapping 746-755
# ex20ins: ins/dups overlapping 762-774 "EGFR exon 20 insertion mutations are heterogeneous at the molecular level but can be characterized as inframe insertions or duplications of between 3 and 21bp (corresponding to 1 to 7 amino acids) clustered between amino acid positions 762 and 774 of the EGFR protein"
# T790M: p.T790M
# uncommon: all the others
ec[(ec$Variant_Classification=="In_Frame_Del") & ( (ec$First_Codon<=755) | (ec$Last_Codon>=746) ),"Class"] = "common"
ec["p.L858R","Class"] = "common"
ec[(ec$Variant_Classification=="In_Frame_Ins") & ( (ec$First_Codon<=774) | (ec$Last_Codon>=762) ),"Class"] = "ex20ins"
ec["p.T790M","Class"] = "T790M"
ec[is.na(ec$Class),"Class"] = "uncommon" # to be checked with OncoKB/BoostDM/Alphamissense
save(ec, file = paste0(OutDir,"ec_annotated.RData"))

####### Pruning OncoKB'ed or BoostDM'ed mutations
### Reading OncoKB annotations
ok = dan.read( paste0(DataDir,"oncokb_egfr_annotation_20Oct2023.txt"),header=F )
ok = data.frame(row.names = paste0("p.",ok[seq(1,nrow(ok),3),]), alteration = paste0("p.",ok[seq(1,nrow(ok),3),]), oncogenic = ok[1+seq(1,nrow(ok),3),], mutation_effect = ok[2+seq(1,nrow(ok),3),], stringsAsFactors=F)
ok = ok[ok$oncogenic %in% c( "Likely Oncogenic","Oncogenic","Resistance" ),]
commonz = intersect(ec$Mutation, rownames(ok))
ec$oncokb_oncogenic_or_resistance = "no"
ec[commonz,"oncokb_oncogenic_or_resistance"] = "yes"
# double-checking (missense are fine, they should have been seen)
ec[ (ec$oncokb_oncogenic_or_resistance=="no") & (!(ec$Variant_Classification=="Missense_Mutation")),"Mutation"]
ec[ (ec$oncokb_oncogenic_or_resistance=="no") & (!(ec$Variant_Classification=="Missense_Mutation")),"Class"]
rownames(ok[(!(rownames(ok) %in% commonz)) & (nchar(ok$alteration)>8),])
# Adding the following:
# "p.Exon 19 in-frame deletions (729_761del)"
# "p.Exon 19 in-frame insertions (729_761ins)"
ec[ ((ec$oncokb_oncogenic_or_resistance=="no") & ((ec$Variant_Classification %in% c("In_Frame_Del","In_Frame_Ins") ))) & ((ec$First_Codon>=729) & (ec$First_Codon<=761)  ) ,"oncokb_oncogenic_or_resistance"] = "yes"
# Adding the following:
# "p.Exon 20 in-frame insertions (762_823ins)"
ec[ ((ec$oncokb_oncogenic_or_resistance=="no") & ((ec$Variant_Classification %in% c("In_Frame_Ins") ))) & ((ec$First_Codon>=762) & (ec$First_Codon<=823)  ) ,"oncokb_oncogenic_or_resistance"] = "yes"
# ec[ (ec$oncokb_oncogenic_or_resistance=="no") & ((ec$Variant_Classification %in% c("In_Frame_Ins") )),"First_Codon"]
# ec[ (ec$oncokb_oncogenic_or_resistance=="no") & ((ec$Variant_Classification %in% c("In_Frame_Ins") )),"Mutation"]

### Reading BoostDM data
# only missense variants, easy
boo = dan.read(paste0(DataDir,"boostDM-EGFR.LUAD-prediction/app/tables/EGFR.LUAD.prediction.tsv"))
boo = boo[boo$boostDM_score>0.5,]
boo = boo[!duplicated(boo$aachange),]
rownames(boo) = paste0( "p.",boo$aachange )
intersect(rownames(boo),rownames(ok))
commonz = intersect(rownames(boo),ec$Mutation)
ec$boostdm_driver = "no"
ec[commonz,"boostdm_driver"] = "yes"
dtable(ec$boostdm_driver,ec$oncokb_oncogenic_or_resistance,ec$Class)
ec[((ec$boostdm_driver=="yes") & (ec$oncokb_oncogenic_or_resistance=="yes")) & (ec$Class=="common"),] # p.L858R, as expected
mean(boo[ec[(ec$boostdm_driver=="yes") & (ec$oncokb_oncogenic_or_resistance=="yes"),"Mutation"],"boostDM_score"])
mean(boo[ec[(ec$boostdm_driver=="yes") & (ec$oncokb_oncogenic_or_resistance=="no"),"Mutation"],"boostDM_score"])
## BoostDM scores tend to be higher for OncoKB'ed variants, as expected

######## Let's check clustering of uncommon (before and after BoostDm's)
load(file = paste0( OutDir,"ec_annotated.RData" ) )
ec = ec[substr(ec$Mutation,1,1)!="c",]
# format for https://www.cbioportal.org/mutation_mapper
cb = data.frame( Hugo_Symbol="EGFR",Sample_ID=ec$Class,Protein_Change=gsub("p.","",ec$Mutation) )
prefix = "all"
dan.write(cb,file=paste0( OutDir,"LollipopPlot_ec_exons18_21_",prefix,"_table.txt" ))
this_cb = cb[cb$Sample_ID=="common",]
prefix = "common"
dan.write(this_cb,file=paste0( OutDir,"LollipopPlot_ec_exons18_21_",prefix,"_table.txt" ))
this_cb = cb[cb$Sample_ID=="ex20ins",] # my calls are all correct
prefix = "ex20ins"
dan.write(this_cb,file=paste0( OutDir,"LollipopPlot_ec_exons18_21_",prefix,"_table.txt" ))
this_cb = cb[cb$Sample_ID=="uncommon",]
prefix = "uncommon"
dan.write(this_cb,file=paste0( OutDir,"LollipopPlot_ec_exons18_21_",prefix,"_table.txt" ))
load(file = paste0( OutDir,"ec_annotated.RData" ) )
ec = ec[substr(ec$Mutation,1,1)!="c",]
cb = data.frame( Hugo_Symbol="EGFR",Sample_ID=ec$Class,Protein_Change=gsub("p.","",ec$Mutation) )
this_cb = cb[cb$Sample_ID=="uncommon",]
prefix = "uncommon"
dan.write(this_cb,file=paste0( OutDir,"LollipopPlot_ec_",prefix,"_table.txt" ))

tcb = cb[cb$Class=="uncommon",c("Hugo_Symbol","Sample_ID","Protein_Change")]
prefix = "uncommon"
dan.write(tcb,file=paste0( OutDir,"LollipopPlot_Run04_",prefix,"_table.txt" ))
head(tcb)

### Tagging mutations with CBioPortal
