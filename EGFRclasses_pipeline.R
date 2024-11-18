
### sourcing helper and analysis functions
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_utils.R")
source( "/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_functions.R" )

### Define pipeline run name/ID and MainDir
MainDir = "/mnt/ndata/daniele/alfredo_egfr/Processed/EGFRclasses_pipeline/Repro/"
dir.create(MainDir)

### Constructing EGFR mutation list
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_MutationListConstruction_step1.R")
# Tagging uncommon variants with cBioPortal: go to https://www.cbioportal.org/mutation_mapper, upload the file in MainDir 'LollipopPlot_ec_exons18_21_uncommon_table.txt' and download to MainDir the annotation of mutations, calling it 'LollipopPlot_ec_exons18_21_uncommon_table_CBioPortal_annotation.tsv'
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_MutationListConstruction_step2.R")

### Preprocessing clinical data
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_Preprocessing.R")

### Step 1 - assigning patients
OutDir = paste0(MainDir,"Preprocessing/")
dir.create(OutDir)
patients_assignment( OutDir, retain = c( "common","uncommon","compound"), formatted_clinical_folder = paste0(MainDir,"Raw_clinical_formatting/") )
compounding_uncommons( OutDir )
vaf_class_analysis( OutDir )

### Step 2 - clinical associations
OutDir = paste0(MainDir,"Clinical/")
dir.create(OutDir)
clinical_associations( OutDir )

### Step 3 - mutational signatures analysis
OutDir = paste0(MainDir,"MutSignatures/") # analysis on targeted sequencing supported by: https://pubmed.ncbi.nlm.nih.gov/33571361/
dir.create(OutDir)
mutationalSignatures_preprocessing( OutDir )
mutationalSignatures_associations_binary( OutDir, threshold = 0.3 )
mutationalSignatures_associations_binary( OutDir, threshold = 0 )
mutationalSignatures_associations_binary( OutDir, threshold = 0.5 )
mutationalSignatures_egfr_muts_contexts( OutDir ) # only common vs uncommon

### Step 4 - genomic associations
OutDir = paste0(MainDir,"Genomic_TmbFga/")
dir.create(OutDir)
TmbFga_associations( OutDir )
OutDir = paste0(MainDir,"Genomic_coMutations/")
dir.create(OutDir)
mutations_associations_alldatasets( OutDir )
kras_further_analyses( OutDir )
# mutations_associations_accounting_for_TMB( OutDir, nbins = 10 )
mutations_associations_accounting_for_TMB( OutDir, nbins = 8 )
# mutations_associations_accounting_for_TMB( OutDir, nbins = 6 )
OutDir = paste0(MainDir,"Genomic_coCNAs/")
dir.create(OutDir)
cna_associations_alldatasets( OutDir, ncopies = 2 )
cna_associations_alldatasets( OutDir, ncopies = 1 )

### Step 5 - gene expression associations
OutDir = paste0(MainDir,"Expression/")
dir.create(OutDir)
DifferentialExpression_pooled( OutDir,comparisons_adj.P.Val_thresh=0.1 )
ComparisonDegs( OutDir,comparisons_adj.P.Val_thresh=0.1 )
ProliferationScoring( OutDir )

### Step 6 - tumor microenvironment associations
OutDir = paste0(MainDir,"Microenvironment/")
dir.create(OutDir)
purity_associations(OutDir)
tme_associations(OutDir)

### Step 7 - therapy associations
OutDir = paste0(MainDir,"Therapy/")
dir.create(OutDir)
therapy_associations( OutDir )
OutDir = paste0(MainDir,"Therapy/Regimens_FirstDrug_PatientLevel/")
dir.create(OutDir)
regimen_associations_FirstDrug_PatientLevel( OutDir )
OutDir = paste0(MainDir,"Therapy/Regimens_FirstDrug_RegimenLevel/") # patients can be repeated
dir.create(OutDir)
regimen_associations_FirstDrug_RegimenLevel( OutDir )
OutDir = paste0(MainDir,"Therapy/Regimens_SingleDrugLevel/") # patients can be repeated
dir.create(OutDir)
regimen_associations_SingleDrugLevel( OutDir )

### Step 8 - suppl tables for paper submission
OutDir = paste0( MainDir,"Paper_SupplTables/" )
dir.create(OutDir)
supplementary_tables_formatting( OutDir )

### Step 9 - paper revision
OutDir = paste0( MainDir,"Revision/L858R_compound_vs_single/" )
dir.create(OutDir, recursive = TRUE)
revision_L858R_compound_vs_single( OutDir )
OutDir = paste0( MainDir,"Revision/ex19del_compound_vs_single/" )
dir.create(OutDir, recursive = TRUE)
revision_ex19del_compound_vs_single( OutDir )
OutDir = paste0( MainDir,"Revision/DoubleUncommon_compound_vs_anyCommon/" )
dir.create(OutDir, recursive = TRUE)
revision_DoubleUncommon_compound_vs_anyCommon( OutDir )
OutDir = paste0( MainDir,"Revision/krasmut_vs_treatment/" )
dir.create(OutDir, recursive = TRUE)
krasmut_vs_treatment( OutDir )











