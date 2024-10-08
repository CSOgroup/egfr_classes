
### sourcing helper and analysis functions
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_utils.R")
source( "/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_functions.R" )

### constructing EGFR mutation list
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_MutationListConstruction.R")

### Preprocessing clinical data
source("/mnt/ndata/daniele/alfredo_egfr/Scripts/EGFRclasses_Preprocessing.R")

load(file = "/mnt/ndata/daniele/alfredo_egfr/Processed/EGFR_classes/ec_exons18_21_oncogenic_alphamissense_filtered.RData" ) # new version, including china

### Define pipeline run name/ID and MainDir
MainDir = "/mnt/ndata/daniele/alfredo_egfr/Processed/EGFRclasses_pipeline/Repro/"
dir.create(MainDir)

### Step 1 - assigning patients
OutDir = paste0(MainDir,"Preprocessing/")
dir.create(OutDir)
patients_assignment( OutDir, retain = c( "common","uncommon","compound"), using_AlphaMissense = TRUE )
compounding_uncommons( OutDir )
vaf_class_analysis( OutDir )
compare_borgeaud_robichaux( ec,OutDir )
# Run mutation_mapper here: https://www.cbioportal.org/mutation_mapper for lollipop plot

### Step 2 - clinical associations
OutDir = paste0(MainDir,"Clinical/")
dir.create(OutDir)
clinical_associations( OutDir )

### Step 3 - mutational signatures analysis
OutDir = paste0(MainDir,"MutSignatures/") # analysis on targeted sequencing supported by: https://pubmed.ncbi.nlm.nih.gov/33571361/
dir.create(OutDir)
mutationalSignatures_preprocessing( OutDir )
# mutationalSignatures_associations( OutDir )
mutationalSignatures_associations_binary( OutDir, threshold = 0.3 )
mutationalSignatures_associations_binary( OutDir, threshold = 0 )
mutationalSignatures_associations_binary( OutDir, threshold = 0.5 )
mutationalSignatures_egfr_muts_contexts( OutDir ) # only common vs uncommon

### Step 4 - genomic associations
OutDir = paste0(MainDir,"Genomic_TmbFga/")
dir.create(OutDir)
TmbFga_associations( OutDir )
tracerx_analyses( OutDir )
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
ViperAnalysis(OutDir)
DifferentialProtein(OutDir)

### Step 6 - DNA methylation associations
OutDir = paste0(MainDir,"Methylation/")
dir.create(OutDir)
contr.matrix <- makeContrasts(
	   CommonVsUncommon = common-uncommon, 
	   CommonVsCompound = common-compound,
	   UncommonVsCompound = uncommon-compound,
	   levels = c( "common","uncommon","compound" ))
methylation_associations(OutDir, contr.matrix)

### Step 7 - tumor microenvironment associations
OutDir = paste0(MainDir,"Microenvironment/")
dir.create(OutDir)
purity_associations(OutDir)
tme_associations(OutDir)

### Step 8 - therapy associations
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

### Various
OutDir = paste0( MainDir,"Paper_SupplTables/" )
dir.create(OutDir)
supplementary_tables_formatting( OutDir )













