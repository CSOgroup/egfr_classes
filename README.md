# egfr_classes

Repository containing the scripts for the analyses in the manuscript _Decoding the Clinical and Molecular Signatures of EGFR Common, Compound, and Uncommon Mutations in Non-Small Cell Lung Cancer_

All the steps are performed following __EGFRclasses_pipeline.R__. The clinical and molecular data for the various datasets was retreived from the sources cited in the manuscript. 

__EGFRclasses_MutationListConstruction_step1.R__, __EGFRclasses_MutationListConstruction_step2.R__:
Construction of the list of EGFR mutations, assignment to the classes, filtering.

__EGFRclasses_Preprocessing.R__:
Collection and harmonization of clinical data, assignment to EGFR classes.

__EGFRclasses_functions.R__:
Functions implementing each analysis.

__EGFRclasses_pipeline.R__:
Pipeline calling the analysis functions.

__EGFRclasses_utils.R__:
Generic helper functions for plotting.
