
library(reshape2)
library(ggpubr)
library(limma)
library(edgeR)
source("/mnt/ndata/daniele/lung_multiregion/sc/cell-state-inference/utils/dan.functions.R")
DataDir = "/mnt/ndata/daniele/alfredo_egfr/Data/"
CommonDataDir = "/mnt/ed2/daniele/Common_Data/"
ordered_classes = c( "common_ex19del","common_L858R", "common","uncommon","compound","T790M","ex20ins" )
colorz_classes = c("deepskyblue4","deepskyblue3","steelblue4","tomato3","mediumpurple4","darksalmon","orange2")
clco = data.frame(row.names = ordered_classes, classes = ordered_classes, colorz = colorz_classes)
size_labels = 6/.pt

compare_borgeaud_robichaux = function( OutDir ){
	load(file = "/mnt/ndata/daniele/alfredo_egfr/Processed/EGFR_classes/ec_exons18_21_oncogenic_alphamissense_filtered.RData" ) # new version, including china
	ec = ec[ec$Class %in% c( "common","uncommon" ),]
	ec = ec[substr(ec$Mutation,1,2)=="p.",] # exclude c. mutations of TRACERx
	ec$Mutation = substr(ec$Mutation,3,nchar(ec$Mutation))
	ec = ec[!duplicated(ec$Mutation),]
	library(eulerr)
	# read robichaux TableS4
	er = dan.read(paste0(DataDir,"Robi_mutations_tableS4.txt"))
	er = er[!grepl(",",er$mutations),]
	er = er[er$four_classes %in% c("classical-like","PACC"),]
	# common vs classical/PACC
	# er also contains also ex19del. Let's include them directly from ec, and remove the duplicates
	er = rbind(er, data.frame( four_classes="classical-like",six_classes="classical-like",mutations=ec[ec$Class=="common","Mutation"],stringsAsFactors=F ))
	er = er[!duplicated(er$mutations),]
	er = er[!(er$mutations=="Ex19del"),]
	classical = sum(er$four_classes=="classical-like")
	pacc = sum(er$four_classes=="PACC")
	common = sum(ec$Class=="common")
	uncommon = sum(ec$Class=="uncommon")
	classical_pacc = 0
	classical_common = length(intersect(er[er$four_classes=="classical-like","mutations"],ec[ec$Class=="common","Mutation"]))
	pacc_common = length(intersect(er[er$four_classes=="PACC","mutations"],ec[ec$Class=="common","Mutation"]))
	classical_uncommon = length(intersect(er[er$four_classes=="classical-like","mutations"],ec[ec$Class=="uncommon","Mutation"]))
	pacc_uncommon = length(intersect(er[er$four_classes=="PACC","mutations"],ec[ec$Class=="uncommon","Mutation"]))
	fit1 <- euler(c("A"=classical-classical_common,"B"=common-classical_common-pacc_common,"C"=pacc-pacc_common,"A&B"=classical_common,"A&C"=0,"B&C"=pacc_common,"A&B&C"=0))
	pdf(paste0(OutDir,"borgeaud_vs_robichaux_common.pdf"),2.5,2.5, useDingbats=F,onefile=F,pointsize=6)
	plot(fit1,labels = c("Classical-like","common","PACC"),fills = c("pink1","steelblue4","dodgerblue3"),edges = FALSE,fontsize=6,quantities = list(fontsize = 6))
	dev.off()

	# uncommon vs classical/PACC
	fit1 <- euler(c("A"=classical-classical_uncommon,"B"=common-classical_uncommon-pacc_uncommon,"C"=pacc-pacc_uncommon,"A&B"=classical_uncommon,"A&C"=0,"B&C"=pacc_uncommon,"A&B&C"=0))
	pdf(paste0(OutDir,"borgeaud_vs_robichaux_uncommon.pdf"),2.5,2.5, useDingbats=F,onefile=F,pointsize=6)
	plot(fit1,labels = c("Classical-like","uncommon","PACC"),fills = c("pink1","tomato3","dodgerblue3"),edges = FALSE,fontsize=6,quantities = list(fontsize = 6))
	dev.off()
	
	xvenn = list(Robichaux_et_al_2021=pr[(xpval<cutoff_fdr) & (x>0),"gene"],Borgeaud_et_al_2024=pr[(ypval<cutoff_fdr) & (y>0),"gene"])
	pdf(paste0(OutDir,"ATG3OE_vs_ATG7OE_up_venn.pdf"),7,7)
	ggvenn( xvenn,fill_color = c("darkorange1", "goldenrod"),stroke_size = 0.5,show_elements=F, set_name_size = 4,label_sep='\n', text_size=3)
	dev.off()
}

patients_assignment = function( OutDir, retain, using_AlphaMissense = TRUE, split_common = FALSE, compound_to_uncommon = FALSE ){
	for (dataset in c( "TCGA","Chen","Genie","Zhang","TRACERx421","Origimed" ) ){
		if ( using_AlphaMissense ){
			load(file = paste0( DataDir,"using_AlphaMissense/Clin_",dataset,"_complete.RData" ))
		} else {
			load(file = paste0( DataDir,"Clin_",dataset,"_complete.RData" ))	
		}
		
		Clin = Clin[Clin$egfr_class_consensus %in% retain,]
		if (split_common){
			double_commons = Clin[(Clin$egfr_class_consensus=="common") & ((Clin$egfr_class_2=="common") %in% c(T) ),]
			double_commons = double_commons[ (double_commons$egfr_mut_1=="p.L858R") | (double_commons$egfr_mut_2=="p.L858R"), ]
			Clin = Clin[!(rownames(Clin) %in% rownames(double_commons)),]
			Clin[Clin$egfr_class_consensus=="common","egfr_class_consensus"] = ifelse( Clin[Clin$egfr_class_consensus=="common","egfr_mut_1"]=="p.L858R","common_L858R","common_ex19del" )
		}
		if (compound_to_uncommon){
			Clin[Clin$egfr_class_consensus=="compound","egfr_class_consensus"] = "uncommon"
		}
		save(Clin, file = paste0( OutDir,"Clin_",dataset,".RData" ))
		if ( using_AlphaMissense ){
			load(file = paste0( DataDir,"using_AlphaMissense/Clin2_",dataset,"_SampleLevel_complete.RData" ))
		} else {
			load(file = paste0( DataDir,"Clin2_",dataset,"_SampleLevel_complete.RData" ))
		}
		Clin2 = Clin2[Clin2$egfr_class_consensus %in% retain,]
		if (split_common){
			double_commons = Clin2[(Clin2$egfr_class_consensus=="common") & ((Clin2$egfr_class_2=="common") %in% c(T) ),]
			double_commons = double_commons[ (double_commons$egfr_mut_1=="p.L858R") | (double_commons$egfr_mut_2=="p.L858R"), ]
			Clin2 = Clin2[!(rownames(Clin2) %in% rownames(double_commons)),]
			Clin2[Clin2$egfr_class_consensus=="common","egfr_class_consensus"] = ifelse( Clin2[Clin2$egfr_class_consensus=="common","egfr_mut_1"]=="p.L858R","common_L858R","common_ex19del" )
		}
		if (compound_to_uncommon){
			Clin2[Clin2$egfr_class_consensus=="compound","egfr_class_consensus"] = "uncommon"
		}
		save(Clin2, file = paste0( OutDir,"Clin2_",dataset,"_SampleLevel.RData" ))
	}
	first = TRUE
	for (dataset in c( "Genie","TCGA","Chen","TRACERx421","Zhang","Origimed" )){
		load(file = paste0( OutDir,"Clin_",dataset,".RData" ))
		maxid = max(as.numeric(gsub("egfr_mut_","",colnames(Clin)[substr(colnames(Clin),1,9)=="egfr_mut_"])))
		this_cb = dan.df( 0,c("Hugo_Symbol","Sample_ID","Protein_Change","Class" ))
		for (id in 1:maxid){
			this_cb = rbind(this_cb,data.frame( Hugo_Symbol="EGFR",Sample_ID=Clin$Patient,Protein_Change=gsub("p.","",Clin[,paste0( "egfr_mut_",id )]),Class=Clin[,paste0( "egfr_class_",id )],stringsAsFactors=F ))
		}
		this_cb = this_cb[!is.na(this_cb$Class),]
		this_cb$Dataset = dataset
		if (first){
			cb = this_cb
			first = FALSE
		} else {
			cb = rbind(cb,this_cb)
		}
	}
	## Lollipop plot of uncommon using patients (to see the most recurrent)
	cba = cb[!duplicated(cb$Protein_Change),c( "Protein_Change","Class" )]
	rownames(cba) = cba$Protein_Change
	tcb = cb[,c("Hugo_Symbol","Sample_ID","Protein_Change")]
	prefix = "all"
	dan.write(tcb,file=paste0( OutDir,"LollipopPlot_",prefix,"_table.txt" ))
	dcat( paste0("Frequency L858R: ",sum(tcb$Protein_Change=="L858R")/nrow(tcb) ))
	dcat( paste0("Frequency E746_A750del: ",sum(tcb$Protein_Change=="E746_A750del")/nrow(tcb) ))
	tabb = sort(dtable(tcb$Protein_Change),decreasing=T)
	tabb = tabb[tabb>=20]
	dff = data.frame(Protein_Change = names(tabb), Occurrences = as.numeric(tabb), Class = cba[names(tabb),"Class"],stringsAsFactors=F)
	colorz_classes = clco[rownames(clco) %in% (dff$Class),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (dff$Class),"classes"]
	fileName = paste0( OutDir,"MutOccurrences_",prefix,"_barplot.pdf" )
	dan.barplot( fileName, factor(dff$Protein_Change,levels = as.character(dff$Protein_Change)), dff$Occurrences, fill = factor(dff$Class,levels=ordered_classes), xlab = "", ylab = "Occurrences", ylimLeft = NULL, ylimRight = NULL, filllab = "Class", fillColors = colorz_classes, plotTitle = "EGFR mutation occurrences (only >= 20 shown), any class", fileWidth = 4, fileHeight = 2 )
	dan.save(dff,fileName)
	tcb = cb[cb$Class=="uncommon",c("Hugo_Symbol","Sample_ID","Protein_Change")]
	dcat( paste0("Frequency G719A: ",sum(tcb$Protein_Change=="G719A")/nrow(tcb)))
	prefix = "uncommon"
	dan.write(tcb,file=paste0( OutDir,"LollipopPlot_",prefix,"_table.txt" ))
	tabb = sort(dtable(tcb$Protein_Change),decreasing=T)
	tabb = tabb[tabb>=20]
	dff = data.frame(Protein_Change = names(tabb), Occurrences = as.numeric(tabb), Class = cba[names(tabb),"Class"],stringsAsFactors=F)
	colorz_classes = clco[rownames(clco) %in% (dff$Class),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (dff$Class),"classes"]
	fileName = paste0( OutDir,"MutOccurrences_",prefix,"_barplot.pdf" )
	dan.barplot( fileName, factor(dff$Protein_Change,levels = as.character(dff$Protein_Change)), dff$Occurrences, fill = factor(dff$Class,levels=ordered_classes), xlab = "", ylab = "Occurrences", ylimLeft = NULL, ylimRight = NULL, filllab = "Class", fillColors = colorz_classes, plotTitle = "EGFR mutation occurrences (only >= 20 shown), uncommon only", fileWidth = 4, fileHeight = 2 )
	dan.save(dff,fileName)
	tcb = cb[cb$Class=="common",c("Hugo_Symbol","Sample_ID","Protein_Change")]
	prefix = "common"
	dan.write(tcb,file=paste0( OutDir,"LollipopPlot_",prefix,"_table.txt" ))
	tabb = sort(dtable(tcb$Protein_Change),decreasing=T)
	tabb = tabb[tabb>=20]
	dff = data.frame(Protein_Change = names(tabb), Occurrences = as.numeric(tabb), Class = cba[names(tabb),"Class"],stringsAsFactors=F)
	colorz_classes = clco[rownames(clco) %in% (dff$Class),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (dff$Class),"classes"]
	fileName = paste0( OutDir,"MutOccurrences_",prefix,"_barplot.pdf" )
	dan.barplot( fileName, factor(dff$Protein_Change,levels = as.character(dff$Protein_Change)), dff$Occurrences, fill = factor(dff$Class,levels=ordered_classes), xlab = "", ylab = "Occurrences", ylimLeft = NULL, ylimRight = NULL, filllab = "Class", fillColors = colorz_classes, plotTitle = "EGFR mutation occurrences (only >= 20 shown), common only", fileWidth = 4, fileHeight = 2 )
	tcb = cb[,c("Hugo_Symbol","Sample_ID","Protein_Change","Class")]
	prefix = "common_aggregated"
	tcb[(tcb$Class=="common") & (tcb$Protein_Change!="L858R"),"Protein_Change"] = "ex19del"
	tabb = sort(dtable(tcb$Protein_Change),decreasing=T)
	tabb = tabb[tabb>=20]
	cba["ex19del","Class"] = "common"
	dff = data.frame(Protein_Change = names(tabb), Occurrences = as.numeric(tabb), Class = cba[names(tabb),"Class"],stringsAsFactors=F)
	colorz_classes = clco[rownames(clco) %in% (dff$Class),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (dff$Class),"classes"]
	fileName = paste0( OutDir,"MutOccurrences_",prefix,"_barplot.pdf" )
	dan.barplot( fileName, factor(dff$Protein_Change,levels = as.character(dff$Protein_Change)), dff$Occurrences, fill = factor(dff$Class,levels=ordered_classes), xlab = "", ylab = "Occurrences", ylimLeft = NULL, ylimRight = NULL, filllab = "Class", fillColors = colorz_classes, plotTitle = "EGFR mutation occurrences (only >= 20 shown)", fileWidth = 4, fileHeight = 2 )
}

compounding_uncommons = function( OutDir ){
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Clin = Clin[Clin$egfr_class_consensus!="common",]
	tc = dan.df(0,c( "EGFR_consensus","EGFR_uncommon_single","Dataset" ))
	for (rn in rownames(Clin)){
		consensus = Clin[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_uncommon_single"))
		for (i in 1:5){
			if ( (Clin[rn,paste0("egfr_class_",i)]=="uncommon") %in% c(T) ) { 
				a = Clin[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ))]
				colnames(a) = c("EGFR_consensus","EGFR_uncommon_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = "Genie"
		tc = rbind(tc,tcc)
	}
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clin = Clin[Clin$egfr_class_consensus!="common",]
	for (rn in rownames(Clin)){
		consensus = Clin[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_uncommon_single"))
		for (i in 1:2){
			if ( (Clin[rn,paste0("egfr_class_",i)]=="uncommon") %in% c(T) ) { 
				a = Clin[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ))]
				colnames(a) = c("EGFR_consensus","EGFR_uncommon_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = "TCGA"
		tc = rbind(tc,tcc)
	}
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clin = Clin[Clin$egfr_class_consensus!="common",]
	for (rn in rownames(Clin)){
		consensus = Clin[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_uncommon_single"))
		for (i in 1:3){
			if ( (Clin[rn,paste0("egfr_class_",i)]=="uncommon") %in% c(T) ) { 
				a = Clin[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ))]
				colnames(a) = c("EGFR_consensus","EGFR_uncommon_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = "Chen"
		tc = rbind(tc,tcc)
	}
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	Clin = Clin[Clin$egfr_class_consensus!="common",]
	for (rn in rownames(Clin)){
		consensus = Clin[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_uncommon_single"))
		for (i in 1:2){
			if ( (Clin[rn,paste0("egfr_class_",i)]=="uncommon") %in% c(T) ) { 
				a = Clin[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ))]
				colnames(a) = c("EGFR_consensus","EGFR_uncommon_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = "TRACERx421"
		tc = rbind(tc,tcc)
	}
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	Clin = Clin[Clin$egfr_class_consensus!="common",]
	for (rn in rownames(Clin)){
		consensus = Clin[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_uncommon_single"))
		for (i in 1:2){
			if ( (Clin[rn,paste0("egfr_class_",i)]=="uncommon") %in% c(T) ) { 
				a = Clin[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ))]
				colnames(a) = c("EGFR_consensus","EGFR_uncommon_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = "Zhang"
		tc = rbind(tc,tcc)
	}
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	Clin = Clin[Clin$egfr_class_consensus!="common",]
	for (rn in rownames(Clin)){
		consensus = Clin[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_uncommon_single"))
		for (i in 1:3){
			if ( (Clin[rn,paste0("egfr_class_",i)]=="uncommon") %in% c(T) ) { 
				a = Clin[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ))]
				colnames(a) = c("EGFR_consensus","EGFR_uncommon_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = "Origimed"
		tc = rbind(tc,tcc)
	}
	save(tc,file=paste0(OutDir,"UncommonVariant_enriched_inCompound_volcanoPlot_tc.RData"))
	fileName = paste0(OutDir,"UncommonVariant_enriched_inCompound_volcanoPlot.pdf")
	cu = tc
	tabb = table(cu[,"EGFR_consensus"],cu$EGFR_uncommon_single)
	# tabb = tabb[,colSums(tabb)>=3]
	cdf = dan.df(colnames(tabb),c( "Variant","Chisq_residual_compound","Chisq_pval" ))
	for (g in colnames(tabb)){
		cdf[g,"Variant"] = g
		tabb2 = table(cu[,"EGFR_consensus"],cu$EGFR_uncommon_single==g)
		if (g %in% c( "p.L861Q","p.S768I" )){ print(tabb2) }
		ch = chisq.test(tabb2)
		cdf[g,"Chisq_pval"] = ch$p.value
		cdf[g,"Chisq_residual_compound"] = ch$residuals["compound","TRUE"]
	}
	cdf = cdf[order(cdf$Chisq_pval),]
	cdf$Chisq_qval = p.adjust(cdf$Chisq_pval,method='BH')
	cdf$fill = ifelse( (cdf$Chisq_qval>=0.1) | (abs(cdf$Chisq_residual_compound)<=1),"Not significant",ifelse( cdf$Chisq_residual_compound<0,"Higher in uncommon","Higher in compound"))
	cdf$fill = factor(cdf$fill,levels=c( "Higher in uncommon","Higher in compound","Not significant" ))
	fillColorz = c( "tomato3","mediumpurple4","gray88" )
	cdf$repel_labelz = rownames(cdf)
	cdf[cdf$fill=="Not significant","repel_labelz"] = ""
	dan.scatterplot( fileName, cdf$Chisq_residual_compound, -log10(cdf$Chisq_pval), fill = cdf$fill, xlab = "Chi-square residual in compound", ylab = "-log10(Chi-square p-value)", filllab = "", fillColors = fillColorz, dotSize = 2, repel_labels=cdf$repel_labelz, coord_fixed = FALSE, fileWidth = 5, fileHeight = 2.5 )
	dan.save(cdf,fileName)
}

vaf_class_analysis = function( OutDir ){
	dataset = "Genie"
	load(file = paste0( OutDir,"../Preprocessing/","Clin_",dataset,".RData" ))
	load(file = paste0( OutDir,"Clin2_",dataset,"_SampleLevel.RData" ))
	tc = dan.df(0,c( "EGFR_consensus","EGFR_mut_single","EGFR_class_single","Dataset","Sample","Patient" ))
	for (rn in rownames(Clin2)){
		consensus = Clin2[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_mut_single","EGFR_class_single"))
		for (i in 1:5){
			if (!is.na(Clin2[rn,paste0("egfr_class_",i)])){
				a = Clin2[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ),paste0( "egfr_class_",i ))]	
				colnames(a) = c("EGFR_consensus","EGFR_mut_single","EGFR_class_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = dataset
		tcc$Sample = rn
		tcc$Patient = Clin2[rn,"Patient"]
		tc = rbind(tc,tcc)
	}
	load(file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
	maf_ss = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
	maf_ss$vaf = as.numeric(maf_ss$t_alt_count)/as.numeric(maf_ss$t_depth)
	for (rn in rownames(tc)){
		s = tc[rn,"Sample"]
		mut = tc[rn,"EGFR_mut_single"]
		tc[rn,"vaf"] = mean(maf_ss[(maf_ss$Tumor_Sample_Barcode %in% s) & (maf_ss$HGVSp_Short==mut),"vaf"])
	}
	tcAll = tc
	dataset = "TCGA"
	load(file = paste0( OutDir,"../Preprocessing/","Clin_",dataset,".RData" ))
	load(file = paste0( OutDir,"Clin2_",dataset,"_SampleLevel.RData" ))
	tc = dan.df(0,c( "EGFR_consensus","EGFR_mut_single","EGFR_class_single","Dataset","Sample","Patient" ))
	for (rn in rownames(Clin2)){
		consensus = Clin2[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_mut_single","EGFR_class_single"))
		for (i in 1:2){
			if (!is.na(Clin2[rn,paste0("egfr_class_",i)])){
				a = Clin2[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ),paste0( "egfr_class_",i ))]	
				colnames(a) = c("EGFR_consensus","EGFR_mut_single","EGFR_class_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = dataset
		tcc$Sample = rn
		tcc$Patient = Clin2[rn,"Patient"]
		tc = rbind(tc,tcc)
	}
	load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
	maf_ss = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
	maf_ss$vaf = as.numeric(maf_ss$t_alt_count)/as.numeric(maf_ss$t_depth)
	for (rn in rownames(tc)){
		s = tc[rn,"Sample"]
		mut = tc[rn,"EGFR_mut_single"]
		tc[rn,"vaf"] = mean(maf_ss[(substr(maf_ss$Tumor_Sample_Barcode,1,16) %in% s) & (maf_ss$HGVSp_Short==mut),"vaf"])
	}
	tcAll = rbind(tcAll,tc)
	dataset = "Chen"
	load(file = paste0( OutDir,"../Preprocessing/","Clin_",dataset,".RData" ))
	load(file = paste0( OutDir,"Clin2_",dataset,"_SampleLevel.RData" ))
	tc = dan.df(0,c( "EGFR_consensus","EGFR_mut_single","EGFR_class_single","Dataset","Sample","Patient" ))
	for (rn in rownames(Clin2)){
		consensus = Clin2[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_mut_single","EGFR_class_single"))
		for (i in 1:3){
			if (!is.na(Clin2[rn,paste0("egfr_class_",i)])){
				a = Clin2[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ),paste0( "egfr_class_",i ))]	
				colnames(a) = c("EGFR_consensus","EGFR_mut_single","EGFR_class_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = dataset
		tcc$Sample = rn
		tcc$Patient = Clin2[rn,"Patient"]
		tc = rbind(tc,tcc)
	}
	maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
	maf$Patient = maf$Tumor_Sample_Barcode
	maf_ss = maf
	maf_ss = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
	maf_ss$vaf = as.numeric(maf_ss$t_alt_count)/(as.numeric(maf_ss$t_ref_count)+as.numeric(maf_ss$t_alt_count))
	for (rn in rownames(tc)){
		s = tc[rn,"Sample"]
		mut = tc[rn,"EGFR_mut_single"]
		tc[rn,"vaf"] = mean(maf_ss[(maf_ss$Tumor_Sample_Barcode %in% s) & (maf_ss$HGVSp_Short==mut),"vaf"])
	}
	tcAll = rbind(tcAll,tc)
	dataset = "TRACERx421"
	load(file = paste0( OutDir,"../Preprocessing/","Clin_",dataset,".RData" ))
	load(file = paste0( OutDir,"Clin2_",dataset,"_SampleLevel.RData" ))
	tc = dan.df(0,c( "EGFR_consensus","EGFR_mut_single","EGFR_class_single","Dataset","Sample","Patient" ))
	for (rn in rownames(Clin2)){
		consensus = Clin2[rn,"egfr_class_consensus"]
		tcc = dan.df(0,c("EGFR_consensus","EGFR_mut_single","EGFR_class_single"))
		for (i in 1:2){
			if (!is.na(Clin2[rn,paste0("egfr_class_",i)])){
				a = Clin2[rn,c("egfr_class_consensus",paste0( "egfr_mut_",i ),paste0( "egfr_class_",i ))]	
				colnames(a) = c("EGFR_consensus","EGFR_mut_single","EGFR_class_single")
				tcc = rbind(tcc,a)
			}
		}
		tcc$Dataset = dataset
		tcc$Sample = rn
		tcc$Patient = Clin2[rn,"Patient"]
		tc = rbind(tc,tcc)
	}
	library(fst)
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	mafr = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221123_TRACERx421_mutation_table_region.fst"))
	maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$tumour_id
	maf_ss$Patient = maf_ss$patient_id
	maf_ss = maf_ss[maf_ss$Hugo_Symbol=="EGFR",]
	for (rn in rownames(maf_ss)){
		thiz = mafr[mafr$mutation_id==maf_ss[rn,"mutation_id"],]
		maf_ss[rn,"vaf"] = mean(as.numeric(thiz$var_count)/(as.numeric(thiz$depth)))
	}
	maf_ss$HGVSp_Short = maf_ss$AAChange
	maf_ss[is.na(maf_ss$HGVSp_Short),"HGVSp_Short"] = maf_ss[is.na(maf_ss$HGVSp_Short),"NucleotideChange"]
	for (rn in rownames(tc)){
		s = tc[rn,"Sample"]
		mut = tc[rn,"EGFR_mut_single"]
		tc[rn,"vaf"] = mean(maf_ss[(maf_ss$patient_id %in% s) & (maf_ss$HGVSp_Short==mut),"vaf"])
	}
	tcAll = rbind(tcAll,tc)
	## no vaf info for Origimed, Zhang
	tcAll_save = tcAll
	tcAll = tcAll_save[(((tcAll_save$vaf<=1) & (tcAll_save$vaf>=0))) %in% c(T),]
	save(tcAll,file=paste0(OutDir,"tcAll_vafs.RData"))
	x = factor(tcAll[,"Dataset"])
	y = tcAll$vaf
	dan.densityPlot( paste0(OutDir,"vaf_across_Datasets.pdf"), y, x, groupinglab = "", xlab = paste0("Variant allele frequency"), ylab = "Density", show_medians = F, plotTitle = "",xlimLeft = 0, xlimRight = 1, groupingColors = c("red","blue","green","yellow" ), fileWidth = 3.5, fileHeight = 2 )
	colorz_classes = clco[rownames(clco) %in% (tc$EGFR_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (tc$EGFR_consensus),"classes"]
	# consensus class
	x = factor(tcAll[,"EGFR_consensus"],levels = ordered_classes)
	plotTitle = paste0( "Kruskal-Wallis test, p-value = ", signif(kruskal.test(y~x)$p.value,2))
	dan.densityPlot( paste0(OutDir,"vaf_across_EGFR_consensus.pdf"), y, x, groupinglab = "", xlab = paste0("Variant allele frequency"), ylab = "Density", show_medians = F, plotTitle = plotTitle,xlimLeft = 0, xlimRight = 1, groupingColors = colorz_classes, fileWidth = 3.5, fileHeight = 2 )
	# single mutations
	x = factor(tcAll[,"EGFR_class_single"],levels = c( "common","uncommon" ))
	plotTitle = paste0( "Wilcoxon test, p-value = ", signif(wilcox.test(y~x)$p.value,2))
	these_colorz = colorz_classes[ordered_classes %in% levels(x)]
	dan.densityPlot( paste0(OutDir,"vaf_across_EGFR_single_mutations.pdf"), y, x, groupinglab = "", xlab = paste0("Variant allele frequency"), ylab = "Density", show_medians = F, plotTitle = plotTitle,xlimLeft = 0, xlimRight = 1, groupingColors = these_colorz, fileWidth = 3.5, fileHeight = 2 )
	# single uncommon - consensus uncommon vs compound
	thiz = tcAll[tcAll$EGFR_class_single=="uncommon",]
	x = factor(thiz[,"EGFR_consensus"],levels = c( "uncommon","compound" ))
	y = thiz$vaf
	plotTitle = paste0( "Wilcoxon test, p-value = ", signif(wilcox.test(y~x)$p.value,2))
	these_colorz = colorz_classes[ordered_classes %in% levels(x)]
	dan.densityPlot( paste0(OutDir,"vaf_across_EGFR_single_uncommon_acrossConsensus.pdf"), y, x, groupinglab = "", xlab = paste0("Variant allele frequency"), ylab = "Density", show_medians = F, plotTitle = plotTitle,xlimLeft = 0, xlimRight = 1, groupingColors = these_colorz, fileWidth = 3.5, fileHeight = 2 )
	# single common - consensus common vs compound
	thiz = tcAll[tcAll$EGFR_class_single=="common",]
	x = factor(thiz[,"EGFR_consensus"],levels = c( "common","compound" ))
	y = thiz$vaf
	plotTitle = paste0( "Wilcoxon test, p-value = ", signif(wilcox.test(y~x)$p.value,2))
	these_colorz = colorz_classes[ordered_classes %in% levels(x)]
	dan.densityPlot( paste0(OutDir,"vaf_across_EGFR_single_common_acrossConsensus.pdf"), y, x, groupinglab = "", xlab = paste0("Variant allele frequency"), ylab = "Density", show_medians = F, plotTitle = plotTitle,xlimLeft = 0, xlimRight = 1, groupingColors = these_colorz, fileWidth = 3.5, fileHeight = 2 )
}

fisher_tests = function( clin, vars ){
	for (var in vars){
		dcat(var)
		tc = clin[!is.na(clin[,var]),]
		ft = fisher.test(table(tc[,var],tc$egfr_class_consensus),simulate.p.value = TRUE)
		print(ft)
	}
}

clinical_associations = function( OutDir ){

	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clint = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TCGA_SampleLevel.RData" ))
	Clin2t = Clin2
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clinc = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Chen_SampleLevel.RData" ))
	Clin2c = Clin2
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	Clinx = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	Clinz = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	Clino = Clin


	colorz_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"classes"]

	## Let's do it at patient level. Survival excluded, for now
	rm(Clin)
	rm(Clin2)
	common_vars = c( "Stage","Sex","Race","Ethnicity","Smoking","Pattern","Age" )
	Cling$Age = as.numeric(Cling$Age)
	Clin = Cling
	save(Clin, file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Cling[,c("Patient","egfr_class_consensus", common_vars )]

	Clint$Stage = NA
	Clint[ Clint$pathologic_stage %in% c( "stage i","stage ia","stage ib" )  ,"Stage" ] = "I"
	Clint[ Clint$pathologic_stage %in% c( "stage iia","stage iib" )  ,"Stage" ] = "II"
	Clint[ Clint$pathologic_stage %in% c( "stage iiia","stage iiib" )  ,"Stage" ] = "III"
	Clint[ Clint$pathologic_stage %in% c( "stage iv" )  ,"Stage" ] = "IV"
	Clint$Sex = Clint$gender
	Clint$Race = Clint$race
	Clint[(Clint$Race=="black or african american") %in% c(T),"Race"] = "black"
	Clint$Ethnicity = tolower(Clint$ethnicity)
	Clint[(Clint$Ethnicity=="not hispanic or latino") %in% c(T),"Ethnicity"] = "non-hispanic"
	Clint$Smoking = NA
	Clint[(Clint$tobacco_smoking_history=="1") %in% c(T),"Smoking"] = "never smoker"
	Clint[(Clint$tobacco_smoking_history %in% c("2","3","4") ) %in% c(T),"Smoking"] = "ever smoker"
	Clint$Age = as.numeric(Clint$age_at_initial_pathologic_diagnosis)
	Clin = Clint
	save(Clin, file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))

	Clint = Clint[,c("Patient","egfr_class_consensus", common_vars )]

	Clinc$Sex = tolower(Clinc$Gender)
	Clinc$Race = "asian"
	Clinc$Ethnicity = "non-hispanic"
	Clinc$Smoking = ifelse( (Clinc$Smoker=="Yes") %in% c(T),"ever smoker", ifelse( (Clinc$Smoker=="No") %in% c(T),"never smoker",NA ))
	Clinc$Age = as.numeric(Clinc$Age)
	Clin = Clinc
	save(Clin, file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))

	Clinc = Clinc[,c("Patient","egfr_class_consensus", common_vars )]

	Clinx = Clinx[,c("Patient","egfr_class_consensus", common_vars )]
	Clinz = Clinz[,c("Patient","egfr_class_consensus", common_vars )]
	Clino = Clino[,c("Patient","egfr_class_consensus", common_vars )]

	## Tests
	dcat( "Genie" )
	fisher_tests( Cling, common_vars )
	dcat( "TCGA" )
	fisher_tests( Clint, c( "Stage","Sex","Race","Smoking","Pattern" ) )
	dcat( "Chen" )
	fisher_tests( Clinc, c( "Stage","Sex","Smoking","Pattern" ) )
	dcat( "TRACERx421" )
	fisher_tests( Clinx, c( "Stage","Sex","Race","Smoking","Pattern" ) )
	dcat( "Origimed" )
	fisher_tests( Clino, c( "Stage","Sex","Smoking" ) )

	Clinall = rbind(rbind(rbind(Cling,Clint),rbind(Clinc,Clinx) ),Clinz)
	Clinall = rbind(Clinall,Clino)

	dcat( "All datasets" )
	fisher_tests( Clinall, c( "Stage","Sex","Race","Smoking","Pattern" ) )

	# Sex: negative association
 	tc = Clinall[!is.na(Clinall[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = 6/.pt) + geom_text(aes(label=round(Percentage,1)), size = 6/.pt, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_Sex_onlyAll_barplot.pdf"),2.8,2.8,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()

	tc = Clinall[!is.na(Clinall[,"Stage"]),]
	tabb = table(tc[,"Stage"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Stage","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Stage=="I") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Stage=="I"))
	mt[(mt$Stage=="II") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Stage=="II"))
	mt[(mt$Stage=="III") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Stage=="III"))
	mt[(mt$Stage=="IV") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Stage=="IV"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=Stage, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("Tumor stage") + ggtitle(paste0("All datasets\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = 6/.pt) + geom_text(aes(label=round(Percentage,1)), size = 6/.pt, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_Stage_onlyAll_barplot.pdf"),3,2.5,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()

	tc = Clinall[!is.na(Clinall[,"Race"]),]
	# tc = tc[!( tc$Race %in% c( "native american","pacific islander" ) ),]
	tc$Race2 = "non-asian"
	tc$Race2[(tc$Race=="asian") %in% c(T)] = "asian"
	tc$Race = tc$Race2
	tabb = table(tc[,"Race"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race","EGFR_class","Percentage" )
	mt$Race = factor( mt$Race, levels = c( "asian","non-asian" ) )
	mt$TotalPatients = NA
	mt[(mt$Race=="asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="asian"))
	mt[(mt$Race=="non-asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="non-asian"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=Race, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = 6/.pt) + geom_text(aes(label=round(Percentage,1)), size = 6/.pt, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_Race_onlyAll_barplot.pdf"),2.8,2.8,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	print(pa)
	dev.off()

	tc = Clinall[!is.na(Clinall[,"Smoking"]),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = 6/.pt) + geom_text(aes(label=round(Percentage,1)), size = 6/.pt, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_Smoking_onlyAll_barplot.pdf"),2.8,2.8,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()

	tc = Clinall[!is.na(Clinall[,"Age"]),]
	tc$Age = as.numeric(tc$Age)
	fileName = paste0(OutDir,"EGFRclasses_Age_onlyAll_boxplot.pdf")
	tc$egfr_class_consensus = factor(tc$egfr_class_consensus,levels=ordered_classes)
	jitterColors = dan.expand_colors( tc$egfr_class_consensus, ordered_classes, colorz_classes )
	dan.boxplots(fileName,x = tc$egfr_class_consensus, y=tc$Age, xlab="", ylab="Age (years)", labelycoo=max(tc$Age), xColors=colorz_classes, jitterColors=jitterColors,jitterDotSize=1, fileWidth = 2.5, fileHeight = 2.5)
 


	kruskal.test(Age~egfr_class_consensus,data=Clinall)
	kruskal.test(Age~egfr_class_consensus,data=Cling)
	kruskal.test(Age~egfr_class_consensus,data=Clinc)
	kruskal.test(Age~egfr_class_consensus,data=Clint)
	kruskal.test(Age~egfr_class_consensus,data=Clinx)
	kruskal.test(Age~egfr_class_consensus,data=Clino)

	# Sex
	tc = Cling[!is.na(Cling[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pg = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Genie\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clint[!is.na(Clint[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pt = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TCGA\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinc[!is.na(Clinc[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pc = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("ChenEAS\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinx[!is.na(Clinx[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	px = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TRACERx421\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinz[!is.na(Clinz[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pz = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Zhang\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clino[!is.na(Clino[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	po = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Origimed\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinall[!is.na(Clinall[,"Sex"]),]
	tabb = table(tc[,"Sex"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Sex=="female") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="female"))
	mt[(mt$Sex=="male") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Sex=="male"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	pdf(paste0(OutDir,"EGFRclasses_Sex_All_barplot.pdf"),21,5,onefile=FALSE)
	plot = ggarrange(pg, pt, pc, px,pz,po, pa, common.legend=TRUE,ncol = 7, nrow = 1)
	print(plot)
	dev.off()

	# Race
	tc = Cling[!is.na(Cling[,"Race"]),]
	tc = tc[!( tc$Race %in% c( "native american","pacific islander" ) ),]
	tabb = table(tc[,"Race"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race","EGFR_class","Percentage" )
	mt$Race = factor( mt$Race, levels = c( "asian","black","white" ) )
	mt$TotalPatients = NA
	mt[(mt$Race=="white") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="white"))
	mt[(mt$Race=="black") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="black"))
	mt[(mt$Race=="asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="asian"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pg = ggplot(data=mt, aes(x=Race, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Genie\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clint[!is.na(Clint[,"Race"]),]
	tabb = table(tc[,"Race"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race","EGFR_class","Percentage" )
	mt$Race = factor( mt$Race, levels = c( "asian","black","white" ) )
	mt$TotalPatients = NA
	mt[(mt$Race=="white") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="white"))
	mt[(mt$Race=="black") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="black"))
	mt[(mt$Race=="asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="asian"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pt = ggplot(data=mt, aes(x=Race, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TCGA\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinc[!is.na(Clinc[,"Race"]),]
	tabb = table(tc[,"Race"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race","EGFR_class","Percentage" )
	mt$Race = factor( mt$Race, levels = c( "asian","black","white" ) )
	mt$TotalPatients = NA
	mt[(mt$Race=="white") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="white"))
	mt[(mt$Race=="black") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="black"))
	mt[(mt$Race=="asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="asian"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pc = ggplot(data=mt, aes(x=Race, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("ChenEAS\n" )) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinx[!is.na(Clinx[,"Race"]),]
	tabb = table(tc[,"Race"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race","EGFR_class","Percentage" )
	mt$Race = factor( mt$Race, levels = c( "asian","black","white" ) )
	mt$TotalPatients = NA
	mt[(mt$Race=="white") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="white"))
	mt[(mt$Race=="black") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="black"))
	mt[(mt$Race=="asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="asian"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	px = ggplot(data=mt, aes(x=Race, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TRACERx421\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clino[!is.na(Clino[,"Race"]),]
	tabb = table(tc[,"Race"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race","EGFR_class","Percentage" )
	mt$Race = factor( mt$Race, levels = c( "asian","black","white" ) )
	mt$TotalPatients = NA
	mt[(mt$Race=="white") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="white"))
	mt[(mt$Race=="black") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="black"))
	mt[(mt$Race=="asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="asian"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	po = ggplot(data=mt, aes(x=Race, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Origimed\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinall[!is.na(Clinall[,"Race"]),]
	tc = tc[!( tc$Race %in% c( "native american","pacific islander" ) ),]
	tabb = table(tc[,"Race"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race","EGFR_class","Percentage" )
	mt$Race = factor( mt$Race, levels = c( "asian","black","white" ) )
	mt$TotalPatients = NA
	mt[(mt$Race=="white") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="white"))
	mt[(mt$Race=="black") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="black"))
	mt[(mt$Race=="asian") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Race=="asian"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=Race, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	pdf(paste0(OutDir,"EGFRclasses_Race_All_barplot.pdf"),18,5,onefile=FALSE)
	plot=ggarrange(pg, pt, pc, px, po, pa, common.legend=TRUE,ncol = 6, nrow = 1)
	print(plot)
	dev.off()

	# Smoking
	tc = Cling[!is.na(Cling[,"Smoking"]),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pg = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Genie\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clint[!is.na(Clint[,"Smoking"]),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pt = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TCGA\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinc[!is.na(Clinc[,"Smoking"]),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pc = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("ChenEAS\n","Chi-square test, p = ",signif(ch,2))) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinx[!is.na(Clinx[,"Smoking"]),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	px = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TRACERx421\n","Chi-square test, p = ",signif(ch,2))) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clino[!is.na(Clino[,"Smoking"]),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	po = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Origimed\n","Chi-square test, p = ",signif(ch,2))) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinz[!is.na(Clinz[,"Smoking"]),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pz = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Zhang\n")) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	tc = Clinall[!is.na(Clinall[,"Smoking"]),]
	tc = tc[!( tc$Smoking %in% c( "native american","pacific islander" ) ),]
	tabb = table(tc[,"Smoking"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Smoking=="ever smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="ever smoker"))
	mt[(mt$Smoking=="never smoker") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Smoking=="never smoker"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))

	pdf(paste0(OutDir,"EGFRclasses_Smoking_All_barplot.pdf"),21,5,onefile=FALSE)
	plot=ggarrange(pg, pt, pc, px,pz,po, pa, common.legend=TRUE,ncol = 7, nrow = 1)
	print(plot)
	dev.off()

	## Are variables independent or not? Using Clinall only. Not really independent..
	tc = Clinall
	tc = tc[!( tc$Race %in% c( "native american","pacific islander" ) ),]
	tabb = table(tc$Sex,tc$Race)
	chisq.test(tabb)$p.value
	t(apply(tabb,1, function(x) x/sum(x)))*100

	tc = Clinall
	tc = tc[!( tc$Race %in% c( "native american","pacific islander" ) ),]
	tabb = table(tc$Smoking,tc$Race)
	chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	tabb

	tc = Clinall
	tabb = table(tc$Smoking,tc$Sex)
	chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	tabb

	# Let's try to account for variables
	tc = Clinall
	# tc = tc[!( tc$Race %in% c( "native american","pacific islander" ) ),]
	tc$Race2 = "non-asian"
	tc$Race2[(tc$Race=="asian") %in% c(T)] = "asian"
	tc$Race2[is.na(tc$Race)] = NA
	dtable(tc$Race2)
	tabb = table(tc$Smoking,tc$egfr_class_consensus,tc$Sex,tc$Race2)
	tabb
	# test if sex difference depends on race difference
	tc1 = tc[(tc$Race2=="asian") %in% c(T),]
	tabb1= table(tc1[,"Sex"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pAsian = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Asians (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pAsian = pAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	tc1 = tc[(tc$Race2=="non-asian") %in% c(T),]
	tabb1= table(tc1[,"Sex"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pNonAsian = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Non-asians (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pNonAsian = pNonAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	pdf(paste0(OutDir,"EGFRclasses_Sex_SplittedByRace_All_barplot.pdf"),4,3,onefile=FALSE)
	plot=ggarrange(pAsian, pNonAsian, common.legend=TRUE,ncol = 2, nrow = 1)
	print(plot)
	dev.off()
	# test if sex difference depends on smoking difference
	tc1 = tc[(tc$Smoking=="ever smoker") %in% c(T),]
	tabb1= table(tc1[,"Sex"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pAsian = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Ever smokers (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	tc1 = tc[(tc$Smoking=="never smoker") %in% c(T),]
	tabb1= table(tc1[,"Sex"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Sex","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pNonAsian = ggplot(data=mt, aes(x=Sex, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Never smokers (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_Sex_SplittedBySmoking_All_barplot.pdf"),7,5,onefile=FALSE)
	plot=ggarrange(pAsian, pNonAsian, common.legend=TRUE,ncol = 2, nrow = 1)
	print(plot)
	dev.off()

	# test if smoking difference depends on race difference
	tc1 = tc[(tc$Race2=="asian") %in% c(T),]
	tabb1= table(tc1[,"Smoking"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pAsian = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Asians (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pAsian = pAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	tc1 = tc[(tc$Race2=="non-asian") %in% c(T),]
	tabb1= table(tc1[,"Smoking"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pNonAsian = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Non-asians (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pNonAsian = pNonAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	pdf(paste0(OutDir,"EGFRclasses_Smoking_SplittedByRace_All_barplot.pdf"),4,3,onefile=FALSE)
	plot=ggarrange(pAsian, pNonAsian, common.legend=TRUE,ncol = 2, nrow = 1)
	print(plot)
	dev.off()
	# test if smoking difference depends on sex difference
	tc1 = tc[(tc$Sex=="female") %in% c(T),]
	tabb1= table(tc1[,"Smoking"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pAsian = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Females (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pAsian = pAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	tc1 = tc[(tc$Sex=="male") %in% c(T),]
	tabb1= table(tc1[,"Smoking"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Smoking","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pNonAsian = ggplot(data=mt, aes(x=Smoking, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Males (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pNonAsian = pNonAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	pdf(paste0(OutDir,"EGFRclasses_Smoking_SplittedBySex_All_barplot.pdf"),4,3,onefile=FALSE)
	plot=ggarrange(pAsian, pNonAsian, common.legend=TRUE,ncol = 2, nrow = 1)
	print(plot)
	dev.off()

	# test if race difference depends on smoking difference
	tc1 = tc[(tc$Smoking=="ever smoker") %in% c(T),]
	tabb1= table(tc1[,"Race2"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race2","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pAsian = ggplot(data=mt, aes(x=Race2, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Ever smokers (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pAsian = pAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	tc1 = tc[(tc$Smoking=="never smoker") %in% c(T),]
	tabb1= table(tc1[,"Race2"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race2","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pNonAsian = ggplot(data=mt, aes(x=Race2, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Never smokers (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pNonAsian = pNonAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	pdf(paste0(OutDir,"EGFRclasses_Race_SplittedBySmoking_All_barplot.pdf"),4,3,onefile=FALSE)
	plot=ggarrange(pAsian, pNonAsian, common.legend=TRUE,ncol = 2, nrow = 1)
	print(plot)
	dev.off()
	# test if race difference depends on sex difference
	tc1 = tc[(tc$Sex=="female") %in% c(T),]
	tabb1= table(tc1[,"Race2"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race2","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pAsian = ggplot(data=mt, aes(x=Race2, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Females (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pAsian = pAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	tc1 = tc[(tc$Sex=="male") %in% c(T),]
	tabb1= table(tc1[,"Race2"],tc1$egfr_class_consensus)
	ch = chisq.test(tabb1)$p.value
	tabb = t(apply(tabb1,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Race2","EGFR_class","Percentage" )
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pNonAsian = ggplot(data=mt, aes(x=Race2, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Males (N = ",sum(tabb1),")\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pNonAsian = pNonAsian + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	pdf(paste0(OutDir,"EGFRclasses_Race_SplittedBySex_All_barplot.pdf"),4,3,onefile=FALSE)
	plot=ggarrange(pAsian, pNonAsian, common.legend=TRUE,ncol = 2, nrow = 1)
	print(plot)
	dev.off()

	## Ok, it seems we can drop Sex as covariate. Let's focus on Race2 and Smoking
	## Multivariate. see https://stats.stackexchange.com/questions/147925/chi2-of-multidimensional-data and https://online.stat.psu.edu/stat504/book/export/html/720

	#   Cochran-Mantel-Haenszel tests

	tc = Clinall
	tc$Race2 = "non-asian"
	tc$Race2[(tc$Race=="asian") %in% c(T)] = "asian"
	tc$Race2[is.na(tc$Race)] = NA
	dtable(tc$Race2)

	print( "Association with Smoking accounting for Race" )
	tabb = table(tc$egfr_class_consensus,tc$Smoking,tc$Race2)
	cmh = mantelhaen.test(tabb) # association between Smoking and egfr_class while controlling for Race2
	print(cmh)
	tc = Clinall
	print( "Association with Smoking accounting for Sex" )
	tabb = table(tc$egfr_class_consensus,tc$Smoking,tc$Sex)
	cmh = mantelhaen.test(tabb) # association between Smoking and egfr_class while controlling for Sex
	print(cmh)

	tc = Clinall
	tc$Race2 = "non-asian"
	tc$Race2[(tc$Race=="asian") %in% c(T)] = "asian"
	tc$Race2[is.na(tc$Race)] = NA
	print( "Association with Race accounting for Smoking" )
	tabb = table(tc$egfr_class_consensus,tc$Race2,tc$Smoking)
	cmh = mantelhaen.test(tabb) # association between Race2 and egfr_class while controlling for Smoking
	print(cmh)

	print( "Association with Race accounting for Sex" )
	tabb = table(tc$egfr_class_consensus,tc$Race2,tc$Sex)
	cmh = mantelhaen.test(tabb) # association between race and egfr_class while controlling for Sex
	print(cmh)

	print( "Association with Sex accounting for Race" )
	tabb = table(tc$egfr_class_consensus,tc$Sex,tc$Race2)
	cmh = mantelhaen.test(tabb) # association between race and egfr_class while controlling for Sex
	print(cmh)

	tc = Clinall
	print( "Association with Sex accounting for Smoking" )
	tabb = table(tc$egfr_class_consensus,tc$Sex,tc$Smoking)
	cmh = mantelhaen.test(tabb) # association between race and egfr_class while controlling for Sex
	print(cmh)

	################ other vars: has_metastasis, where_metastasis, age

	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clint = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TCGA_SampleLevel.RData" ))
	Clin2t = Clin2
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clinc = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Chen_SampleLevel.RData" ))
	Clin2c = Clin2

	## has_metastasis
	met_tcga = dan.read("/mnt/ndata/daniele/various/Processed/Arvind/TCGA_allTT_MetastasisPresence.txt")
	rownames(met_tcga) = met_tcga$Patient
	commonz = intersect(rownames(met_tcga),rownames(Clint))
	Clint = cbind(Clint[commonz,],met_tcga[commonz,c( "Lymphnode_metastasis","OtherOrgans_metastasis","Metastasis_presence" )])
	# overall metastasis presence
	tc = Clint
	tabb = table(tc[,"Metastasis_presence"],tc$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Metastasis_presence","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$Metastasis_presence=="yes") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Metastasis_presence=="yes",na.rm=T))
	mt[(mt$Metastasis_presence=="no") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc$Metastasis_presence=="no",na.rm=T))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pdf(paste0(OutDir,"EGFRclasses_MetastasisPresence_TCGA_barplot.pdf"),2.5,2.5)
	pt = ggplot(data=mt, aes(x=Metastasis_presence, y=Percentage, fill=EGFR_class)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("Presence of metastasis") + ggtitle(paste0("TCGA\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = size_labels,) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pt = pt + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt')) 
	print(pt)
	dev.off()

	################# survival
	library(survival)
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clint = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clinc = Clin

	### TCGA
	this_Clin = Clint
	a = as.numeric(as.character(this_Clin$days_to_death))
	b = as.numeric(as.character(this_Clin$days_to_last_followup))
	vital_status_num <- vector(mode="numeric", length=length(this_Clin$vital_status))
	times <- vector(mode="numeric", length=length(vital_status_num))
	for (v in 1:length(vital_status_num))
	{
	 if (this_Clin$vital_status[v]=="alive")
	 {
	    vital_status_num[v] <- 0
	    times[v] <- b[v]
	 }
	 else
	 {
	    vital_status_num[v] <- 1
	    times[v] <- a[v]
	 }
	}
	this_Clin$Times <- times
	this_Clin$vital_status_num <- vital_status_num
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "TCGA_KM_SurvivalBy_EGFRclass.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("TCGA - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()

	### Chen
	this_Clin = Clinc
	this_Clin$Times = as.numeric(this_Clin$OS.Month)
	this_Clin$vital_status_num = NA
	this_Clin[this_Clin$OS.Status=="Dead","vital_status_num"] = 1
	this_Clin[this_Clin$OS.Status=="Alive","vital_status_num"] = 0
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "ChenEAS_KM_SurvivalBy_EGFRclass.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("ChenEAS - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()

	### Genie 
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "Genie_KM_SurvivalBy_EGFRclass.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Genie - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()
}

tracerx_analyses = function( OutDir ){
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	em = dan.read(paste0(DataDir,"TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/20221110_TRACERx421_evolutionary_metrics.tsv")) # contains loooots of metrics - can be used
	rownames(em) = em$tumour_id
	varz = c( "clonal_mixing_index","num_branches","branching_index","perc_subclonal","TMB","num_subclones","wGII","wFLOH","ploidy","purity","genom_frac_clonal_event","genom_frac_subclonal_event", "genom_frac_event","SBS4_clonal","SBS4_subclonal","mean_smoking_sig","mean_apobec_sig")
	Clin = cbind(Clin,em[rownames(Clin),varz])
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	pdf(paste0(OutDir,"TRACERx421_evo_metrics.pdf"),6,5)
	for (ct in varz){
		x = factor(Clin$egfr_class_consensus,levels=ordered_classes)
		y = as.numeric(Clin[,ct])
	  	plot=dan.boxplots.multipages(x=x,y=y,signifTest="kruskal",xlab="",labelycoo = max(y,na.rm=T)+0.01,ylab=ct,xColors = colorz_classes, jitterColors = Clin$colorz,includeJitters = T)
	  	print(plot)
	}
	dev.off()
}

mutations_associations_alldatasets = function( OutDir ){

	hotspots = dan.read(paste0(DataDir,"hotspots_mskcc_chang2016.txt"))
	hotspots$mut = paste0("p.",hotspots$Reference_Amino_Acid,hotspots$Amino_Acid_Position,hotspots$Variant_Amino_Acid )
	genez = dan.read(paste0( DataDir,"mina_natgen_2020_genes.txt" ))
	rownames(genez) = genez$gene
	###### Genie
	dcat( "Processing Genie",1 )
	###### Genie
	gi = dan.read(paste0(DataDir,"genie_14.0/genomic_information.txt"))
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2
	colorz_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"classes"]	
	load(file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
	gene_universe = unique(maf_ss$Hugo_Symbol)
	save(gene_universe, file = paste0( OutDir,"Genie_gene_universe.RData" ))
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	maf_ss = maf_ss[maf_ss$Sample %in% Clin2$Sample,]
	maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	these_patients = unique(maf_ss$Patient)
	maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation","Nonsense_Mutation","Splice_Site")),]

	tmaf_kras = dan.df(0,c( "Patient","KRAS_variant","KRAS_variant_site","KRAS_vaf" ))
	
	tClin2 = Clin2g
	tclin = Cling
	tclin = Cling[rownames(Cling) %in% tClin2$Patient,]
	tmaf = maf_ss[(maf_ss$Patient %in% tclin$Patient),]
	gam = dan.df(these_patients,rownames(genez), data=NA)
	assays = sort(unique(Clin2g$SEQ_ASSAY_ID))
	for (a in assays){
		gia = gi[gi$SEQ_ASSAY_ID==a,]
		ug = intersect(unique(gia$Hugo_Symbol[nchar(gia$Hugo_Symbol)>0]),colnames(gam))
		patients = intersect(unique(tClin2[Clin2g$SEQ_ASSAY_ID==a,"Patient"]),rownames(gam))
		gam[patients,ug] = 0
	}
	for (cn in colnames(gam)){
		if (cn=="KRAS"){ 
			ttmaf = tmaf[(tmaf$Hugo_Symbol==cn) & ((tmaf$Variant_Classification=="Missense_Mutation") & (tmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"])),]
			if (nrow(ttmaf)>0) {tmaf_kras = rbind(tmaf_kras,data.frame(Patient=ttmaf$Patient,KRAS_variant=ttmaf$HGVSp_Short,KRAS_variant_site=paste0(substr(ttmaf$HGVSp_Short,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$HGVSp_Short)),KRAS_vaf=as.numeric(ttmaf$t_alt_count)/as.numeric(ttmaf$t_depth),stringsAsFactors=F))}
		} else {
			these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
			if (cn %in% rownames(genez)){
				if (genez[cn,"role"]=="og"){
				these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
				} else {
					these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
				}
			}
			ttmaf = tmaf[(tmaf$Variant_Classification %in% these_var_class) & (tmaf$Hugo_Symbol==cn),]
			ttmaf = ttmaf[(ttmaf$Variant_Classification!="Missense_Mutation") | (ttmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"]),]
		}
		if (nrow(ttmaf)==0){ next }
		gam[unique(ttmaf$Patient),cn] = 1
	}
	tclin$Dataset = "Genie"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = gam
	clinAll = tclin

	dcat( "Processing TCGA",1 )
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clint = Clin
	load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
	gene_universe = unique(maf_ss$Hugo_Symbol)
	save(gene_universe, file = paste0( OutDir,"TCGA_gene_universe.RData" ))
	maf_ss$Patient = substr(maf_ss$Tumor_Sample_Barcode,1,12)
	maf_ss$Sample = substr(maf_ss$Tumor_Sample_Barcode,1,16)
	tclin = Clint[Clint$Patient %in% maf_ss$Patient,]
	maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	these_patients = unique(maf_ss$Patient)
	maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation","Nonsense_Mutation","Splice_Site"),]
	tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	for (cn in colnames(gam)){
		if (cn=="KRAS"){ 
			ttmaf = tmaf[(tmaf$Hugo_Symbol==cn) & ((tmaf$Variant_Classification=="Missense_Mutation") & (tmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"])),]
			if (nrow(ttmaf)>0) {tmaf_kras = rbind(tmaf_kras,data.frame(Patient=ttmaf$Patient,KRAS_variant=ttmaf$HGVSp_Short,KRAS_variant_site=paste0(substr(ttmaf$HGVSp_Short,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$HGVSp_Short)),KRAS_vaf=as.numeric(ttmaf$t_alt_count)/as.numeric(ttmaf$t_depth),stringsAsFactors=F))}
		} else {
			these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
			if (cn %in% rownames(genez)){
				if (genez[cn,"role"]=="og"){
				these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
				} else {
					these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
				}
			}
			ttmaf = tmaf[(tmaf$Variant_Classification %in% these_var_class) & (tmaf$Hugo_Symbol==cn),]
			ttmaf = ttmaf[(ttmaf$Variant_Classification!="Missense_Mutation") | (ttmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"]),]
		}
		if (nrow(ttmaf)==0){ next }
		gam[unique(ttmaf$Patient),cn] = 1
	}
	tclin$Dataset = "TCGA"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	dcat( "Processing ChenEAS",1 )
	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clinc = Clin
	maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
	gene_universe = unique(maf$Hugo_Symbol)
	save(gene_universe, file = paste0( OutDir,"ChenEAS_gene_universe.RData" ))
	maf$Patient = maf$Tumor_Sample_Barcode
	maf_ss = maf
	tclin = Clinc[Clinc$Patient %in% maf_ss$Patient,]
	maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	these_patients = unique(maf_ss$Patient)
	maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins", "missense","nonsense_mutation","Splice_Site")),]
	tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	for (cn in colnames(gam)){
		if (cn=="KRAS"){ 
			ttmaf = tmaf[(tmaf$Hugo_Symbol==cn) & ((tmaf$Variant_Classification=="missense") & (tmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"])),]
			if (nrow(ttmaf)>0) {tmaf_kras = rbind(tmaf_kras,data.frame(Patient=ttmaf$Patient,KRAS_variant=ttmaf$HGVSp_Short,KRAS_variant_site=paste0(substr(ttmaf$HGVSp_Short,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$HGVSp_Short)),KRAS_vaf=as.numeric(ttmaf$t_alt_count)/(as.numeric(ttmaf$t_ref_count)+as.numeric(ttmaf$t_alt_count)),stringsAsFactors=F))}
			} else {
			these_var_class = c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins", "missense","nonsense_mutation","Splice_Site")
			if (cn %in% rownames(genez)){
				if (genez[cn,"role"]=="og"){
				these_var_class = c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins", "missense")
				} else {
					these_var_class = c("frame_shift_del","frame_shift_ins", "in_frame_del", "in_frame_ins","nonsense_mutation","Splice_Site","missense")
				}
			}
			ttmaf = tmaf[(tmaf$Variant_Classification %in% these_var_class) & (tmaf$Hugo_Symbol==cn),]
			ttmaf = ttmaf[(ttmaf$Variant_Classification!="missense") | (ttmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"]),]
		}
		if (nrow(ttmaf)==0){ next }
		gam[unique(ttmaf$Patient),cn] = 1
	}
	tclin$Dataset = "ChenEAS"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	dcat( "Processing Zhang",1 )
	###### Zhang
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Zhang_SampleLevel.RData" ))
	Clinz = Clin
	Clin2z = Clin2
	maf = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_mutations.txt"))
	gene_universe = unique(maf$Hugo_Symbol)
	save(gene_universe, file = paste0( OutDir,"Zhang_gene_universe.RData" ))
	maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	these_patients = unique(maf_ss$Patient)
	maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation","Nonsense_Mutation","Splice_Site"),]
	tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	for (cn in colnames(gam)){
		if (cn=="KRAS"){ 
			ttmaf = tmaf[(tmaf$Hugo_Symbol==cn) & ((tmaf$Variant_Classification=="Missense_Mutation") & (tmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"])),]
			if (nrow(ttmaf)>0) {tmaf_kras = rbind(tmaf_kras,data.frame(Patient=ttmaf$Patient,KRAS_variant=ttmaf$HGVSp_Short,KRAS_variant_site=paste0(substr(ttmaf$HGVSp_Short,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$HGVSp_Short)),KRAS_vaf=NA,stringsAsFactors=F))}
			} else {
			these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
			if (cn %in% rownames(genez)){
				if (genez[cn,"role"]=="og"){
				these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
				} else {
					these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
				}
			}
			ttmaf = tmaf[(tmaf$Variant_Classification %in% these_var_class) & (tmaf$Hugo_Symbol==cn),]
			ttmaf = ttmaf[(ttmaf$Variant_Classification!="Missense_Mutation") | (ttmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"]),]
		}
		if (nrow(ttmaf)==0){ next }
		gam[unique(ttmaf$Patient),cn] = 1
	}
	tclin$Dataset = "Zhang"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	dcat( "Processing TRACERx421",1 )
	###### TRACERx421
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TRACERx421_SampleLevel.RData" ))
	Clinz = Clin
	Clin2z = Clin2
	library(fst)
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	gene_universe = unique(maf$Hugo_Symbol)
	save(gene_universe, file = paste0( OutDir,"TRACERx421_gene_universe.RData" ))
	maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$tumour_id
	maf_ss$Patient = maf_ss$patient_id
	mafe = maf_ss[!(maf_ss$Hugo_Symbol=="EGFR"),]
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	these_patients = unique(maf_ss$Patient)
	maf_ss$mut = maf_ss$AAChange
	maf_ss[is.na(maf_ss$mut),"mut"] = paste0("p.",maf_ss[is.na(maf_ss$mut),"ref"],maf_ss[is.na(maf_ss$mut),"start"],maf_ss[is.na(maf_ss$mut),"var"])
	tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$mut)),]
	gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	for (cn in colnames(gam)){
		ttmaf = tmaf[(tmaf$Hugo_Symbol==cn) & (tmaf$DriverMut),]
		if (nrow(ttmaf)==0){ next }
		gam[unique(ttmaf$Patient),cn] = 1
		if (cn=="KRAS"){ 
			tmaf_kras = rbind(tmaf_kras,data.frame(Patient=ttmaf$Patient,KRAS_variant=ttmaf$mut,KRAS_variant_site=paste0(substr(ttmaf$mut,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$mut)),KRAS_vaf=NA,stringsAsFactors=F))
		}
	}
	tclin$Dataset = "TRACERx421"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	###### Origimed
	dcat( "Processing Origimed",1 )
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Origimed_SampleLevel.RData" ))
	Clinz = Clin
	Clin2z = Clin2
	colorz_classes = clco[rownames(clco) %in% (Clin2z$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin2z$egfr_class_consensus),"classes"]
	maf = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_mutations_extended.txt"))
	gene_universe = unique(maf$Hugo_Symbol)
	save(gene_universe, file = paste0( OutDir,"Origimed_gene_universe.RData" ))
	maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	rownames(Clin2) = Clin2$Sample
	maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	these_patients = unique(maf_ss$Patient)
	maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation","Nonsense_Mutation","Splice_Site"),]
	tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	for (cn in colnames(gam)){
		if (cn=="KRAS"){ 
			ttmaf = tmaf[(tmaf$Hugo_Symbol==cn) & ((tmaf$Variant_Classification=="Missense_Mutation") & (tmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"])),]
			if (nrow(ttmaf)>0) {tmaf_kras = rbind(tmaf_kras,data.frame(Patient=ttmaf$Patient,KRAS_variant=ttmaf$HGVSp_Short,KRAS_variant_site=paste0(substr(ttmaf$HGVSp_Short,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$HGVSp_Short)),KRAS_vaf=NA,stringsAsFactors=F))}
			} else {
			these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
			if (cn %in% rownames(genez)){
				if (genez[cn,"role"]=="og"){
				these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation")
				} else {
					these_var_class = c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins","Nonsense_Mutation","Splice_Site","Missense_Mutation")
				}
			}
			ttmaf = tmaf[(tmaf$Variant_Classification %in% these_var_class) & (tmaf$Hugo_Symbol==cn),]
			ttmaf = ttmaf[(ttmaf$Variant_Classification!="Missense_Mutation") | (ttmaf$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol==cn,"mut"]),]

		}
		if (nrow(ttmaf)==0){ next }
		gam[unique(ttmaf$Patient),cn] = 1
	}
	tclin$Dataset = "Origimed"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	tclin = clinAll
	save(gamAll,file=paste0( OutDir,"gamAll_unfiltered.RData" ))
	save(tclin,file=paste0( OutDir,"tclin_gamAll_unfiltered.RData" ))

	gamAll = gamAll[,colnames(gamAll)!="EGFR"]
	gamAll = gamAll[rownames(tclin),]
	all(rownames(gamAll) %in% rownames(tclin))
	cols_common = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="common",]),],na.rm=T)
	cols_uncommon = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="uncommon",]),],na.rm=T)
	cols_compound = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="compound",]),],na.rm=T)
	gamAll = gamAll[,((cols_common>=3) | (cols_uncommon>=3)) | (cols_compound>=3)]
	tclin = tclin[rownames(gamAll),]
	cdf = dan.df(0,c( "Gene","EGFR_class","nMut","Percentage","Chisq_pval","winner" ))
	for (g in colnames(gamAll)){
		this_tclin = tclin
		this_tclin$this_gene = gamAll[,g]
		this_tclin = this_tclin[!is.na(this_tclin$this_gene),]
		this_tclin[this_tclin$this_gene==0,"this_gene"] = "wt"
		this_tclin[this_tclin$this_gene==1,"this_gene"] = "mut"
		this_tclin$this_gene = factor(this_tclin$this_gene,levels=c("mut","wt"))
		this_tclin$egfr_class_consensus = factor(this_tclin$egfr_class_consensus,levels=ordered_classes)
		tabb = table(this_tclin$this_gene,this_tclin$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tcdf = data.frame(row.names = ordered_classes, Gene=rep(g,length(ordered_classes) ),EGFR_class=ordered_classes,Percentage=NA,Chisq_pval=NA,winner="no",stringsAsFactors=F)
		for (cl in ordered_classes){ tcdf[cl,"nMut"] = tabb["mut",cl] }
		for (cl in ordered_classes){ tcdf[cl,"Percentage"] = 100*tabb["mut",cl]/sum(tabb[,cl]) }
		tcdf[which.max( tcdf$Percentage),"winner"] = "yes"
		tcdf[tcdf$winner=="yes","Chisq_pval"] = ch
		cdf = rbind(cdf,tcdf)
	}
	cdf$Chisq_qval = NA
	cdf[cdf$winner=="yes","Chisq_qval"] = p.adjust(cdf[cdf$winner=="yes","Chisq_pval"],method="BH")
	cdf$EGFR_class = factor(cdf$EGFR_class, levels = ordered_classes)
	aa=cdf[order(cdf$Chisq_pval),]
	cdf_signif_genes = unique(cdf[(cdf$Chisq_qval<0.1) %in% c(T),"Gene"])
	cdf = cdf[cdf$Gene %in% cdf_signif_genes,]
	save(cdf,file=paste0( OutDir,"AllDatasets_signifComutated_table.RData" ))
	dcat( paste0("Uncommon percentage: ", 100*sum(cdf[(cdf$winner=="yes"),"EGFR_class"]=="uncommon")/nrow(cdf[(cdf$winner=="yes"),]) ))

	tcdf = cdf[cdf$Gene=="TP53",]
	tcdf$Chisq_pval = signif(tcdf$Chisq_pval,2)
	tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)] = paste0("p-value = ",tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)])
	tcdf$Gene = "TP53 mutations\n(putative loss-of-function)"
	pdf(paste0( OutDir,"EGFRclasses_vs_TP53mutant_AllDatasets.pdf" ),2.5,2.5)
	pl = ggplot(data=tcdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, mutant patients = ",sum(tcdf$nMut)," out of ",sum(!is.na(gamAll$TP53))) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 0,size=12)) + geom_text(aes(label=Chisq_pval,y=Percentage),vjust=-0.2,size = size_labels) + geom_text(aes(label=round(Percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	pl = pl + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pl)
	dev.off()

	tcdf = cdf[cdf$Gene=="KRAS",]
	tcdf$Chisq_pval = signif(tcdf$Chisq_pval,2)
	tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)] = paste0("p-value = ",tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)])
	tcdf$Gene = "KRAS mutations\n(missense hotspots)"
	pdf(paste0( OutDir,"EGFRclasses_vs_KRASmutant_AllDatasets.pdf" ),2.5,2.5)
	pl = ggplot(data=tcdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, mutant patients = ",sum(tcdf$nMut)," out of ",sum(!is.na(gamAll$KRAS))) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 0,size=12)) + geom_text(aes(label=Chisq_pval,y=Percentage),vjust=-0.2,size = size_labels) + geom_text(aes(label=round(Percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	pl = pl + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pl)
	dev.off()

	cdf = cdf[order(cdf$Chisq_pval),]
	genes_keep = cdf[1:20,"Gene"]
	cdf = cdf[cdf$Gene %in% genes_keep,]
	cdf$Gene = factor(cdf$Gene, levels = unique(cdf$Gene))
	cdf$Chisq_pval = signif(cdf$Chisq_pval,2)
	cdf[order(cdf$Gene),]

	pdf(paste0( OutDir,"AllDatasets_comutated_barplot_top20.pdf" ),8,2.5)
	pl = ggplot(data=cdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Co-mutated genes, Chi-square test adjusted p-values (top 20 shown)") ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=Chisq_pval,y=Percentage),vjust=-0.2,size = size_labels) + geom_text(aes(label=round(Percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	pl = pl + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pl)
	dev.off()

	save(gamAll,file=paste0( OutDir,"gamAll.RData" ))
	save(tclin,file=paste0( OutDir,"tclin_gamAll.RData" ))
	save(tmaf_kras,file=paste0( OutDir,"tmaf_kras.RData" ))
}

mutations_associations_accounting_for_TMB = function( OutDir, nbins = 8 ){
	# # computing coverage for each Genie assay (bedtools)
	# gi = dan.read(paste0(DataDir,"genie_14.0/genomic_information.txt"))
	# dir.create(paste0(DataDir,"genie_14.0/assays_coverages/"))
	# for (assay in unique(gi$SEQ_ASSAY_ID)){
	# 	gia = gi[ gi$SEQ_ASSAY_ID==assay,]
	# 	# removing duplicated entries
	# 	gia = gia[!duplicated( paste0(gia$Chromosome,"_",gia$Start_Position,"_",gia$End_Position )),c("Chromosome","Start_Position","End_Position")]
	# 	gia = gia[order(gia$Chromosome,gia$Start_Position),]
	# 	write.table(gia, file = paste0(DataDir,"genie_14.0/assays_coverages/assay_",assay,".bed"), row.names = F, col.names = F, sep = "\t", quote = F)
	# }
	# ### execute in bash
	# # cd "/mnt/ndata/daniele/alfredo_egfr/Data/genie_14.0/assays_coverages/"
	# # for file in assay_*.bed
	# # do
	# #     /mnt/ndata/daniele/lung_multiregion/Scripts/Tools/bedtools2/bin/mergeBed -i $file > $file"_merged.txt"
	# # done
	# assay_coverages = data.frame(row.names=unique(gi$SEQ_ASSAY_ID), assay=unique(gi$SEQ_ASSAY_ID), coverage = NA, stringsAsFactors=F)
	# for (assay in unique(gi$SEQ_ASSAY_ID)){
	# 	cov = read.table(file = paste0(DataDir,"genie_14.0/assays_coverages/assay_",assay,".bed_merged.txt"), header = F, sep = "\t", quote = '',stringsAsFactors=F)
	# 	assay_coverages[assay,"coverage"] = sum(cov$V3-cov$V2+1)
	# }
	# assay_coverages[order(assay_coverages$coverage,decreasing=T),]
	# # the "WAKE" assays have gene intervals of 1bp. excluding them
	# assay_coverages[assay_coverages$coverage<10000,"coverage"] = NA
	# # estimate TMB for Genie from Nmut and assay coverages. Compare with the reported one (available for much fewer cases)
	# ### I tried, but TMB estimated in this way varies too much. It is possible that the assay coverages are not real.

	# ### Doing the same analysis but on synonymous variants
	# genez = dan.read(paste0( DataDir,"mina_natgen_2020_genes.txt" ))
	# rownames(genez) = genez$gene
	# ###### Genie
	# dcat( "Processing Genie",1 )
	# ###### Genie
	# gi = dan.read(paste0(DataDir,"genie_14.0/genomic_information.txt"))
	# load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	# Cling = Clin
	# load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	# Clin2g = Clin2
	# colorz_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"colorz"]
	# ordered_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"classes"]	
	# load(file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
	# gene_universe = unique(maf_ss$Hugo_Symbol)
	# save(gene_universe, file = paste0( OutDir,"Genie_gene_universe.RData" ))
	# maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	# maf_ss = maf_ss[maf_ss$Sample %in% Clin2$Sample,]
	# maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	# these_patients = unique(maf_ss$Patient)
	# maf_ss = maf_ss[(maf_ss$Variant_Classification %in% c("Silent","3'Flank","3'UTR","5'Flank","5'UTR","Intron")),]
	# tClin2 = Clin2g
	# tclin = Cling
	# tclin = Cling[rownames(Cling) %in% tClin2$Patient,]
	# tmaf = maf_ss[(maf_ss$Patient %in% tclin$Patient),]
	# gam = dan.df(these_patients,rownames(genez), data=NA)
	# assays = sort(unique(Clin2g$SEQ_ASSAY_ID))
	# for (a in assays){
	# 	gia = gi[gi$SEQ_ASSAY_ID==a,]
	# 	ug = intersect(unique(gia$Hugo_Symbol[nchar(gia$Hugo_Symbol)>0]),colnames(gam))
	# 	patients = intersect(unique(tClin2[Clin2g$SEQ_ASSAY_ID==a,"Patient"]),rownames(gam))
	# 	gam[patients,ug] = 0
	# }
	# for (cn in colnames(gam)){
	# 	ttmaf = tmaf[(tmaf$Hugo_Symbol==cn),]
	# 	if (nrow(ttmaf)==0){ next }
	# 	gam[unique(ttmaf$Patient),cn] = 1
	# }
	# tclin$Dataset = "Genie"
	# gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	# gamAll = gam
	# clinAll = tclin

	# dcat( "Processing TCGA",1 )
	# ###### TCGA
	# load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	# Clint = Clin
	# load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
	# gene_universe = unique(maf_ss$Hugo_Symbol)
	# save(gene_universe, file = paste0( OutDir,"TCGA_gene_universe.RData" ))
	# maf_ss$Patient = substr(maf_ss$Tumor_Sample_Barcode,1,12)
	# maf_ss$Sample = substr(maf_ss$Tumor_Sample_Barcode,1,16)
	# tclin = Clint[Clint$Patient %in% maf_ss$Patient,]
	# maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	# these_patients = unique(maf_ss$Patient)
	# maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Silent","3'Flank","3'UTR","5'Flank","5'UTR","Intron"),]
	# tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	# gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	# for (cn in colnames(gam)){
	# 	ttmaf = tmaf[(tmaf$Hugo_Symbol==cn),]
	# 	if (nrow(ttmaf)==0){ next }
	# 	gam[unique(ttmaf$Patient),cn] = 1
	# }
	# tclin$Dataset = "TCGA"
	# gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	# gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	# clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	# dcat( "Processing ChenEAS",1 )
	# ###### ChenEAS
	# load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	# Clinc = Clin
	# maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
	# gene_universe = unique(maf$Hugo_Symbol)
	# save(gene_universe, file = paste0( OutDir,"ChenEAS_gene_universe.RData" ))
	# maf$Patient = maf$Tumor_Sample_Barcode
	# maf_ss = maf
	# tclin = Clinc[Clinc$Patient %in% maf_ss$Patient,]
	# maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	# these_patients = unique(maf_ss$Patient)
	# maf_ss = maf[(maf$Variant_Classification %in% c("silent mutation","3'Flank","3'UTR","5'Flank","5'UTR","Intron")),]
	# tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	# gam = dan.df(unique(tmaf$Patient),intersect(gene_universe,rownames(genez)), data=0)
	# for (cn in colnames(gam)){
	# 	ttmaf = tmaf[(tmaf$Hugo_Symbol==cn),]
	# 	if (nrow(ttmaf)==0){ next }
	# 	gam[unique(ttmaf$Patient),cn] = 1
	# }
	# tclin$Dataset = "ChenEAS"
	# gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	# gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	# clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	# dcat( "Processing Zhang",1 )
	# ###### Zhang
	# load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	# load(file = paste0( OutDir,"../Preprocessing/","Clin2_Zhang_SampleLevel.RData" ))
	# Clinz = Clin
	# Clin2z = Clin2
	# maf = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_mutations.txt"))
	# gene_universe = unique(maf$Hugo_Symbol)
	# save(gene_universe, file = paste0( OutDir,"Zhang_gene_universe.RData" ))
	# maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
	# maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	# maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	# tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	# maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	# these_patients = unique(maf_ss$Patient)
	# maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Silent","Intron","5'UTR","3'UTR","5'Flank","3'Flank"),]
	# tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	# gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	# for (cn in colnames(gam)){
	# 	ttmaf = tmaf[(tmaf$Hugo_Symbol==cn),]
	# 	if (nrow(ttmaf)==0){ next }
	# 	gam[unique(ttmaf$Patient),cn] = 1
	# }
	# tclin$Dataset = "Zhang"
	# gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	# gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	# clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	# dcat( "Processing TRACERx421",1 )
	# ###### TRACERx421
	# load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	# load(file = paste0( OutDir,"../Preprocessing/","Clin2_TRACERx421_SampleLevel.RData" ))
	# Clinz = Clin
	# Clin2z = Clin2
	# library(fst)
	# maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	# gene_universe = unique(maf$Hugo_Symbol)
	# save(gene_universe, file = paste0( OutDir,"TRACERx421_gene_universe.RData" ))
	# maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
	# maf_ss$Sample = maf_ss$tumour_id
	# maf_ss$Patient = maf_ss$patient_id
	# tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	# maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	# these_patients = unique(maf_ss$Patient)
	# maf_ss = maf_ss[((maf_ss$exonic.func %in% c("synonymous SNV")) %in% c(T)),]
	# maf_ss$mut = maf_ss$AAChange
	# maf_ss[is.na(maf_ss$mut),"mut"] = paste0("p.",maf_ss[is.na(maf_ss$mut),"ref"],maf_ss[is.na(maf_ss$mut),"start"],maf_ss[is.na(maf_ss$mut),"var"])
	# tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$mut)),]
	# gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	# for (cn in colnames(gam)){
	# 	ttmaf = tmaf[(tmaf$Hugo_Symbol==cn),]
	# 	if (nrow(ttmaf)==0){ next }
	# 	gam[unique(ttmaf$Patient),cn] = 1
	# }
	# tclin$Dataset = "TRACERx421"
	# gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	# gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	# clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	# ###### Origimed
	# dcat( "Processing Origimed",1 )
	# load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	# load(file = paste0( OutDir,"../Preprocessing/","Clin2_Origimed_SampleLevel.RData" ))
	# Clinz = Clin
	# Clin2z = Clin2
	# colorz_classes = clco[rownames(clco) %in% (Clin2z$egfr_class_consensus),"colorz"]
	# ordered_classes = clco[rownames(clco) %in% (Clin2z$egfr_class_consensus),"classes"]
	# maf = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_mutations_extended.txt"))
	# gene_universe = unique(maf$Hugo_Symbol)
	# save(gene_universe, file = paste0( OutDir,"Origimed_gene_universe.RData" ))
	# maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
	# maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	# rownames(Clin2) = Clin2$Sample
	# maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	# tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	# maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	# these_patients = unique(maf_ss$Patient)
	# maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Silent","3'Flank","3'UTR","5'Flank","5'UTR","Intron"),]
	# tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$HGVSp_Short)),]
	# gam = dan.df(these_patients,intersect(gene_universe,rownames(genez)), data=0)
	# for (cn in colnames(gam)){
	# 	ttmaf = tmaf[(tmaf$Hugo_Symbol==cn),]
	# 	if (nrow(ttmaf)==0){ next }
	# 	gam[unique(ttmaf$Patient),cn] = 1
	# }
	# tclin$Dataset = "Origimed"
	# gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	# gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	# clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	# gamAll = gamAll[,colnames(gamAll)!="EGFR"]
	# gamAll = gamAll[rownames(clinAll),]
	# tclin = clinAll
	# all(rownames(gamAll) %in% rownames(tclin))
	# gamAll = gamAll[,colSums(gamAll,na.rm=T)>=3]
	# tclin = tclin[rownames(gamAll),]
	# cdf = dan.df(0,c( "Gene","EGFR_class","nMut","Percentage","Chisq_pval","winner" ))
	# for (g in colnames(gamAll)){
	# 	this_tclin = tclin
	# 	this_tclin$this_gene = gamAll[,g]
	# 	this_tclin = this_tclin[!is.na(this_tclin$this_gene),]
	# 	this_tclin[this_tclin$this_gene==0,"this_gene"] = "wt"
	# 	this_tclin[this_tclin$this_gene==1,"this_gene"] = "mut"
	# 	this_tclin$this_gene = factor(this_tclin$this_gene,levels=c("mut","wt"))
	# 	this_tclin$egfr_class_consensus = factor(this_tclin$egfr_class_consensus,levels=ordered_classes)
	# 	tabb = table(this_tclin$this_gene,this_tclin$egfr_class_consensus)
	# 	ch = chisq.test(tabb)$p.value
	# 	tcdf = data.frame(row.names = ordered_classes, Gene=rep(g,length(ordered_classes) ),EGFR_class=ordered_classes,Percentage=NA,Chisq_pval=NA,winner="no",stringsAsFactors=F)
	# 	for (cl in ordered_classes){ tcdf[cl,"nMut"] = tabb["mut",cl] }
	# 	for (cl in ordered_classes){ tcdf[cl,"Percentage"] = 100*tabb["mut",cl]/sum(tabb[,cl]) }
	# 	tcdf[which.max( tcdf$Percentage),"winner"] = "yes"
	# 	tcdf[tcdf$winner=="yes","Chisq_pval"] = ch
	# 	cdf = rbind(cdf,tcdf)
	# }
	# cdf$Chisq_qval = NA
	# cdf[cdf$winner=="yes","Chisq_qval"] = p.adjust(cdf[cdf$winner=="yes","Chisq_pval"],method="BH")
	# cdf$EGFR_class = factor(cdf$EGFR_class, levels = ordered_classes)
	# aa=cdf[order(cdf$Chisq_pval),]
	# cdf_signif_genes = unique(cdf[(cdf$Chisq_qval<0.1) %in% c(T),"Gene"])
	# scdf = cdf[cdf$Gene %in% cdf_signif_genes,]
	# dcat( paste0("Uncommon percentage: ", 100*sum(scdf[(scdf$winner=="yes"),"EGFR_class"]=="uncommon")/nrow(scdf[(scdf$winner=="yes"),]) ))

	# tcdf = cdf[cdf$Gene=="TP53",]
	# tcdf$Chisq_pval = signif(tcdf$Chisq_pval,2)
	# tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)] = paste0("Chi-square p-value = ",tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)])
	# pdf(paste0( OutDir,"EGFRclasses_vs_TP53mutant_AllDatasets_SynOnly.pdf" ),5,5)
	# pl = ggplot(data=tcdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\nTP53-mutant patients = ",sum(tcdf$nMut)," out of ",sum(!is.na(gamAll$TP53))) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=Chisq_pval,y=Percentage),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	# print(pl)
	# dev.off()

	# tcdf = cdf[cdf$Gene=="KRAS",]
	# tcdf$Chisq_pval = signif(tcdf$Chisq_pval,2)
	# tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)] = paste0("Chi-square p-value = ",tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)])
	# pdf(paste0( OutDir,"EGFRclasses_vs_KRASmutant_AllDatasets_SynOnly.pdf" ),5,5)
	# pl = ggplot(data=tcdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\nTP53-mutant patients = ",sum(tcdf$nMut)," out of ",sum(!is.na(gamAll$KRAS))) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=Chisq_pval,y=Percentage),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	# print(pl)
	# dev.off()

	# cdf = cdf[order(cdf$Chisq_pval),]
	# genes_keep = cdf[1:20,"Gene"]
	# cdf = cdf[cdf$Gene %in% genes_keep,]
	# cdf$Gene = factor(cdf$Gene, levels = unique(cdf$Gene))
	# cdf$Chisq_pval = signif(cdf$Chisq_pval,2)
	# cdf[order(cdf$Gene),]

	# pdf(paste0( OutDir,"AllDatasets_comutated_barplot_top20_SynOnly.pdf" ),nrow(cdf)/3,6)
	# pl = ggplot(data=cdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\n","Co-mutated genes, Chi-square test adjusted p-values (top 20 shown)") ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=Chisq_pval,y=Percentage),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	# print(pl)
	# dev.off()

	# save(gamAll,file=paste0( OutDir,"gamAll_SynOnly.RData" ))
	# save(tclin,file=paste0( OutDir,"tclin_gamAll_SynOnly.RData" ))

	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2
	colorz_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"classes"]	

	######## trying to match TMB
	load(file=paste0( OutDir,"gamAll.RData" ))
	load(file=paste0( OutDir,"tclin_gamAll.RData" ))
	load(file=paste0( OutDir,"../Genomic_TmbFga/TMB_nonsyn_Clinall.RData" ))
	gamAll = gamAll[rownames(gamAll) %in% rownames(Clinall),]
	tclin = tclin[rownames(tclin) %in% rownames(tclin),]
	cols_common = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="common",]),],na.rm=T)
	cols_uncommon = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="uncommon",]),],na.rm=T)
	cols_compound = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="compound",]),],na.rm=T)
	gamAll = gamAll[,((cols_common>=3) | (cols_uncommon>=3)) | (cols_compound>=3)]
	tclin = tclin[rownames(gamAll),]
	cdf = dan.df(0,c( "Gene","EGFR_class","nMut","Percentage","Chisq_pval","winner" ))
	for (g in colnames(gamAll)){
		this_tclin = tclin
		this_tclin$this_gene = gamAll[,g]
		this_tclin = this_tclin[!is.na(this_tclin$this_gene),]
		this_tclin[this_tclin$this_gene==0,"this_gene"] = "wt"
		this_tclin[this_tclin$this_gene==1,"this_gene"] = "mut"
		this_tclin$this_gene = factor(this_tclin$this_gene,levels=c("mut","wt"))
		this_tclin$egfr_class_consensus = factor(this_tclin$egfr_class_consensus,levels=ordered_classes)
		tabb = table(this_tclin$this_gene,this_tclin$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tcdf = data.frame(row.names = ordered_classes, Gene=rep(g,length(ordered_classes) ),EGFR_class=ordered_classes,Percentage=NA,Chisq_pval=NA,winner="no",stringsAsFactors=F)
		for (cl in ordered_classes){ tcdf[cl,"nMut"] = tabb["mut",cl] }
		for (cl in ordered_classes){ tcdf[cl,"Percentage"] = 100*tabb["mut",cl]/sum(tabb[,cl]) }
		tcdf[which.max( tcdf$Percentage),"winner"] = "yes"
		tcdf[tcdf$winner=="yes","Chisq_pval"] = ch
		cdf = rbind(cdf,tcdf)
	}
	cdf$Chisq_qval = NA
	cdf[cdf$winner=="yes","Chisq_qval"] = p.adjust(cdf[cdf$winner=="yes","Chisq_pval"],method="BH")
	cdf$EGFR_class = factor(cdf$EGFR_class, levels = ordered_classes)
	aa=cdf[order(cdf$Chisq_pval),]
	cdf_signif_genes = unique(cdf[(cdf$Chisq_qval<0.1) %in% c(T),"Gene"])
	scdf = cdf[cdf$Gene %in% cdf_signif_genes,]
	dcat( paste0("Uncommon percentage: ", 100*sum(scdf[(scdf$winner=="yes"),"EGFR_class"]=="uncommon")/nrow(scdf[(scdf$winner=="yes"),]) ))

	tclin_matched = tclin
	tclin_matched$TMB_nonsynonymous = Clinall[rownames(tclin_matched),"TMB_nonsynonymous"]
	dcat( "Unmatched TMB, quantiles (common, uncommon, compound): " )
	print(quantile(tclin_matched[tclin_matched$egfr_class_consensus=="common","TMB_nonsynonymous"]) )
	print(quantile(tclin_matched[tclin_matched$egfr_class_consensus=="uncommon","TMB_nonsynonymous"]) )
	print(quantile(tclin_matched[tclin_matched$egfr_class_consensus=="compound","TMB_nonsynonymous"]) )
	# matching TMB
	
	qq = quantile(tclin_matched[,"TMB_nonsynonymous"], seq(0,1,1/nbins))
	qq[1] = 0
	qq[length(qq)] = Inf # just to fix extremes
	tclin_matched$tmb_bin = NA
	for (class in c( "common","uncommon","compound" )){
		tclin_matched[tclin_matched$egfr_class_consensus==class,"tmb_bin"] = NA
		for (i in 1:(length(qq)-1) ){
			tclin_matched[(tclin_matched$egfr_class_consensus==class) & ((tclin_matched$TMB_nonsynonymous>=qq[i]) & (tclin_matched$TMB_nonsynonymous<qq[i+1])),"tmb_bin"] = i
		}	
	}
	dtable(tclin_matched$tmb_bin,tclin_matched$Dataset)
	matching_table = dtable(tclin_matched$tmb_bin,tclin_matched$egfr_class_consensus)
	# uncommon is the least numerous. We will downsample compound and common, trying to keep as many samples as possible
	matching_table_ratios = t(t(matching_table)/rowSums(t(matching_table)))
	matching_table_counts = matching_table
	for (class in c( "common","compound" )){
		matching_table_counts[,class] = floor((min(matching_table[,class]/matching_table_ratios[,"uncommon"]))*matching_table_ratios[,"uncommon"])
	}	
	set.seed(123)
	scdf_list = list()
	for (i in c(1:1000)){
		if (i==10){dcat(i,1)}
		tclin = dan.df(0,colnames(tclin_matched))
		for (class in c( "common","compound" )){
			for (bin in rownames(matching_table_counts)){
				thiz = tclin_matched[(tclin_matched$egfr_class_consensus==class) & (tclin_matched$tmb_bin==bin),]
				tclin = rbind(tclin,thiz[sample(rownames(thiz),matching_table_counts[bin,class]),])
			}
		}
		tclin = rbind(tclin,tclin_matched[tclin_matched$egfr_class_consensus=="uncommon",])
		# dcat( "Matched TMB, quantiles (common, uncommon, compound): " )
		# print(quantile(tclin[tclin$egfr_class_consensus=="common","TMB_nonsynonymous"]) )
		# print(quantile(tclin[tclin$egfr_class_consensus=="uncommon","TMB_nonsynonymous"]) )
		# print(quantile(tclin[tclin$egfr_class_consensus=="compound","TMB_nonsynonymous"]) )
		cdf = dan.df(0,c( "Gene","EGFR_class","nMut","Percentage","Chisq_pval","winner" ))
		for (g in colnames(gamAll)){
			this_tclin = tclin
			this_tclin$this_gene = gamAll[rownames(this_tclin),g]
			this_tclin = this_tclin[!is.na(this_tclin$this_gene),]
			this_tclin[this_tclin$this_gene==0,"this_gene"] = "wt"
			this_tclin[this_tclin$this_gene==1,"this_gene"] = "mut"
			this_tclin$this_gene = factor(this_tclin$this_gene,levels=c("mut","wt"))
			this_tclin$egfr_class_consensus = factor(this_tclin$egfr_class_consensus,levels=ordered_classes)
			tabb = table(this_tclin$this_gene,this_tclin$egfr_class_consensus)
			ch = chisq.test(tabb)$p.value
			tcdf = data.frame(row.names = ordered_classes, Gene=rep(g,length(ordered_classes) ),EGFR_class=ordered_classes,Percentage=NA,Chisq_pval=NA,winner="no",stringsAsFactors=F)
			for (cl in ordered_classes){ tcdf[cl,"nMut"] = tabb["mut",cl] }
			for (cl in ordered_classes){ tcdf[cl,"Percentage"] = 100*tabb["mut",cl]/sum(tabb[,cl]) }
			tcdf[which.max( tcdf$Percentage),"winner"] = "yes"
			tcdf[tcdf$winner=="yes","Chisq_pval"] = ch
			cdf = rbind(cdf,tcdf)
		}
		cdf$Chisq_qval = NA
		cdf[cdf$winner=="yes","Chisq_qval"] = p.adjust(cdf[cdf$winner=="yes","Chisq_pval"],method="BH")
		cdf$EGFR_class = factor(cdf$EGFR_class, levels = ordered_classes)
		scdf_list[[i]] = cdf[cdf$Gene %in% cdf_signif_genes,]
	}
	### mean ratios in TMB-matched scenarios for TP53 and KRAS
	tcdf = scdf[scdf$Gene %in% c("TP53","KRAS"),]
	for (class in c( "common","uncommon","compound" )){
		tcdf[(tcdf$Gene=="TP53") & (tcdf$EGFR_class==class), "mean_percentage"] = mean(sapply(scdf_list,function(x) x[(x$Gene=="TP53") & (x$EGFR_class==class), "Percentage"]))
		tcdf[(tcdf$Gene=="TP53") & (tcdf$EGFR_class==class), "sd_percentage"] = sd(sapply(scdf_list,function(x) x[(x$Gene=="TP53") & (x$EGFR_class==class), "Percentage"]))	
		tcdf[(tcdf$Gene=="KRAS") & (tcdf$EGFR_class==class), "mean_percentage"] = mean(sapply(scdf_list,function(x) x[(x$Gene=="KRAS") & (x$EGFR_class==class), "Percentage"]))
		tcdf[(tcdf$Gene=="KRAS") & (tcdf$EGFR_class==class), "sd_percentage"] = sd(sapply(scdf_list,function(x) x[(x$Gene=="KRAS") & (x$EGFR_class==class), "Percentage"]))
	}
	tcdf[(tcdf$Gene=="TP53") & (tcdf$EGFR_class=="uncommon"), "percent_confirmed"] = 100*sum(sapply(scdf_list,function(x) ( (x[(x$Gene=="TP53") & (x$EGFR_class=="uncommon"), "winner"]=="yes") %in% c(T)) & ( (x[(x$Gene=="TP53") & (x$EGFR_class=="uncommon"), "Chisq_qval"]<0.1) %in% c(T))  ))/length(scdf_list)
	tcdf[(tcdf$Gene=="KRAS") & (tcdf$EGFR_class=="uncommon"), "percent_confirmed"] = 100*sum(sapply(scdf_list,function(x) ( (x[(x$Gene=="KRAS") & (x$EGFR_class=="uncommon"), "winner"]=="yes") %in% c(T)) & ( (x[(x$Gene=="KRAS") & (x$EGFR_class=="uncommon"), "Chisq_qval"]<0.1) %in% c(T))  ))/length(scdf_list)
	
	### mean ratios in TMB-matched scenarios for all genes in scdf
	tcdf = scdf
	for (class in c( "common","uncommon","compound" )){
		for (g in unique(tcdf$Gene)){
			tcdf[(tcdf$Gene==g) & (tcdf$EGFR_class==class), "mean_percentage"] = mean(sapply(scdf_list,function(x) x[(x$Gene==g) & (x$EGFR_class==class), "Percentage"]))
			tcdf[(tcdf$Gene==g) & (tcdf$EGFR_class==class), "sd_percentage"] = sd(sapply(scdf_list,function(x) x[(x$Gene==g) & (x$EGFR_class==class), "Percentage"]))
		}
	}
	for (g in unique(tcdf$Gene)){
		winner_class = as.character(tcdf[(tcdf$Gene==g) & (!is.na(tcdf$Chisq_pval)),"EGFR_class"])
		tcdf[(tcdf$Gene==g) & (tcdf$EGFR_class==winner_class), "percent_confirmed"] = 100*sum(sapply(scdf_list,function(x) ( (x[(x$Gene==g) & (x$EGFR_class==winner_class), "winner"]=="yes") %in% c(T)) & ( (x[(x$Gene==g) & (x$EGFR_class==winner_class), "Chisq_qval"]<0.1) %in% c(T))  ))/length(scdf_list)
	}
	save(tcdf,file=paste0( OutDir,"EGFRclasses_vs_DriverMuts_AllDatasets_MatchingTmb.RData"))
	
	tcdf = tcdf[tcdf$Gene %in% c("TP53","KRAS"),]
	tcdf$percent_confirmed[!is.na(tcdf$percent_confirmed)] = paste0("Significant in ",tcdf$percent_confirmed[!is.na(tcdf$percent_confirmed)],"% of trials" )
	pdf(paste0( OutDir,"EGFRclasses_vs_TP53mutKRASmut_AllDatasets_MatchingTmb.pdf" ),4.5,3)
	pl = ggplot(data=tcdf, aes(x=Gene, y=mean_percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + geom_linerange( aes(ymin=mean_percentage-sd_percentage, ymax=mean_percentage+sd_percentage),stat="identity", position=position_dodge(width = 0.9),linewidth = 0.4) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0( "TMB-matched subsampled cases, ",length(scdf_list)," trials" )) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=percent_confirmed,y=mean_percentage),vjust=-0.2, size = size_labels) + geom_text(aes(label=round(mean_percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	pl=pl+theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pl)
	dev.off()	
}

cna_associations_alldatasets = function( OutDir, ncopies = 2 ){ 
	genez = dan.read(paste0( DataDir,"mina_natgen_2020_genes.txt" ))
	rownames(genez) = genez$gene
	###### Genie
	dcat( "Processing Genie",1 )
	gi = dan.read(paste0(DataDir,"genie_14.0/genomic_information.txt"))
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2
	colorz_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin2g$egfr_class_consensus),"classes"]	
	load(file = paste0( DataDir,"genie_14.0/cna_Genie_EGFRmutants.RData" ))
	dtable(unlist(cna))
	tClin2 = Clin2g[Clin2g$Sample %in% colnames(cna),]
	tclin = Cling
	tclin = Cling[rownames(Cling) %in% tClin2$Patient,]
	cna = cna[,intersect(rownames(Clin2), colnames(cna))]
	gam = dan.df(unique(tclin$Patient),rownames(genez), data=NA)
	for (cn in colnames(gam)){
		dcat(cn)
		tcna = cna[,!is.na(cna[cn,])]
		if (sum(!is.na(cna[cn,]))<2) {next}
		if (genez[cn,"role"]=="og"){
				samples_alt = colnames(cna)[which(tcna[cn,]>=ncopies)]
				samples_wt = colnames(cna)[which(tcna[cn,]<ncopies)]
			} else {
				samples_alt = colnames(cna)[which(tcna[cn,]<=(-ncopies))]
				samples_wt = colnames(cna)[which(tcna[cn,]>(-ncopies))]
			}
		gam[unique(tClin2[samples_wt,"Patient"]),cn] = 0
		gam[unique(tClin2[samples_alt,"Patient"]),cn] = 1
	}
	tclin$Dataset = "Genie"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = gam
	clinAll = tclin

	dcat( "Processing TCGA",1 )
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TCGA_SampleLevel.RData" ))
	load(file = paste0( CommonDataDir,"cna_LUAD_all.RData" ))
	cna = t(cna_ss)
	cna = cna[,colnames(cna) %in% rownames(Clin)]
	tclin = Clin[colnames(cna),]
	gam = dan.df(unique(tclin$Patient),intersect(rownames(cna),rownames(genez)), data=0)
	for (cn in colnames(gam)){
		tcna = cna[,!is.na(cna[cn,])]
		if (sum(!is.na(cna[cn,]))<2) {next}
		if (genez[cn,"role"]=="og"){
				samples_alt = colnames(cna)[which(tcna[cn,]>=ncopies)]
				samples_wt = colnames(cna)[which(tcna[cn,]<ncopies)]
			} else {
				samples_alt = colnames(cna)[which(tcna[cn,]<=(-ncopies))]
				samples_wt = colnames(cna)[which(tcna[cn,]>(-ncopies))]
			}
		gam[samples_wt,cn] = 0
		gam[samples_alt,cn] = 1
	}
	tclin$Dataset = "TCGA"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	dcat( "Processing ChenEAS",1 )
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	cna = dan.read( "/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/cnv.tsv" )
	cna = cna[!duplicated(cna$Hugo_Symbol),]
	rownames(cna) = cna$Hugo_Symbol
	cna$Hugo_Symbol = NULL
	cna$Entrez_Gene_Id = NULL
	colnames(cna) = gsub("\\.","-",colnames(cna))
	cna = cna[,colnames(cna) %in% rownames(Clin)]
	tclin = Clin[colnames(cna),]
	gam = dan.df(unique(tclin$Patient),intersect(rownames(cna),rownames(genez)), data=0)
	for (cn in colnames(gam)){
		tcna = cna[,!is.na(cna[cn,])]
		if (sum(!is.na(cna[cn,]))<2) {next}
		if (genez[cn,"role"]=="og"){
				samples_alt = colnames(cna)[which(tcna[cn,]>=ncopies)]
				samples_wt = colnames(cna)[which(tcna[cn,]<ncopies)]
			} else {
				samples_alt = colnames(cna)[which(tcna[cn,]<=(-ncopies))]
				samples_wt = colnames(cna)[which(tcna[cn,]>(-ncopies))]
			}
		gam[samples_wt,cn] = 0
		gam[samples_alt,cn] = 1
	}
	tclin$Dataset = "ChenEAS"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	dcat( "Processing Zhang",1 )
	dcat( "CNA not available",2 )

	dcat( "Processing TRACERx421",1 )
	dcat( "CNA not readily available",2 )

	###### Origimed
	dcat( "Processing Origimed",1 )
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Origimed_SampleLevel.RData" ))
	cna = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_cna.txt"))
	colnames(cna) = gsub("\\.","-",colnames(cna))
	rownames(cna) = cna$Hugo_Symbol
	cna$Hugo_Symbol = NULL
	cna$Entrez_Gene_Id = NULL	
	tClin2 = Clin2[Clin2$Sample %in% colnames(cna),]
	tclin = Clin[rownames(Clin) %in% tClin2$Patient,]
	cna = cna[,intersect(tClin2$Sample, colnames(cna))]
	rownames(tClin2) = tClin2$Sample
	gam = dan.df(unique(tclin$Patient),intersect(rownames(cna),rownames(genez)), data=0)
	for (cn in colnames(gam)){
		dcat(cn)
		tcna = cna[,!is.na(cna[cn,])]
		if (sum(!is.na(cna[cn,]))<2) {next}
		if (genez[cn,"role"]=="og"){
				samples_alt = colnames(cna)[which(tcna[cn,]>=ncopies)]
				samples_wt = colnames(cna)[which(tcna[cn,]<ncopies)]
			} else {
				samples_alt = colnames(cna)[which(tcna[cn,]<=(-ncopies))]
				samples_wt = colnames(cna)[which(tcna[cn,]>(-ncopies))]
			}
		gam[unique(tClin2[samples_wt,"Patient"]),cn] = 0
		gam[unique(tClin2[samples_alt,"Patient"]),cn] = 1
	}
	tclin$Dataset = "Origimed"
	gam = cbind(gam,dan.df(rownames(gam),rownames(genez)[!(rownames(genez) %in% colnames(gam))],NA ))
	gamAll = rbind(gamAll[,intersect(colnames(gamAll),colnames(gam))],gam[,intersect(colnames(gamAll),colnames(gam))])
	clinAll = rbind(clinAll[,intersect(colnames(clinAll),colnames(tclin))], tclin[,intersect(colnames(clinAll),colnames(tclin))])

	# gamAll = gamAll[,colnames(gamAll)!="EGFR"]
	gamAll = gamAll[rownames(clinAll),]
	tclin = clinAll
	all(rownames(gamAll) %in% rownames(tclin))
	cols_common = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="common",]),],na.rm=T)
	cols_uncommon = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="uncommon",]),],na.rm=T)
	cols_compound = colSums(gamAll[rownames(tclin[tclin$egfr_class_consensus=="compound",]),],na.rm=T)
	gamAll = gamAll[,((cols_common>=3) | (cols_uncommon>=3)) | (cols_compound>=3)]
	tclin = tclin[rownames(gamAll),]
	cdf = dan.df(0,c( "Gene","EGFR_class","nMut","Percentage","Chisq_pval","winner" ))
	for (g in colnames(gamAll)){
		this_tclin = tclin
		this_tclin$this_gene = gamAll[,g]
		this_tclin = this_tclin[!is.na(this_tclin$this_gene),]
		this_tclin[this_tclin$this_gene==0,"this_gene"] = "wt"
		this_tclin[this_tclin$this_gene==1,"this_gene"] = "mut"
		this_tclin$this_gene = factor(this_tclin$this_gene,levels=c("mut","wt"))
		this_tclin$egfr_class_consensus = factor(this_tclin$egfr_class_consensus,levels=ordered_classes)
		tabb = table(this_tclin$this_gene,this_tclin$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tcdf = data.frame(row.names = ordered_classes, Gene=rep(g,length(ordered_classes) ),EGFR_class=ordered_classes,Percentage=NA,Chisq_pval=NA,winner="no",stringsAsFactors=F)
		for (cl in ordered_classes){ tcdf[cl,"nMut"] = tabb["mut",cl] }
		for (cl in ordered_classes){ tcdf[cl,"Percentage"] = 100*tabb["mut",cl]/sum(tabb[,cl]) }
		tcdf[which.max( tcdf$Percentage),"winner"] = "yes"
		tcdf[tcdf$winner=="yes","Chisq_pval"] = ch
		cdf = rbind(cdf,tcdf)
	}
	cdf$Chisq_qval = NA
	cdf[cdf$winner=="yes","Chisq_qval"] = p.adjust(cdf[cdf$winner=="yes","Chisq_pval"],method="BH")
	cdf$EGFR_class = factor(cdf$EGFR_class, levels = ordered_classes)
	cdf=cdf[order(cdf$Chisq_pval),]
	cdf_signif_genes = unique(cdf[(cdf$Chisq_qval<0.1) %in% c(T),"Gene"])
	cdf = cdf[cdf$Gene %in% cdf_signif_genes,]
	save(cdf,file=paste0( OutDir,"AllDatasets_signifCoCNAs_table.RData" ))
	dcat( paste0("Uncommon percentage: ", 100*sum(cdf[(cdf$winner=="yes"),"EGFR_class"]=="uncommon")/nrow(cdf[(cdf$winner=="yes"),]) ))

	tcdf = cdf[cdf$Gene=="EGFR",]
	tcdf$Chisq_pval = signif(tcdf$Chisq_pval,2)
	tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)] = paste0("p-value = ",tcdf$Chisq_pval[!is.na(tcdf$Chisq_pval)])
	tcdf$Gene = paste0("EGFR amplifications\n(GISTIC>=",ncopies,")")
	pdf(paste0( OutDir,"cnas_EGFRclasses_vs_EGFRamp_AllDatasets_ncopies",ncopies,".pdf" ),2.5,2.5)
	pl = ggplot(data=tcdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets\nEGFR-amplified patients = ",sum(tcdf$nMut)," out of ",sum(!is.na(gamAll$EGFR))) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 0,size=12)) + geom_text(aes(label=Chisq_pval,y=Percentage),vjust=-0.2,size = size_labels) + geom_text(aes(label=round(Percentage,1),size = size_labels), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	pl = pl + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pl)
	dev.off()

	cdf = cdf[order(cdf$Chisq_pval),]
	genes_keep = cdf[1:20,"Gene"]
	cdf = cdf[cdf$Gene %in% genes_keep,]
	cdf$Gene = factor(cdf$Gene, levels = unique(cdf$Gene))
	cdf$Chisq_pval = signif(cdf$Chisq_pval,2)
	cdf$Chisq_qval = signif(cdf$Chisq_qval,2)
	cdf[order(cdf$Gene),]

	pdf(paste0( OutDir,"cnas_AllDatasets_coaltered_barplot_top20_ncopies",ncopies,".pdf" ),4,2.5)
	pl = ggplot(data=cdf, aes(x=Gene, y=Percentage, fill=EGFR_class)) + geom_bar(stat="identity", position=position_dodge(),colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, co-altered genes\nChi-square test adjusted p-values") ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) + geom_text(aes(label=Chisq_qval,y=Percentage),vjust=-0.2,size = size_labels) + geom_text(aes(label=round(Percentage,1)), size = size_labels, vjust = 1.5, position = position_dodge(.9))
	pl = pl + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pl)
	dev.off()

	save(gamAll,file=paste0( OutDir,"cnas_gamAll.RData" ))
	save(tclin,file=paste0( OutDir,"cnas_tclin_gamAll.RData" ))
}

kras_further_analyses = function( OutDir ){	
	hotspots = dan.read(paste0(DataDir,"hotspots_mskcc_chang2016.txt"))
	hotspots$mut = paste0("p.",hotspots$Reference_Amino_Acid,hotspots$Amino_Acid_Position,hotspots$Variant_Amino_Acid )

	load(file=paste0( OutDir,"../Preprocessing/tcAll_vafs.RData" ))
	load(file=paste0( OutDir,"tmaf_kras.RData" ))
	load(file=paste0( OutDir,"tclin_gamAll.RData" ))
	for (p in tmaf_kras$Patient){
		tmaf_kras[tmaf_kras$Patient==p,"egfr_class_consensus"] = tclin[p,"egfr_class_consensus"]
		tmaf_kras[tmaf_kras$Patient==p,"Dataset"] = tclin[p,"Dataset"]
	}
	ordered_classes = clco[rownames(clco) %in% (tmaf_kras$egfr_class_consensus),"classes"]

	tabb = table(tmaf_kras$egfr_class_consensus,tmaf_kras[,"KRAS_variant_site"])
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","KRAS_variant_site","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="compound"))
	mt[(mt$EGFR_class=="T790M") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="T790M"))
	mt[(mt$EGFR_class=="ex20ins") & (mt$KRAS_variant_site=="p.Q61"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="ex20ins"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$Percentage_label = ifelse((mt$KRAS_variant_site=="p.G12") & (mt$Percentage>0),paste0("G12: ",round(mt$Percentage),"%"),"")
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=KRAS_variant_site)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("EGFR class vs KRAS hotspot variant sites\n","Chi-square test, p = ",signif(ch,2)) ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2, size = size_labels) + geom_text(aes(label=Percentage_label), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_KRASvariantsite.pdf"),3,2.5,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()

	tabb = table(tmaf_kras$egfr_class_consensus,tmaf_kras[,"KRAS_variant"])
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","KRAS_variant","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$KRAS_variant=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$KRAS_variant=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$KRAS_variant=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="compound"))
	mt[(mt$EGFR_class=="T790M") & (mt$KRAS_variant_site=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="T790M"))
	mt[(mt$EGFR_class=="ex20ins") & (mt$KRAS_variant_site=="p.Q61L"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="ex20ins"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$Percentage_label = ifelse((mt$KRAS_variant=="p.G12C") & (mt$Percentage>0),paste0("G12C: ",round(mt$Percentage),"%"),"")
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=KRAS_variant)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("EGFR class vs KRAS hotspot variants\n","Chi-square test, p = ",signif(ch,2)) ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2, size = size_labels) + geom_text(aes(label=Percentage_label), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_KRASvariant.pdf"),3,2.5,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()

	dtable(tmaf_kras$Dataset) # Genie, Origimed
	load(file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	maf_ss = merge(maf_ss,Clin2[,c( "Patient","Sample" )],by.x="Tumor_Sample_Barcode",by.y="Sample",all.x=T)
	maf_ss[is.na(maf_ss$Patient),"Patient"] = "others"
	ttmaf = maf_ss[(maf_ss$Hugo_Symbol=="KRAS") & ((maf_ss$Variant_Classification=="Missense_Mutation") & (maf_ss$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol=="KRAS","mut"])),]
	ttmaf = ttmaf[!(ttmaf$Patient %in% tmaf_kras$Patient),]
	ttmaf_kras = data.frame(Patient=ttmaf$Tumor_Sample_Barcode,KRAS_variant=ttmaf$HGVSp_Short,KRAS_variant_site=paste0(substr(ttmaf$HGVSp_Short,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$HGVSp_Short)),KRAS_vaf=as.numeric(ttmaf$t_alt_count)/as.numeric(ttmaf$t_depth),egfr_class_consensus="EGFR-wt",Dataset="Genie",stringsAsFactors=F)
	tmaf_kras = rbind(tmaf_kras,ttmaf_kras)

	maf_ss = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_mutations_extended.txt"))
	Clin2 = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_clinical_sample.txt"))
	maf_ss = merge(maf_ss,Clin2[,c( "PATIENT_ID","SAMPLE_ID" )],by.x="Tumor_Sample_Barcode",by.y="SAMPLE_ID",all.x=T)
	maf_ss$Patient = maf_ss$PATIENT_ID
	maf_ss[is.na(maf_ss$Patient),"Patient"] = "others"
	ttmaf = maf_ss[(maf_ss$Hugo_Symbol=="KRAS") & ((maf_ss$Variant_Classification=="Missense_Mutation") & (maf_ss$HGVSp_Short %in% hotspots[hotspots$Hugo_Symbol=="KRAS","mut"])),]
	ttmaf = ttmaf[!(ttmaf$Patient %in% tmaf_kras$Patient),]
	ttmaf_kras = data.frame(Patient=ttmaf$Tumor_Sample_Barcode,KRAS_variant=ttmaf$HGVSp_Short,KRAS_variant_site=paste0(substr(ttmaf$HGVSp_Short,1,3),gsub(".*?([0-9]+).*", "\\1", ttmaf$HGVSp_Short)),KRAS_vaf=NA,egfr_class_consensus="EGFR-wt",Dataset="Genie",stringsAsFactors=F)
	tmaf_kras = rbind(tmaf_kras,ttmaf_kras)

	tabb = table(tmaf_kras$egfr_class_consensus,tmaf_kras[,"KRAS_variant_site"])
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","KRAS_variant_site","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="EGFR-wt") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="EGFR-wt"))
	mt[(mt$EGFR_class=="common") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="compound"))
	mt[(mt$EGFR_class=="T790M") & (mt$KRAS_variant_site=="p.G12"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="T790M"))
	mt[(mt$EGFR_class=="ex20ins") & (mt$KRAS_variant_site=="p.Q61"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="ex20ins"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=c("EGFR-wt",ordered_classes))
	mt$Percentage_label = ifelse((mt$KRAS_variant_site=="p.G12") & (mt$Percentage>0),paste0("G12: ",round(mt$Percentage),"%"),"")
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=KRAS_variant_site)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("EGFR class vs KRAS hotspot variant sites\n","Chi-square test, p = ",signif(ch,2)) ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2, size = size_labels) + geom_text(aes(label=Percentage_label), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_KRASvariantsite_withEGFRwt.pdf"),4,2.5,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()

	tabb = table(tmaf_kras$egfr_class_consensus,tmaf_kras[,"KRAS_variant"])
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","KRAS_variant","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="EGFR-wt") & (mt$KRAS_variant=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="EGFR-wt"))
	mt[(mt$EGFR_class=="common") & (mt$KRAS_variant=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$KRAS_variant=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$KRAS_variant=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="compound"))
	mt[(mt$EGFR_class=="T790M") & (mt$KRAS_variant_site=="p.G12C"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="T790M"))
	mt[(mt$EGFR_class=="ex20ins") & (mt$KRAS_variant_site=="p.Q61L"),"TotalPatients"] = paste0( "N = ",sum(tmaf_kras$egfr_class_consensus=="ex20ins"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=c("EGFR-wt",ordered_classes))
	mt$Percentage_label = ifelse((mt$KRAS_variant=="p.G12C") & (mt$Percentage>0),paste0("G12C: ",round(mt$Percentage),"%"),"")
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=KRAS_variant)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("EGFR class vs KRAS hotspot variants\n","Chi-square test, p = ",signif(ch,2)) ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2, size = size_labels) + geom_text(aes(label=Percentage_label), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(paste0(OutDir,"EGFRclasses_KRASvariant_withEGFRwt.pdf"),4,2.5,onefile=FALSE)
	pa = pa + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()

	### vafs analyses
	tmaf_kras = tmaf_kras[!is.na(tmaf_kras$KRAS_vaf),]
	colorz_classes = c("black", clco[rownames(clco) %in% (tmaf_kras$egfr_class_consensus),"colorz"])
	ordered_classes = c("EGFR-wt", clco[rownames(clco) %in% (tmaf_kras$egfr_class_consensus),"classes"])
	# consensus class
	y = tmaf_kras$KRAS_vaf
	x = factor(tmaf_kras[,"egfr_class_consensus"],levels = ordered_classes)
	plotTitle = paste0( "Kruskal-Wallis test, p-value = ", signif(kruskal.test(y~x)$p.value,2))
	dan.densityPlot( paste0(OutDir,"vaf_across_EGFR_consensus_withEGFRwt.pdf"), y, x, groupinglab = "", xlab = paste0("KRAS mutation VAF"), ylab = "Density", show_medians = F, plotTitle = plotTitle,xlimLeft = 0, xlimRight = 1, groupingColors = colorz_classes, fileWidth = 3, fileHeight = 2 )

	load(file=paste0( OutDir,"../Preprocessing/tcAll_vafs.RData" ))
	load(file=paste0( OutDir,"tmaf_kras.RData" ))
	load(file=paste0( OutDir,"tclin_gamAll.RData" ))
	for (p in tmaf_kras$Patient){
		tmaf_kras[tmaf_kras$Patient==p,"egfr_class_consensus"] = tclin[p,"egfr_class_consensus"]
		tmaf_kras[tmaf_kras$Patient==p,"Dataset"] = tclin[p,"Dataset"]
	}
	ordered_classes = clco[rownames(clco) %in% (tmaf_kras$egfr_class_consensus),"classes"]

	tcAll$KRAS_status = "wt"
	tcAll[tcAll$Patient %in% tmaf_kras$Patient,"KRAS_status"] = "mut"
	thiz = tcAll[tcAll$EGFR_consensus=="uncommon",]
	y = thiz$vaf
	x = factor(thiz[,"KRAS_status"],levels = c( "wt","mut" ))
	plotTitle = paste0( "Wilcoxon test, p-value = ", signif(wilcox.test(y~x)$p.value,2))
	dan.densityPlot( paste0(OutDir,"vaf_inUncommons_acrossKRASstatus.pdf"), y, x, groupinglab = "KRAS", xlab = paste0("EGFR mutation VAF"), ylab = "Density", show_medians = F, plotTitle = plotTitle,xlimLeft = 0, xlimRight = 1, groupingColors = c("tomato3","black"), fileWidth = 2.5, fileHeight = 2 )

	tcAll$KRAS_status = "wt"
	tcAll[tcAll$Patient %in% tmaf_kras$Patient,"KRAS_status"] = "mut"
	thiz = tcAll[tcAll$EGFR_consensus=="common",]
	y = thiz$vaf
	x = factor(thiz[,"KRAS_status"],levels = c( "wt","mut" ))
	plotTitle = paste0( "Wilcoxon test, p-value = ", signif(wilcox.test(y~x)$p.value,2))
	dan.densityPlot( paste0(OutDir,"vaf_inCommons_acrossKRASstatus.pdf"), y, x, groupinglab = "KRAS", xlab = paste0("EGFR mutation VAF"), ylab = "Density", show_medians = F, plotTitle = plotTitle,xlimLeft = 0, xlimRight = 1, groupingColors = c("steelblue4","black"), fileWidth = 2.5, fileHeight = 2 )

	thiz_kras = aggregate(KRAS_vaf~Patient,data=tmaf_kras,FUN='mean')
	rownames(thiz_kras) = thiz_kras$Patient
	thiz = tcAll[tcAll$Patient %in% rownames(thiz_kras),]
	thiz = aggregate(vaf~Patient,data=thiz,FUN='mean')
	rownames(thiz) = thiz$Patient
	thiz_kras$EGFR_vaf = thiz[rownames(thiz_kras),"vaf"]
	cor(thiz_kras$EGFR_vaf,thiz_kras$KRAS_vaf)
	thiz_kras$egfr_class_consensus = tclin[rownames(thiz_kras),"egfr_class_consensus"]
	plotTitle = paste0( "Spearman R = ",signif(cor(thiz_kras$KRAS_vaf,thiz_kras$EGFR_vaf,method='spearman'),2),", p-value = ",signif(cor.test(thiz_kras$KRAS_vaf,thiz_kras$EGFR_vaf,method='spearman')$p.value,2) )
	dan.scatterplot(paste0(OutDir,"vaf_EGFR_KRAS_matched_scatterplot.pdf"),x=thiz_kras$KRAS_vaf,xlab='mean KRAS mutation VAF',ylab='mean EGFR mutation VAF',y=thiz_kras$EGFR_vaf,fill=thiz_kras$egfr_class_consensus,filllab="EGFR mutation class",fillColors=c("steelblue4","tomato3"),plotTitle=plotTitle,plotFitLine=T,FitLineColor='gray22',fileWidth=2.5,fileHeight=2,coord_fixed=T )

	thiz_kras_common = thiz_kras[thiz_kras$egfr_class_consensus=="common",]
	plotTitle = paste0( "Spearman R = ",signif(cor(thiz_kras_common$KRAS_vaf,thiz_kras_common$EGFR_vaf,method='spearman'),2),", p-value = ",signif(cor.test(thiz_kras_common$KRAS_vaf,thiz_kras_common$EGFR_vaf,method='spearman')$p.value,2) )
	dan.scatterplot(paste0(OutDir,"vaf_EGFR_KRAS_matched_scatterplot_common.pdf"),x=thiz_kras_common$KRAS_vaf,xlab='mean KRAS mutation VAF',ylab='mean EGFR mutation VAF',y=thiz_kras_common$EGFR_vaf,fill=thiz_kras_common$egfr_class_consensus,filllab="EGFR mutation class",fillColors=c("steelblue4"),plotTitle=plotTitle,plotFitLine=T,FitLineColor='gray22',fileWidth=2.5,fileHeight=2,coord_fixed=T )

	thiz_kras_common = thiz_kras[thiz_kras$egfr_class_consensus=="uncommon",]
	plotTitle = paste0( "Spearman R = ",signif(cor(thiz_kras_common$KRAS_vaf,thiz_kras_common$EGFR_vaf,method='spearman'),2),", p-value = ",signif(cor.test(thiz_kras_common$KRAS_vaf,thiz_kras_common$EGFR_vaf,method='spearman')$p.value,2) )
	dan.scatterplot(paste0(OutDir,"vaf_EGFR_KRAS_matched_scatterplot_uncommon.pdf"),x=thiz_kras_common$KRAS_vaf,xlab='mean KRAS mutation VAF',ylab='mean EGFR mutation VAF',y=thiz_kras_common$EGFR_vaf,fill=thiz_kras_common$egfr_class_consensus,filllab="EGFR mutation class",fillColors=c("tomato3"),plotTitle=plotTitle,plotFitLine=T,FitLineColor='gray22',fileWidth=2.5,fileHeight=2,coord_fixed=T )


	### did they arise post-treatment?
	load(file=paste0( OutDir,"../Preprocessing/tcAll_vafs.RData" ))
	load(file=paste0( OutDir,"tmaf_kras.RData" ))
	load(file=paste0( OutDir,"tclin_gamAll.RData" ))
	for (p in tmaf_kras$Patient){
		tmaf_kras[tmaf_kras$Patient==p,"egfr_class_consensus"] = tclin[p,"egfr_class_consensus"]
		tmaf_kras[tmaf_kras$Patient==p,"Dataset"] = tclin[p,"Dataset"]
	}
	ordered_classes = clco[rownames(clco) %in% (tmaf_kras$egfr_class_consensus),"classes"]

	Clin = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_clinical_patient.txt"))
	Clin = Clin[Clin$PATIENT_ID %in% tmaf_kras$Patient,] # Origimed: 1 "Other_Treatments" (not targeted therapy), 2 "Treatment-naive"
	load(file=paste0(OutDir,"../Therapy/clinall_treatments.RData"))
	clinall = clinall[clinall$Patient %in% tmaf_kras$Patient,] # 4 received targeted treatment, 3 didn't, the others are NA's
	reg = read.csv(paste0(DataDir,"genie_BPC_NSCLC/regimen_cancer_level_dataset.csv"),stringsAsFactors=F) # here, unclear whether treatment came before or after sequencing
}

mutationalSignatures_preprocessing = function( OutDir ){

	library(deconstructSigs)

	### Zhang
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Zhang_SampleLevel.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clinz = Clin
	maf = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_mutations.txt"))
	maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	rownames(Clin2) = Clin2$Patient
	Clin$Sample = Clin2[rownames(Clin),"Sample"]
	rownames(Clin) = Clin$Sample
	maf_ss$Patient = Clin[maf_ss$Tumor_Sample_Barcode,"Patient"]
	rownames(Clin) = Clin$Patient
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	ts_df = data.frame( Patient = maf_snp$Patient, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	save(ts_df, file = paste0(OutDir, "Zhang_ts_df.RData"))
	load(file = paste0(OutDir, "Zhang_ts_df.RData")) # ts_df
	sigs.input = mut.to.sigs.input( mut.ref = ts_df, sample.id = "Patient", chr = "chr", pos = "pos", ref = "ref", alt = "alt" )
	save(sigs.input, file = paste0(OutDir, "Zhang_sigs_input.RData"))
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = ts_df$Patient[1], 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "Zhang_Sign_weights_unnormalized.RData" ))

	### TRACERx421
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TRACERx421_SampleLevel.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clinz = Clin
	library(fst)
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$tumour_id
	maf_ss$Patient = maf_ss$patient_id
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$ref) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$var) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	ts_df = data.frame( Patient = maf_snp$Patient, chr = paste0(maf_snp$chr), pos = as.numeric(maf_snp$start), ref = maf_snp$ref, alt = maf_snp$var) # for all of them, end-start=0
	save(ts_df, file = paste0(OutDir, "TRACERx421_ts_df.RData"))
	load(file = paste0(OutDir, "TRACERx421_ts_df.RData")) # ts_df
	sigs.input = mut.to.sigs.input( mut.ref = ts_df, sample.id = "Patient", chr = "chr", pos = "pos", ref = "ref", alt = "alt" )
	save(sigs.input, file = paste0(OutDir, "TRACERx421_sigs_input.RData"))
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = ts_df$Patient[1], 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "TRACERx421_Sign_weights_unnormalized.RData" ))

	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clint = Clin
	load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
	maf_ss$Patient = substr(maf_ss$Tumor_Sample_Barcode,1,12)
	maf_ss$Sample = substr(maf_ss$Tumor_Sample_Barcode,1,16)
	tclin = Clint[Clint$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Match_Norm_Seq_Allele1) %in% c("A", "T", "C", "G")) & ((maf$Match_Norm_Seq_Allele2) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	ts_df = data.frame( Patient = substr(maf_snp$Tumor_Sample_Barcode,1,12), chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Match_Norm_Seq_Allele1, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	save(ts_df, file = paste0(OutDir, "TCGA_ts_df.RData"))
	load(file = paste0(OutDir, "TCGA_ts_df.RData")) # ts_df
	sigs.input = mut.to.sigs.input( mut.ref = ts_df, sample.id = "Patient", chr = "chr", pos = "pos", ref = "ref", alt = "alt" )
	save(sigs.input, file = paste0(OutDir, "TCGA_sigs_input.RData"))
	# try both exome2genome and default (no normalization)
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = "TCGA-05-4402", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'exome2genome')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
		prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'exome2genome')
		Sign_weights_df[rn,] = prova$weights
	  	index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "TCGA_Sign_weights_exome2genome.RData" ))
	# try both exome2genome and default (no normalization)
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = "TCGA-05-4402", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
		prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
		Sign_weights_df[rn,] = prova$weights
	  	index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "TCGA_Sign_weights_unnormalized.RData" ))
	## compare the two
	load( file = paste0(OutDir, "TCGA_Sign_weights_exome2genome.RData" ))
	aa = Sign_weights_df
	load( file = paste0(OutDir, "TCGA_Sign_weights_unnormalized.RData" ))
	bb = Sign_weights_df
	ma = melt(aa)
	mb = melt(bb)
	fileName = paste0(OutDir, "TCGA_exome2genome_vs_unnormalized.pdf")
	x = ma$value
	y = mb$value
	dan.scatterplot(fileName, x=x, y=y, fill=ma$variable, fillColors=dan.colors(ma$variable), plotFitLine = T, plotBisector=T, dotSize = 2, plotTitle=paste0("Pearson r = ",signif(cor(x,y),2)), filllab="Signature", xlab="exome2genome weigths", ylab="unnormalized weigths", fileWidth = 8 )

	### Origimed
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clinz = Clin
	maf = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_mutations_extended.txt"))
	maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin$Sample,]
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	rownames(Clin) = Clin$Sample
	maf_ss$Patient = Clin[maf_ss$Tumor_Sample_Barcode,"Patient"]
	rownames(Clin) = Clin$Patient
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	ts_df = data.frame( Patient = maf_snp$Patient, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	save(ts_df, file = paste0(OutDir, "Origimed_ts_df.RData"))
	load(file = paste0(OutDir, "Origimed_ts_df.RData")) # ts_df
	sigs.input = mut.to.sigs.input( mut.ref = ts_df, sample.id = "Patient", chr = "chr", pos = "pos", ref = "ref", alt = "alt" )
	save(sigs.input, file = paste0(OutDir, "Origimed_sigs_input.RData"))
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = ts_df$Patient[1], 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "Origimed_Sign_weights_unnormalized.RData" ))

	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clinc = Clin
	maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
	maf$Patient = maf$Tumor_Sample_Barcode
	tclin = Clinc[Clinc$Patient %in% maf$Patient,]
	maf = maf[maf$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	ts_df = data.frame( Patient = maf_snp$Patient, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	save(ts_df, file = paste0(OutDir, "ChenEAS_ts_df.RData"))
	load(file = paste0(OutDir, "ChenEAS_ts_df.RData")) # ts_df
	sigs.input = mut.to.sigs.input( mut.ref = ts_df, sample.id = "Patient", chr = "chr", pos = "pos", ref = "ref", alt = "alt" )
	save(sigs.input, file = paste0(OutDir, "ChenEAS_sigs_input.RData"))
	# try both exome2genome and default (no normalization)
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = "BGI-EX02", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'exome2genome')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'exome2genome')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "ChenEAS_Sign_weights_exome2genome.RData" ))
	# try both exome2genome and default (no normalization)
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = "BGI-EX02", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "ChenEAS_Sign_weights_unnormalized.RData" ))
	## compare the two
	load( file = paste0(OutDir, "ChenEAS_Sign_weights_exome2genome.RData" ))
	aa = Sign_weights_df
	load( file = paste0(OutDir, "ChenEAS_Sign_weights_unnormalized.RData" ))
	bb = Sign_weights_df
	ma = melt(aa)
	mb = melt(bb)
	fileName = paste0(OutDir, "ChenEAS_exome2genome_vs_unnormalized.pdf")
	x = ma$value
	y = mb$value
	dan.scatterplot(fileName, x=x, y=y, fill=ma$variable, fillColors=dan.colors(ma$variable), plotFitLine = T, plotBisector=T, dotSize = 2, plotTitle=paste0("Pearson r = ",signif(cor(x,y),2)), filllab="Signature", xlab="exome2genome weigths", ylab="unnormalized weigths", fileWidth = 8 )

	###### Genie
	gi = dan.read(paste0(DataDir,"genie_14.0/genomic_information.txt"))
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2
	load(file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	### let's find a gene universe that we can test
	assays = sort(unique(Clin2g$SEQ_ASSAY_ID))
	asdf = data.frame(row.names = assays, assay = assays, ngenes = NA, nsamples = NA, npatients = NA, stringsAsFactors=FALSE)
	for (a in assays){
	   gia = gi[gi$SEQ_ASSAY_ID==a,]
	   ug = unique(gia$Hugo_Symbol[nchar(gia$Hugo_Symbol)>0])
	   asdf[a,"ngenes"] = length(ug)
	   asdf[a,"nsamples"] = length(Clin2g[Clin2g$SEQ_ASSAY_ID==a,"Sample"])
	   asdf[a,"npatients"] = length(unique(Clin2g[Clin2g$SEQ_ASSAY_ID==a,"Patient"]))
	}
	asdf = asdf[order(asdf$ngenes,decreasing=T),]
	# Choice: remove panels with less than 100 genes
	asdf = asdf[asdf$ngenes>=100,]
	# asdf = asdf[asdf$npatients>=100,]
	asdf = asdf[order(asdf$nsamples,decreasing=T),]
	first = TRUE
	for (a in rownames(asdf)){
	   gia = gi[gi$SEQ_ASSAY_ID==a,]
	   ug = unique(gia$Hugo_Symbol[nchar(gia$Hugo_Symbol)>0])
	   if (first){
	      common_genes = ug
	      first = FALSE
	   } else {
	      common_genes = intersect(common_genes,ug)
	   }
	   dcat(a)
	   dcat(length(common_genes),1)
	}
	tClin2 = Clin2g[Clin2g$SEQ_ASSAY_ID %in% rownames(asdf),]
	tclin = Cling[rownames(Cling) %in% tClin2$Patient,]
	maf = maf_ss[(maf_ss$Patient %in% tclin$Patient) & (maf_ss$Hugo_Symbol %in% common_genes),]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	ts_df = data.frame( Patient = maf_snp$Patient, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	save(ts_df, file = paste0(OutDir, "Genie_ts_df.RData"))
	load(file = paste0(OutDir, "Genie_ts_df.RData")) # ts_df
	sigs.input = mut.to.sigs.input( mut.ref = ts_df, sample.id = "Patient", chr = "chr", pos = "pos", ref = "ref", alt = "alt" )
	save(sigs.input, file = paste0(OutDir, "Genie_sigs_input.RData"))
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = "GENIE-DFCI-038224", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "Genie_Sign_weights_unnormalized.RData" ))


	### Now, same but for Cosmic signatures

	###### TCGA
	load(file = paste0(OutDir, "TCGA_sigs_input.RData"))
	# try both exome2genome and default (no normalization)
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = "TCGA-05-4402", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	 prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	 Sign_weights_df[rn,] = prova$weights
	    index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "Cosmic_TCGA_Sign_weights_unnormalized.RData" ))


	###### ChenEAS
	load(file = paste0(OutDir, "ChenEAS_sigs_input.RData"))
	# try both exome2genome and default (no normalization)
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = "BGI-EX02", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "Cosmic_ChenEAS_Sign_weights_unnormalized.RData" ))


	###### Genie
	load( file = paste0(OutDir, "Genie_sigs_input.RData"))
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = "GENIE-DFCI-038224", 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "Cosmic_Genie_Sign_weights_unnormalized.RData" ))


	###### Origimed
	load( file = paste0(OutDir, "Origimed_sigs_input.RData"))
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = rownames(sigs.input)[1], 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.cosmic, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(OutDir, "Cosmic_Origimed_Sign_weights_unnormalized.RData" ))
}

mutationalSignatures_associations = function( OutDir ){

	#### Downstream analyses
	this_OutDir = paste0(OutDir,"All_signatures_associations/")
	dir.create(this_OutDir)
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load(file = paste0(OutDir, "TCGA_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[substr(rownames(Sign_weights_df),1,12) %in% rownames(Clin),]
	rownames(tw) = substr(rownames(tw),1,12)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"TCGA_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	load(file = paste0(OutDir, "ChenEAS_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"ChenEAS_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	load(file = paste0(OutDir, "Genie_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"Genie_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### Origimed
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	load(file = paste0(OutDir, "Origimed_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"Origimed_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}

	#### Cosmic, Downstream analyses
	this_OutDir = paste0(OutDir,"All_signatures_associations_Cosmic/")
	dir.create(this_OutDir)
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	load(file = paste0(OutDir, "Cosmic_TCGA_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[substr(rownames(Sign_weights_df),1,12) %in% rownames(Clin),]
	rownames(tw) = substr(rownames(tw),1,12)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"Cosmic_TCGA_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	load(file = paste0(OutDir, "Cosmic_ChenEAS_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"Cosmic_ChenEAS_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	load(file = paste0(OutDir, "Cosmic_Genie_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"Cosmic_Genie_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	load(file = paste0(OutDir, "Cosmic_Origimed_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in colnames(tw)){
	   fileName = paste0(this_OutDir,"Cosmic_Origimed_EGFRclass_vs_MutSign",sign,".pdf")
	   x = Clin$egfr_class_consensus
	   y = Clin[,sign]
	   dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign," weights"), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = Clin$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
}

mutationalSignatures_egfr_muts_contexts = function( OutDir ){
	library(deconstructSigs)
	this_OutDir = paste0(OutDir, "egfr_muts_contexts/" )
	dir.create(this_OutDir,showWarnings=F)
	###### Genie
	gi = dan.read(paste0(DataDir,"genie_14.0/genomic_information.txt"))
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2
	load(file = paste0(DataDir,'genie_14.0/maf_Genie_LUAD_all.RData'))
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	maf_ss$Patient = Clin2[maf_ss$Tumor_Sample_Barcode,"Patient"]
	### let's find a gene universe that we can test
	assays = sort(unique(Clin2g$SEQ_ASSAY_ID))
	asdf = data.frame(row.names = assays, assay = assays, ngenes = NA, nsamples = NA, npatients = NA, stringsAsFactors=FALSE)
	for (a in assays){
	   gia = gi[gi$SEQ_ASSAY_ID==a,]
	   ug = unique(gia$Hugo_Symbol[nchar(gia$Hugo_Symbol)>0])
	   asdf[a,"ngenes"] = length(ug)
	   asdf[a,"nsamples"] = length(Clin2g[Clin2g$SEQ_ASSAY_ID==a,"Sample"])
	   asdf[a,"npatients"] = length(unique(Clin2g[Clin2g$SEQ_ASSAY_ID==a,"Patient"]))
	}
	asdf = asdf[order(asdf$ngenes,decreasing=T),]
	# Choice: remove panels with less than 100 genes
	asdf = asdf[asdf$ngenes>=100,]
	# asdf = asdf[asdf$npatients>=100,]
	asdf = asdf[order(asdf$nsamples,decreasing=T),]
	first = TRUE
	for (a in rownames(asdf)){
	   gia = gi[gi$SEQ_ASSAY_ID==a,]
	   ug = unique(gia$Hugo_Symbol[nchar(gia$Hugo_Symbol)>0])
	   if (first){
	      common_genes = ug
	      first = FALSE
	   } else {
	      common_genes = intersect(common_genes,ug)
	   }
	   dcat(a)
	   dcat(length(common_genes),1)
	}
	tClin2 = Clin2g[Clin2g$SEQ_ASSAY_ID %in% rownames(asdf),]
	tclin = Cling[rownames(Cling) %in% tClin2$Patient,]
	maf = maf_ss[(maf_ss$Patient %in% tclin$Patient) & (maf_ss$Hugo_Symbol %in% common_genes),]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	maf_snp = maf_snp[maf_snp$Hugo_Symbol=="EGFR",]
	cnn = colnames(tclin)[grepl( "egfr_mut_",colnames(tclin) )]
	all_muts = unlist(tclin[,cnn])
	all_muts = as.character(all_muts[!is.na(all_muts)])
	cnn = colnames(tclin)[(grepl( "egfr_class_",colnames(tclin) )) & (!grepl( "consensus",colnames(tclin) ))]
	all_classes = unlist(tclin[,cnn])
	all_classes = as.character(all_classes[!is.na(all_classes)])
	maf_snp = maf_snp[maf_snp$HGVSp_Short %in% all_muts,]
	all_classes = all_classes[all_muts %in% maf_snp$HGVSp_Short]
	all_muts = all_muts[all_muts %in% maf_snp$HGVSp_Short]
	maf_snp$egfr_class = "uncommon"
	maf_snp[maf_snp$HGVSp_Short=="p.L858R","egfr_class"] = "common"
	ts_df = data.frame( variant = maf_snp$HGVSp_Short,Patient = maf_snp$egfr_class, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	# ts_df = ts_df[!duplicated(paste0(ts_df$pos,"_",ts_df$ref,"_",ts_df$alt)),]
	ts_df_all = ts_df

	### Zhang
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Zhang_SampleLevel.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clinz = Clin
	maf = dan.read(paste0(DataDir,"lung_cancer_never_smokers_nci_2022/data_mutations.txt"))
	maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	rownames(Clin2) = Clin2$Patient
	Clin$Sample = Clin2[rownames(Clin),"Sample"]
	rownames(Clin) = Clin$Sample
	maf_ss$Patient = Clin[maf_ss$Tumor_Sample_Barcode,"Patient"]
	rownames(Clin) = Clin$Patient
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	maf_snp = maf_snp[maf_snp$Hugo_Symbol=="EGFR",]
	cnn = colnames(tclin)[grepl( "egfr_mut_",colnames(tclin) )]
	all_muts = unlist(tclin[,cnn])
	all_muts = as.character(all_muts[!is.na(all_muts)])
	cnn = colnames(tclin)[(grepl( "egfr_class_",colnames(tclin) )) & (!grepl( "consensus",colnames(tclin) ))]
	all_classes = unlist(tclin[,cnn])
	all_classes = as.character(all_classes[!is.na(all_classes)])
	maf_snp = maf_snp[maf_snp$HGVSp_Short %in% all_muts,]
	all_classes = all_classes[all_muts %in% maf_snp$HGVSp_Short]
	all_muts = all_muts[all_muts %in% maf_snp$HGVSp_Short]
	maf_snp$egfr_class = "uncommon"
	maf_snp[maf_snp$HGVSp_Short=="p.L858R","egfr_class"] = "common"
	ts_df = data.frame( variant = maf_snp$HGVSp_Short, Patient = maf_snp$egfr_class, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	# ts_df = ts_df[!duplicated(paste0(ts_df$pos,"_",ts_df$ref,"_",ts_df$alt)),]
	ts_df_all = rbind(ts_df_all,ts_df)

	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TRACERx421_SampleLevel.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clinz = Clin
	library(fst)
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$tumour_id
	maf_ss$Patient = maf_ss$patient_id
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$ref) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$var) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	maf_snp = maf_snp[maf_snp$Hugo_Symbol=="EGFR",]
	cnn = colnames(tclin)[grepl( "egfr_mut_",colnames(tclin) )]
	all_muts = unlist(tclin[,cnn])
	all_muts = as.character(all_muts[!is.na(all_muts)])
	cnn = colnames(tclin)[(grepl( "egfr_class_",colnames(tclin) )) & (!grepl( "consensus",colnames(tclin) ))]
	all_classes = unlist(tclin[,cnn])
	all_classes = as.character(all_classes[!is.na(all_classes)])
	maf_snp = maf_snp[maf_snp$AAChange %in% all_muts,]
	all_classes = all_classes[all_muts %in% maf_snp$AAChange]
	all_muts = all_muts[all_muts %in% maf_snp$AAChange]
	maf_snp$egfr_class = "uncommon"
	maf_snp[maf_snp$AAChange=="p.L858R","egfr_class"] = "common"
	ts_df = data.frame( variant = maf_snp$AAChange,Patient = maf_snp$egfr_class, chr = paste0(maf_snp$chr), pos = as.numeric(maf_snp$start), ref = maf_snp$ref, alt = maf_snp$var) # for all of them, end-start=0
	# ts_df = ts_df[!duplicated(paste0(ts_df$pos,"_",ts_df$ref,"_",ts_df$alt)),]
	ts_df_all = rbind(ts_df_all,ts_df)

	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clint = Clin
	load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
	maf_ss$Patient = substr(maf_ss$Tumor_Sample_Barcode,1,12)
	maf_ss$Sample = substr(maf_ss$Tumor_Sample_Barcode,1,16)
	tclin = Clint[Clint$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Match_Norm_Seq_Allele1) %in% c("A", "T", "C", "G")) & ((maf$Match_Norm_Seq_Allele2) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	maf_snp = maf_snp[maf_snp$Hugo_Symbol=="EGFR",]
	cnn = colnames(tclin)[grepl( "egfr_mut_",colnames(tclin) )]
	all_muts = unlist(tclin[,cnn])
	all_muts = as.character(all_muts[!is.na(all_muts)])
	cnn = colnames(tclin)[(grepl( "egfr_class_",colnames(tclin) )) & (!grepl( "consensus",colnames(tclin) ))]
	all_classes = unlist(tclin[,cnn])
	all_classes = as.character(all_classes[!is.na(all_classes)])
	maf_snp = maf_snp[maf_snp$HGVSp_Short %in% all_muts,]
	all_classes = all_classes[all_muts %in% maf_snp$HGVSp_Short]
	all_muts = all_muts[all_muts %in% maf_snp$HGVSp_Short]
	maf_snp$egfr_class = "uncommon"
	maf_snp[maf_snp$HGVSp_Short=="p.L858R","egfr_class"] = "common"
	ts_df = data.frame( variant = maf_snp$HGVSp_Short,Patient = maf_snp$egfr_class, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Match_Norm_Seq_Allele1, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	# ts_df = ts_df[!duplicated(paste0(ts_df$pos,"_",ts_df$ref,"_",ts_df$alt)),]
	ts_df_all = rbind(ts_df_all,ts_df)


	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clinc = Clin
	maf = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/snv_indel.maf"))
	maf$Patient = maf$Tumor_Sample_Barcode
	tclin = Clinc[Clinc$Patient %in% maf$Patient,]
	maf = maf[maf$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	maf_snp = maf_snp[maf_snp$Hugo_Symbol=="EGFR",]
	cnn = colnames(tclin)[grepl( "egfr_mut_",colnames(tclin) )]
	all_muts = unlist(tclin[,cnn])
	all_muts = as.character(all_muts[!is.na(all_muts)])
	cnn = colnames(tclin)[(grepl( "egfr_class_",colnames(tclin) )) & (!grepl( "consensus",colnames(tclin) ))]
	all_classes = unlist(tclin[,cnn])
	all_classes = as.character(all_classes[!is.na(all_classes)])
	maf_snp = maf_snp[maf_snp$HGVSp_Short %in% all_muts,]
	all_classes = all_classes[all_muts %in% maf_snp$HGVSp_Short]
	all_muts = all_muts[all_muts %in% maf_snp$HGVSp_Short]
	maf_snp$egfr_class = "uncommon"
	maf_snp[maf_snp$HGVSp_Short=="p.L858R","egfr_class"] = "common"
	ts_df = data.frame( variant = maf_snp$HGVSp_Short,Patient = maf_snp$egfr_class, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	# ts_df = ts_df[!duplicated(paste0(ts_df$pos,"_",ts_df$ref,"_",ts_df$alt)),]
	ts_df_all = rbind(ts_df_all,ts_df)

	### Origimed
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Clinz = Clin
	maf = dan.read(paste0(DataDir,"china_pan_origimed_2020/data_mutations_extended.txt"))
	maf_ss = maf[maf$Tumor_Sample_Barcode %in% Clin$Sample,]
	maf_ss$Sample = maf_ss$Tumor_Sample_Barcode
	rownames(Clin) = Clin$Sample
	maf_ss$Patient = Clin[maf_ss$Tumor_Sample_Barcode,"Patient"]
	rownames(Clin) = Clin$Patient
	tclin = Clinz[Clinz$Patient %in% maf_ss$Patient,]
	maf = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	condition_normal = ((maf$Reference_Allele) %in% c("A", "T", "C", "G"))
	condition_tumor = ((maf$Tumor_Seq_Allele2) %in% c("A", "T", "C", "G")) & ((maf$Tumor_Seq_Allele1) %in% c("A", "T", "C", "G"))
	maf_snp = maf[condition_normal & condition_tumor,]
	maf_snp = maf_snp[maf_snp$Hugo_Symbol=="EGFR",]
	cnn = colnames(tclin)[grepl( "egfr_mut_",colnames(tclin) )]
	all_muts = unlist(tclin[,cnn])
	all_muts = as.character(all_muts[!is.na(all_muts)])
	cnn = colnames(tclin)[(grepl( "egfr_class_",colnames(tclin) )) & (!grepl( "consensus",colnames(tclin) ))]
	all_classes = unlist(tclin[,cnn])
	all_classes = as.character(all_classes[!is.na(all_classes)])
	maf_snp = maf_snp[maf_snp$HGVSp_Short %in% all_muts,]
	all_classes = all_classes[all_muts %in% maf_snp$HGVSp_Short]
	all_muts = all_muts[all_muts %in% maf_snp$HGVSp_Short]
	maf_snp$egfr_class = "uncommon"
	maf_snp[maf_snp$HGVSp_Short=="p.L858R","egfr_class"] = "common"
	ts_df = data.frame( variant = maf_snp$HGVSp_Short,Patient = maf_snp$egfr_class, chr = paste0("chr" ,maf_snp$Chromosome), pos = maf_snp$Start_Position, ref = maf_snp$Reference_Allele, alt = maf_snp$Tumor_Seq_Allele2) # for all of them, end-start=0
	# ts_df = ts_df[!duplicated(paste0(ts_df$pos,"_",ts_df$ref,"_",ts_df$alt)),]
	ts_df_all = rbind(ts_df_all,ts_df)

	dim(ts_df_all)
	# ts_df_all = ts_df_all[!duplicated(paste0(ts_df_all$pos,"_",ts_df_all$ref,"_",ts_df_all$alt)),]

	save(ts_df_all, file = paste0(this_OutDir, "All_ts_df.RData"))

	load(file = paste0(this_OutDir, "All_ts_df.RData")) # ts_df
	sigs.input = mut.to.sigs.input( mut.ref = ts_df_all, sample.id = "Patient", chr = "chr", pos = "pos", ref = "ref", alt = "alt" )
	save(sigs.input, file = paste0(this_OutDir, "All_sigs_input.RData"))
	prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = ts_df_all$Patient[1], 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	Sign_weights_df = data.frame(matrix(nrow = nrow(sigs.input), ncol = length(prova$weights)))
	rownames(Sign_weights_df) = rownames(sigs.input)
	colnames(Sign_weights_df) = names(prova$weights)
	tot = length(rownames(Sign_weights_df))
	index = 1
	for (rn in rownames(Sign_weights_df) )
	{
	   if (index %in% c( round((1:9)*0.1*tot)))
	   {
	      cat("\n","...",index,"\n")
	   }
	   prova = whichSignatures(tumor.ref = sigs.input, 
	                           signatures.ref = signatures.nature2013, 
	                           sample.id = rn, 
	                           contexts.needed = TRUE,
	                           tri.counts.method = 'default')
	   Sign_weights_df[rn,] = prova$weights
	   index = index + 1
	}
	save(Sign_weights_df, file = paste0(this_OutDir, "All_Sign_weights_unnormalized.RData" ))
	load( file = paste0(this_OutDir, "All_Sign_weights_unnormalized.RData" ))
	Sign_weights_df = Sign_weights_df/rowSums(Sign_weights_df)
	ww = Sign_weights_df[,colSums(Sign_weights_df)>0]
	
	mt = data.frame( Signature=c( rep(colnames(ww),2) ), EGFR_class = rep(c( "common","uncommon" ),each=ncol(ww)), Percentage = NA, N_mutations = NA, stringsAsFactors=F  )
	mt[1,"N_mutations"] = paste0( "N = ",rowSums(sigs.input)["common"])
	mt[1+ncol(ww),"N_mutations"] = paste0( "N = ",rowSums(sigs.input)["uncommon"])
	for (sub in colnames(ww)){
		mt[(mt$Signature==sub) & (mt$EGFR_class=="common"),"Percentage"] = sum(ww["common",sub] )
		mt[(mt$Signature==sub) & (mt$EGFR_class=="uncommon"),"Percentage"] = sum(ww["uncommon",sub] )
	}
	mt[mt$EGFR_class=="common","Percentage"] = 100*mt[mt$EGFR_class=="common","Percentage"]/sum(mt[mt$EGFR_class=="common","Percentage"])
	mt[mt$EGFR_class=="uncommon","Percentage"] = 100*mt[mt$EGFR_class=="uncommon","Percentage"]/sum(mt[mt$EGFR_class=="uncommon","Percentage"])
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=Signature)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Normalized signature weight") + xlab("") + ggtitle("All datasets combined" ) + labs(fill = "Signature") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=N_mutations,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(paste0( this_OutDir,"common_vs_uncommon_signatures.pdf" ),5,5)
	print(pa)
	dev.off()

	pie = mt[mt$EGFR_class=="uncommon",]
	library(ggplot2)
	library(dplyr)
	pie$percent_rounded = round(pie$Percentage,1)
	pie$label = paste0( pie$percent_rounded,"%" )
	# Basic piechart
	plot = ggplot(pie, aes(x="", y=Percentage, fill=Signature)) +
	  geom_bar(stat="identity", width=1, color="white") +
	  coord_polar("y", start=0) +
	  theme_void() + 
	  geom_text(aes(x=1.25,label = label),
	            position = position_stack(vjust = 0.5),size=size_labels) 
	pdf( paste0(this_OutDir,"Pie_uncommonIndividual_signatures.pdf"),2.5,2.5 )
	plot=plot+theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()

	mt = data.frame( Substitution=c( "C>A","C>G","C>T","T>A","T>C","T>G","C>A","C>G","C>T","T>A","T>C","T>G" ), EGFR_class = rep(c( "common","uncommon" ),each=6), Percentage = NA, N_mutations = NA, stringsAsFactors=F  )
	mt[1,"N_mutations"] = paste0( "N = ",rowSums(sigs.input)["common"])
	mt[7,"N_mutations"] = paste0( "N = ",rowSums(sigs.input)["uncommon"])
	for (sub in c( "C>A","C>G","C>T","T>A","T>C","T>G" )){
		mt[(mt$Substitution==sub) & (mt$EGFR_class=="common"),"Percentage"] = sum(sigs.input["common",grepl( sub,colnames(sigs.input) )] )
		mt[(mt$Substitution==sub) & (mt$EGFR_class=="uncommon"),"Percentage"] = sum(sigs.input["uncommon",grepl( sub,colnames(sigs.input) )] )
	}
	mt[mt$EGFR_class=="common","Percentage"] = 100*mt[mt$EGFR_class=="common","Percentage"]/sum(mt[mt$EGFR_class=="common","Percentage"])
	mt[mt$EGFR_class=="uncommon","Percentage"] = 100*mt[mt$EGFR_class=="uncommon","Percentage"]/sum(mt[mt$EGFR_class=="uncommon","Percentage"])
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=Substitution)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of substitutions") + xlab("") + ggtitle("All datasets combined" ) + labs(fill = "Substitution") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=N_mutations,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(paste0( this_OutDir,"common_vs_uncommon_substitutions.pdf" ),5,5)
	print(pa)
	dev.off()

	pie = mt[mt$EGFR_class=="uncommon",]
	library(ggplot2)
	library(dplyr)
	pie$percent_rounded = round(pie$Percentage,1)
	pie$label = paste0( pie$percent_rounded,"%" )
	pie$Substitution = factor(pie$Substitution,levels = c( "T>A","T>C","T>G","C>A","C>G","C>T" )) # to match colors
	# Basic piechart
	plot = ggplot(pie, aes(x="", y=Percentage, fill=Substitution)) +
	  geom_bar(stat="identity", width=1, color="white") +
	  coord_polar("y", start=0) +
	  theme_void() + 
	  geom_text(aes(x=1.25,label = label),
	            position = position_stack(vjust = 0.5),size=size_labels) 
	pdf( paste0(this_OutDir,"Pie_uncommonIndividual_Substitutions.pdf"),2.5,2.5 )
	plot=plot+theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(plot)
	dev.off()
}

mutationalSignatures_associations_binary = function( OutDir, threshold = 0 ){

	#### Downstream analyses
	this_OutDir = paste0(OutDir,"All_signatures_associations_binary_threshold",threshold,"/")
	dir.create(this_OutDir)
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load(file = paste0(OutDir, "TCGA_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[substr(rownames(Sign_weights_df),1,12) %in% rownames(Clin),]
	rownames(tw) = substr(rownames(tw),1,12)
	tw = tw/rowSums(tw)
	tw = 1*(tw>threshold)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	tc = Clin
	Clin$Dataset = "TCGA"
	Clin_all = Clin
	tcga_pa = list()
	reversed_tcga_pa = list()
	for (sign in colnames(tw)){
		fileName = paste0(this_OutDir,"TCGA_EGFRclass_vs_MutSign",sign,".pdf")
		x = Clin$egfr_class_consensus
		y = Clin[,sign]
		tabb = table(tc[,sign],tc$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==0))
		mt[(mt[,"Signature"]==1) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==1))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=EGFR_class)) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TCGA, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(fileName,5,5)
		print(pa)
		dev.off()
		tcga_pa[[sign]] = pa

		x = Clin[,sign]
		y = Clin$egfr_class_consensus
		tabb = table(x,y)
		ch = chisq.test(tabb)$p.value
		tabb = (apply(tabb,2, function(x) x/sum(x)))*100
		mt = data.frame(Signature=rep(c(0,1),length(unique(Clin[,"egfr_class_consensus"]))),EGFR_class=rep( unique(Clin[,"egfr_class_consensus"]),each=2 ),Percentage=NA,stringsAsFactors=F)
		for (sval in c(0,1)){
			for (tclass in unique( mt$EGFR_class )){
				mt[(mt$Signature==sval) & ( mt$EGFR_class==tclass ),"Percentage"] = 100*sum((Clin[,sign]==sval) & (Clin[,"egfr_class_consensus"]==tclass))/sum(Clin[,"egfr_class_consensus"]==tclass)
			}
		}
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="common"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="uncommon"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="compound"))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=as.character(Signature))) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0(unique(Clin$Dataset),", ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) + 
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		reversed_tcga_pa[[sign]] = pall

		if (sign=="Signature.4"){
			Clin = Clin[(Clin$Smoking %in% c( "ever smoker","never smoker" )) %in% c(T),]
			fileName = paste0(this_OutDir,"TCGA_smoking_vs_MutSign",sign,".pdf")
			x = Clin$Smoking
			y = Clin[,sign]
			tabb = table(Clin[,sign],Clin$Smoking)
			ch = chisq.test(tabb)$p.value
			tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
			mt = melt(tabb)
			colnames(mt) = c( "Signature","Smoking","Percentage" )
			mt$TotalPatients = NA
			mt[(mt[,"Signature"]==0) & (mt$Smoking=="ever smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==0))
			mt[(mt[,"Signature"]==1) & (mt$Smoking=="never smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==1))
			pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=Smoking)) +
				geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TCGA, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + labs(fill = "Smoking") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
				geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
			pdf(fileName,5,5)
			print(pa)
			dev.off()
		}

	}
	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	load(file = paste0(OutDir, "ChenEAS_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	tw = tw/rowSums(tw)
	tw = 1*(tw>threshold)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	## Barplot: association between egfr_class_consensus and signature N yes/no
	tc = Clin
	Clin$Dataset = "ChenEAS"
	commonz = intersect(colnames(Clin_all),colnames(Clin))
	Clin_all = rbind(Clin_all[,commonz],Clin[,commonz])
	chen_pa = list()
	reversed_chen_pa = list()
	for (sign in colnames(tw)){
		fileName = paste0(this_OutDir,"ChenEAS_EGFRclass_vs_MutSign",sign,".pdf")
		x = Clin$egfr_class_consensus
		y = Clin[,sign]
		tabb = table(tc[,sign],tc$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==0))
		mt[(mt[,"Signature"]==1) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==1))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=EGFR_class)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("ChenEAS, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(fileName,5,5)
		print(pa)
		dev.off()
		chen_pa[[sign]] = pa

		x = Clin[,sign]
		y = Clin$egfr_class_consensus
		tabb = table(x,y)
		ch = chisq.test(tabb)$p.value
		tabb = (apply(tabb,2, function(x) x/sum(x)))*100
		mt = data.frame(Signature=rep(c(0,1),length(unique(Clin[,"egfr_class_consensus"]))),EGFR_class=rep( unique(Clin[,"egfr_class_consensus"]),each=2 ),Percentage=NA,stringsAsFactors=F)
		for (sval in c(0,1)){
			for (tclass in unique( mt$EGFR_class )){
				mt[(mt$Signature==sval) & ( mt$EGFR_class==tclass ),"Percentage"] = 100*sum((Clin[,sign]==sval) & (Clin[,"egfr_class_consensus"]==tclass))/sum(Clin[,"egfr_class_consensus"]==tclass)
			}
		}
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="common"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="uncommon"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="compound"))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=as.character(Signature))) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0(unique(Clin$Dataset),", ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) + 
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		reversed_chen_pa[[sign]] = pall

		if (sign=="Signature.4"){
			Clin = Clin[(Clin$Smoking %in% c( "ever smoker","never smoker" )) %in% c(T),]
			fileName = paste0(this_OutDir,"ChenEAS_smoking_vs_MutSign",sign,".pdf")
			x = Clin$Smoking
			y = Clin[,sign]
			tabb = table(Clin[,sign],Clin$Smoking)
			ch = chisq.test(tabb)$p.value
			tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
			mt = melt(tabb)
			colnames(mt) = c( "Signature","Smoking","Percentage" )
			mt$TotalPatients = NA
			mt[(mt[,"Signature"]==0) & (mt$Smoking=="ever smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==0))
			mt[(mt[,"Signature"]==1) & (mt$Smoking=="never smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==1))
			pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=Smoking)) +
				geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("ChenEAS, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + labs(fill = "Smoking") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
				geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
			pdf(fileName,5,5)
			print(pa)
			dev.off()
		}
	}
	###### Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	load(file = paste0(OutDir, "Genie_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	tw = tw/rowSums(tw)
	tw = 1*(tw>threshold)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	## Barplot: association between egfr_class_consensus and signature N yes/no
	tc = Clin
	Clin$Dataset = "Genie"
	commonz = intersect(colnames(Clin_all),colnames(Clin))
	Clin_all = rbind(Clin_all[,commonz],Clin[,commonz])
	genie_pa = list()
	reversed_genie_pa = list()
	for (sign in colnames(tw)){
		fileName = paste0(this_OutDir,"Genie_EGFRclass_vs_MutSign",sign,".pdf")
		x = Clin$egfr_class_consensus
		y = Clin[,sign]
		tabb = table(tc[,sign],tc$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==0))
		mt[(mt[,"Signature"]==1) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==1))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=EGFR_class)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Genie, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(fileName,5,5)
		print(pa)
		dev.off()
		genie_pa[[sign]] = pa

		x = Clin[,sign]
		y = Clin$egfr_class_consensus
		tabb = table(x,y)
		ch = chisq.test(tabb)$p.value
		tabb = (apply(tabb,2, function(x) x/sum(x)))*100
		mt = data.frame(Signature=rep(c(0,1),length(unique(Clin[,"egfr_class_consensus"]))),EGFR_class=rep( unique(Clin[,"egfr_class_consensus"]),each=2 ),Percentage=NA,stringsAsFactors=F)
		for (sval in c(0,1)){
			for (tclass in unique( mt$EGFR_class )){
				mt[(mt$Signature==sval) & ( mt$EGFR_class==tclass ),"Percentage"] = 100*sum((Clin[,sign]==sval) & (Clin[,"egfr_class_consensus"]==tclass))/sum(Clin[,"egfr_class_consensus"]==tclass)
			}
		}
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="common"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="uncommon"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="compound"))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=as.character(Signature))) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0(unique(Clin$Dataset),", ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) + 
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		reversed_genie_pa[[sign]] = pall

		if (sign=="Signature.4"){
			Clin = Clin[(Clin$Smoking %in% c( "ever smoker","never smoker" )) %in% c(T),]
			fileName = paste0(this_OutDir,"Genie_smoking_vs_MutSign",sign,".pdf")
			x = Clin$Smoking
			y = Clin[,sign]
			tabb = table(Clin[,sign],Clin$Smoking)
			ch = chisq.test(tabb)$p.value
			tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
			mt = melt(tabb)
			colnames(mt) = c( "Signature","Smoking","Percentage" )
			mt$TotalPatients = NA
			mt[(mt[,"Signature"]==0) & (mt$Smoking=="ever smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==0))
			mt[(mt[,"Signature"]==1) & (mt$Smoking=="never smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==1))
			pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=Smoking)) +
				geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Genie, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + labs(fill = "Smoking") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
				geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
			pdf(fileName,5,5)
			print(pa)
			dev.off()
		}
	}
	###### Origimed
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	load(file = paste0(OutDir, "Origimed_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	tw = tw/rowSums(tw)
	tw = 1*(tw>threshold)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	## Barplot: association between egfr_class_consensus and signature N yes/no
	tc = Clin
	Clin$Dataset = "Origimed"
	commonz = intersect(colnames(Clin_all),colnames(Clin))
	Clin_all = rbind(Clin_all[,commonz],Clin[,commonz])
	origimed_pa = list()
	reversed_origimed_pa = list()
	for (sign in colnames(tw)){
		fileName = paste0(this_OutDir,"Origimed_EGFRclass_vs_MutSign",sign,".pdf")
		x = Clin$egfr_class_consensus
		y = Clin[,sign]
		tabb = table(tc[,sign],tc$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==0))
		mt[(mt[,"Signature"]==1) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==1))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=EGFR_class)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Origimed, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(fileName,5,5)
		print(pa)
		dev.off()
		origimed_pa[[sign]] = pa

		x = Clin[,sign]
		y = Clin$egfr_class_consensus
		tabb = table(x,y)
		ch = chisq.test(tabb)$p.value
		tabb = (apply(tabb,2, function(x) x/sum(x)))*100
		mt = data.frame(Signature=rep(c(0,1),length(unique(Clin[,"egfr_class_consensus"]))),EGFR_class=rep( unique(Clin[,"egfr_class_consensus"]),each=2 ),Percentage=NA,stringsAsFactors=F)
		for (sval in c(0,1)){
			for (tclass in unique( mt$EGFR_class )){
				mt[(mt$Signature==sval) & ( mt$EGFR_class==tclass ),"Percentage"] = 100*sum((Clin[,sign]==sval) & (Clin[,"egfr_class_consensus"]==tclass))/sum(Clin[,"egfr_class_consensus"]==tclass)
			}
		}
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="common"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="uncommon"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="compound"))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=as.character(Signature))) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0(unique(Clin$Dataset),", ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) + 
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		reversed_origimed_pa[[sign]] = pall

		if (sign=="Signature.4"){
			Clin = Clin[(Clin$Smoking %in% c( "ever smoker","never smoker" )) %in% c(T),]
			fileName = paste0(this_OutDir,"Origimed_smoking_vs_MutSign",sign,".pdf")
			x = Clin$Smoking
			y = Clin[,sign]
			tabb = table(Clin[,sign],Clin$Smoking)
			ch = chisq.test(tabb)$p.value
			tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
			mt = melt(tabb)
			colnames(mt) = c( "Signature","Smoking","Percentage" )
			mt$TotalPatients = NA
			mt[(mt[,"Signature"]==0) & (mt$Smoking=="ever smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==0))
			mt[(mt[,"Signature"]==1) & (mt$Smoking=="never smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==1))
			pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=Smoking)) +
				geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Origimed, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + labs(fill = "Smoking") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
				geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
			pdf(fileName,5,5)
			print(pa)
			dev.off()
		}
	}
	###### Zhang
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	load(file = paste0(OutDir, "Zhang_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	tw = tw/rowSums(tw)
	tw = 1*(tw>threshold)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	## Barplot: association between egfr_class_consensus and signature N yes/no
	tc = Clin
	Clin$Dataset = "Zhang"
	commonz = intersect(colnames(Clin_all),colnames(Clin))
	Clin_all = rbind(Clin_all[,commonz],Clin[,commonz])
	zhang_pa = list()
	reversed_zhang_pa = list()
	for (sign in colnames(tw)){
		fileName = paste0(this_OutDir,"Zhang_EGFRclass_vs_MutSign",sign,".pdf")
		x = Clin$egfr_class_consensus
		y = Clin[,sign]
		tabb = table(tc[,sign],tc$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==0))
		mt[(mt[,"Signature"]==1) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==1))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=EGFR_class)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("Zhang, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(fileName,5,5)
		print(pa)
		dev.off()
		zhang_pa[[sign]] = pa

		x = Clin[,sign]
		y = Clin$egfr_class_consensus
		tabb = table(x,y)
		ch = chisq.test(tabb)$p.value
		tabb = (apply(tabb,2, function(x) x/sum(x)))*100
		mt = data.frame(Signature=rep(c(0,1),length(unique(Clin[,"egfr_class_consensus"]))),EGFR_class=rep( unique(Clin[,"egfr_class_consensus"]),each=2 ),Percentage=NA,stringsAsFactors=F)
		for (sval in c(0,1)){
			for (tclass in unique( mt$EGFR_class )){
				mt[(mt$Signature==sval) & ( mt$EGFR_class==tclass ),"Percentage"] = 100*sum((Clin[,sign]==sval) & (Clin[,"egfr_class_consensus"]==tclass))/sum(Clin[,"egfr_class_consensus"]==tclass)
			}
		}
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="common"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="uncommon"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="compound"))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=as.character(Signature))) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0(unique(Clin$Dataset),", ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) + 
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		reversed_zhang_pa[[sign]] = pall
	}
	###### TRACERx421
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	load(file = paste0(OutDir, "TRACERx421_Sign_weights_unnormalized.RData" ))
	Clin = Clin[rownames(Clin) %in% rownames(Sign_weights_df),]
	tw = Sign_weights_df[rownames(Sign_weights_df) %in% rownames(Clin),]
	tw = tw/rowSums(tw)
	tw = 1*(tw>threshold)
	Clin = cbind(Clin,tw[rownames(Clin),])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	## Barplot: association between egfr_class_consensus and signature N yes/no
	tc = Clin
	Clin$Dataset = "TRACERx421"
	commonz = intersect(colnames(Clin_all),colnames(Clin))
	Clin_all = rbind(Clin_all[,commonz],Clin[,commonz])
	tracerx_pa = list()
	reversed_tracerx_pa = list()
	for (sign in colnames(tw)){
		fileName = paste0(this_OutDir,"TRACERx421_EGFRclass_vs_MutSign",sign,".pdf")
		x = Clin$egfr_class_consensus
		y = Clin[,sign]
		tabb = table(tc[,sign],tc$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==0))
		mt[(mt[,"Signature"]==1) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(tc[,sign]==1))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=EGFR_class)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TRACERx421, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(fileName,5,5)
		print(pa)
		dev.off()
		tracerx_pa[[sign]] = pa

		x = Clin[,sign]
		y = Clin$egfr_class_consensus
		tabb = table(x,y)
		ch = chisq.test(tabb)$p.value
		tabb = (apply(tabb,2, function(x) x/sum(x)))*100
		mt = data.frame(Signature=rep(c(0,1),length(unique(Clin[,"egfr_class_consensus"]))),EGFR_class=rep( unique(Clin[,"egfr_class_consensus"]),each=2 ),Percentage=NA,stringsAsFactors=F)
		for (sval in c(0,1)){
			for (tclass in unique( mt$EGFR_class )){
				mt[(mt$Signature==sval) & ( mt$EGFR_class==tclass ),"Percentage"] = 100*sum((Clin[,sign]==sval) & (Clin[,"egfr_class_consensus"]==tclass))/sum(Clin[,"egfr_class_consensus"]==tclass)
			}
		}
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="common"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="uncommon"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin[,"egfr_class_consensus"]=="compound"))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=as.character(Signature))) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0(unique(Clin$Dataset),", ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) + 
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		reversed_tracerx_pa[[sign]] = pall

		if (sign=="Signature.4"){
			Clin = Clin[(Clin$Smoking %in% c( "ever smoker","never smoker" )) %in% c(T),]
			fileName = paste0(this_OutDir,"TRACERx421_smoking_vs_MutSign",sign,".pdf")
			x = Clin$Smoking
			y = Clin[,sign]
			tabb = table(Clin[,sign],Clin$Smoking)
			ch = chisq.test(tabb)$p.value
			tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
			mt = melt(tabb)
			colnames(mt) = c( "Signature","Smoking","Percentage" )
			mt$TotalPatients = NA
			mt[(mt[,"Signature"]==0) & (mt$Smoking=="ever smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==0))
			mt[(mt[,"Signature"]==1) & (mt$Smoking=="never smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin[,sign]==1))
			pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=Smoking)) +
				geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("TRACERx421, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + labs(fill = "Smoking") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
				geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
			pdf(fileName,5,5)
			print(pa)
			dev.off()
		}
	}
	#### across datasets
	for (sign in colnames(tw)){
		x = Clin_all$egfr_class_consensus
		y = Clin_all[,sign]
		tabb = table(Clin_all[,sign],Clin_all$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,sign]==0))
		mt[(mt[,"Signature"]==1) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,sign]==1))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=EGFR_class)) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = size_labels) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(paste0(this_OutDir,"All_EGFRclass_vs_MutSign_",sign,".pdf"),2.5,2.5,onefile=FALSE)
		pall=pall+theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		# plot=ggarrange(tcga_pa[[sign]], chen_pa[[sign]], genie_pa[[sign]], origimed_pa[[sign]], zhang_pa[[sign]], tracerx_pa[[sign]], pall, common.legend=TRUE,ncol = 7, nrow = 1)
		print(pall)
		dev.off()
	}
	#### across datasets, but the other way round
	for (sign in colnames(tw)){
		x = Clin_all[,sign]
		y = Clin_all$egfr_class_consensus
		tabb = table(x,y)
		ch = chisq.test(tabb)$p.value
		tabb = (apply(tabb,2, function(x) x/sum(x)))*100
		mt = melt(tabb)
		colnames(mt) = c( "Signature","EGFR_class","Percentage" )
		mt$TotalPatients = NA
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,"egfr_class_consensus"]=="common"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,"egfr_class_consensus"]=="uncommon"))
		mt[(mt[,"Signature"]==0) & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,"egfr_class_consensus"]=="compound"))
		mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
		pall = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=as.character(Signature))) +
			geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) + 
			geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
		pdf(paste0(this_OutDir,"reversed_All_EGFRclass_vs_MutSign_",sign,".pdf"),22,5,onefile=FALSE)
		plot=ggarrange(reversed_tcga_pa[[sign]], reversed_chen_pa[[sign]], reversed_genie_pa[[sign]], reversed_origimed_pa[[sign]], reversed_zhang_pa[[sign]], reversed_tracerx_pa[[sign]], pall, common.legend=TRUE,ncol = 7, nrow = 1)
		print(plot)
		dev.off()
	}

	fileName = paste0(this_OutDir,"UncommonVariant_enriched_inSignature4_volcanoPlot.pdf")
	cu = Clin_all[Clin_all$egfr_class_consensus=="uncommon",]
	tabb = table(cu[,"Signature.4"],cu$egfr_mut_1)
	tabb = tabb[,colSums(tabb)>=3]
	cdf = dan.df(colnames(tabb),c( "Variant","Chisq_residual_smoker","Chisq_pval" ))
	for (g in colnames(tabb)){
		cdf[g,"Variant"] = g
		tabb2 = table(cu[,"Signature.4"],cu$egfr_mut_1==g)
		ch = chisq.test(tabb2)
		cdf[g,"Chisq_pval"] = ch$p.value
		cdf[g,"Chisq_residual_smoker"] = ch$residuals["1","TRUE"]
	}
	cdf = cdf[order(cdf$Chisq_pval),]
	cdf$Chisq_qval = p.adjust(cdf$Chisq_pval,method="BH")
	cdf$fill = ifelse( (cdf$Chisq_qval>=0.1) | (abs(cdf$Chisq_residual_smoker)<=1),"Not significant",ifelse( cdf$Chisq_residual_smoker<0,"Higher in non-smokers","Higher in smokers"))
	cdf$fill = factor(cdf$fill,levels=c( "Higher in non-smokers","Higher in smokers","Not significant" ))
	fillColorz = c( "gray55","orange","gray88" )
	cdf$repel_labelz = rownames(cdf)
	cdf[cdf$fill=="Not significant","repel_labelz"] = ""
	dan.scatterplot( fileName, cdf$Chisq_residual_smoker, -log10(cdf$Chisq_pval), fill = cdf$fill, xlab = "Chi-square residual in smokers", ylab = "-log10(Chi-square p-value)", filllab = "", fillColors = fillColorz, dotSize = 3, repel_labels=cdf$repel_labelz, coord_fixed = FALSE, fileWidth = 7, fileHeight = 4 )

	save(Clin_all,file=paste0(this_OutDir,"Clin_all.RData"))
	load(file=paste0(this_OutDir,"Clin_all.RData"))

	sign = "Signature.4"
	Clin_all = Clin_all[(Clin_all$Smoking %in% c( "ever smoker","never smoker" )) %in% c(T),]
	fileName = paste0(this_OutDir,"AllDatasets_smoking_vs_MutSign",sign,".pdf")
	x = Clin_all$Smoking
	y = Clin_all[,sign]
	tabb = table(Clin_all[,sign],Clin_all$Smoking)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Signature","Smoking","Percentage" )
	mt$TotalPatients = NA
	mt[(mt[,"Signature"]==0) & (mt$Smoking=="ever smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,sign]==0))
	mt[(mt[,"Signature"]==1) & (mt$Smoking=="never smoker"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,sign]==1))
	pa = ggplot(data=mt, aes(x=as.character(Signature), y=Percentage, fill=Smoking)) +
		geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + labs(fill = "Smoking") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
		geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(fileName,5,5)
	print(pa)
	dev.off()

	# Clin_all$Confirmed_smoker = NA
	# Clin_all[(Clin_all$Smoking=="ever smoker") & (Clin_all[,sign]==1),"Confirmed_smoker"] = "Ever smoker and\n active signature 4"
	# Clin_all[(Clin_all$Smoking=="never smoker") & (Clin_all[,sign]==0),"Confirmed_smoker"] = "Never smoker and\n inactive signature 4"
	# Clin_all = Clin_all[!is.na(Clin_all$Confirmed_smoker),]

	# fileName = paste0(this_OutDir,"AllDatasets_EGFRclass_vs_ConfirmedSmoking.pdf")
	# x = Clin_all$egfr_class_consensus
	# y = Clin_all$Confirmed_smoker
	# tabb = table(Clin_all$Confirmed_smoker,Clin_all$egfr_class_consensus)
	# ch = chisq.test(tabb)$p.value
	# tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	# mt = melt(tabb)
	# colnames(mt) = c( "Signature","EGFR_class","Percentage" )
	# mt$TotalPatients = NA
	# mt[(mt[,"Signature"]=="Never smoker and\n inactive signature 4") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all$Confirmed_smoker=="Never smoker and\n inactive signature 4"))
	# mt[(mt[,"Signature"]=="Ever smoker and\n active signature 4") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all$Confirmed_smoker=="Ever smoker and\n active signature 4"))
	# mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	# mt$Signature = factor(mt$Signature,levels=c( "Never smoker and\n inactive signature 4","Ever smoker and\n active signature 4" ))
	# pa = ggplot(data=mt, aes(x=Signature, y=Percentage, fill=EGFR_class)) +
	# geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	# geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	# pdf(fileName,5,5)
	# print(pa)
	# dev.off()

	# fileName = paste0(this_OutDir,"reversed_AllDatasets_EGFRclass_vs_ConfirmedSmoking.pdf")
	# tabb = table(Clin_all$Confirmed_smoker,Clin_all$egfr_class_consensus)
	# ch = chisq.test(tabb)$p.value
	# tabb = (apply(tabb,2, function(x) x/sum(x)))*100
	# mt = melt(tabb)
	# colnames(mt) = c( "Signature","EGFR_class","Percentage" )
	# mt$TotalPatients = NA
	# mt[(mt[,"Signature"]=="Ever smoker and\n active signature 4") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,"egfr_class_consensus"]=="common"))
	# mt[(mt[,"Signature"]=="Ever smoker and\n active signature 4") & (mt$EGFR_class=="uncommon"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,"egfr_class_consensus"]=="uncommon"))
	# mt[(mt[,"Signature"]=="Ever smoker and\n active signature 4") & (mt$EGFR_class=="compound"),"TotalPatients"] = paste0( "N = ",sum(Clin_all[,"egfr_class_consensus"]=="compound"))
	# mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	# mt$Signature = factor(mt$Signature,levels=c( "Never smoker and\n inactive signature 4","Ever smoker and\n active signature 4" ))
	# pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=Signature)) +
	# geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "Signature",values=c( "gray55","orange" )) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	# geom_text(aes(label=TotalPatients,y=100),vjust=-0.2) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	# pdf(fileName,5,5)
	# print(pa)
	# dev.off()

	# fileName = paste0(this_OutDir,"UncommonVariant_enriched_inConfirmedSmokers_volcanoPlot.pdf")
	# cu = Clin_all[Clin_all$egfr_class_consensus=="uncommon",]
	# tabb = table(cu[,"Confirmed_smoker"],cu$egfr_mut_1)
	# tabb = tabb[,colSums(tabb)>=3]
	# cdf = dan.df(colnames(tabb),c( "Variant","Chisq_residual_smoker","Chisq_pval" ))
	# for (g in colnames(tabb)){
	# 	cdf[g,"Variant"] = g
	# 	tabb2 = table(cu[,"Confirmed_smoker"],cu$egfr_mut_1==g)
	# 	ch = chisq.test(tabb2)
	# 	cdf[g,"Chisq_pval"] = ch$p.value
	# 	cdf[g,"Chisq_residual_smoker"] = ch$residuals["Ever smoker and\n active signature 4","TRUE"]
	# }
	# cdf = cdf[order(cdf$Chisq_pval),]
	# cdf$fill = ifelse( (cdf$Chisq_pval>=0.05) | (abs(cdf$Chisq_residual_smoker)<=1),"Not significant",ifelse( cdf$Chisq_residual_smoker<0,"Higher in non-smokers","Higher in smokers"))
	# cdf$fill = factor(cdf$fill,levels=c( "Higher in non-smokers","Higher in smokers","Not significant" ))
	# fillColorz = c( "gray55","orange","gray88" )
	# cdf$repel_labelz = rownames(cdf)
	# cdf[cdf$fill=="Not significant","repel_labelz"] = ""
	# dan.scatterplot( fileName, cdf$Chisq_residual_smoker, -log10(cdf$Chisq_pval), fill = cdf$fill, xlab = "Chi-square residual in smokers", ylab = "-log10(Chi-square p-value)", filllab = "", fillColors = fillColorz, dotSize = 3, repel_labels=cdf$repel_labelz, coord_fixed = FALSE, fileWidth = 7, fileHeight = 7 )

	Clin_all$Confirmed_smoker = NA
	Clin_all[(Clin_all$Smoking=="ever smoker") & (Clin_all[,sign]==1),"Confirmed_smoker"] = "Ever smoker and\n active signature 4"
	Clin_all[(Clin_all$Smoking=="ever smoker") & (Clin_all[,sign]==0),"Confirmed_smoker"] = "Ever smoker but\n inactive signature 4"
	Clin_all[(Clin_all$Smoking=="never smoker") & (Clin_all[,sign]==1),"Confirmed_smoker"] = "Never smoker but\n active signature 4"
	Clin_all[(Clin_all$Smoking=="never smoker") & (Clin_all[,sign]==0),"Confirmed_smoker"] = "Never smoker and\n inactive signature 4"
	Clin_all = Clin_all[!is.na(Clin_all$Confirmed_smoker),]

	fileName = paste0(this_OutDir,"AllDatasets_EGFRclass_vs_ConfirmedSmoking_extended.pdf")
	x = Clin_all$egfr_class_consensus
	y = Clin_all$Confirmed_smoker
	tabb = table(Clin_all$Confirmed_smoker,Clin_all$egfr_class_consensus)
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "Signature","EGFR_class","Percentage" )
	mt$TotalPatients = NA
	mt[(mt[,"Signature"]=="Never smoker and\n inactive signature 4") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all$Confirmed_smoker=="Never smoker and\n inactive signature 4"))
	mt[(mt[,"Signature"]=="Ever smoker and\n active signature 4") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all$Confirmed_smoker=="Ever smoker and\n active signature 4"))
	mt[(mt[,"Signature"]=="Never smoker but\n active signature 4") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all$Confirmed_smoker=="Never smoker but\n active signature 4"))
	mt[(mt[,"Signature"]=="Ever smoker but\n inactive signature 4") & (mt$EGFR_class=="common"),"TotalPatients"] = paste0( "N = ",sum(Clin_all$Confirmed_smoker=="Ever smoker but\n inactive signature 4"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$Signature = factor(mt$Signature,levels=c( "Never smoker and\n inactive signature 4","Ever smoker but\n inactive signature 4","Never smoker but\n active signature 4","Ever smoker and\n active signature 4" ))
	pa = ggplot(data=mt, aes(x=Signature, y=Percentage, fill=EGFR_class)) +
	geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("All datasets, ",sign,"\n","Chi-square test, p = ",signif(ch,2)) ) + scale_fill_manual(name = "EGFR mutation class",values=colorz_classes) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	geom_text(aes(label=TotalPatients,y=100),vjust=-0.2,size = size_labels) + geom_text(aes(label=round(Percentage,1)), size = size_labels, position = position_stack(vjust = 0.5))
	pdf(fileName,3,2.8)
	pa=pa+theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(pa)
	dev.off()
}

CnaSv_associations = function( OutDir ){

	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	Cling = Clin
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	Clin2g = Clin2

	load(file = paste0( DataDir,"genie_14.0/cna_Genie_EGFRmutants.RData" ))
	dtable(unlist(cna))
	cna = cna[ rowSums((cna==(-2)) | (cna==(2)), na.rm=T )>=5, ] # test only genes altered in at least 5 patients
	tclin = Clin2[intersect(rownames(Clin2), colnames(cna)),]
	cna = cna[,intersect(rownames(Clin2), colnames(cna))]
	tclin = cbind(tclin, t(cna) )
	cdf = dan.df(rownames(cna),c( "Gene","Chisq_pval" ))
	for (g in rownames(cdf) ){
		tclin[(tclin[,g] %in% c( -1.5,-1,0,1 )) %in% c(T),g ] = 0
		tclin$egfr_class_consensus = factor(tclin$egfr_class_consensus,levels=ordered_classes)
		tabb = table(tclin[,g],tclin$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		cdf[g,"Gene"] = g
		cdf[g,"Chisq_pval"] = ch
	}
	cdf$Chisq_qval = p.adjust(cdf$Chisq_pval,method="BH")
	cdf = cdf[order(cdf$Chisq_qval),]
	save(cdf, file = paste0( OutDir,"Genie_coCNA.RData" ))
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TCGA_SampleLevel.RData" ))
	load(file = paste0( CommonDataDir,"cna_LUAD_all.RData" ))
	cna = t(cna_ss)
	cna = cna[,colnames(cna) %in% rownames(Clin)]
	cna = cna[ rowSums((cna==(-2)) | (cna==(2)), na.rm=T )>=5, ] # test only genes altered in at least 5 patients
	tclin = Clin[rownames(Clin) %in% colnames(cna),]
	tclin = cbind(tclin, t(cna) )
	cdf = dan.df(rownames(cna),c( "Gene","Chisq_pval" ))
	for (g in rownames(cdf) ){
		tclin[(tclin[,g] %in% c( -1.5,-1,0,1 )) %in% c(T),g ] = 0
		tclin$egfr_class_consensus = factor(tclin$egfr_class_consensus,levels=ordered_classes)
		tabb = table(tclin[,g],tclin$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		cdf[g,"Gene"] = g
		cdf[g,"Chisq_pval"] = ch
	}
	cdf$Chisq_qval = p.adjust(cdf$Chisq_pval,method="BH")
	cdf = cdf[order(cdf$Chisq_qval),]
	save(cdf, file = paste0( OutDir,"TCGA_coCNA.RData" ))
	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	cna = dan.read( "/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/cnv.tsv" )
	cna = cna[!duplicated(cna$Hugo_Symbol),]
	rownames(cna) = cna$Hugo_Symbol
	cna$Hugo_Symbol = NULL
	cna$Entrez_Gene_Id = NULL
	colnames(cna) = gsub("\\.","-",colnames(cna))
	cna = cna[,colnames(cna) %in% rownames(Clin)]
	cna = cna[ rowSums((cna==(-2)) | (cna==(2)), na.rm=T )>=5, ] # test only genes altered in at least 5 patients
	tclin = Clin[rownames(Clin) %in% colnames(cna),]
	tclin = cbind(tclin, t(cna) )
	cdf = dan.df(rownames(cna),c( "Gene","Chisq_pval" ))
	for (g in rownames(cdf) ){
		tclin[(tclin[,g] %in% c( -1.5,-1,0,1 )) %in% c(T),g ] = 0
		tclin$egfr_class_consensus = factor(tclin$egfr_class_consensus,levels=ordered_classes)
		tabb = table(tclin[,g],tclin$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		cdf[g,"Gene"] = g
		cdf[g,"Chisq_pval"] = ch
	}
	cdf$Chisq_qval = p.adjust(cdf$Chisq_pval,method="BH")
	cdf = cdf[order(cdf$Chisq_qval),]
	save(cdf, file = paste0( OutDir,"Chen_coCNA.RData" )) # the arm with several zinc-finger protein is amplified in 4 patient with common EGFR variants

	################# SV
	###### Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_Genie_SampleLevel.RData" ))
	sv = dan.read(paste0(DataDir,"genie_14.0/data_sv.txt"))
	tclin = Clin2[rownames(Clin2) %in% sv$Sample_Id,]
	sv = sv[sv$Sample_Id %in% tclin$Sample,]
	sv = sv[nchar(sv$Site1_Hugo_Symbol)>0,]
	testable = names(dtable(sv$Site1_Hugo_Symbol))[dtable(sv$Site1_Hugo_Symbol)>=5]
	sv = sv[sv$Site1_Hugo_Symbol %in% testable,]
	cdf = dan.df(testable,c( "Gene","Chisq_pval" ))
	for (g in rownames(cdf) ){
		tclin$this_gene = "wt"
		ttmaf = sv[sv$Site1_Hugo_Symbol==g,]	
		tclin[unique(ttmaf$Sample_Id),"this_gene"] = "sv"
		tclin$this_gene = factor(tclin$this_gene,levels=c("wt","sv"))
		tclin$egfr_class_consensus = factor(tclin$egfr_class_consensus,levels=ordered_classes)
		tabb = table(tclin$this_gene,tclin$egfr_class_consensus)
		ch = chisq.test(tabb)$p.value
		cdf[g,"Gene"] = g
		cdf[g,"Chisq_pval"] = ch
	}
	cdf$Chisq_qval = p.adjust(cdf$Chisq_pval,method="BH")
	cdf = cdf[order(cdf$Chisq_qval),]
	save(cdf, file = paste0( OutDir,"Genie_coSV.RData" ))
}

TmbFga_associations = function( OutDir ){

	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]

	###### Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = rank(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsynonymous"])
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = (Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))/(max(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"])-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))
	for (sign in c( "TMB","TMB_nonsynonymous","FGA","TMB_nonsyn_rankNorm" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"Genie_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = (as.numeric(tc[,sign]))
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	Clin$Dataset = "Genie"
	Clinall = Clin[,c( "Dataset","Patient","egfr_class_consensus","TMB_nonsynonymous","TMB_nonsyn_rankNorm" )]

	###### Origimed
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = rank(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsynonymous"])
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = (Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))/(max(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"])-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))
	for (sign in c( "TMB_nonsynonymous","TMB_nonsyn_rankNorm" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"Origimed_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = (as.numeric(tc[,sign]))
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	Clin$Dataset = "Origimed"
	Clinall = rbind(Clinall,Clin[,c( "Dataset","Patient","egfr_class_consensus","TMB_nonsynonymous","TMB_nonsyn_rankNorm" )])

	###### Zhang
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Zhang.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = rank(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsynonymous"])
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = (Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))/(max(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"])-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))
	for (sign in c( "TMB_nonsynonymous","TMB_nonsyn_rankNorm" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"Zhang_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,sign])
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	Clin$Dataset = "Zhang"
	Clinall = rbind(Clinall,Clin[,c( "Dataset","Patient","egfr_class_consensus","TMB_nonsynonymous","TMB_nonsyn_rankNorm" )])

	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	load(paste0(CommonDataDir,'maf_LUAD_all.RData'))
	maf_ss = maf_ss[maf_ss$Variant_Classification %in% c("Frame_Shift_Del","Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation","Nonsense_Mutation","Splice_Site"),]
	maf_ss$Patient = substr(maf_ss$Tumor_Sample_Barcode,1,12)
	maf_ss = maf_ss[maf_ss$Patient %in% rownames(Clin),]
	Clin = Clin[unique(maf_ss$Patient),]
	for (rn in rownames(Clin)){
		Clin[rn,"TMB_nonsynonymous"] = sum(maf_ss$Patient==rn)/(50) # from Nature paper TCGA-LUAD methods
	}
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	seg = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/","TCGA-LUAD.cnv.tsv")) # from Xena
	seg$Patient = substr(seg$sample, 1,12)
	seg = seg[(seg$Patient %in% Clin$Patient),]
	# Some repeated patients. Keep the one with lower lexicographic ID
	temp_seg = seg[order(seg$sample),]
	temp_seg = seg[!duplicated(seg$Patient),]
	seg = seg[seg$sample %in% temp_seg$sample,]
	Clin$FGA = NA
	for (rn in rownames(Clin))
	{
	  this_seg = seg[seg$Patient==rn,]
	  total_len = sum(this_seg$End-this_seg$Start)
	  altered_len = sum(this_seg[abs(this_seg$value)>0.2,"End"]-this_seg[abs(this_seg$value)>0.2,"Start"])
	  Clin[rn,"FGA"] = altered_len/total_len
	}
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = rank(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsynonymous"])
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = (Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))/(max(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"])-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))
	for (sign in c( "TMB_nonsynonymous","FGA","TMB_nonsyn_rankNorm" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"TCGA_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,sign])
		if (sign=="TMB_nonsynonymous") { y = log10(y+1) }
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	Clin$Dataset = "TCGA"
	Clinall = rbind(Clinall,Clin[,c( "Dataset","Patient","egfr_class_consensus","TMB_nonsynonymous","TMB_nonsyn_rankNorm" )])

	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	tmb = dan.read( file = paste0( DataDir,"Chen_TMB_TableS3.txt" ) )
	rownames(tmb) = tmb$Patient.ID
	commonz = intersect(rownames(Clin),rownames(tmb))
	tmb$TMB = tmb$Total.Mutation.Burden.MB
	tmb$TMB_nonsynonymous = tmb$Total.Nonsynonymous.Mutation.Burden.MB
	Clin = cbind(Clin[commonz,],tmb[commonz,c( "TMB","TMB_nonsynonymous" )])
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	seg = dan.read(file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/","Chen2020/seg")) # from Xena
	Clin$FGA = NA
	for (rn in rownames(Clin))
	{
	  this_seg = seg[seg$ID==rn,]
  	  total_len = sum(this_seg$loc.end-this_seg$loc.start)
      altered_len = sum(this_seg[abs(this_seg$seg.mean)>0.2,"loc.end"]-this_seg[abs(this_seg$seg.mean)>0.2,"loc.start"])
	  Clin[rn,"FGA"] = altered_len/total_len
	}
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = rank(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsynonymous"])
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = (Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))/(max(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"])-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))
	for (sign in c( "TMB","TMB_nonsynonymous","FGA","TMB_nonsyn_rankNorm" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"Chen_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,sign])
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	Clin$Dataset = "ChenEAS"
	Clinall = rbind(Clinall,Clin[,c( "Dataset","Patient","egfr_class_consensus","TMB_nonsynonymous","TMB_nonsyn_rankNorm" )])

	###### TRACERx421
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TRACERx421_SampleLevel.RData" ))
	library(fst)
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	maf_ss = maf[maf$tumour_id %in% Clin2$Sample,]
	maf_ss$Sample = maf_ss$tumour_id
	maf_ss$Patient = maf_ss$patient_id
	mafe = maf_ss
	maf_ss = mafe[((mafe$exonic.func %in% c("frameshift substitution","frameshift insertion", "nonframeshift insertion", "nonframeshift substitution", "nonsynonymous","nonsynonymous SNV")) %in% c(T)),]
	tclin = Clin[Clin$Patient %in% maf_ss$Patient,]
	maf_ss = maf_ss[maf_ss$Patient %in% tclin$Patient,]
	maf_ss$mut = maf_ss$AAChange
	maf_ss[is.na(maf_ss$mut),"mut"] = maf_ss[is.na(maf_ss$mut),"NucleotideChange"]
	tmaf = maf_ss[!duplicated(paste0(maf_ss$Hugo_Symbol, maf_ss$Patient, maf_ss$mut)),]
	for (rn in rownames(Clin)){
		Clin[rn,"TMB_nonsynonymous"] = sum(tmaf$Patient==rn)/(50) # SureSelect V5. from https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=exomeProbesets 
	}
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = rank(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsynonymous"])
	Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"] = (Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))/(max(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"])-min(Clin[!is.na(Clin$TMB_nonsynonymous),"TMB_nonsyn_rankNorm"]))
	for (sign in c( "TMB_nonsynonymous","TMB_nonsyn_rankNorm" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"TRACERx_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = (as.numeric(tc[,sign]))
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	Clin$Dataset = "TRACERx421"
	Clinall = rbind(Clinall,Clin[,c( "Dataset","Patient","egfr_class_consensus","TMB_nonsynonymous","TMB_nonsyn_rankNorm" )])

	Clinall$TMB_nonsynonymous = as.numeric(Clinall$TMB_nonsynonymous)
	Clinall = Clinall[!is.na(Clinall$TMB_nonsynonymous),]
	ccsave = Clinall
	qup = quantile(Clinall$TMB_nonsynonymous,seq(0,1,0.01))['99%']
	qdown = quantile(Clinall$TMB_nonsynonymous,seq(0,1,0.01))['1%']
	Clinall = Clinall[(Clinall$TMB_nonsynonymous<qup) & (Clinall$TMB_nonsynonymous>qdown),]
	Clinall$egfr_class_consensus = factor(Clinall$egfr_class_consensus, levels = ordered_classes)
	Clinall$colorz = "steelblue4"
	Clinall$colorz[Clinall$egfr_class_consensus=="uncommon"] = "tomato3"
	Clinall$colorz[Clinall$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clinall$colorz[Clinall$egfr_class_consensus=="T790M"] = "darksalmon"
	Clinall$colorz[Clinall$egfr_class_consensus=="ex20ins"] = "orange2"
	Clinall$ranked_TMB_nonsynonymous = rank(Clinall$TMB_nonsynonymous)
	
	for (sign in c( "TMB_nonsynonymous","TMB_nonsyn_rankNorm" )){
		tc = Clinall[!is.na(Clinall[,sign]),]
		x = tc$egfr_class_consensus
		y = (as.numeric(tc[,sign]))
		fileName = paste0(OutDir,"AllDatasets_Datasets_vs_",sign,"_densities.pdf")
		dan.densityPlot( fileName, y, tc$Dataset, groupinglab = "Dataset", xlab = paste0(sign), ylab = "Density", show_medians = F, plotTitle = "", groupingColors = dan.colors(length(unique(tc$Dataset))), fileWidth = 4, fileHeight = 2.5 )
		fileName = paste0(OutDir,"AllDatasets_Datasets_vs_",sign,"_boxplots.pdf")
		pdf(fileName,3,2)
		ylimLeft = min(sapply( unique(tc$Dataset), function(lev) {quantile(y[tc$Dataset==lev],seq(0,1,0.25))["25%"]-(quantile(y[tc$Dataset==lev],seq(0,1,0.25))["75%"]-quantile(y[tc$Dataset==lev],seq(0,1,0.25))["25%"])*1.5} ))
		ylimRight = max(sapply( unique(tc$Dataset), function(lev) {quantile(y[tc$Dataset==lev],seq(0,1,0.25))["75%"]+(quantile(y[tc$Dataset==lev],seq(0,1,0.25))["75%"]-quantile(y[tc$Dataset==lev],seq(0,1,0.25))["25%"])*1.5} ))
		plot=dan.boxplots.multipages( tc$Dataset, y, fill = NULL, xlab = "Dataset", ylab = paste0(sign), filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = ylimLeft, ylimRight = ylimRight,comparisons = NULL, labelycoo = max(y)+0.01, xColors = dan.colors(length(unique(tc$Dataset))), labelJitteredPoints = NULL, includeJitters = F )
		print(plot)
		dev.off()
		fileName = paste0(OutDir,"AllDatasets_Datasets_vs_",sign,"_violins.pdf")
		pdf(fileName,3,2)
		plot=dan.violinplots.multipages( tc$Dataset, y, fill = NULL, xlab = "Dataset", ylab = paste0(sign), filllab = "", plotTitle = "", signifTest = NULL, ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y)+0.01, xColors = dan.colors(length(unique(tc$Dataset))), fillColors = "default", jitterColors = NULL, labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
		print(plot)
		dev.off()
		fileName = paste0(OutDir,"AllDatasets_EGFRclass_vs_",sign,"_boxplots.pdf")
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 1, fileWidth = 3, fileHeight = 2.5)
		fileName = paste0(OutDir,"AllDatasets_EGFRclass_vs_",sign,"_boxplots_byDataset.pdf")
		dan.boxplots( fileName, x, y, fill = tc$Dataset, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01,jitterDotSize = 1, fileWidth = 3, fileHeight = 2.5)
		fileName = paste0(OutDir,"AllDatasets_EGFRclass_vs_",sign,"_violins.pdf")
		pdf(fileName,3,2.6)
		plot=dan.violinplots.multipages( x, y, fill = NULL, xlab = "", ylab = paste0(sign), filllab = "", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, fillColors = "default", jitterColors = tc$colorz, labelJitteredPoints = NULL, includeJitters = F, hlines_coo = NULL, hlines_labels = NULL )
		print(plot)
		dev.off()
		fileName = paste0(OutDir,"AllDatasets_EGFRclass_vs_",sign,"_densities.pdf")
		plotTitle = paste0( "Kruskal-Wallis test, p-value = ", signif(kruskal.test(y~x)$p.value,2))
		dan.densityPlot( fileName, y, x, groupinglab = "", xlab = paste0(sign), ylab = "Density", show_medians = F, plotTitle = plotTitle, groupingColors = colorz_classes, fileWidth = 4, fileHeight = 2.5 )
	}
	save(Clinall,file=paste0( OutDir,"TMB_nonsyn_Clinall.RData" ))
}

DifferentialExpression_pooled = function( OutDir, comparisons_adj.P.Val_thresh ){
	###### TCGA
	RnaDataDir = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/"
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load( file = paste0(RnaDataDir, "TCGA_expression_gdc_htseq_counts/TCGA_LUAD_ge_rawcounts_GeneSymbols.RData") )
	rc = ge
	cn = as.character(colnames(rc))
	cn = cn[substr(cn,14,15) %in% c("01")] # only primary
	rc = rc[,cn]
	cn_short = substr(cn,1,12)
	colnames(rc) = cn_short
	commonz = intersect(rownames(Clin),colnames(rc))
	rc = rc[,commonz]
	Clin_expr = Clin[commonz,]
	samples = rownames(Clin_expr)
	rc = rc[,samples]
	Clin_expr$Dataset = "TCGA"
	clinAll = Clin_expr[,c( "Dataset","egfr_class_consensus" )]
	rcAll = rc
	###### Chen
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load( file = paste0("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/ge_Chen2020_expcounts_GeneSymbols.RData") )
	rc = ge
	commonz = intersect(rownames(Clin),colnames(rc))
	rc = rc[,commonz]
	Clin_expr = Clin[commonz,]
	samples = rownames(Clin_expr)
	rc = rc[,samples]
	Clin_expr$Dataset = "ChenEAS"
	clinAll = rbind(clinAll,Clin_expr[,c( "Dataset","egfr_class_consensus" )])
	rcAll = cbind(rcAll[intersect(rownames(rcAll),rownames(rc)),],rc[intersect(rownames(rcAll),rownames(rc)),])
	###### TRACERx421
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	mafr = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221123_TRACERx421_mutation_table_region.fst"))
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	maf = maf[(maf$Hugo_Symbol=="EGFR") & (maf$patient_id %in% rownames(Clin)),]
	mafr = mafr[mafr$mutation_id %in% maf$mutation_id,]
	region_level_clin = dan.df(0,c( "Patient","Region","egfr_class_consensus" ))
	# common
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="common", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="common",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	# uncommon
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="uncommon", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="uncommon",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	# compound - actually, all regions of the two compound cases have both compounded mutations (yeah!)
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="compound", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="compound",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	region_level_clin$Region = gsub("\\.","-",gsub(":","_",region_level_clin$Region))
	ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_counts_mat.fst"))
	rownames(ge) = ge$gene_id
	ge$gene_id = NULL
	region_level_clin = region_level_clin[region_level_clin$Region %in% colnames(ge),]
	region_level_clin = region_level_clin[!duplicated(region_level_clin$Region),]
	rc = ge[,region_level_clin$Region]
	region_level_clin$Dataset = "TRACERx421"
	rownames(region_level_clin) = region_level_clin$Region
	clinAll = rbind(clinAll,region_level_clin[,c( "Dataset","egfr_class_consensus" )])
	rcAll = cbind(rcAll[intersect(rownames(rcAll),rownames(rc)),],rc[intersect(rownames(rcAll),rownames(rc)),])

	y = DGEList(counts = rcAll, group = clinAll$egfr_class_consensus, remove.zeros = TRUE , genes = rownames(rcAll))

	isexpr = rep(FALSE,nrow(y))
	for (this_class in ordered_classes){
		c1 = rownames(clinAll[clinAll$egfr_class_consensus==this_class,])
		c1 = rowSums(cpm(y[,c1])>1) >= 3
		isexpr = isexpr | c1
	}
	y = y[isexpr,]
	y = calcNormFactors(y)
	Group = clinAll$egfr_class_consensus
	Dataset = clinAll$Dataset
	design = model.matrix(~0+Group+Dataset)
	colnames(design) = gsub("Group","",colnames(design) )
	mc = voom(y, design, plot=F)
	fit = lmFit(mc,design)
	contr.matrix_tcga = makeContrasts(CommonVsUncommon = common-uncommon, CommonVsCompound = common-compound, UncommonVsCompound = uncommon-compound, levels = colnames(design))
	fit_pair = contrasts.fit(fit, contrasts = contr.matrix_tcga[colnames(design),])
	fit_pair = eBayes(fit_pair)
	pp = as.vector(fit_pair$p.value)
	ppadjust = p.adjust(pp,method = 'BH')
	ppadj_mat = matrix(ppadjust, nrow = nrow(fit_pair$p.value),ncol = ncol(fit_pair$p.value))
	rownames(ppadj_mat) = rownames(fit_pair$p.value)
	colnames(ppadj_mat) = colnames(fit_pair$p.value)
	# sum(abs(dtTab)!=1*(ppadj_mat<0.1)) # This method gets adj p-values also correcting for across-genes hypotheses (i.e. *exactly* as decideTests with method = 'global')
	final_table = as.data.frame(ppadj_mat)
	colnames(final_table) = paste0("adjPval_",colnames(ppadj_mat))
	top = topTable(fit_pair, adjust="BH", num=Inf) # coef=2,
	dt1 = decideTests(fit_pair, method = "global", adjust.method="BH", p.value=comparisons_adj.P.Val_thresh)
	summary(dt1)
	dtTab = data.frame(dt1@.Data, stringsAsFactors=F)
	dtTab_signif = dtTab[rowSums(dtTab!=0)>0,]
	final = cbind(top,final_table[rownames(top),])
	fit_pair = data.frame(fit_pair$p.value)
	colnames(fit_pair) = paste0("nominalPval_",colnames(fit_pair))
	final = cbind(final,fit_pair[rownames(final),])
	if (length(colnames(dt1))!=0){
			for (cn in colnames(dt1)){
			final[rownames(final),paste0("signif_",cn)] = dtTab[rownames(final),cn]
		}
	}
	save(final, file = paste0(OutDir,"Limma_table_AllDatasets_AllGenes_AllVsAll_comparisonsSignLevel",comparisons_adj.P.Val_thresh,".RData"))
}

ProliferationScoring = function( OutDir ){

	library(singscore)
	# cp = read.table(file=paste0("processed/for_Amaia_7Dec2023/","LocardPaulet_proliferationSignature.txt"),header=T)
	cp = read.table(file=paste0(DataDir,"GOBP_CELL_CYCLE.v2024.1.Hs.grp"),header=T)
	colnames(cp) = "gene"
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))

	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/TCGA_expression_gdc_TPM/TCGA_LUAD_ge_TPM_GeneSymbols.RData")
	# load( file = paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Data/TCGA_expression_gdc_htseq_counts/TCGA_LUAD_ge_rawcounts_GeneSymbols.RData") )
	Clin = Clin[Clin$egfr_class_consensus %in% c( "common","uncommon" ),]
	cn = as.character(colnames(ge))
	cn = cn[substr(cn,14,15) %in% c("01")] # only primary
	ge = ge[,cn]
	cn_short = substr(cn,1,12)
	colnames(ge) = cn_short
	commonz = intersect(rownames(Clin),colnames(ge))
	ge = ge[,commonz]
	logtpm = log2(ge+1)	
	Clin = Clin[commonz,]
	rankData = rankGenes(logtpm)
	scoredf = simpleScore(rankData, upSet = cp$gene)
	Clin[rownames(scoredf),"LP_proliferation"] = scoredf$TotalScore
	Clin$Dataset = "TCGA"
	Clin$Sample = Clin$Patient
	fileName = paste0(OutDir,"Proliferation_acrossClass_TCGA.pdf")
	x = factor(Clin$egfr_class_consensus,levels=c( "common","uncommon" ))
	y = Clin$LP_proliferation
	xcolorz = colorz_classes[ordered_classes %in% levels(x)]
	jittercolorz = dan.expand_colors(as.character(x),levels(x),xcolorz)
	dan.boxplots( fileName, x, y, xlab = "", ylab = "GO:BP Cell Cycle score", filllab = "EGFR class", plotTitle = "", signifTest = "wilcox.test", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y), xColors = xcolorz, fillColors = "default", jitterColors = jittercolorz, labelJitteredPoints = NULL, jitterDotSize = 2, fileWidth = 3, fileHeight = 2.5, hlines_coo = NULL, hlines_labels = NULL )
	clinall = Clin[,c( "Dataset","Sample","egfr_class_consensus","LP_proliferation" )]

	library(fst)
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	mafr = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221123_TRACERx421_mutation_table_region.fst"))
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	maf = maf[(maf$Hugo_Symbol=="EGFR") & (maf$patient_id %in% rownames(Clin)),]
	mafr = mafr[mafr$mutation_id %in% maf$mutation_id,]
	region_level_clin = dan.df(0,c( "Patient","Region","egfr_class_consensus" ))
	# common
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="common", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="common",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	# uncommon
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="uncommon", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="uncommon",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	# compound - actually, all regions of the two compound cases have both compounded mutations (yeah!)
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="compound", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="compound",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	region_level_clin$Region = gsub("\\.","-",gsub(":","_",region_level_clin$Region))
	# ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_counts_mat.fst"))
	ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"))
	rownames(ge) = ge$gene_id
	ge$gene_id = NULL
	region_level_clin = region_level_clin[region_level_clin$Region %in% colnames(ge),]
	region_level_clin = region_level_clin[!duplicated(region_level_clin$Region),]
	rc = ge[,region_level_clin$Region]
	region_level_clin$Dataset = "TRACERx421"
	rownames(region_level_clin) = region_level_clin$Region
	logtpm = log2(ge[,region_level_clin$Region]+1)
	rankData = rankGenes(logtpm)
	scoredf = simpleScore(rankData, upSet = cp$gene)
	region_level_clin[rownames(scoredf),"LP_proliferation"] = scoredf$TotalScore
	region_level_clin$Sample = region_level_clin$Region
	region_level_clin = region_level_clin[region_level_clin$egfr_class_consensus %in% c( "common","uncommon" ),]
	fileName = paste0(OutDir,"Proliferation_acrossClass_TRACERx.pdf")
	x = factor(region_level_clin$egfr_class_consensus,levels=c( "common","uncommon" ))
	y = region_level_clin$LP_proliferation
	xcolorz = colorz_classes[ordered_classes %in% levels(x)]
	jittercolorz = dan.expand_colors(as.character(x),levels(x),xcolorz)
	dan.boxplots( fileName, x, y, xlab = "", ylab = "GO:BP Cell Cycle score", filllab = "EGFR class", plotTitle = "", signifTest = "wilcox.test", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y), xColors = xcolorz, fillColors = "default", jitterColors = jittercolorz, labelJitteredPoints = NULL, jitterDotSize = 2, fileWidth = 3, fileHeight = 2.5, hlines_coo = NULL, hlines_labels = NULL )


	# ###### Chen
	# load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	# colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	# ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	# load("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/ge_Chen2020_expcounts_GeneSymbols.RData") # ge_Chen2020_expcounts_GeneSymbols.RData ge_Chen2020_normalized_GeneSymbols.RData
	# # load("/mnt/ndata/daniele/lung_multiregion/Data/Chen2020/ge_Chen2020_normalized_GeneSymbols.RData") # ge_Chen2020_expcounts_GeneSymbols.RData ge_Chen2020_normalized_GeneSymbols.RData
	# rc = log2(ge+1)
	# commonz = intersect(rownames(Clin),colnames(rc))
	# logtpm = rc[,commonz]
	# Clin = Clin[commonz,]
	# rankData = rankGenes(logtpm)
	# scoredf = simpleScore(rankData, upSet = cp$gene)
	# Clin[rownames(scoredf),"LP_proliferation"] = scoredf$TotalScore
	# Clin$Sample = Clin$Patient
	# Clin$Dataset = "ChenEAS"
	# clinall = rbind(clinall,Clin[,c( "Dataset","Sample","egfr_class_consensus","LP_proliferation" )])

	fileName = paste0(OutDir,"Proliferation_acrossClass_acrossDatasets.pdf")
	x = clinall$Dataset
	y = clinall$LP_proliferation
	dan.boxplots( fileName, x, y, fill = clinall$egfr_class_consensus, xlab = "", ylab = "GO:BP Cell Cycle score", filllab = "EGFR class", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = max(y), xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, jitterDotSize = 1, fileWidth = 4, fileHeight = 2.5, hlines_coo = NULL, hlines_labels = NULL )
}

ComparisonDegs = function( OutDir, comparisons_adj.P.Val_thresh ){
	load(paste0( OutDir,"Limma_table_AllDatasets_AllGenes_AllVsAll_comparisonsSignLevel",comparisons_adj.P.Val_thresh,".RData" ))
	top = final
	dtable(top$signif_CommonVsUncommon)
	dtable(top$signif_CommonVsCompound)
	dtable(top$signif_UncommonVsCompound)

	markers_common1 = top[(top$signif_CommonVsUncommon==1) & (top$CommonVsUncommon>1),"genes"]
	markers_common2 = top[(top$signif_CommonVsCompound==1) & (top$CommonVsCompound>1),"genes"]

	markers_uncommon1 = top[(top$signif_CommonVsUncommon==(-1) ) & (top$CommonVsUncommon<(-1) ),"genes"]
	markers_uncommon2 = top[(top$signif_UncommonVsCompound==1) & (top$UncommonVsCompound>1),"genes"]

	markers_compound1 = top[(top$signif_CommonVsCompound==(-1) ) & (top$CommonVsCompound<(-1) ),"genes"]
	markers_compound2 = top[(top$signif_UncommonVsCompound==(-1) ) & (top$UncommonVsCompound<(-1) ),"genes"]

	gene_universe = rownames(top)
	msig_df_H = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
	msig_df_GO_BP = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
	msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP ))
	msig_list = split(x=msig_df$gene_symbol, f=msig_df$gs_name)
	fgRes = fgsea::fora(pathways = msig_list,
	                       genes = markers_common1,
	                       universe = gene_universe)
			fgRes = fgRes[fgRes$padj <= 0.05,]
			a = as.data.frame(fgRes)
			a = apply(a,2,as.character)
			a = as.data.frame(fgRes)
			a$overlapGenes = NULL
			a$pval = NULL
			# a$padj = signif(a$padj,2)
	# nothing
	# common_gsea = a[c(1,5,6,9,14,15,23),]
	# common_gsea$log_qval = -log10(common_gsea$padj)
	# common_gsea = common_gsea[order(common_gsea$log_qval),]
	# common_gsea$pathway = factor(common_gsea$pathway,levels=as.character(common_gsea$pathway))
	# pdf(paste0(OutDir,"shared_common_barplot.pdf"),7,3)
	# ggplot(common_gsea, aes(x=pathway, y=log_qval)) + geom_bar(stat='identity',fill="steelblue4") + xlab("" )+ylab( "-log10(adjusted p-value)" ) +scale_y_reverse()+ coord_flip() + theme_classic()
	# dev.off()

	gene_universe = rownames(top)
	msig_df_H = msigdbr::msigdbr(species = "Homo sapiens", category = "H")
	msig_df_GO_BP = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
	msig_df = dplyr::bind_rows(list( msig_df_H,msig_df_GO_BP ))
	msig_list = split(x=msig_df$gene_symbol, f=msig_df$gs_name)
	fgRes = fgsea::fora(pathways = msig_list,
	                       genes = markers_uncommon1,
	                       universe = gene_universe)
			fgRes = fgRes[fgRes$padj <= 0.05,]
			a = as.data.frame(fgRes)
			a = apply(a,2,as.character)
			a = as.data.frame(fgRes)
	a = a[order(a$padj),]
	a$overlapGenes = NULL
	dan.write(a,paste0(OutDir,"../Paper_SupplTables/TableS4.txt" ))			
			a$pval = NULL
			# a$padj = signif(a$padj,2)

	uncommon_gsea = a[1:20,]
	uncommon_gsea$log_qval = -log10(uncommon_gsea$padj)
	uncommon_gsea = uncommon_gsea[order(uncommon_gsea$log_qval),]
	uncommon_gsea$pathway = factor(uncommon_gsea$pathway,levels=as.character(uncommon_gsea$pathway))

	pdf(paste0(OutDir,"AllDatasets_uncommon_barplot.pdf"),3.8,3)
	plot=(ggplot(uncommon_gsea, aes(x=pathway, y=log_qval)) + geom_bar(stat='identity',fill="tomato3") + xlab("" )+ylab( "-log10(adjusted p-value)" ) + coord_flip() + theme_classic())
	plot=plot+theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6))
	print(plot)
	dev.off()

	ymax = max(c(-log10(top$nominalPval_CommonVsUncommon),-log10(top$nominalPval_CommonVsCompound),-log10(top$nominalPval_UncommonVsCompound)))
	# three volcanos
	top$Marker = "none"
	top[top$signif_CommonVsUncommon==(1),"Marker"] = "common"
	top[top$signif_CommonVsUncommon==(-1),"Marker"] = "uncommon"
	selected_genes1 = top[top$signif_CommonVsUncommon==1,]
	selected_genes1 = selected_genes1[order(abs(selected_genes1$nominalPval_CommonVsUncommon),decreasing=F),"genes"][1:10]
	selected_genes2 = top[top$signif_CommonVsUncommon==(-1),]
	selected_genes2 = selected_genes2[order(abs(selected_genes2$nominalPval_CommonVsUncommon),decreasing=F),"genes"][1:10]
	top$repel_labelz = top$genes
	top[!(top$genes %in% c( selected_genes1,selected_genes2 )),"repel_labelz"] = ""
	dan.scatterplot(paste0(OutDir,"volcano_common_vs_uncommon.pdf"),x=top$CommonVsUncommon, y=-log10(top$nominalPval_CommonVsUncommon), fill=factor(top$Marker,levels=c( "common","none","uncommon" )),fillColors=c( "steelblue4","gray","tomato3" ),repel_labels=top$repel_labelz,
		xlab="log2(fold change)",ylab="-log10(p-value)",filllab="Marker of",dotSize=1,fileWidth=4.5,fileHeight=3,ylimLeft=0,ylimRight=ymax) # in cell cycle: PSRC1, SASS6, SKP2, TEX15

	top$Marker = "none"
	top[top$signif_CommonVsCompound==(1),"Marker"] = "common"
	top[top$signif_CommonVsCompound==(-1),"Marker"] = "compound"
	selected_genes1 = top[top$signif_CommonVsCompound==1,]
	selected_genes1 = selected_genes1[order(abs(selected_genes1$nominalPval_CommonVsCompound),decreasing=F),"genes"][1:10]
	selected_genes2 = top[top$signif_CommonVsCompound==(-1),]
	selected_genes2 = selected_genes2[order(abs(selected_genes2$nominalPval_CommonVsCompound),decreasing=F),"genes"][1:10]
	top$repel_labelz = top$genes
	top[!(top$genes %in% c( selected_genes1,selected_genes2 )),"repel_labelz"] = ""
	dan.scatterplot(paste0(OutDir,"volcano_common_vs_compound.pdf"),x=top$CommonVsCompound, y=-log10(top$nominalPval_CommonVsCompound), fill=factor(top$Marker,levels=c( "common","none","compound" )),fillColors=c( "steelblue4","gray","mediumpurple4" ),repel_labels=top$repel_labelz,
		xlab="log2(fold change)",ylab="-log10(p-value)",filllab="Marker of",dotSize=1,fileWidth=4.2,fileHeight=3,ylimLeft=0,ylimRight=ymax) # in cell cycle: PSRC1, SASS6, SKP2, TEX15

	top$Marker = "none"
	top[top$signif_UncommonVsCompound==(1),"Marker"] = "uncommon"
	top[top$signif_UncommonVsCompound==(-1),"Marker"] = "compound"
	selected_genes1 = top[top$signif_UncommonVsCompound==1,]
	selected_genes1 = selected_genes1[order(abs(selected_genes1$nominalPval_UncommonVsCompound),decreasing=F),"genes"][1:10]
	selected_genes2 = top[top$signif_UncommonVsCompound==(-1),]
	selected_genes2 = selected_genes2[order(abs(selected_genes2$nominalPval_UncommonVsCompound),decreasing=F),"genes"][1:10]
	top$repel_labelz = top$genes
	top[!(top$genes %in% c( selected_genes1,selected_genes2 )),"repel_labelz"] = ""
	dan.scatterplot(paste0(OutDir,"volcano_uncommon_vs_compound.pdf"),x=top$UncommonVsCompound, y=-log10(top$nominalPval_UncommonVsCompound), fill=factor(top$Marker,levels=c( "uncommon","none","compound" )),fillColors=c( "tomato3","gray","mediumpurple4" ),repel_labels=top$repel_labelz,
		xlab="log2(fold change)",ylab="-log10(p-value)",filllab="Marker of",dotSize=1,fileWidth=4.2,fileHeight=3,ylimLeft=0,ylimRight=ymax) # in cell cycle: PSRC1, SASS6, SKP2, TEX15
}

DifferentialProtein = function( OutDir ){
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	load("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/TCGA_RPPA/rppa_tcga_luad.RData")
	commonz = intersect(rownames(Clin),rownames(ge))
	Clin = Clin[commonz,]
	ge = ge[commonz,]
	ge$Pattern = NULL
	vdf = dan.df( colnames(ge),c( "Protein","pval","qval" ) )
	for (v in rownames(vdf)){
		vdf[v,"Protein"] = v
		Clin$this_tr = as.numeric(ge[rownames(Clin),v])
		vdf[v,"pval"] = kruskal.test(this_tr~egfr_class_consensus,data=Clin)$p.value
	}
	vdf$qval = p.adjust(vdf$pval,method="BH")
	vdf = vdf[order(vdf$pval),]
	save(vdf, file = paste0(OutDir,"TCGA_RPPA.RData"))
}

ViperAnalysis = function( OutDir ){
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	viper = dan.read( file = paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoSsignature_VIPER/ClassicalSignature_29Oct_3PatientsSplitted_qval0.1_logfc1/tcga_viper_SingleSample_GeneSymbols.txt"), row.names = 1)
	colnames(viper) = gsub("\\.","-",colnames(viper))
	commonz = intersect(colnames(viper),rownames(Clin))
	Clin = Clin[commonz,]
	viper = viper[,commonz]

	vdf = dan.df( rownames(viper),c( "TR","pval","qval" ) )
	viper = t(viper)
	for (v in rownames(vdf)){
		vdf[v,"TR"] = v
		Clin$this_tr = viper[rownames(Clin),v]
		vdf[v,"pval"] = kruskal.test(this_tr~egfr_class_consensus,data=Clin)$p.value
	}

	vdf$qval = p.adjust(vdf$pval,method="BH")
	vdf = vdf[order(vdf$pval),]
	save(vdf, file = paste0(OutDir,"TCGA_viper.RData"))
}

purity_associations = function( OutDir ){
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	###### Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in c( "Purity" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"Genie_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,sign])
		if (sign=="TMB_nonsynonymous") { y = log10(y+1) }
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### Origimed
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Origimed.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	for (sign in c( "TUMOR_PURTITY" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"Origimed_EGFRclass_vs_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,sign])
		if (sign=="TMB_nonsynonymous") { y = log10(y+1) }
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	pur = dan.read( paste0(DataDir,"aran_purities.txt") )
	load(file = paste0( OutDir,"../Preprocessing/","Clin2_TCGA_SampleLevel.RData" ))
	pur = pur[pur$Sample.ID %in% rownames(Clin2),]
	rownames(pur) = substr(pur$Sample.ID,1,12)
	commonz = intersect(rownames(pur),rownames(Clin))
	these_vars = c("ESTIMATE","ABSOLUTE","LUMP","IHC","CPE")
	Clin = cbind(Clin[commonz,],pur[commonz,these_vars])
	for (sign in c( these_vars )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"TCGA_EGFRclass_vs_Purity_",sign,".pdf")
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,sign])
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y,na.rm=T)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	pur = dan.read( paste0(DataDir,"Chen_GIS031_clinical_data.tsv") )
	rownames(pur) = pur$Patient.ID
	commonz = intersect(rownames(pur),rownames(Clin))
	these_vars = c("Purity","Stage")
	Clin = cbind(Clin[commonz,],pur[commonz,these_vars])
	for (sign in c( "Purity" )){
		tc = Clin[!is.na(Clin[,sign]),]
		fileName = paste0(OutDir,"Chen_EGFRclass_vs_Purity.pdf")
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,sign])
		dan.boxplots( fileName, x, y, fill = NULL, xlab = "", ylab = paste0(sign), plotTitle = "", signifTest = "kruskal", comparisons = NULL, labelycoo = max(y,na.rm=T)+0.01, xColors = colorz_classes, jitterColors = tc$colorz, jitterDotSize = 2.5, fileWidth = 6, fileHeight = 5)
	}
}

tme_associations = function( OutDir ){
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	dataset="TCGA"
	load(paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternMicroEnv/ConsensusTME/All_datasets/ctme_",dataset,".RData"))
	ctme = t(ctme)
	commonz = intersect(rownames(ctme),rownames(Clin))
	colnames(ctme) = paste0("ctme_",colnames(ctme) )
	tc = cbind(Clin[commonz,],ctme[commonz,])
	pdf(paste0(OutDir,"TCGA_EGFRclass_vs_ctme.pdf"),6,5)
	for (ct in colnames(ctme)){
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,ct])
	  	plot=dan.boxplots.multipages(x=x,y=y,signifTest="kruskal",xlab="",labelycoo = max(y,na.rm=T)+0.01,ylab=ct,xColors = colorz_classes, jitterColors = tc$colorz,includeJitters = T)
	  	print(plot)
	}
	dev.off()
	tc_tcga = tc

	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	tme = dan.read(file=paste0(DataDir,"TCGA_HE_TME_Szczurek2022.txt"))
	cts = c( "m_STROMA","m_MIXED","m_IMMUNE","m_VESSEL","m_BRONCHI","m_NECROSIS","m_LUNG","t_TUMOR","t_STROMA","t_MIXED","t_IMMUNE","t_VESSEL","t_BRONCHI","t_NECROSIS","t_LUNG","ITLR","Shannon","Simpson" )
	tc = Clin[Clin$Patient %in% tme$patient_id, ]
	for (rn in rownames(tc)){
		for (ct in cts){
			this_tme = tme[tme$patient_id==rn,]
			tc[rn,ct] = mean(this_tme[,ct],na.rm=T)
		}
	}
	pdf(paste0(OutDir,"TCGA_EGFRclass_vs_HEtme.pdf"),6,5)
	for (ct in cts){
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,ct])
	  	plot=dan.boxplots.multipages(x=x,y=y,signifTest="kruskal",xlab="",labelycoo = max(y,na.rm=T)+0.01,ylab=ct,xColors = colorz_classes, jitterColors = tc$colorz,includeJitters = T)
	  	print(plot)
	}
	dev.off()

	###### ChenEAS
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	Clin$egfr_class_consensus = factor(Clin$egfr_class_consensus, levels = ordered_classes)
	Clin$colorz = "steelblue4"
	Clin$colorz[Clin$egfr_class_consensus=="uncommon"] = "tomato3"
	Clin$colorz[Clin$egfr_class_consensus=="compound"] = "mediumpurple4"
	Clin$colorz[Clin$egfr_class_consensus=="T790M"] = "darksalmon"
	Clin$colorz[Clin$egfr_class_consensus=="ex20ins"] = "orange2"
	dataset="Chen"
	load(paste0("/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternMicroEnv/ConsensusTME/All_datasets/ctme_",dataset,".RData"))
	ctme = t(ctme)
	commonz = intersect(rownames(ctme),rownames(Clin))
	colnames(ctme) = paste0("ctme_",colnames(ctme) )
	tc = cbind(Clin[commonz,],ctme[commonz,])
	colorz_classes = clco[rownames(clco) %in% (tc$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (tc$egfr_class_consensus),"classes"]
	pdf(paste0(OutDir,"Chen_EGFRclass_vs_ctme.pdf"),6,5)
	for (ct in colnames(ctme)){
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,ct])
	  	plot=dan.boxplots.multipages(x=x,y=y,signifTest="kruskal",xlab="",labelycoo = max(y,na.rm=T)+0.01,ylab=ct,xColors = colorz_classes, jitterColors = tc$colorz,includeJitters = T)
	  	print(plot)
	}
	dev.off()
	tc_chen = tc

	library(ConsensusTME)
	###### TRACERx421
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TRACERx421.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	mafr = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221123_TRACERx421_mutation_table_region.fst"))
	maf = read_fst(paste0(DataDir,"TRACERx_421/figurecode/data/20221109_TRACERx421_mutation_table.fst"))
	maf = maf[(maf$Hugo_Symbol=="EGFR") & (maf$patient_id %in% rownames(Clin)),]
	mafr = mafr[mafr$mutation_id %in% maf$mutation_id,]
	region_level_clin = dan.df(0,c( "Patient","Region","egfr_class_consensus" ))
	# common
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="common", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="common",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	# uncommon
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="uncommon", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="uncommon",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	# compound - actually, all regions of the two compound cases have both compounded mutations (yeah!)
	for (rn in rownames(Clin[Clin$egfr_class_consensus=="compound", ]) ){
		mm = maf[maf$patient_id==rn,"mutation_id"]
		print(length(mm))
		mmr = mafr[mafr$mutation_id %in% mm,]
		thiz = data.frame(Patient=rn,Region=mmr$RegionID,egfr_class_consensus="compound",stringsAsFactors=F)
		region_level_clin = rbind(region_level_clin,thiz)
	}
	region_level_clin$Region = gsub("\\.","-",gsub(":","_",region_level_clin$Region))
	ge = read_fst(paste0("/mnt/ndata/daniele/alfredo_egfr/Data/TRACERx_421/transcriptomics_scripts_data_updated/20221014_transcriptomic_DATA/2022-10-17_rsem_tpm_mat.fst"))
	rownames(ge) = ge$gene_id
	ge$gene_id = NULL
	region_level_clin = region_level_clin[region_level_clin$Region %in% colnames(ge),]
	region_level_clin = region_level_clin[!duplicated(region_level_clin$Region),]
	rownames(region_level_clin) = region_level_clin$Region
	logtpm = log2(ge+0.0001)
	logtpm = logtpm[,substr(colnames(logtpm),1,8) %in% rownames(Clin)]
	logtpm = logtpm[,!grepl("_N",colnames(logtpm))]
	ctme = consensusTMEAnalysis(as.matrix(logtpm), cancer = "LUAD", statMethod = "singScore")
	ctme = ctme$Scores
	rownames(ctme) = paste0("ctme_",rownames(ctme))
	ctme = data.frame(t(ctme))
	# ctme$Patient = substr(rownames(ctme),1,8)
	# ctme = aggregate(.~Patient,data=ctme,FUN='mean')
	# rownames(ctme) = ctme$Patient	
	# ctme$Patient = NULL
	commonz = intersect(rownames(ctme),rownames(region_level_clin))
	tc = cbind(region_level_clin[commonz,],ctme[commonz,])
	colorz_classes = clco[rownames(clco) %in% (tc$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (tc$egfr_class_consensus),"classes"]

	pdf(paste0(OutDir,"TRACERx_EGFRclass_vs_ctme.pdf"),6,5)
	for (ct in colnames(ctme)){
		x = tc$egfr_class_consensus
		y = as.numeric(tc[,ct])
	  	plot=dan.boxplots.multipages(x=x,y=y,signifTest="kruskal",xlab="",labelycoo = max(y,na.rm=T)+0.01,ylab=ct,xColors = colorz_classes, jitterColors = "black",includeJitters = T)
	  	print(plot)
	}
	dev.off()
	tc_tracerx = tc
	tcAll = rbind(tc_tcga[,c( "Patient","egfr_class_consensus", colnames(ctme))],rbind(tc_chen[,c( "Patient","egfr_class_consensus", colnames(ctme))],tc_tracerx[,c( "Patient","egfr_class_consensus", colnames(ctme))]))
	library(ComplexHeatmap)
	pdf(paste0(OutDir, "tcAll.pdf"),8,3)
	print(ComplexHeatmap::Heatmap(t(tcAll[,colnames(ctme)]),cluster_rows=TRUE, cluster_columns=FALSE, column_names_rot = 45,use_raster = FALSE))
	dev.off()

	pdf(paste0(OutDir,"All_EGFRclass_vs_ctme.pdf"),6,5)
	for (ct in colnames(ctme)){
		x = tcAll$egfr_class_consensus
		y = as.numeric(tcAll[,ct])
	  	plot=dan.boxplots.multipages(x=x,y=y,signifTest="kruskal",xlab="",labelycoo = max(y,na.rm=T)+0.01,ylab=ct,xColors = colorz_classes, jitterColors = "black",includeJitters = T)
	  	print(plot)
	}
	dev.off()
	dan.write(tcAll,file=paste0(OutDir,"All_EGFRclass_vs_ctme.txt"))
}

methylation_associations = function( OutDir, contr.matrix ){
	###### TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	colorz_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Clin$egfr_class_consensus),"classes"]
	load("/mnt/ndata/daniele/lung_multiregion/methylation/Data/TCGA_LUAD_all/tcga.luad.all.RData") 
	# Only primary
	meth = methylome.tumor[,substr(colnames(methylome.tumor),14,15)=="01"]
	cn = colnames(meth)
	cn_format = gsub("\\.","-",substr(cn,1,12))
	colnames(meth) = cn_format
	meth = meth[rowSums(is.na(meth))==0,]
	pseudo = 0.0000001
	meth[meth==1] = meth[meth==1]-pseudo
	meth[meth==0] = meth[meth==0]+pseudo
	mc = log2(meth/(1-meth))
	commonz = intersect(rownames(Clin),colnames(mc))
	mc = mc[,commonz]
	Clin_meth = Clin[commonz,]
	rownames(Clin_meth) = Clin_meth$Patient
	mc = mc[,rownames(Clin_meth)]
	Group = Clin_meth$egfr_class_consensus
	design = model.matrix(~0+Group)
	colnames(design) = gsub("Group","",colnames(design) )
	fit = lmFit(mc,design)
	fit_pair = contrasts.fit(fit, contrasts = contr.matrix[colnames(design),])
	fit_pair = eBayes(fit_pair)
	top = topTable(fit_pair, adjust="BH", num=Inf) # coef=2,
	comparisons_adj.P.Val_thresh = 0.25
	dt1 = decideTests(fit_pair, method = "global", adjust.method="BH", p.value=comparisons_adj.P.Val_thresh)
	summary(dt1)
	dtTab = data.frame(dt1@.Data, stringsAsFactors=F)
	dtTab_signif = dtTab[rowSums(dtTab!=0)>0,]
	save(top, file = paste0(OutDir,"TopTable_DiffMethacrossEgfrClasses_TCGA.RData"))
	final = top
	if (length(colnames(dt1))!=0){
			for (cn in colnames(dt1)){
			final[rownames(final),paste0("signif_",cn)] = dtTab[rownames(final),cn]
		}
	}
	save(final, file = paste0(OutDir,"Limma_table_TCGA_AllMethProbes_AllVsAll_comparisonsSignLevel",comparisons_adj.P.Val_thresh,".RData"))
}

therapy_associations = function( OutDir ){
	library(survival)
	### Filling in Genie
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Genie.RData" ))
	######### MSK organotropism (also..)
	clin = dan.read( paste0(DataDir,"luad_mskcc_2023_met_organotropism/data_clinical_sample.txt"))
	clin$Patient = paste0("GENIE-MSK-",clin$PATIENT_ID)
	commonz = intersect(rownames(Clin),clin$Patient)
	clin = clin[clin$Patient %in% commonz,]
	clin = clin[!duplicated(clin$Patient),]
	rownames(clin) = clin$Patient
	varz = c( "treatment_chemo","treatment_immuno","treatment_targeted" )
	for (v in varz){ Clin[,v] = NA }
	Clin[rownames(clin),"treatment_chemo"] = tolower(clin$PRE_SAMPLE_CHEMOTHERAPY)
	Clin[clin[clin$PRE_SAMPLE_CHEMOTHERAPY!="Yes","Patient"],"treatment_chemo"] = tolower(clin[clin[clin$PRE_SAMPLE_CHEMOTHERAPY!="Yes","Patient"],"POST_SAMPLE_CHEMOTHERAPY"])
	Clin[rownames(clin),"treatment_immuno"] = tolower(clin$PRE_SAMPLE_IMMUNOTHERAPY)
	Clin[clin[clin$PRE_SAMPLE_IMMUNOTHERAPY!="Yes","Patient"],"treatment_immuno"] = tolower(clin[clin[clin$PRE_SAMPLE_IMMUNOTHERAPY!="Yes","Patient"],"POST_SAMPLE_IMMUNOTHERAPY"])
	Clin[rownames(clin),"treatment_targeted"] = tolower(clin$PRE_SAMPLE_TARGETED)
	Clin[clin[clin$PRE_SAMPLE_TARGETED!="Yes","Patient"],"treatment_targeted"] = tolower(clin[clin[clin$PRE_SAMPLE_TARGETED!="Yes","Patient"],"POST_SAMPLE_TARGETED"])
	for (v in varz){ Clin[(nchar(Clin[,v])<2) %in% c(T),v] = NA }
	for (v in varz){ print(dtable(Clin[,v])) }
	################### MSK 2017 ###################
	msk = dan.read( paste0(DataDir,"lung_msk_2017_clinical_data.tsv") )
	msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
	commonz = intersect(msk$Patient,Clin$Patient)
	msk = msk[msk$Patient %in% commonz,]
	msk = msk[!duplicated(msk$Patient),]
	rownames(msk) = msk$Patient
	Clin[msk[(msk$Immunotherapy=="YES"),"Patient"],"treatment_immuno" ] = "yes"
	Clin[ intersect( Clin[!((Clin$treatment_immuno=="yes") %in% c(T)),"Patient" ],msk[(msk$Immunotherapy=="NO"),"Patient"] ),"treatment_immuno" ] = "no"
	Clin[msk[(msk$Target.Therapy=="YES"),"Patient"],"treatment_targeted" ] = "yes"
	Clin[ intersect( Clin[!((Clin$treatment_targeted=="yes") %in% c(T)),"Patient" ],msk[(msk$Target.Therapy=="NO"),"Patient"] ),"treatment_targeted" ] = "no"
	Clin[msk[(msk$Chemotherapy=="YES"),"Patient"],"treatment_chemo" ] = "yes"
	Clin[ intersect( Clin[!((Clin$treatment_chemo=="yes") %in% c(T)),"Patient" ],msk[(msk$Chemotherapy=="NO"),"Patient"] ),"treatment_chemo" ] = "no"
	for (v in varz){ print(dtable(Clin[,v])) }
	################### MSK Mind Oct 2022 ###################
	msk = dan.read( paste0(DataDir,"lung_msk_mind_2020_clinical_data.tsv") )
	msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
	commonz = intersect(msk$Patient,Clin$Patient)
	msk = msk[msk$Patient %in% commonz,]
	msk = msk[!duplicated(msk$Patient),]
	rownames(msk) = msk$Patient
	Clin[rownames(msk),"treatment_immuno"] = "yes"
	for (v in varz){ print(dtable(Clin[,v])) }
	################### MSK NSCLC BM Aug 2023 ###################
	msk = dan.read( paste0(DataDir,"bm_nsclc_mskcc_2023_clinical_data.tsv") )
	msk$Patient = paste0("GENIE-MSK-",msk$Patient.ID)
	commonz = intersect(msk$Patient,Clin$Patient)
	msk = msk[msk$Patient %in% commonz,]
	msk = msk[!duplicated(msk$Patient),]
	rownames(msk) = msk$Patient
	Clin[ intersect( Clin[!((Clin$treatment_chemo=="yes") %in% c(T)),"Patient" ],msk[(msk$Treatment.prior.to.BM.resection=="No treatment"),"Patient"] ),"treatment_chemo" ] = "no"
	Clin[ intersect( Clin[!((Clin$treatment_targeted=="yes") %in% c(T)),"Patient" ],msk[(msk$Treatment.prior.to.BM.resection=="No treatment"),"Patient"] ),"treatment_targeted" ] = "no"
	Clin[ intersect( Clin[!((Clin$treatment_immuno=="yes") %in% c(T)),"Patient" ],msk[(msk$Treatment.prior.to.BM.resection=="No treatment"),"Patient"] ),"treatment_immuno" ] = "no"
	Clin[msk[(msk$Prior.TKI.at.any.time=="Yes") %in% c(T),"Patient"],"treatment_targeted" ] = "yes"
	for (v in varz){ print(dtable(Clin[,v])) }
	Cling = Clin

	### Filling in TCGA
	load(file = paste0( OutDir,"../Preprocessing/","Clin_TCGA.RData" ))
	Clinical_extended = as.data.frame(t(read.table(paste0(CommonDataDir,"Clinical_extended/LUAD.Clinical.txt"), quote = '', sep = "\t", header = FALSE, row.names = 1)), stringsAsFactors = FALSE)
	Clinical_extended[,"Patient"] <- toupper(as.character(Clinical_extended[,"bcr_patient_barcode"]))
	rownames(Clinical_extended) = Clinical_extended[,"Patient"]
	commonz = intersect(Clinical_extended$Patient,Clin$Patient)
	Clinical_extended = Clinical_extended[commonz,]	
	dtable(Clinical_extended$radiation_therapy) 
	dtable(Clinical_extended$targeted_molecular_therapy)
	Clin$treatment_targeted = NA
	Clin[rownames(Clinical_extended),"treatment_targeted"] = Clinical_extended$targeted_molecular_therapy
	dtable(Clin$treatment_targeted)
	Clint = Clin
	
	### Filling in Chen
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	chen_indir = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoSsignature_onChen/"
	load(file = paste0(chen_indir,"ClinChen.RData"))
	ClinChen = ClinChen[rownames(Clin),]
	dtable(ClinChen$TKI.treatment)
	dtable(ClinChen$Chemo.treatment)
	Clin$treatment_chemo = tolower(ClinChen$Chemo.treatment)
	Clin$treatment_targeted = tolower(ClinChen$TKI.treatment)
	Clinc = Clin

	### Filling in Chen
	load(file = paste0( OutDir,"../Preprocessing/","Clin_Chen.RData" ))
	chen_indir = "/mnt/ndata/daniele/lung_multiregion/rna-seq/Processed/PatternSignatures/LtoSsignature_onChen/"
	load(file = paste0(chen_indir,"ClinChen.RData"))
	ClinChen = ClinChen[rownames(Clin),]
	dtable(ClinChen$TKI.treatment)
	dtable(ClinChen$Chemo.treatment)
	Clin$treatment_chemo = tolower(ClinChen$Chemo.treatment)
	Clin$treatment_targeted = tolower(ClinChen$TKI.treatment)
	Clinc = Clin

	Cling$Dataset = "Genie"
	Cling$survival_data = !is.na(Cling$OS_status)
	Clint$Dataset = "TCGA"
	Clint$survival_data = !is.na(Clint$vital_status)
	Clinc$Dataset = "ChenEAS"
	Clinc$survival_data = !is.na(Clinc$OS.Status)

	### how many, across datasets and across EGFR variant type, have targeted treatment data and survival data
	colorz_classes = clco[rownames(clco) %in% (Cling$egfr_class_consensus),"colorz"]
	ordered_classes = clco[rownames(clco) %in% (Cling$egfr_class_consensus),"classes"]
	clinall = rbind(Cling[,c( "Dataset","Patient","egfr_class_consensus","treatment_targeted","survival_data" )],
					Clint[,c( "Dataset","Patient","egfr_class_consensus","treatment_targeted","survival_data" )],
					Clinc[,c( "Dataset","Patient","egfr_class_consensus","treatment_targeted","survival_data" )])
	save(clinall,file=paste0(OutDir,"clinall_treatments.RData"))
	clinall[is.na(clinall$treatment_targeted),"treatment_targeted"] = "unavailable"


	## overall 
	tc = clinall[clinall$Dataset=="Genie",]
	tc2 = tc[tc$treatment_targeted!="unavailable",]
	tabb = table(tc2[,"egfr_class_consensus"],tc2$treatment_targeted)
	ch = chisq.test(tabb)$p.value # no difference
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no","unavailable" ))
	colorz_groups = c( "forestgreen","orange","gray" )
	pg = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "Genie" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	tc = clinall[clinall$Dataset=="TCGA",]
	tc2 = tc[tc$treatment_targeted!="unavailable",]
	tabb = table(tc2[,"egfr_class_consensus"],tc2$treatment_targeted)
	ch = chisq.test(tabb)$p.value # no difference
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no","unavailable" ))
	colorz_groups = c( "forestgreen","orange","gray" )
	pt = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "TCGA" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	tc = clinall[clinall$Dataset=="ChenEAS",]
	tc2 = tc[tc$treatment_targeted!="unavailable",]
	tabb = table(tc2[,"egfr_class_consensus"],tc2$treatment_targeted)
	ch = chisq.test(tabb)$p.value # no difference
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no","unavailable" ))
	colorz_groups = c( "forestgreen","orange","gray" )
	pc = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "ChenEAS" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	tc = clinall
	tc2 = tc[tc$treatment_targeted!="unavailable",]
	tabb = table(tc2[,"egfr_class_consensus"],tc2$treatment_targeted)
	ch = chisq.test(tabb)$p.value # no difference
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no","unavailable" ))
	colorz_groups = c( "forestgreen","orange","gray" )
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "All datasets" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	pdf(paste0(OutDir,"TargetedTreatment_EGFRclasses_All_barplot.pdf"),14,5,onefile=FALSE)
	plot=ggarrange(pg, pt, pc, pa, common.legend=TRUE,ncol = 4, nrow = 1)
	print(plot)
	dev.off()

	## same but including only those with survival data and available therapy
	clinall = clinall[!(clinall$treatment_targeted=="unavailable") & ( clinall$survival_data ), ]
	tc = clinall[clinall$Dataset=="Genie",]
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no" ))
	colorz_groups = c( "forestgreen","orange" )
	pg = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "Genie" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	tc = clinall[clinall$Dataset=="TCGA",]
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no" ))
	colorz_groups = c( "forestgreen","orange" )
	pt = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "TCGA" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	tc = clinall[clinall$Dataset=="ChenEAS",]
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no" ))
	colorz_groups = c( "forestgreen","orange" )
	pc = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "ChenEAS" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	tc = clinall
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_targeted)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_targeted","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_targeted=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_targeted = factor(mt$treatment_targeted, levels=c( "yes","no" ))
	colorz_groups = c( "forestgreen","orange" )
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_targeted)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "All datasets" ) + scale_fill_manual(name = "Targeted treatment",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)

	pdf(paste0(OutDir,"TargetedTreatment_EGFRclasses_OnlySurvivalDataAvailable_barplot.pdf"),14,5,onefile=FALSE)
	plot=ggarrange(pg, pt, pc, pa, common.legend=TRUE,ncol = 4, nrow = 1)
	print(plot)
	dev.off()


	### same for immunotherapy (only genie)
	tc = Cling
	tc[is.na(tc$treatment_immuno),"treatment_immuno"] = "unavailable"
	tc2 = tc[tc$treatment_immuno!="unavailable",]
	tabb = table(tc2[,"egfr_class_consensus"],tc2$treatment_immuno)
	ch = chisq.test(tabb)$p.value # no difference
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_immuno)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_immuno","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_immuno = factor(mt$treatment_immuno, levels=c( "yes","no","unavailable" ))
	colorz_groups = c( "forestgreen","orange","gray" )
	pg = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_immuno)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "Genie" ) + scale_fill_manual(name = "Immunotherapy",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)
	pdf(paste0(OutDir,"Immunotherapy_EGFRclasses_Genie_barplot.pdf"),6,5,onefile=FALSE)
	print(pg)
	dev.off()

	### same for immunotherapy (only genie), no unavailable
	tc = Cling[!is.na(Cling$treatment_immuno),]
	tc2 = tc[tc$treatment_immuno!="unavailable",]
	tabb = table(tc2[,"egfr_class_consensus"],tc2$treatment_immuno)
	ch = chisq.test(tabb)$p.value # no difference
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_immuno)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_immuno","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_immuno = factor(mt$treatment_immuno, levels=c( "yes","no" ))
	colorz_groups = c( "forestgreen","orange" )
	pg = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_immuno)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( paste0("Genie, Chi-square p-val = ",signif(ch,2)) ) + scale_fill_manual(name = "Immunotherapy",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)
	pdf(paste0(OutDir,"Immunotherapy_EGFRclasses_Genie_NoUnavailable_barplot.pdf"),6,5,onefile=FALSE)
	print(pg)
	dev.off()

	### same for immunotherapy (only genie), no unavailable, only having survival data
	tc = Cling[(!is.na(Cling$treatment_immuno)) & (Cling$survival_data),]
	tc2 = tc[tc$treatment_immuno!="unavailable",]
	tabb = table(tc2[,"egfr_class_consensus"],tc2$treatment_immuno)
	ch = chisq.test(tabb)$p.value # no difference
	tabb = table(tc[,"egfr_class_consensus"],tc$treatment_immuno)
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","treatment_immuno","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$treatment_immuno=="no"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	mt$treatment_immuno = factor(mt$treatment_immuno, levels=c( "yes","no" ))
	colorz_groups = c( "forestgreen","orange" )
	pg = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=treatment_immuno)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle( "Genie" ) + scale_fill_manual(name = "Immunotherapy",values=colorz_groups) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)
	pdf(paste0(OutDir,"Immunotherapy_EGFRclasses_Genie_NoUnavailableWithSurvivalData_barplot.pdf"),6,5,onefile=FALSE)
	print(pg)
	dev.off()


	### Survivals

	################# survival
	### TCGA
	this_Clin = Clint[!is.na(Clint$treatment_targeted),]
	a = as.numeric(as.character(this_Clin$days_to_death))
	b = as.numeric(as.character(this_Clin$days_to_last_followup))
	vital_status_num <- vector(mode="numeric", length=length(this_Clin$vital_status))
	times <- vector(mode="numeric", length=length(vital_status_num))
	for (v in 1:length(vital_status_num))
	{
	 if (this_Clin$vital_status[v]=="alive")
	 {
	    vital_status_num[v] <- 0
	    times[v] <- b[v]
	 }
	 else
	 {
	    vital_status_num[v] <- 1
	    times[v] <- a[v]
	 }
	}
	this_Clin$Times <- times
	this_Clin$vital_status_num <- vital_status_num
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByTargetedTreatment_overall_TCGA.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("TCGA\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()
	this_Clin = this_Clin[this_Clin$egfr_class_consensus=="common",]
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByTargetedTreatment_CommonOnly_TCGA.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("TCGA\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()

	### Chen
	this_Clin = Clinc[!is.na(Clinc$treatment_targeted),]
	this_Clin$Times = as.numeric(this_Clin$OS.Month)
	this_Clin$vital_status_num = NA
	this_Clin[this_Clin$OS.Status=="Dead","vital_status_num"] = 1
	this_Clin[this_Clin$OS.Status=="Alive","vital_status_num"] = 0
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByTargetedTreatment_overall_ChenEAS.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("ChenEAS\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()
	this_Clin = this_Clin[this_Clin$egfr_class_consensus=="common",]
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByTargetedTreatment_CommonOnly_ChenEAS.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("ChenEAS\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()

	### Genie 
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[!is.na(this_Clin$treatment_targeted),]
	# this_Clin = this_Clin[!(this_Clin$Patient %in% patients_treated[1:32]),]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByTargetedTreatment_overall_Genie.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()
	this_Clin_save = this_Clin
	for (cl in c( "common","uncommon","compound" )){
		this_Clin = this_Clin_save[this_Clin_save$egfr_class_consensus==cl,]
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByTargetedTreatment_",cl,"Only_Genie.pdf")
		pdf(fileName,6,6,useDingbats=F)
		plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
		legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
		dev.off()
	}

	## of course, we need to account for stage. Not possible on TCGA, Chen (too few). Let's try stages III and IV in Genie.
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[(!is.na(this_Clin$treatment_targeted)) & ((this_Clin$Stage=="III") %in% c(T) ),]
	# this_Clin = this_Clin[!(this_Clin$Patient %in% patients_treated[1:32]),]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByTargetedTreatment_overall_GenieStageIII.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()
	this_Clin_save = this_Clin
	for (cl in c( "common" )){
		this_Clin = this_Clin_save[this_Clin_save$egfr_class_consensus==cl,]
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByTargetedTreatment_",cl,"Only_GenieStageIII.pdf")
		pdf(fileName,6,6,useDingbats=F)
		plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
		legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
		dev.off()
	}
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[(!is.na(this_Clin$treatment_targeted)) & ((this_Clin$Stage=="IV") %in% c(T) ),]
	# this_Clin = this_Clin[!(this_Clin$Patient %in% patients_treated[1:32]),]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByTargetedTreatment_overall_GenieStageIV.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()
	this_Clin_save = this_Clin
	for (cl in c( "common","uncommon","compound" )){
		this_Clin = this_Clin_save[this_Clin_save$egfr_class_consensus==cl,]
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_targeted = factor(this_Clin$treatment_targeted,levels=c("no","yes" )) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~treatment_targeted, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~treatment_targeted, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByTargetedTreatment_",cl,"Only_GenieStageIV.pdf")
		pdf(fileName,6,6,useDingbats=F)
		plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
		legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
		dev.off()
	}

	## Let's do it differently. Let's take all treated=yes and split into classes
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[!is.na(this_Clin$treatment_targeted),]
	this_Clin = this_Clin[this_Clin$treatment_targeted=="yes",]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByClasses_OnlyTargetedTreated_Genie.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Genie - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[!is.na(this_Clin$treatment_targeted),]
	this_Clin = this_Clin[this_Clin$treatment_targeted=="no",]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByClasses_OnlyTargetedUntreated_Genie.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Genie - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()

	## let's check immuno treatment. Only stage IV
	### Genie 
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[(!is.na(this_Clin$treatment_immuno)) & ((this_Clin$Stage=="IV")  ),]
	this_Clin = this_Clin[(!is.na(this_Clin$treatment_immuno)) & ((this_Clin$Stage=="IV") %in% c(T) ),]
	# this_Clin = this_Clin[!(this_Clin$Patient %in% patients_treated[1:32]),]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_immuno = factor(this_Clin$treatment_immuno,levels=c("no","yes" )) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~treatment_immuno, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~treatment_immuno, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByImmunotherapy_overall_GenieStageIV.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by immunotherapy, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
	dev.off()
	this_Clin_save = this_Clin
	for (cl in c( "common","uncommon","compound" )){
		this_Clin = this_Clin_save[this_Clin_save$egfr_class_consensus==cl,]
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, treatment_immuno = factor(this_Clin$treatment_immuno,levels=c("no","yes" )) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~treatment_immuno, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~treatment_immuno, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByImmunotherapy_",cl,"Only_GenieStageIV.pdf")
		pdf(fileName,6,6,useDingbats=F)
		plot(km_gs,mark.time=T, col=c( "orange","forestgreen" ), main = paste0("Genie\nsurvival by targeted treatment, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
		legend(x = "topright", legend = c( "no","yes" ) , lty=c(1,1,1), col=c( "orange","forestgreen" ))
		dev.off()
	}

	## Let's do it differently. Let's take all treated=yes and split into classes
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[(!is.na(this_Clin$treatment_immuno)) & ((this_Clin$Stage=="IV") %in% c(T) ),]
	this_Clin = this_Clin[this_Clin$treatment_immuno=="yes",]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByClasses_OnlyImmunoTreated_Genie.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Genie - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	this_Clin = this_Clin[(!is.na(this_Clin$treatment_immuno)) & (this_Clin$Stage=="IV"),]
	this_Clin = this_Clin[this_Clin$treatment_immuno=="no",]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "SurvivalByClasses_OnlyImmunoUntreated_Genie.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Genie - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()

	### Genie 
	this_Clin = Cling[!(is.na(Cling$OS_status)),]
	# this_Clin = this_Clin[((this_Clin$Stage=="IV") %in% c(T) ),]
	this_Clin$Times = as.numeric(this_Clin$OS_months)
	this_Clin$vital_status_num = as.numeric(substr(this_Clin$OS_status,1,1))
	Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
	Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
	km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
	km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
	p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
	fileName = paste0(OutDir, "Genie_KM_SurvivalBy_EGFRclass.pdf")
	pdf(fileName,6,6,useDingbats=F)
	plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Genie - survival by EGFR class, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
	legend(x = "topright", legend = ordered_classes, lty=c(1,1,1), col=colorz_classes)
	dev.off()
}

assign_drug_class = function( drug ){
	drug_class_list = list(
		chemo = c( "Carboplatin","Cisplatin","Cyclophosphamide","Dacarbazine","Docetaxel","Doxorubicin Hydrochloride","Etoposide","Fludarabine Phosphate","Fluorouracil","Gemcitabine Hydrochloride","Mitomycin","Nabpaclitaxel","Oxaliplatin","Paclitaxel","Paclitaxel Poliglumex","Pemetrexed Disodium","Temozolomide","Vinblastine Sulfate","Vinorelbine Tartrate" ),
		egfri = c( "Afatinib Dimaleate","Erlotinib Hydrochloride","Gefitinib","Osimertinib","Rociletinib" ),
		egfr_ab = c( "Cetuximab","Necitumumab" ),
		brafi = c( "Dabrafenib","Trametinib" ),
		alki = c( "Alectinib","Crizotinib" ),
		anti_angio = c( "Bevacizumab","Ramucirumab" ),
		hormonal = c( "Degarelix","Anastrozole","Lanreotide Acetate","Letrozole","Tamoxifen Citrate" ),
		immuno = c( "Atezolizumab","Avelumab","Durvalumab","Ipilimumab","Nivolumab","Pembrolizumab","Talimogene Laherparepvec" ),
		others_targeted = c( "Trastuzumab","Rituximab","Aldesleukin" ),
		others_unclear = c( "Investigational Drug","Other antineoplastic","Other NOS" )
		)
	tablez = dan.df(names(drug_class_list),c( "drug_class","drug" ))
	for (rn in rownames(tablez)){
		tablez[rn,"drug_class"] = rn
		tablez[rn,"drug"] = paste(drug_class_list[[rn]],collapse=", ")
	}
	dan.write(tablez,file=paste0(OutDir,"drug_classes.txt"),row.names=F)
	assigned_class = c()
	for (class in names(drug_class_list)){
		if (any(sapply(drug_class_list[[ class ]], function(x) grepl(x,drug)))){ assigned_class = c(assigned_class, class) }
	}
	if (length(assigned_class)>1){ 
		dcat( "Weird, more than one assigned class" ) 
		print(assigned_class)
	}
	if (length(assigned_class)==0){ assigned_class = NA }
	return(assigned_class)
}

assign_egfri_generation = function( drug ){
	drug_class_list = list(
		egfri_1stGen = c( "Erlotinib Hydrochloride","Gefitinib" ),
		egfri_2ndGen = c( "Afatinib Dimaleate"),
		egfri_3rdGen = c( "Osimertinib","Rociletinib" )
		)
	assigned_class = c()
	for (class in names(drug_class_list)){
		if (any(sapply(drug_class_list[[ class ]], function(x) grepl(x,drug)))){ assigned_class = c(assigned_class, class) }
	}
	if (length(assigned_class)>1){ 
		dcat( "Weird, more than one assigned class" ) 
		print(assigned_class)
	}
	if (length(assigned_class)==0){ assigned_class = NA }
	return(assigned_class)
}

regimen_associations_FirstDrug_PatientLevel = function( OutDir ){
	library(survival)
	### Filling in Genie
	load(file = paste0( OutDir,"../../Preprocessing/","Clin_Genie.RData" ))
	reg = read.csv(paste0(DataDir,"genie_BPC_NSCLC/regimen_cancer_level_dataset.csv"),stringsAsFactors=F)
	colz = c("record_id","regimen_number_within_cancer","drugs_dc_ynu","regimen_drugs","drugs_drug_1","drugs_drug_2","drugs_drug_3","drugs_drug_4","drugs_drug_5",
		"dx_drug_start_int_1","dx_drug_start_int_2","dx_drug_start_int_3","dx_drug_start_int_4","dx_drug_start_int_5","dx_drug_end_or_lastadm_int_1","dx_drug_end_or_lastadm_int_2","dx_drug_end_or_lastadm_int_3","dx_drug_end_or_lastadm_int_4","dx_drug_end_or_lastadm_int_5",
		"os_d_status","tt_os_d1_days","tt_os_d2_days","tt_os_d3_days","tt_os_d4_days","tt_os_d5_days","os_g_status","tt_os_g_days","pfs_i_or_m_g_status","tt_pfs_i_or_m_g_days")
	reg = reg[,colz]
	reg = reg[reg$record_id %in% rownames(Clin),]
	nn=unlist(strsplit(reg$regimen_drugs,split=", "))
	ass = c()
	for (n in nn){
		ass=c(ass,assign_drug_class(n))
	}
	dtable(ass)
	# Classify drug
	for (drug_id in c(1:5)){
		for (rn in rownames(reg)){
			reg[rn,paste0( "class_drug_",drug_id )] = ifelse((is.na(reg[rn,paste0( "drugs_drug_",drug_id )])),NA,assign_drug_class( reg[rn,paste0( "drugs_drug_",drug_id )] ))
		}
	}
	dtable(c(reg$class_drug_1,reg$class_drug_2,reg$class_drug_3,reg$class_drug_4,reg$class_drug_5) ) # UTTERLY PERFECT
	# barplots
	tc = Clin[(Clin$Patient %in% reg$record_id),]
	for (p in rownames(tc)){
		regp = reg[(reg$record_id==p),]
		all_drugs_classes = sort(unique(c(regp$class_drug_1,regp$class_drug_2,regp$class_drug_3,regp$class_drug_4,regp$class_drug_5)))
		all_drugs_classes = all_drugs_classes[!(all_drugs_classes=="others_unclear")]
		all_drugs_classes = paste(all_drugs_classes,collapse=",")
		tc[p,"all_drugs_classes"] = all_drugs_classes
	}
	tc = tc[grepl("egfri",tc$all_drugs_classes ),]
	for (p in rownames(tc)){
		regp = reg[(reg$record_id==p),]
		all_drugs_egfri = unlist(strsplit(sort(unlist(regp$regimen_drugs)),split=", " ))
		all_drugs_egfri = unique(sort(all_drugs_egfri[all_drugs_egfri %in% c( "Afatinib Dimaleate","Erlotinib Hydrochloride","Gefitinib","Osimertinib","Rociletinib" )]))
		all_drugs_egfri = sub(" .*", "", all_drugs_egfri)
		all_drugs_egfri = paste(all_drugs_egfri,collapse=",")
		tc[p,"all_drugs_egfri"] = all_drugs_egfri
	}
	tabz = dtable(tc$all_drugs_egfri,tc$egfr_class_consensus)
	keepz = names(which(rowSums(tabz>5)>0))
	tc[!(tc$all_drugs_egfri %in% keepz),"all_drugs_egfri" ] = "other_combinations"
	tabb = t(dtable(tc$all_drugs_egfri,tc$egfr_class_consensus))
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","EGFRi_combination","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=EGFRi_combination)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("EGFRi combination vs EGFR class\n","Chi-square test, p = ",signif(ch,2)) ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)
	pdf(paste0(OutDir,"EGFRclasses_EGFRicombination_Genie.pdf"),6,5,onefile=FALSE)
	print(pa)
	dev.off()

	theze = c( "uncommon","compound" )
	cox_MultiVar_df = dan.df( c("chemo","egfri","immuno","egfri_1stGen","egfri_2ndGen","egfri_3rdGen"),c( "Regimen","Npatients",paste0("HR_",theze), paste0("HR_lower95_",theze), paste0("HR_upper95_",theze),paste0("pval_",theze) ) )
	cox_UniVar_df = dan.df( c("chemo","egfri","immuno","egfri_1stGen","egfri_2ndGen","egfri_3rdGen"),c( "Regimen","Npatients",paste0("HR_",theze), paste0("HR_lower95_",theze), paste0("HR_upper95_",theze),paste0("pval_",theze) ) )

	for (this_treatmentclass in c( "chemo","egfri","immuno" ) ){ # anti_angio
		this_treatmentclass_alias = this_treatmentclass
		if (this_treatmentclass=="immuno") { this_treatmentclass_alias="immunotherapy" }
		if (this_treatmentclass=="egfri") { this_treatmentclass_alias="EGFR inhibitor (any)" }
		if (this_treatmentclass=="chemo") { this_treatmentclass_alias="Chemotherapy" }
		dcat( this_treatmentclass )
		reg$earlier = NA
		reg$earlier_os = NA
		reg$earlier_pfs = NA
		reg$earlier_tut = NA
		for (rn in rownames(reg)){
			i = 1
			while ((is.na(reg[rn,"earlier"])) & (i<=5) ){
				if ((reg[rn,paste0("class_drug_",i)]==this_treatmentclass) %in% c(T) ){
					reg[rn,"earlier"] = reg[rn,paste0("drugs_drug_",i)]
					reg[rn,"earlier_os"] = reg[rn,paste0("tt_os_d",i,"_days")]
					reg[rn,"earlier_pfs"] = reg[rn,paste0("tt_pfs_i_or_m_g_days")]
					reg[rn,"earlier_tut"] = reg[rn,paste0("dx_drug_end_or_lastadm_int_",i)]-reg[rn,paste0("dx_drug_start_int_",i)]
				}
				i=i+1
			}
		}
		reg$Patient = reg$record_id	
		for (rn in rownames(reg)){
			reg[rn,"egfr_class_consensus"] = Clin[reg[rn,"Patient"],"egfr_class_consensus"]
			reg[rn,"Stage"] = Clin[reg[rn,"Patient"],"Stage"]
			reg[rn,"Sex"] = Clin[reg[rn,"Patient"],"Sex"]
			reg[rn,"Age"] = Clin[reg[rn,"Patient"],"Age"]
		}
		this_Clin = reg[!is.na(reg$earlier),]
		this_Clin = this_Clin[order(this_Clin$Patient,this_Clin$regimen_number_within_cancer),]
		this_Clin = this_Clin[!(duplicated(this_Clin$Patient)),]
		colorz_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"colorz"]
		ordered_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"classes"]
		# OS
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin$vital_status_num = as.numeric(this_Clin$os_d_status)
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_OS.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F,pointsize=6)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Overall Survival by EGFR class\nGenie, p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Overall survival")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_OS_CoxRegression_CorrectingStage.txt"))
		# PFS
		this_Clin$Times = as.numeric(this_Clin$earlier_pfs)
		this_Clin$vital_status_num = as.numeric(this_Clin$pfs_i_or_m_g_status)
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df = Surv_df[!is.na(Surv_df$vital_status),]
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_PFS.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F,pointsize=6)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Progression-Free Survival by EGFR class\nGenie, p-val = ",signif(p.val,3)), xlab = "Time (days)", ylab = "Progression-Free Survival")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_PFS_CoxRegression_CorrectingStage.txt"))
		# Time under treatment (TUT) - event if treatment is discontinued or patient dies
		print(dtable(this_Clin$drugs_dc_ynu,is.na(this_Clin$earlier_tut)))
		this_Clin = this_Clin[(this_Clin$drugs_dc_ynu) %in% c("Yes","No"),]
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","Times"] = as.numeric(this_Clin[this_Clin$drugs_dc_ynu=="Yes","earlier_tut"])
		this_Clin$vital_status_num = 0
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","vital_status_num"] = 1
		this_Clin[this_Clin$os_d_status==1,"vital_status_num"] = 1
		this_Clin$egfr_class_consensus = factor(this_Clin$egfr_class_consensus,levels=ordered_classes)
		fileName = paste0(OutDir,"tut_boxplots_",this_treatmentclass,".pdf")
		dan.boxplots( fileName, this_Clin$egfr_class_consensus, this_Clin$Times, xlab = "", ylab = "Time Under Treatment", signifTest = "kruskal", labelycoo = max(this_Clin$Times), xColors = colorz_classes, jitterColors = dan.expand_colors(this_Clin$egfr_class_consensus,ordered_classes,colorz_classes), labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 3, fileHeight = 2.3, hlines_coo = NULL, hlines_labels = NULL )
		
		print(aggregate(Times~egfr_class_consensus,data=this_Clin,FUN='median'))

		Surv_df = data.frame(Patient = this_Clin$Patient, drugs_dc_ynu = this_Clin$drugs_dc_ynu, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage, Sex = this_Clin$Sex, Age = this_Clin$Age, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df = Surv_df[(Surv_df$drugs_dc_ynu) %in% c("Yes","No"),]
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		print(km_gs)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F,pointsize=6)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Time Under Treatment with ",this_treatmentclass_alias," \nby EGFR class, Genie, p-val = ",signif(p.val,2)), xlab = "Time (days)", ylab = "Treatment continuation probability")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage + Age + Sex")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_MultiVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_MultiVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_MultiVar_df[this_treatmentclass,"Npatients"] = a$n
	   	coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_UniVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_UniVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_UniVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_UniVar_df[this_treatmentclass,"Npatients"] = a$n
	}

	# Classify egfri generation
	for (drug_id in c(1:5)){
		for (rn in rownames(reg)){
			reg[rn,paste0( "class_drug_",drug_id )] = ifelse((is.na(reg[rn,paste0( "drugs_drug_",drug_id )])),NA,assign_egfri_generation( reg[rn,paste0( "drugs_drug_",drug_id )] ))
		}
	}
	reg = reg[!is.na(reg$class_drug_1),]

	for (this_treatmentclass in c( "egfri_1stGen","egfri_2ndGen","egfri_3rdGen" ) ){
		if (this_treatmentclass=="egfri_1stGen"){ this_treatmentclass_alias = "EGFR inhibitor (1st gen.)" }
		if (this_treatmentclass=="egfri_2ndGen"){ this_treatmentclass_alias = "EGFR inhibitor (2nd gen.)" }
		if (this_treatmentclass=="egfri_3rdGen"){ this_treatmentclass_alias = "EGFR inhibitor (3rd gen.)" }
		dcat( this_treatmentclass )
		reg$earlier = NA
		reg$earlier_os = NA
		reg$earlier_pfs = NA
		reg$earlier_tut = NA
		for (rn in rownames(reg)){
			i = 1
			while ((is.na(reg[rn,"earlier"])) & (i<=5) ){
				if ((reg[rn,paste0("class_drug_",i)]==this_treatmentclass) %in% c(T) ){
					reg[rn,"earlier"] = reg[rn,paste0("drugs_drug_",i)]
					reg[rn,"earlier_os"] = reg[rn,paste0("tt_os_d",i,"_days")]
					reg[rn,"earlier_pfs"] = reg[rn,paste0("tt_pfs_i_or_m_g_days")]
					reg[rn,"earlier_tut"] = reg[rn,paste0("dx_drug_end_or_lastadm_int_",i)]-reg[rn,paste0("dx_drug_start_int_",i)]
				}
				i=i+1
			}
		}
		reg$Patient = reg$record_id	
		for (rn in rownames(reg)){
			reg[rn,"egfr_class_consensus"] = Clin[reg[rn,"Patient"],"egfr_class_consensus"]
			reg[rn,"Stage"] = Clin[reg[rn,"Patient"],"Stage"]
		}
		this_Clin = reg[!is.na(reg$earlier),]
		this_Clin = this_Clin[order(this_Clin$Patient,this_Clin$regimen_number_within_cancer),]
		this_Clin = this_Clin[!(duplicated(this_Clin$Patient)),]
		colorz_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"colorz"]
		ordered_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"classes"]
		# OS
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin$vital_status_num = as.numeric(this_Clin$os_d_status)
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		print(km_gs)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_OS.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F,pointsize=6)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Overall Survival by EGFR class\nGenie, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_OS_CoxRegression_CorrectingStage.txt"))
		# PFS
		this_Clin$Times = as.numeric(this_Clin$earlier_pfs)
		this_Clin$vital_status_num = as.numeric(this_Clin$pfs_i_or_m_g_status)
		Surv_df = data.frame(Patient = this_Clin$Patient, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df = Surv_df[!is.na(Surv_df$vital_status),]
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_PFS.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F,pointsize=6)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Progression-Free Survival by EGFR class\nGenie, p-val = ",signif(p.val,3)), xlab = "Time (months)", ylab = "Survival")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_PFS_CoxRegression_CorrectingStage.txt"))
		# Time under treatment (TUT) - event if treatment is discontinued or patient dies
		print(dtable(this_Clin$drugs_dc_ynu,is.na(this_Clin$earlier_tut)))
		this_Clin = this_Clin[(this_Clin$drugs_dc_ynu) %in% c("Yes","No"),]
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","Times"] = as.numeric(this_Clin[this_Clin$drugs_dc_ynu=="Yes","earlier_tut"])
		this_Clin$vital_status_num = 0
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","vital_status_num"] = 1
		this_Clin[this_Clin$os_d_status==1,"vital_status_num"] = 1
		Surv_df = data.frame(Patient = this_Clin$Patient, drugs_dc_ynu = this_Clin$drugs_dc_ynu, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage,Stage = this_Clin$Stage, Sex = this_Clin$Sex, Age = this_Clin$Age, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F,pointsize=6)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Time Under Treatment with ",this_treatmentclass_alias," \nby EGFR class, Genie, p-val = ",signif(p.val,2)), xlab = "Time (months)", ylab = "Treatment continuation probability")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage + Age + Sex")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_MultiVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_MultiVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_MultiVar_df[this_treatmentclass,"Npatients"] = a$n
	   	coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_UniVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_UniVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_UniVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_UniVar_df[this_treatmentclass,"Npatients"] = a$n
	}
	library(forestplot)
	# clip = c(0.2,10)
	tabletext = cbind(
	     c("Regimen", cox_MultiVar_df$Regimen, "Summary"),
	     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
	     c("HR\n(uncommon)", signif(cox_MultiVar_df$HR_EGFR_classuncommon,3), signif(weighted.mean(x = cox_MultiVar_df$HR_EGFR_classuncommon, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(uncommon)", formatC(signif(cox_MultiVar_df$pval_EGFR_classuncommon,2)), NA),
	     c("HR\n(compound)", signif(cox_MultiVar_df$HR_EGFR_classcompound,3), signif(weighted.mean(x = cox_MultiVar_df$HR_EGFR_classcompound, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(compound)", formatC(signif(cox_MultiVar_df$pval_EGFR_classcompound,2)), NA),
	     rep("      ",length(c("Regimen", rownames(cox_MultiVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllRegimens.pdf"),8,nrow(cox_MultiVar_df)/1.5, useDingbats=F,onefile=F,pointsize=6)
	   fp=forestplot(tabletext, 
	            legend = c("uncommon","compound"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(5,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_MultiVar_df$HR_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_EGFR_classcompound,NA)), 
	            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_lower95_EGFR_classcompound,NA)), 
	            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_upper95_EGFR_classcompound,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("tomato3","mediumpurple4"),line=c("tomato3","mediumpurple4")))
	   print(fp)
	   dev.off()

	tabletext = cbind(
	     c("Regimen", cox_UniVar_df$Regimen, "Summary"),
	     c("# patients", cox_UniVar_df$Npatients, sum(cox_UniVar_df$Npatients)),
	     c("HR\n(uncommon)", signif(cox_UniVar_df$HR_EGFR_classuncommon,3), signif(weighted.mean(x = cox_UniVar_df$HR_EGFR_classuncommon, w = cox_UniVar_df$Npatients),3)),
	     c("p-value\n(uncommon)", formatC(signif(cox_UniVar_df$pval_EGFR_classuncommon,2)), NA),
	     c("HR\n(compound)", signif(cox_UniVar_df$HR_EGFR_classcompound,3), signif(weighted.mean(x = cox_UniVar_df$HR_EGFR_classcompound, w = cox_UniVar_df$Npatients),3)),
	     c("p-value\n(compound)", formatC(signif(cox_UniVar_df$pval_EGFR_classcompound,2)), NA),
	     rep("      ",length(c("Regimen", rownames(cox_UniVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_UniVar_df_ForestPlot_AllRegimens.pdf"),8,nrow(cox_UniVar_df)/1.5, useDingbats=F,onefile=F,pointsize=6)
	   fp=forestplot(tabletext, 
	            legend = c("uncommon","compound"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(5,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_UniVar_df$Npatients/max(cox_UniVar_df$Npatients),NA)/4,c(NA,cox_UniVar_df$Npatients/max(cox_UniVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_UniVar_df$HR_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_EGFR_classcompound,NA)), 
	            lower = cbind(c(NA,cox_UniVar_df$HR_lower95_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_lower95_EGFR_classcompound,NA)), 
	            upper = cbind(c(NA,cox_UniVar_df$HR_upper95_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_upper95_EGFR_classcompound,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_UniVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("tomato3","mediumpurple4"),line=c("tomato3","mediumpurple4")))
	   print(fp)
	   dev.off()
}

regimen_associations_FirstDrug_RegimenLevel = function( OutDir ){
	library(survival)
	### Filling in Genie
	load(file = paste0( OutDir,"../../Preprocessing/","Clin_Genie.RData" ))
	reg = read.csv(paste0(DataDir,"genie_BPC_NSCLC/regimen_cancer_level_dataset.csv"),stringsAsFactors=F)
	colz = c("record_id","regimen_number_within_cancer","drugs_dc_ynu","regimen_drugs","drugs_drug_1","drugs_drug_2","drugs_drug_3","drugs_drug_4","drugs_drug_5",
		"dx_drug_start_int_1","dx_drug_start_int_2","dx_drug_start_int_3","dx_drug_start_int_4","dx_drug_start_int_5","dx_drug_end_or_lastadm_int_1","dx_drug_end_or_lastadm_int_2","dx_drug_end_or_lastadm_int_3","dx_drug_end_or_lastadm_int_4","dx_drug_end_or_lastadm_int_5",
		"os_d_status","tt_os_d1_days","tt_os_d2_days","tt_os_d3_days","tt_os_d4_days","tt_os_d5_days","os_g_status","tt_os_g_days","pfs_i_or_m_g_status","tt_pfs_i_or_m_g_days")
	reg = reg[,colz]
	reg = reg[reg$record_id %in% rownames(Clin),]
	nn=unlist(strsplit(reg$regimen_drugs,split=", "))
	ass = c()
	for (n in nn){
		ass=c(ass,assign_drug_class(n))
	}
	dtable(ass)
	# Classify drug
	for (drug_id in c(1:5)){
		for (rn in rownames(reg)){
			reg[rn,paste0( "class_drug_",drug_id )] = ifelse((is.na(reg[rn,paste0( "drugs_drug_",drug_id )])),NA,assign_drug_class( reg[rn,paste0( "drugs_drug_",drug_id )] ))
		}
	}
	dtable(c(reg$class_drug_1,reg$class_drug_2,reg$class_drug_3,reg$class_drug_4,reg$class_drug_5) ) # UTTERLY PERFECT
	# barplots
	tc = Clin[(Clin$Patient %in% reg$record_id),]
	for (p in rownames(tc)){
		regp = reg[(reg$record_id==p),]
		all_drugs_classes = sort(unique(c(regp$class_drug_1,regp$class_drug_2,regp$class_drug_3,regp$class_drug_4,regp$class_drug_5)))
		all_drugs_classes = all_drugs_classes[!(all_drugs_classes=="others_unclear")]
		all_drugs_classes = paste(all_drugs_classes,collapse=",")
		tc[p,"all_drugs_classes"] = all_drugs_classes
	}
	tc = tc[grepl("egfri",tc$all_drugs_classes ),]
	for (p in rownames(tc)){
		regp = reg[(reg$record_id==p),]
		all_drugs_egfri = unlist(strsplit(sort(unlist(regp$regimen_drugs)),split=", " ))
		all_drugs_egfri = unique(sort(all_drugs_egfri[all_drugs_egfri %in% c( "Afatinib Dimaleate","Erlotinib Hydrochloride","Gefitinib","Osimertinib","Rociletinib" )]))
		all_drugs_egfri = sub(" .*", "", all_drugs_egfri)
		all_drugs_egfri = paste(all_drugs_egfri,collapse=",")
		tc[p,"all_drugs_egfri"] = all_drugs_egfri
	}
	tabz = dtable(tc$all_drugs_egfri,tc$egfr_class_consensus)
	keepz = names(which(rowSums(tabz>5)>0))
	tc[!(tc$all_drugs_egfri %in% keepz),"all_drugs_egfri" ] = "other_combinations"
	tabb = t(dtable(tc$all_drugs_egfri,tc$egfr_class_consensus))
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","EGFRi_combination","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=EGFRi_combination)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("EGFRi combination vs EGFR class\n","Chi-square test, p = ",signif(ch,2)) ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)
	pdf(paste0(OutDir,"EGFRclasses_EGFRicombination_Genie.pdf"),6,5,onefile=FALSE)
	print(pa)
	dev.off()

	theze = c( "uncommon","compound" )
	cox_MultiVar_df = dan.df( c("chemo","egfri","immuno","egfri_1stGen","egfri_2ndGen","egfri_3rdGen"),c( "Regimen","Npatients",paste0("HR_",theze), paste0("HR_lower95_",theze), paste0("HR_upper95_",theze),paste0("pval_",theze) ) )
	cox_UniVar_df = dan.df( c("chemo","egfri","immuno","egfri_1stGen","egfri_2ndGen","egfri_3rdGen"),c( "Regimen","Npatients",paste0("HR_",theze), paste0("HR_lower95_",theze), paste0("HR_upper95_",theze),paste0("pval_",theze) ) )

	for (this_treatmentclass in c( "chemo","egfri","immuno" ) ){ # anti_angio
		this_treatmentclass_alias = this_treatmentclass
		if (this_treatmentclass=="immuno") { this_treatmentclass_alias="immunotherapy" }
		if (this_treatmentclass=="egfri") { this_treatmentclass_alias="EGFR inhibitor (any)" }
		if (this_treatmentclass=="chemo") { this_treatmentclass_alias="Chemotherapy" }
		dcat( this_treatmentclass )
		reg$earlier = NA
		reg$earlier_os = NA
		reg$earlier_pfs = NA
		reg$earlier_tut = NA
		for (rn in rownames(reg)){
			i = 1
			while ((is.na(reg[rn,"earlier"])) & (i<=5) ){
				if ((reg[rn,paste0("class_drug_",i)]==this_treatmentclass) %in% c(T) ){
					reg[rn,"earlier"] = reg[rn,paste0("drugs_drug_",i)]
					reg[rn,"earlier_os"] = reg[rn,paste0("tt_os_d",i,"_days")]
					reg[rn,"earlier_pfs"] = reg[rn,paste0("tt_pfs_i_or_m_g_days")]
					reg[rn,"earlier_tut"] = reg[rn,paste0("dx_drug_end_or_lastadm_int_",i)]-reg[rn,paste0("dx_drug_start_int_",i)]
				}
				i=i+1
			}
		}
		reg$Patient = reg$record_id	
		for (rn in rownames(reg)){
			reg[rn,"egfr_class_consensus"] = Clin[reg[rn,"Patient"],"egfr_class_consensus"]
			reg[rn,"Stage"] = Clin[reg[rn,"Patient"],"Stage"]
			reg[rn,"Sex"] = Clin[reg[rn,"Patient"],"Sex"]
			reg[rn,"Age"] = Clin[reg[rn,"Patient"],"Age"]
		}
		this_Clin = reg[!is.na(reg$earlier),]
		this_Clin = this_Clin[order(this_Clin$Patient,this_Clin$regimen_number_within_cancer),]
		# this_Clin = this_Clin[!(duplicated(this_Clin$Patient)),]
		colorz_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"colorz"]
		ordered_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"classes"]
		# Time under treatment (TUT) - event if treatment is discontinued or patient dies
		print(dtable(this_Clin$drugs_dc_ynu,is.na(this_Clin$earlier_tut)))
		this_Clin = this_Clin[(this_Clin$drugs_dc_ynu) %in% c("Yes","No"),]
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","Times"] = as.numeric(this_Clin[this_Clin$drugs_dc_ynu=="Yes","earlier_tut"])
		this_Clin$vital_status_num = 0
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","vital_status_num"] = 1
		this_Clin[this_Clin$os_d_status==1,"vital_status_num"] = 1
		this_Clin$egfr_class_consensus = factor(this_Clin$egfr_class_consensus,levels=ordered_classes)
		fileName = paste0(OutDir,"tut_boxplots_",this_treatmentclass,".pdf")
		dan.boxplots( fileName, this_Clin$egfr_class_consensus, this_Clin$Times, xlab = "", ylab = "Time Under Treatment", signifTest = "kruskal", labelycoo = max(this_Clin$Times), xColors = colorz_classes, jitterColors = dan.expand_colors(this_Clin$egfr_class_consensus,ordered_classes,colorz_classes), labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 3, fileHeight = 2.3, hlines_coo = NULL, hlines_labels = NULL )

		Surv_df = data.frame(Patient = this_Clin$Patient, drugs_dc_ynu = this_Clin$drugs_dc_ynu, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage, Sex = this_Clin$Sex, Age = this_Clin$Age, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df = Surv_df[(Surv_df$drugs_dc_ynu) %in% c("Yes","No"),]
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Time Under Treatment with ",this_treatmentclass_alias," \nby EGFR class, Genie, p-val = ",signif(p.val,2)), xlab = "Time (days)", ylab = "Treatment continuation probability")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage + Age + Sex")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_MultiVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_MultiVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_MultiVar_df[this_treatmentclass,"Npatients"] = a$n
	   	coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_UniVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_UniVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_UniVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_UniVar_df[this_treatmentclass,"Npatients"] = a$n
	}

	# Classify egfri generation
	for (drug_id in c(1:5)){
		for (rn in rownames(reg)){
			reg[rn,paste0( "class_drug_",drug_id )] = ifelse((is.na(reg[rn,paste0( "drugs_drug_",drug_id )])),NA,assign_egfri_generation( reg[rn,paste0( "drugs_drug_",drug_id )] ))
		}
	}
	reg = reg[!is.na(reg$class_drug_1),]

	for (this_treatmentclass in c( "egfri_1stGen","egfri_2ndGen","egfri_3rdGen" ) ){
		if (this_treatmentclass=="egfri_1stGen"){ this_treatmentclass_alias = "EGFR inhibitor (1st gen.)" }
		if (this_treatmentclass=="egfri_2ndGen"){ this_treatmentclass_alias = "EGFR inhibitor (2nd gen.)" }
		if (this_treatmentclass=="egfri_3rdGen"){ this_treatmentclass_alias = "EGFR inhibitor (3rd gen.)" }
		dcat( this_treatmentclass )
		reg$earlier = NA
		reg$earlier_os = NA
		reg$earlier_pfs = NA
		reg$earlier_tut = NA
		for (rn in rownames(reg)){
			i = 1
			while ((is.na(reg[rn,"earlier"])) & (i<=5) ){
				if ((reg[rn,paste0("class_drug_",i)]==this_treatmentclass) %in% c(T) ){
					reg[rn,"earlier"] = reg[rn,paste0("drugs_drug_",i)]
					reg[rn,"earlier_os"] = reg[rn,paste0("tt_os_d",i,"_days")]
					reg[rn,"earlier_pfs"] = reg[rn,paste0("tt_pfs_i_or_m_g_days")]
					reg[rn,"earlier_tut"] = reg[rn,paste0("dx_drug_end_or_lastadm_int_",i)]-reg[rn,paste0("dx_drug_start_int_",i)]
				}
				i=i+1
			}
		}
		reg$Patient = reg$record_id	
		for (rn in rownames(reg)){
			reg[rn,"egfr_class_consensus"] = Clin[reg[rn,"Patient"],"egfr_class_consensus"]
			reg[rn,"Stage"] = Clin[reg[rn,"Patient"],"Stage"]
		}
		this_Clin = reg[!is.na(reg$earlier),]
		this_Clin = this_Clin[order(this_Clin$Patient,this_Clin$regimen_number_within_cancer),]
		# this_Clin = this_Clin[!(duplicated(this_Clin$Patient)),]
		colorz_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"colorz"]
		ordered_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"classes"]
		# Time under treatment (TUT) - event if treatment is discontinued or patient dies
		print(dtable(this_Clin$drugs_dc_ynu,is.na(this_Clin$earlier_tut)))
		this_Clin = this_Clin[(this_Clin$drugs_dc_ynu) %in% c("Yes","No"),]
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","Times"] = as.numeric(this_Clin[this_Clin$drugs_dc_ynu=="Yes","earlier_tut"])
		this_Clin$vital_status_num = 0
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","vital_status_num"] = 1
		this_Clin[this_Clin$os_d_status==1,"vital_status_num"] = 1
		Surv_df = data.frame(Patient = this_Clin$Patient, drugs_dc_ynu = this_Clin$drugs_dc_ynu, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage,Stage = this_Clin$Stage, Sex = this_Clin$Sex, Age = this_Clin$Age, RegimenNumber = this_Clin$regimen_number_within_cancer, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment.pdf")
		pdf(fileName,2.8,2.8,useDingbats=F)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Time Under Treatment with ",this_treatmentclass_alias," \nby EGFR class, Genie, p-val = ",signif(p.val,2)), xlab = "Time (months)", ylab = "Treatment continuation probability")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage + Age + Sex")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_MultiVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_MultiVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_MultiVar_df[this_treatmentclass,"Npatients"] = a$n
	   	coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_UniVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_UniVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_UniVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_UniVar_df[this_treatmentclass,"Npatients"] = a$n
	}
	library(forestplot)
	# clip = c(0.2,10)
	tabletext = cbind(
	     c("Regimen", cox_MultiVar_df$Regimen, "Summary"),
	     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
	     c("HR\n(uncommon)", signif(cox_MultiVar_df$HR_EGFR_classuncommon,3), signif(weighted.mean(x = cox_MultiVar_df$HR_EGFR_classuncommon, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(uncommon)", formatC(signif(cox_MultiVar_df$pval_EGFR_classuncommon,2)), NA),
	     c("HR\n(compound)", signif(cox_MultiVar_df$HR_EGFR_classcompound,3), signif(weighted.mean(x = cox_MultiVar_df$HR_EGFR_classcompound, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(compound)", formatC(signif(cox_MultiVar_df$pval_EGFR_classcompound,2)), NA),
	     rep("      ",length(c("Regimen", rownames(cox_MultiVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllRegimens.pdf"),6,nrow(cox_MultiVar_df)/2, useDingbats=F,onefile=F)
	   fp=forestplot(tabletext, 
	            legend = c("uncommon","compound"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(5,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_MultiVar_df$HR_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_EGFR_classcompound,NA)), 
	            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_lower95_EGFR_classcompound,NA)), 
	            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_upper95_EGFR_classcompound,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("tomato3","mediumpurple4"),line=c("tomato3","mediumpurple4")))
	   print(fp)
	   dev.off()

	tabletext = cbind(
	     c("Regimen", cox_UniVar_df$Regimen, "Summary"),
	     c("# patients", cox_UniVar_df$Npatients, sum(cox_UniVar_df$Npatients)),
	     c("HR\n(uncommon)", signif(cox_UniVar_df$HR_EGFR_classuncommon,3), signif(weighted.mean(x = cox_UniVar_df$HR_EGFR_classuncommon, w = cox_UniVar_df$Npatients),3)),
	     c("p-value\n(uncommon)", formatC(signif(cox_UniVar_df$pval_EGFR_classuncommon,2)), NA),
	     c("HR\n(compound)", signif(cox_UniVar_df$HR_EGFR_classcompound,3), signif(weighted.mean(x = cox_UniVar_df$HR_EGFR_classcompound, w = cox_UniVar_df$Npatients),3)),
	     c("p-value\n(compound)", formatC(signif(cox_UniVar_df$pval_EGFR_classcompound,2)), NA),
	     rep("      ",length(c("Regimen", rownames(cox_UniVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_UniVar_df_ForestPlot_AllRegimens.pdf"),6,nrow(cox_UniVar_df)/2, useDingbats=F,onefile=F)
	   fp=forestplot(tabletext, 
	            legend = c("uncommon","compound"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(5,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_UniVar_df$Npatients/max(cox_UniVar_df$Npatients),NA)/4,c(NA,cox_UniVar_df$Npatients/max(cox_UniVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_UniVar_df$HR_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_EGFR_classcompound,NA)), 
	            lower = cbind(c(NA,cox_UniVar_df$HR_lower95_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_lower95_EGFR_classcompound,NA)), 
	            upper = cbind(c(NA,cox_UniVar_df$HR_upper95_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_upper95_EGFR_classcompound,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_UniVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("tomato3","mediumpurple4"),line=c("tomato3","mediumpurple4")))
	   print(fp)
	   dev.off()
}

regimen_associations_SingleDrugLevel = function( OutDir ){
	library(survival)
	### Filling in Genie
	load(file = paste0( OutDir,"../../Preprocessing/","Clin_Genie.RData" ))
	reg = read.csv(paste0(DataDir,"genie_BPC_NSCLC/regimen_cancer_level_dataset.csv"),stringsAsFactors=F)
	colz = c("record_id","regimen_number_within_cancer","drugs_dc_ynu","regimen_drugs","drugs_drug_1","drugs_drug_2","drugs_drug_3","drugs_drug_4","drugs_drug_5",
		"dx_drug_start_int_1","dx_drug_start_int_2","dx_drug_start_int_3","dx_drug_start_int_4","dx_drug_start_int_5","dx_drug_end_or_lastadm_int_1","dx_drug_end_or_lastadm_int_2","dx_drug_end_or_lastadm_int_3","dx_drug_end_or_lastadm_int_4","dx_drug_end_or_lastadm_int_5",
		"os_d_status","tt_os_d1_days","tt_os_d2_days","tt_os_d3_days","tt_os_d4_days","tt_os_d5_days","os_g_status","tt_os_g_days","pfs_i_or_m_g_status","tt_pfs_i_or_m_g_days")
	reg = reg[,colz]
	reg = reg[reg$record_id %in% rownames(Clin),]
	nn=unlist(strsplit(reg$regimen_drugs,split=", "))
	ass = c()
	for (n in nn){
		ass=c(ass,assign_drug_class(n))
	}
	dtable(ass)
	# Classify drug
	for (drug_id in c(1:5)){
		for (rn in rownames(reg)){
			reg[rn,paste0( "class_drug_",drug_id )] = ifelse((is.na(reg[rn,paste0( "drugs_drug_",drug_id )])),NA,assign_drug_class( reg[rn,paste0( "drugs_drug_",drug_id )] ))
		}
	}
	dtable(c(reg$class_drug_1,reg$class_drug_2,reg$class_drug_3,reg$class_drug_4,reg$class_drug_5) ) # UTTERLY PERFECT
	# barplots
	tc = Clin[(Clin$Patient %in% reg$record_id),]
	for (p in rownames(tc)){
		regp = reg[(reg$record_id==p),]
		all_drugs_classes = sort(unique(c(regp$class_drug_1,regp$class_drug_2,regp$class_drug_3,regp$class_drug_4,regp$class_drug_5)))
		all_drugs_classes = all_drugs_classes[!(all_drugs_classes=="others_unclear")]
		all_drugs_classes = paste(all_drugs_classes,collapse=",")
		tc[p,"all_drugs_classes"] = all_drugs_classes
	}
	tc = tc[grepl("egfri",tc$all_drugs_classes ),]
	for (p in rownames(tc)){
		regp = reg[(reg$record_id==p),]
		all_drugs_egfri = unlist(strsplit(sort(unlist(regp$regimen_drugs)),split=", " ))
		all_drugs_egfri = unique(sort(all_drugs_egfri[all_drugs_egfri %in% c( "Afatinib Dimaleate","Erlotinib Hydrochloride","Gefitinib","Osimertinib","Rociletinib" )]))
		all_drugs_egfri = sub(" .*", "", all_drugs_egfri)
		all_drugs_egfri = paste(all_drugs_egfri,collapse=",")
		tc[p,"all_drugs_egfri"] = all_drugs_egfri
	}
	tabz = dtable(tc$all_drugs_egfri,tc$egfr_class_consensus)
	keepz = names(which(rowSums(tabz>5)>0))
	tc[!(tc$all_drugs_egfri %in% keepz),"all_drugs_egfri" ] = "other_combinations"
	tabb = t(dtable(tc$all_drugs_egfri,tc$egfr_class_consensus))
	ch = chisq.test(tabb)$p.value
	tabb = t(apply(tabb,1, function(x) x/sum(x)))*100
	mt = melt(tabb)
	colnames(mt) = c( "EGFR_class","EGFRi_combination","Percentage" )
	mt$TotalPatients = NA
	mt[(mt$EGFR_class=="common") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="common"))
	mt[(mt$EGFR_class=="uncommon") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="uncommon"))
	mt[(mt$EGFR_class=="compound") & (mt$EGFRi_combination=="Erlotinib"),"TotalPatients"] = paste0( "N = ",sum(tc$egfr_class_consensus=="compound"))
	mt$EGFR_class = factor(mt$EGFR_class, levels=ordered_classes)
	pa = ggplot(data=mt, aes(x=EGFR_class, y=Percentage, fill=EGFRi_combination)) +
	  geom_bar(stat="identity", colour="black", linewidth = 0.1) + ylab("Percentage of patients") + xlab("") + ggtitle(paste0("EGFRi combination vs EGFR class\n","Chi-square test, p = ",signif(ch,2)) ) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12)) +# scale_x_discrete(labels=rownames(ldf)) +
	  geom_text(aes(label=TotalPatients,y=100),vjust=-0.2)
	pdf(paste0(OutDir,"EGFRclasses_EGFRicombination_Genie.pdf"),6,5,onefile=FALSE)
	print(pa)
	dev.off()

	theze = c( "uncommon","compound" )
	cox_MultiVar_df = dan.df( c("chemo","egfri","immuno","egfri_1stGen","egfri_2ndGen","egfri_3rdGen"),c( "Regimen","Npatients",paste0("HR_",theze), paste0("HR_lower95_",theze), paste0("HR_upper95_",theze),paste0("pval_",theze) ) )
	cox_UniVar_df = dan.df( c("chemo","egfri","immuno","egfri_1stGen","egfri_2ndGen","egfri_3rdGen"),c( "Regimen","Npatients",paste0("HR_",theze), paste0("HR_lower95_",theze), paste0("HR_upper95_",theze),paste0("pval_",theze) ) )

	for (this_treatmentclass in c( "chemo","egfri","immuno" ) ){ # anti_angio
		this_treatmentclass_alias = this_treatmentclass
		if (this_treatmentclass=="immuno") { this_treatmentclass_alias="immunotherapy" }
		if (this_treatmentclass=="egfri") { this_treatmentclass_alias="EGFR inhibitor (any)" }
		if (this_treatmentclass=="chemo") { this_treatmentclass_alias="Chemotherapy" }
		dcat( this_treatmentclass )
		reg_expanded = dan.df(0,c("record_id","drugs_dc_ynu","earlier","earlier_os","earlier_tut" ))
		for (rn in rownames(reg)){
			i = 1
			while ((i<=5)){
				if ((reg[rn,paste0("class_drug_",i)]==this_treatmentclass) %in% c(T) ){
					thisreg = data.frame(record_id=reg[rn,"record_id"],drugs_dc_ynu=reg[rn,"drugs_dc_ynu"],earlier=reg[rn,paste0("drugs_drug_",i)],
						earlier_os=reg[rn,paste0("tt_os_d",i,"_days")],earlier_tut=reg[rn,paste0("dx_drug_end_or_lastadm_int_",i)]-reg[rn,paste0("dx_drug_start_int_",i)],stringsAsFactors=F)
					reg_expanded = rbind(reg_expanded,thisreg)
				}
				i=i+1
			}
		}
		reg_expanded$Patient = reg_expanded$record_id
		for (rn in rownames(reg_expanded)){
			reg_expanded[rn,"egfr_class_consensus"] = Clin[reg_expanded[rn,"Patient"],"egfr_class_consensus"]
			reg_expanded[rn,"Stage"] = Clin[reg_expanded[rn,"Patient"],"Stage"]
			reg_expanded[rn,"Sex"] = Clin[reg_expanded[rn,"Patient"],"Sex"]
			reg_expanded[rn,"Age"] = Clin[reg_expanded[rn,"Patient"],"Age"]
		}
		this_Clin = reg_expanded[!is.na(reg_expanded$earlier),]
		# this_Clin = this_Clin[!(duplicated(this_Clin$Patient)),]
		colorz_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"colorz"]
		ordered_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"classes"]
		# Time under treatment (TUT) - event if treatment is discontinued or patient dies
		print(dtable(this_Clin$drugs_dc_ynu,is.na(this_Clin$earlier_tut)))
		this_Clin = this_Clin[(this_Clin$drugs_dc_ynu) %in% c("Yes","No"),]
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","Times"] = as.numeric(this_Clin[this_Clin$drugs_dc_ynu=="Yes","earlier_tut"])
		this_Clin$vital_status_num = 0
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","vital_status_num"] = 1
		this_Clin[this_Clin$os_d_status==1,"vital_status_num"] = 1
		this_Clin$egfr_class_consensus = factor(this_Clin$egfr_class_consensus,levels=ordered_classes)
		fileName = paste0(OutDir,"tut_boxplots_",this_treatmentclass,".pdf")
		dan.boxplots( fileName, this_Clin$egfr_class_consensus, this_Clin$Times, xlab = "", ylab = "Time Under Treatment", signifTest = "kruskal", labelycoo = max(this_Clin$Times), xColors = colorz_classes, jitterColors = dan.expand_colors(this_Clin$egfr_class_consensus,ordered_classes,colorz_classes), labelJitteredPoints = NULL, jitterDotSize = 1.5, fileWidth = 3, fileHeight = 2.3, hlines_coo = NULL, hlines_labels = NULL )

		Surv_df = data.frame(Patient = this_Clin$Patient, drugs_dc_ynu = this_Clin$drugs_dc_ynu, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage, Sex = this_Clin$Sex, Age = this_Clin$Age, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df = Surv_df[(Surv_df$drugs_dc_ynu) %in% c("Yes","No"),]
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment.pdf")
		pdf(fileName,5.3,5.3,useDingbats=F)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Time Under Treatment with ",this_treatmentclass_alias," \nby EGFR class, Genie, p-val = ",signif(p.val,2)), xlab = "Time (days)", ylab = "Treatment continuation probability")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage + Age + Sex")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_MultiVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_MultiVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_MultiVar_df[this_treatmentclass,"Npatients"] = a$n
	   	coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_UniVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_UniVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_UniVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_UniVar_df[this_treatmentclass,"Npatients"] = a$n
	}

	# Classify egfri generation
	for (drug_id in c(1:5)){
		for (rn in rownames(reg)){
			reg[rn,paste0( "class_drug_",drug_id )] = ifelse((is.na(reg[rn,paste0( "drugs_drug_",drug_id )])),NA,assign_egfri_generation( reg[rn,paste0( "drugs_drug_",drug_id )] ))
		}
	}
	reg = reg[!is.na(reg$class_drug_1),]

	for (this_treatmentclass in c( "egfri_1stGen","egfri_2ndGen","egfri_3rdGen" ) ){
		if (this_treatmentclass=="egfri_1stGen"){ this_treatmentclass_alias = "EGFR inhibitor (1st gen.)" }
		if (this_treatmentclass=="egfri_2ndGen"){ this_treatmentclass_alias = "EGFR inhibitor (2nd gen.)" }
		if (this_treatmentclass=="egfri_3rdGen"){ this_treatmentclass_alias = "EGFR inhibitor (3rd gen.)" }
		dcat( this_treatmentclass )
		reg_expanded = dan.df(0,c("record_id","drugs_dc_ynu","earlier","earlier_os","earlier_tut" ))
		for (rn in rownames(reg)){
			i = 1
			while ((i<=5)){
				if ((reg[rn,paste0("class_drug_",i)]==this_treatmentclass) %in% c(T) ){
					thisreg = data.frame(record_id=reg[rn,"record_id"],drugs_dc_ynu=reg[rn,"drugs_dc_ynu"],earlier=reg[rn,paste0("drugs_drug_",i)],
						earlier_os=reg[rn,paste0("tt_os_d",i,"_days")],earlier_tut=reg[rn,paste0("dx_drug_end_or_lastadm_int_",i)]-reg[rn,paste0("dx_drug_start_int_",i)],stringsAsFactors=F)
					reg_expanded = rbind(reg_expanded,thisreg)
				}
				i=i+1
			}
		}
		reg_expanded$Patient = reg_expanded$record_id
		for (rn in rownames(reg_expanded)){
			reg_expanded[rn,"egfr_class_consensus"] = Clin[reg_expanded[rn,"Patient"],"egfr_class_consensus"]
			reg_expanded[rn,"Stage"] = Clin[reg_expanded[rn,"Patient"],"Stage"]
			reg_expanded[rn,"Sex"] = Clin[reg_expanded[rn,"Patient"],"Sex"]
			reg_expanded[rn,"Age"] = Clin[reg_expanded[rn,"Patient"],"Age"]
		}
		this_Clin = reg_expanded[!is.na(reg_expanded$earlier),]
		# this_Clin = this_Clin[!(duplicated(this_Clin$Patient)),]
		colorz_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"colorz"]
		ordered_classes = clco[rownames(clco) %in% (this_Clin$egfr_class_consensus),"classes"]
		# Time under treatment (TUT) - event if treatment is discontinued or patient dies
		print(dtable(this_Clin$drugs_dc_ynu,is.na(this_Clin$earlier_tut)))
		this_Clin = this_Clin[(this_Clin$drugs_dc_ynu) %in% c("Yes","No"),]
		this_Clin$Times = as.numeric(this_Clin$earlier_os)
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","Times"] = as.numeric(this_Clin[this_Clin$drugs_dc_ynu=="Yes","earlier_tut"])
		this_Clin$vital_status_num = 0
		this_Clin[this_Clin$drugs_dc_ynu=="Yes","vital_status_num"] = 1
		this_Clin[this_Clin$os_d_status==1,"vital_status_num"] = 1
		Surv_df = data.frame(Patient = this_Clin$Patient, drugs_dc_ynu = this_Clin$drugs_dc_ynu, vital_status = as.numeric(this_Clin$vital_status_num), Times = this_Clin$Times, Stage = this_Clin$Stage,Stage = this_Clin$Stage, Sex = this_Clin$Sex, Age = this_Clin$Age, EGFR_class = factor(this_Clin$egfr_class_consensus,levels=ordered_classes) )
		Surv_df$SurvObj = with(Surv_df, Surv(Times, as.numeric(as.character(vital_status))))
		km_gs = survfit(SurvObj~EGFR_class, data = Surv_df)
		km_gs_dif = survdiff(SurvObj~EGFR_class, data = Surv_df, rho = 0)
		p.val = 1 - pchisq(km_gs_dif$chisq, length(km_gs_dif$n) - 1)
		fileName = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment.pdf")
		pdf(fileName,5.3,5.3,useDingbats=F)
		plot(km_gs,mark.time=T, col=colorz_classes, main = paste0("Time Under Treatment with ",this_treatmentclass_alias," \nby EGFR class, Genie, p-val = ",signif(p.val,2)), xlab = "Time (months)", ylab = "Treatment continuation probability")
		legend(x = "topright", legend = paste0(ordered_classes," (N=",as.numeric(dtable(Surv_df$EGFR_class)[ ordered_classes ]),")"), lty=c(1,1,1), col=colorz_classes)
		dev.off()
		coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class + Stage + Age + Sex")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_MultiVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_MultiVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_MultiVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_MultiVar_df[this_treatmentclass,"Npatients"] = a$n
	   	coxr = coxph(as.formula(paste0("SurvObj ~ EGFR_class")), data = Surv_df) #  Age + Sex +
	   	a = (summary(coxr))
	   	capture.output(a, file = paste0(OutDir, "SurvivalByClasses_OnlyTreatedWith_",this_treatmentclass,"_Genie_TimeUnderTreatment_CoxRegression_CorrectingStage.txt"))
	   	for (feature in paste0("EGFR_class",theze)){
	   		cox_UniVar_df[this_treatmentclass,paste0("HR_",feature)] = a$coefficients[feature,"exp(coef)"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_lower95_",feature)] = a$conf.int[feature,"lower .95"]
			cox_UniVar_df[this_treatmentclass,paste0("HR_upper95_",feature)] = a$conf.int[feature,"upper .95"]
			cox_UniVar_df[this_treatmentclass,paste0("pval_",feature)] = a$coefficients[feature,"Pr(>|z|)"]
	   	}
	   	cox_UniVar_df[this_treatmentclass,"Regimen"] = this_treatmentclass_alias
	   	cox_UniVar_df[this_treatmentclass,"Npatients"] = a$n
	}
	library(forestplot)
	# clip = c(0.2,10)
	tabletext = cbind(
	     c("Regimen", cox_MultiVar_df$Regimen, "Summary"),
	     c("# patients", cox_MultiVar_df$Npatients, sum(cox_MultiVar_df$Npatients)),
	     c("HR\n(uncommon)", signif(cox_MultiVar_df$HR_EGFR_classuncommon,3), signif(weighted.mean(x = cox_MultiVar_df$HR_EGFR_classuncommon, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(uncommon)", formatC(signif(cox_MultiVar_df$pval_EGFR_classuncommon,2)), NA),
	     c("HR\n(compound)", signif(cox_MultiVar_df$HR_EGFR_classcompound,3), signif(weighted.mean(x = cox_MultiVar_df$HR_EGFR_classcompound, w = cox_MultiVar_df$Npatients),3)),
	     c("p-value\n(compound)", formatC(signif(cox_MultiVar_df$pval_EGFR_classcompound,2)), NA),
	     rep("      ",length(c("Regimen", rownames(cox_MultiVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_MultiVar_df_ForestPlot_AllRegimens.pdf"),12,nrow(cox_MultiVar_df), useDingbats=F,onefile=F)
	   fp=forestplot(tabletext, 
	            legend = c("uncommon","compound"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(5,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4,c(NA,cox_MultiVar_df$Npatients/max(cox_MultiVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_MultiVar_df$HR_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_EGFR_classcompound,NA)), 
	            lower = cbind(c(NA,cox_MultiVar_df$HR_lower95_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_lower95_EGFR_classcompound,NA)), 
	            upper = cbind(c(NA,cox_MultiVar_df$HR_upper95_EGFR_classuncommon,NA),c(NA,cox_MultiVar_df$HR_upper95_EGFR_classcompound,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_MultiVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("tomato3","mediumpurple4"),line=c("tomato3","mediumpurple4")))
	   print(fp)
	   dev.off()

	tabletext = cbind(
	     c("Regimen", cox_UniVar_df$Regimen, "Summary"),
	     c("# patients", cox_UniVar_df$Npatients, sum(cox_UniVar_df$Npatients)),
	     c("HR\n(uncommon)", signif(cox_UniVar_df$HR_EGFR_classuncommon,3), signif(weighted.mean(x = cox_UniVar_df$HR_EGFR_classuncommon, w = cox_UniVar_df$Npatients),3)),
	     c("p-value\n(uncommon)", formatC(signif(cox_UniVar_df$pval_EGFR_classuncommon,2)), NA),
	     c("HR\n(compound)", signif(cox_UniVar_df$HR_EGFR_classcompound,3), signif(weighted.mean(x = cox_UniVar_df$HR_EGFR_classcompound, w = cox_UniVar_df$Npatients),3)),
	     c("p-value\n(compound)", formatC(signif(cox_UniVar_df$pval_EGFR_classcompound,2)), NA),
	     rep("      ",length(c("Regimen", rownames(cox_UniVar_df), "Summary"))))
	   pdf(paste0(OutDir,"cox_UniVar_df_ForestPlot_AllRegimens.pdf"),12,nrow(cox_UniVar_df), useDingbats=F,onefile=F)
	   fp=forestplot(tabletext, 
	            legend = c("uncommon","compound"),
	            align = c("l","l"),
	            fn.ci_norm = c(fpDrawCircleCI,fpDrawCircleCI),
	           hrzl_lines = gpar(col="#444444"),
	           colgap = unit(5,"mm"),
	            boxsize = 1.5*cbind(c(NA,cox_UniVar_df$Npatients/max(cox_UniVar_df$Npatients),NA)/4,c(NA,cox_UniVar_df$Npatients/max(cox_UniVar_df$Npatients),NA)/4),
	            mean = cbind(c(NA,cox_UniVar_df$HR_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_EGFR_classcompound,NA)), 
	            lower = cbind(c(NA,cox_UniVar_df$HR_lower95_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_lower95_EGFR_classcompound,NA)), 
	            upper = cbind(c(NA,cox_UniVar_df$HR_upper95_EGFR_classuncommon,NA),c(NA,cox_UniVar_df$HR_upper95_EGFR_classcompound,NA)),
	           is.summary=c(TRUE,rep(FALSE,length(cox_UniVar_df$Npatients)),TRUE),
	           xlog=T,
	           # clip = clip,
	           xlab = "Hazard ratio",
	           col=fpColors(box=c("tomato3","mediumpurple4"),line=c("tomato3","mediumpurple4")))
	   print(fp)
	   dev.off()
}

DepMap_preprocessing = function( OutDir ){
	# DepMap Public 23Q4, downloaded from https://depmap.org/portal/download/all/ on 28/04/2024
	dep = read.csv(paste0(DataDir,"DepMap_23Q4/Model.csv"),stringsAsFactors=F)
	dep = dep[dep$OncotreeCode=="LUAD",]
	muts = read.csv(paste0(DataDir,"DepMap_23Q4/OmicsSomaticMutations.csv"),stringsAsFactors=F)
	muts = muts[muts$ModelID %in% dep$ModelID,]
	muts$egfr_class_consensus = "wt/passenger"
	rownames(muts) = paste0("u",rownames(muts) )
	dep = dep[dep$ModelID %in% muts$ModelID,]
	# let's classify EGFR variants
	mutse = muts[((muts$HugoSymbol=="EGFR") %in% c(T)) & (nchar(muts$ProteinChange)>1),]
	stringz = gsub("\\_.*","",mutse$ProteinChange )
	stringz = gsub("\\*.*","",stringz)
	stringz = as.numeric(gsub("\\D", "", stringz))
	mutse$First_Codon = stringz
	aa = sapply(1:length(mutse$ProteinChange),function(x) gsub(stringz[x],"",mutse$ProteinChange[x]))
	stringz = as.numeric(gsub("\\D", "", aa))
	mutse$Last_Codon = stringz
	mutse[is.na(mutse$Last_Codon) | ( (mutse$Last_Codon<=mutse$First_Codon) %in% c(T) ),"Last_Codon" ] = mutse[is.na(mutse$Last_Codon) | ( (mutse$Last_Codon<=mutse$First_Codon) %in% c(T) ),"First_Codon" ]
	mutse = mutse[order(mutse$First_Codon),]
	mutse = mutse[(mutse$First_Codon>=688) & (mutse$First_Codon<=875), ]
	mutse[(mutse$VariantInfo=="inframe_deletion") & ( (mutse$First_Codon<=755) | (mutse$Last_Codon>=746) ),"egfr_class_consensus"] = "common"
	mutse[mutse$ProteinChange=="p.L858R","egfr_class_consensus"] = "common"
	mutse[mutse$ProteinChange=="p.T790M","egfr_class_consensus"] = "T790M"
	mutse[(mutse$VariantInfo=="inframe_insertion") & ( (mutse$First_Codon<=774) | (mutse$Last_Codon>=762) ),"egfr_class_consensus"] = "ex20ins"
	mutse[mutse$egfr_class_consensus=="wt/passenger","egfr_class_consensus"] = "uncommon" # the only uncommon, p.A750P, is a true uncommon (found in ec)
	rownames(dep) = dep$ModelID
	for (rn in rownames(dep)){
		if (!(rn %in% mutse$ModelID)){ 
			dep[rn,"egfr_class_consensus"] = "wt/passenger"
		} else {
			all_muts = mutse[mutse$ModelID==rn,"egfr_class_consensus"]
			all_muts_unf = all_muts
			all_muts = all_muts[all_muts %in% c( "common","uncommon" )]
			if (length(all_muts)>0){
				if ((length(unique(all_muts))==1) & (unique(all_muts)=="common")){ 
				dep[rn,"egfr_class_consensus"] = "common" 
				} else {
					dep[rn,"egfr_class_consensus"] = "compound" 
				}
				if ((length(all_muts)==1) & (unique(all_muts)=="uncommon")){ dep[rn,"egfr_class_consensus"] = "uncommon" }
			}
			if (any(tolower(all_muts_unf)=="t790m")) { dep[rn,"egfr_class_consensus"] = "T790M" }
			if (any(all_muts_unf=="ex20ins")) { dep[rn,"egfr_class_consensus"] = "ex20ins" }
		}
	}
	cri = read.csv(paste0(DataDir,"DepMap_23Q4/CRISPRGeneEffect.csv"),stringsAsFactors=F)
	rownames(cri) = cri$X
	cri = cri[rownames(cri) %in% dep$ModelID,]
	dep$egfr_crispr_effect = NA
	dep[rownames(cri),"egfr_crispr_effect"] = cri[,grepl("EGFR",colnames(cri))]
	save(dep,file = paste0(OutDir,"dep.RData"))
}

DepMap_analysis = function( OutDir ){
	load(file = paste0(OutDir,"dep.RData"))
	dep = dep[!is.na(dep$egfr_crispr_effect),]
	x = dep$egfr_class_consensus=="wt/passenger"
	x2 = x
	x2[x] = "EGFR wild-type or \npassenger mutation"
	x2[!x] = "EGFR mutant"
	x2 = factor(x2, levels = c( "EGFR wild-type or \npassenger mutation","EGFR mutant" ))
	y = dep$egfr_crispr_effect
	colorz = rep("gray",nrow(dep))
	colorz[dep$egfr_class_consensus=="common"] = "steelblue4"
	colorz[dep$egfr_class_consensus=="T790M"] = "darksalmon"
	dan.boxplots(paste0( OutDir,"EGFRclass_vs_EGFRdependency.pdf" ),x2,y,xlab="",ylab="Gene dependency score",signifTest="wilcox",labelycoo=max(y),jitterColors=colorz,fileWidth = 5, fileHeight = 5)
}

supplementary_tables_formatting = function( OutDir ){

	load(file=paste0( OutDir,"../Genomic_coMutations/AllDatasets_signifComutated_table.RData" ))
	dan.write(cdf,file = paste0(OutDir, "TableS1.txt"))
	load(file=paste0( OutDir,"../Genomic_coMutations/EGFRclasses_vs_DriverMuts_AllDatasets_MatchingTmb.RData" ))
	dan.write(tcdf,file = paste0(OutDir, "TableS2.txt"))
	load(file = paste0(OutDir,"../Expression/Limma_table_AllDatasets_AllGenes_AllVsAll_comparisonsSignLevel0.1.RData"))
	dan.write(final,file = paste0(OutDir, "TableS3.txt"))
}








