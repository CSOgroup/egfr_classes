
library(ggplot2)
library(ggpubr)
library(gplots)
library(plyr)
library(ggrepel)

dan.save = function( object, file )
{
	formatted_filename = paste0(gsub("\\..*","",file),".RData")
	save(object, file = formatted_filename)
	formatted_filename = paste0(gsub("\\..*","",file),".txt")
	dan.write(object,formatted_filename)
}

dan.read = function( file, row.names = NULL, header = T )
{
	rt = read.table(file = file, header = header, row.names = row.names, stringsAsFactors = F, sep = "\t", quote = '')
	return( rt )
}

dan.write = function( table, file, row.names = F )
{
	write.table(table, file = file, row.names = row.names, col.names = T, sep = "\t", quote = F)
	return()
}

peek = function( table, nrow = 5, ncol = 5 )
{
	nrow = min(c(nrow,nrow(table)))
	ncol = min(c(ncol,ncol(table)))
	print(table[1:nrow,1:ncol])
	print(paste0(nrow(table)," x ",ncol(table) ))
}

dan.colors = function( n, ggplot2_stile = T ){
	if (length(n)>1) { n = length(unique(n)) }
	if (ggplot2_stile){
		hues = seq(15, 375, length = n + 1)
		hcl(h = hues, l = 65, c = 100)[1:n]	
	}
}

dcat = function( string, tabb = 0 )
{
	cat("\n", paste0(rep(".",tabb*5),collapse=""), string, "\n" )
}

ddup = function( df, column )
{
	dupl = unique(df[,column][duplicated(df[,column])] )
	df = df[df[,column] %in% dupl,]
	df = df[order(df[,column]),]
	return(df)
}

dtable = function(...){
   return(table(...,useNA = "ifany"))
}

dan.expand_colors = function( vec, levelz, colorz ){
	colorz_vec = rep(NA,length(vec))
	for (l in levelz){
		colorz_vec[vec==l] = colorz[levelz==l]
	}
	return(colorz_vec)
}

dan.rowMedians = function( table, na.rm = T )
{
	return( apply(table,1,median, na.rm = na.rm) )
}

dan.colMedians = function( table, na.rm = T )
{
	return( apply(table,2,median, na.rm = na.rm) )
}

dan.boxplots = function( fileName, x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, jitterDotSize = 4.5, fileWidth = 4, fileHeight = 3, hlines_coo = NULL, hlines_labels = NULL, legend_position = NULL )
{
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(fill))
	{
		if (!("default" %in% fillColors) )
		{
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,fatten=3) + scale_fill_manual(values = fillColors) + geom_point(pch = 16, size = jitterDotSize, position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=6),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		} else {
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,fatten=3) + geom_point(pch = 16, size = jitterDotSize, position = position_jitterdodge(jitter.width = 0.2,jitter.height = 0)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=6),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		}		
	} else {
		# cat(xlab,ylab) label = "p.format",
		p = ggplot(mapping = aes(y = y, x = x)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,fatten=3) +
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
			theme_classic(base_size=6) + theme(text = element_text(size=6),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
		
		if (!is.null(labelJitteredPoints))
		{
			p = p + geom_text(aes(label=labelJitteredPoints),hjust=0, vjust=0, size = 6/.pt, fontface = "bold", position = position_jitter(width = 0.2, height = 0), color = jitterColors)
		} else {
			p = p + geom_jitter(width=.2,height=0, size = jitterDotSize, alpha = 0.5, color = jitterColors) + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
		}
	}
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo, size = 6/.pt)}
	if (!is.null(hlines_coo))
	{
		p = p + geom_hline( yintercept = hlines_coo, linetype="dashed", color = "gray44" )# + geom_text( aes(0.5, hlines_coo, label = hlines_labels, vjust = -1, hjust = -1), color = "gray44")
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(legend_position)){
		p = p + theme(legend.position=legend_position)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.boxplots.multipages = function( x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	if (!is.null(fill))
	{
		if (!("default" %in% fillColors) )
		{
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,linewidth = 0.2) + scale_fill_manual(values = fillColors) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		} else {
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,linewidth = 0.2) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		}
		if (includeJitters){
				p = p + geom_point(position=position_jitterdodge(jitter.width = 0.1,jitter.height = 0),size=0.2)
			}
	} else {
		# cat(xlab,ylab) label = "p.format",
		p = ggplot(mapping = aes(y = y, x = x)) + geom_boxplot(size=0.3,color = xColors, outlier.shape = NA,linewidth = 0.2) +
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
			theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
		
		if (!is.null(labelJitteredPoints))
		{
			p = p + geom_text(aes(label=labelJitteredPoints),hjust=0, vjust=0, size = 6/.pt, position = position_jitter(width = 0.2, height = 0), color = jitterColors)
		} else {
			if (includeJitters){
				p = p + geom_jitter(width=.1, color = jitterColors,size=0.2)	
			}
		}
	}
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo, label = "p.format", size = 6/.pt)}
	if (!is.null(hlines_coo))
	{
		p = p + geom_hline( yintercept = hlines_coo, linetype="dashed", color = "gray44" ) + geom_text( aes(0.5, hlines_coo, label = hlines_labels, vjust = -1), color = "gray44",size = 6/.pt)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(legend_position)){
		p = p + theme(legend.position=legend_position)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}

dan.violinplots.multipages = function( x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", signifTest = "kruskal", ylimLeft = NULL, ylimRight = NULL,comparisons = NULL, labelycoo = 1, xColors = "black", fillColors = "default", jitterColors = "black", labelJitteredPoints = NULL, includeJitters = T, hlines_coo = NULL, hlines_labels = NULL,legend_position=NULL )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	if (!is.null(fill))
	{
		if (!("default" %in% fillColors) )
		{
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_violin(linewidth=0.2, outlier.shape = NA) + scale_fill_manual(values = fillColors) + scale_color_manual(values=xColors) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		} else {
			p = ggplot(mapping = aes(y = y, x = x, fill = fill)) + geom_violin(linewidth=0.2, outlier.shape = NA) + scale_color_manual(values=xColors) +# geom_point(pch = 16, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
				ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) +
				theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  #+ theme(axis.title.x=element_blank(), axis.text.y = element_text(size = 12))
		}
		if (includeJitters){
				p = p + geom_point(position=position_jitterdodge(jitter.width = 0.2,jitter.height = 0))
			}
	} else {
		# cat(xlab,ylab) label = "p.format",
		p = ggplot(mapping = aes(y = y, x = x, color=x, show.legend = FALSE)) + geom_violin(linewidth=0.2, outlier.shape = NA) + scale_color_manual(values=xColors) +
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + 
			theme_classic(base_size=6) + theme(text = element_text(size=16),axis.text.x = element_text(angle = 45, hjust = 1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
		
		if (!is.null(labelJitteredPoints))
		{
			p = p + geom_text(aes(label=labelJitteredPoints),hjust=0, vjust=0, size = 6/.pt, position = position_jitter(width = 0.2, height = 0), color = jitterColors)
		} else {
			if (includeJitters){
				p = p + geom_jitter(width=.1, color = jitterColors)	
			}
		}
	}
	if (!(is.null(signifTest))) { p = p + stat_compare_means(comparisons = comparisons, method = signifTest, label.y = labelycoo, label = "p.format")}
	if (!is.null(hlines_coo))
	{
		p = p + geom_hline( yintercept = hlines_coo, linetype="dashed", color = "gray44" ) + geom_text( aes(0.5, hlines_coo, label = hlines_labels, vjust = -1), color = "gray44",size = 6/.pt)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	p = p + stat_summary(fun=median, geom="point", size=.1, color=xColors,position = position_dodge(0.9) )
	if (!is.null(legend_position)){
		p = p + theme(legend.position=legend_position)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}

dan.barplots.multipages = function( x, y, fill = NULL, sd = NULL, xlab = "default", ylab = "default", ylimLeft = NULL, ylimRight = NULL, filllab = "default", fillColors = NULL, labelBarsBottom = NULL, plotTitle = "")
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.text.x = element_text(angle = 90, hjust = 1))
		# p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
	} else {
		p = ggplot(mapping = aes(y = y, x = factor(x))) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45)) + theme_classic(base_size=6)
	}
	if (!is.null(fillColors))
	{
		p = p + scale_fill_manual(values = fillColors)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(sd))
	{
		p = p + geom_errorbar(aes(ymin=y-sd, ymax=y+sd), size = 0.2, alpha = 0.5, stat="identity", position=position_dodge2(width = 0.9, preserve = "single"))+ theme(text = element_text(size=20))
	}
	if (!is.null(labelBarsBottom))
	{
		p = p + geom_text(label = labelBarsBottom, y = y + 1*sign(y),size = 6/.pt)#, position = position_dodge2(width=0.9), size=2.7) # (0-(max(y)-min(y))/10)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}

dan.scatterplots.multipages = function( x, y, fill = NULL, xlab = "default", ylab = "default", filllab = "default", plotTitle = "", dotSize = 1, fillColors = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T )
{
	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	} else {
		p = ggplot(mapping = aes(y = y, x = x)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + theme_classic(base_size=6)
	}
	if (!is.null(fillColors))
	{
		p = p + scale_color_manual( values=fillColors,drop = FALSE )
	}
	p = p + geom_point(size = dotSize, alpha = 0.7)
	if (!is.null(plotFitLine)){
		p = p + geom_smooth(method=FitLineMethod, color = FitLineColor, se = plotFitLine_se)	
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	return(p)
}

dan.densityPlot = function( fileName, x, grouping, groupinglab = "default", xlab = "default", ylab = "default", show_means = F, show_medians = F, plotTitle = "",xlimLeft = NULL, xlimRight = NULL, groupingColors = "firebrick", fileWidth = 4, fileHeight = 3 )
{

	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = "Density" }
	if (groupinglab=="default") { groupinglab = deparse(substitute(grouping)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	p = ggplot(mapping = aes(x = x, color = grouping)) + geom_density() + scale_color_manual(values = groupingColors, name = groupinglab ) + 
			ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) +
			theme_classic(base_size=6) + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
	if (show_means)
	{
		mu = data.frame(grp = levels(grouping), grp.mean = NA)
		for (g in mu$grp)
		{
			mu[mu$grp==g,"grp.mean"] = mean(x[as.character(grouping)==g])
		}
		p = p + geom_vline(data=mu, aes(xintercept=grp.mean, color=grp), linetype="dashed")
	}
	if (show_medians)
	{
		mu = data.frame(grp = levels(grouping), grp.mean = NA)
		for (g in mu$grp)
		{
			mu[mu$grp==g,"grp.median"] = median(x[as.character(grouping)==g])
		}
		p = p + geom_vline(data=mu, aes(xintercept=grp.median, color=grp), linetype="dashed")
	}
	if (!is.null(xlimLeft))
	{
		p = p + xlim(xlimLeft, xlimRight)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.scatterplot = function( fileName, x, y, fill = NULL, xlab = "default", ylab = "default", xlimLeft = NULL, xlimRight = NULL, ylimLeft = NULL, ylimRight = NULL, filllab = "default", fillColors = NULL, fillColors_continuous = NULL, plotTitle = "", dotLabels = NULL, dotSize = 1, plotVline = NULL, plotHline = NULL, plotBisector = NULL, plotFitLine = NULL, FitLineMethod = "lm", FitLineColor = "firebrick", plotFitLine_se = T, repel_labels = NULL, coord_fixed = FALSE, coord_flipped = FALSE, fileWidth = 4, fileHeight = 3 )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = x, color = fill)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = filllab ) + theme_classic(base_size=6)
	} else {
		p = ggplot(mapping = aes(y = y, x = x)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + theme_classic(base_size=6)
	}
	if (!is.null(dotLabels))
	{
		p = p + geom_point(shape=NA) + geom_text( label = dotLabels, size = 6/.pt, show.legend = F ) + guides(color = guide_legend(override.aes = list(shape = 19, size = 5)))
	} else {
		p = p + geom_point(size = dotSize, alpha = 0.7)
	}
	if (!is.null(plotBisector))
	{
		if (plotBisector) { p = p + geom_abline( slope=1,intercept=0 ) }
	}
	if (!is.null(plotVline))
	{
		p = p + geom_vline( xintercept = plotVline,linetype="dashed",size=0.2 )
	}
	if (!is.null(plotHline))
	{
		p = p + geom_hline( yintercept = plotHline,linetype="dashed",size=0.2 )
	}
	if (!is.null(plotFitLine))
	{
		if (plotFitLine) { p = p + geom_smooth(method=FitLineMethod, color = FitLineColor, se = plotFitLine_se) }
	}
	if (!is.null(fillColors))
	{
		p = p + scale_color_manual( values=fillColors,drop = FALSE )
	}
	if (!is.null(fillColors_continuous))
	{
		p = p + scale_colour_gradientn( colours=fillColors_continuous )
	}
	if (!is.null(xlimLeft))
	{
		p = p + xlim(xlimLeft, xlimRight)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(repel_labels))
	{
		p = p + geom_label_repel(aes(label = repel_labels), min.segment.length = 0, force=15,size = 6/.pt,max.overlaps=20000,show.legend=FALSE )
	}
	p = p+ theme(text = element_text(size=14))
	if (coord_fixed){
		p = p + coord_fixed()
	}
	if (coord_flipped){
		p = p + coord_flip()
	}
	# p = p + guides(color=guide_legend(ncol=1,override.aes = list(size=2)))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(6, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.lineplot = function( fileName, x, y, group = NULL, xlab = "default", ylab = "default", xlimLeft = NULL, xlimRight = NULL, ylimLeft = NULL, ylimRight = NULL, grouplab = "default", groupColors = NULL, lineShowLegend = T, dotColors = NULL, plotTitle = "", fileWidth = 4, fileHeight = 3 )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (grouplab=="default") { grouplab = deparse(substitute(group)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(group))
	{
		if (!is.null(dotColors))
		{
			p = ggplot(mapping = aes(y = y, x = x, group = group)) + geom_line(aes(color=group), show.legend = lineShowLegend) + geom_point(color=dotColors, size = 5) + geom_text(aes(label=group),color = "gray77",hjust=-0.5, vjust=-0.5,size = 6/.pt) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = grouplab ) + theme_classic(base_size=6)
		} else {
			p = ggplot(mapping = aes(y = y, x = x, group = group)) + geom_line(aes(color=group), show.legend = lineShowLegend) + geom_point(aes(color=group)) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = grouplab ) + theme_classic(base_size=6)
		}
		
	} else {
		p = ggplot(mapping = aes(y = y, x = x, group = 1)) + geom_line() + geom_point() + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( color = grouplab ) + theme_classic(base_size=6)
	}
	if (!is.null(groupColors))
	{
		p = p + scale_color_manual( values=groupColors )
	}
	if (!is.null(xlimLeft))
	{
		p = p + xlim(xlimLeft, xlimRight)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	p = p+ theme(text = element_text(size=14))
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.barplot = function( fileName, x, y, fill = NULL, sd = NULL, xlab = "default", ylab = "default", ylimLeft = NULL, ylimRight = NULL, filllab = "default", fillColors = NULL, labelBarsBottom = NULL, textOnTop = NULL, noxtext=FALSE, plotTitle = "", fileWidth = 4, fileHeight = 3 )
{
	
	if (xlab=="default") { xlab = deparse(substitute(x)) }
	if (ylab=="default") { ylab = deparse(substitute(y)) }
	if (filllab=="default") { filllab = deparse(substitute(fill)) }

	pdf( fileName, width = fileWidth, height = fileHeight , useDingbats = F )
	if (!is.null(fill))
	{
		p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=14), axis.text.x = element_text(angle = 45, hjust = 1))
		if (noxtext){ p = ggplot(mapping = aes(y = y, x = factor(x), fill = fill)) + geom_bar(stat="identity", position="identity") + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=12), axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) }
	} else {
		p = ggplot(mapping = aes(y = y, x = factor(x))) + geom_bar(stat="identity", position=position_dodge2(width = 0.9, preserve = "single")) + ggtitle( plotTitle ) + xlab( xlab ) + ylab( ylab ) + labs( fill = filllab ) + theme_classic(base_size=6) + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1))
	}
	if (!is.null(fillColors))
	{
		p = p + scale_fill_manual(values = fillColors)
	}
	if (!is.null(ylimLeft))
	{
		p = p + ylim(ylimLeft, ylimRight)
	}
	if (!is.null(sd))
	{
		p = p + geom_errorbar(aes(ymin=y-sd, ymax=y+sd), size = 0.2, alpha = 0.5, stat="identity", position=position_dodge2(width = 0.9, preserve = "single"))+ theme(text = element_text(size=20))
	}
	if (!is.null(labelBarsBottom))
	{
		p = p + geom_text(label = labelBarsBottom, y = y + 1*sign(y),size = 6/.pt)#, position = position_dodge2(width=0.9), size=2.7) # (0-(max(y)-min(y))/10)
	}
	p = p + theme(text = element_text(size=6), legend.text=element_text(size=6), legend.title = element_text(size = 6), plot.title = element_text(size = 6), axis.title = element_text(size = 6), axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), legend.key.size = unit(10, 'pt'))
	print(p)
	dev.off()
	return()
}

dan.df = function(rownames, colnames, data = NA, as_df = T){
	if ((length(rownames)==1)) { mat = matrix(nrow = 0, ncol = length(colnames), dimnames = list(NULL,colnames)) }
	if ((length(rownames)!=1)) { mat = matrix(nrow = length(rownames), ncol = length(colnames), dimnames = list(rownames,colnames), data = data) }
	mat = data.frame(mat, stringsAsFactors = F)
	colnames(mat) = colnames
	return( mat )
}
