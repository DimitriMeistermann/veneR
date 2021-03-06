heatmap.DM<-function(matrix,clustering_distance_rows=corrDist,clustering_distance_columns="euclidean",clustering_method_columns="ward.D2",
										 clustering_method_rows="ward.D2",autoFontSizeRow=TRUE,autoFontSizeColumn=TRUE,scale=TRUE,center=TRUE,returnHeatmap=FALSE,name="Z-scored\nexpression",
										 additionnalRowNamesGpar=list(fontface="italic"),additionnalColNamesGpar=list(),colorScale=NULL,cluster_rows=NULL,cluster_columns=NULL,...){
	require("ComplexHeatmap")
	require("circlize")
	args<-list()
	if(is.null(cluster_rows)){
		args$clustering_method_rows<-clustering_method_rows
		args$clustering_distance_rows<-clustering_distance_rows
	}else{
		args$cluster_rows<-cluster_rows
	}
	if(is.null(cluster_columns)){
		args$clustering_method_columns<-clustering_method_columns
		args$clustering_distance_columns<-clustering_distance_columns
	}else{
		args$cluster_columns<-cluster_columns
	}
	matrix<-as.matrix(matrix)
	if(min(apply(matrix,1,sd))==0){
		warning("some row have a 0 sd. sd-based method (correlation distance, scaling) will be desactivated or switched.")
		scale=FALSE
		if(identical(corrDist,clustering_distance_rows)){
			args$clustering_distance_rows<-"euclidean"
		}
	}
	
	if(scale | center) matrix<-rowScale(matrix,scaled=scale,center=center)
	
	if(is.null(colorScale)){
		quantiles<-quantile(unlist(matrix),probs=c(0.05,0.25,0.5,0.75,0.95))
		colorScale<-colorRamp2(breaks=c(quantiles),colors=c("blue","cyan","white","yellow","red"))
	}
	args$col<-colorScale
	args$name=name
	args<-c(args,list(...))
	
	if(autoFontSizeRow) args$row_names_gp=do.call("autoGparFontSizeMatrix",c(list(nrow(matrix)),additionnalRowNamesGpar))
	if(autoFontSizeColumn) args$column_names_gp=do.call("autoGparFontSizeMatrix",c(list(ncol(matrix)),additionnalColNamesGpar))
	
	args$matrix<-matrix
	ht<-do.call("Heatmap",args)
	if(returnHeatmap){
		return(ht)
	}else{
		print(ht)
	}
}


heatmap.DM2<-function(matrix,clustering_distance_rows=corrDist,clustering_distance_columns="euclidean",clustering_method_columns="ward.D2",
											clustering_method_rows="ward.D2",autoFontSizeRow=TRUE,autoFontSizeColumn=TRUE,scale=FALSE,center=FALSE,returnHeatmap=FALSE,name="Z-scored\nexpression",
											additionnalRowNamesGpar=list(),additionnalColNamesGpar=list(),colorScale=NULL,cluster_rows=NULL,cluster_columns=NULL,...){
	require("ComplexHeatmap")
	require("circlize")
	args<-list()
	if(is.null(cluster_rows)){
		args$clustering_method_rows<-clustering_method_rows
		args$clustering_distance_rows<-clustering_distance_rows
	}else{
		args$cluster_rows<-cluster_rows
	}
	if(is.null(cluster_columns)){
		args$clustering_method_columns<-clustering_method_columns
		args$clustering_distance_columns<-clustering_distance_columns
	}else{
		args$cluster_columns<-cluster_columns
	}
	matrix<-as.matrix(matrix)
	if(min(apply(matrix,1,sd))==0 & (scale | identical(corrDist,clustering_distance_rows)) ){
		warning("some row have a 0 sd. sd-based method (correlation distance, scaling) will be desactivated or switched.")
		scale=FALSE
		if(identical(corrDist,clustering_distance_rows)){
			args$clustering_distance_rows<-"euclidean"
		}
	}
	
	if(scale | center) matrix<-rowScale(matrix,scaled=scale,center=center)
	
	if(is.null(colorScale)){
		quantiles<-quantile(unlist(matrix),probs=c(0.05,0.25,0.5,0.75,0.95))
		colorScale<-colorRamp2(breaks=c(quantiles),colors=c("blue","cyan","white","yellow","red"))
	}
	args$col<-colorScale
	args$name=name
	args<-c(args,list(...))
	
	if(autoFontSizeRow) args$row_names_gp=do.call("autoGparFontSizeMatrix",c(list(nrow(matrix)),additionnalRowNamesGpar))
	if(autoFontSizeColumn) args$column_names_gp=do.call("autoGparFontSizeMatrix",c(list(ncol(matrix)),additionnalColNamesGpar))
	
	args$matrix<-matrix
	ht<-do.call("Heatmap",args)
	if(returnHeatmap){
		return(ht)
	}else{
		print(ht)
	}
}


make.umap<-function(data,PCAfirst=FALSE,nDimPCA=NULL,transpose=TRUE,init="spectral",n_neighbors=11, n_components=2, n_epochs=10000, min_dist=0.1, 
										set_op_mix_ratio=1, local_connectivity=1, bandwidth=1, alpha=1, gamma=1, 
										negative_sample_rate=5, spread=1, random_state=NA, transform_state=NA, knn_repeats=1, verbose=F){
	
	#n_neighbors: integer; number of nearest neighbors
	#n_components: integer; dimension of target (output) space
	#metric: character or function; determines how distances between data points are computed. When using a string, available metrics are: euclidean, manhattan. Other availble generalized metrics are: cosine, pearson, pearson2. Note the triangle inequality may not be satisfied by some generalized metrics, hence knn search may not be optimal. When using metric.function as a function, the signature must be function(matrix, origin, target) and should compute a distance between the origin column and the target columns
	#n_epochs: integer; number of iterations performed during layout optimization
	#input: character, use either "data" or "dist"; determines whether the primary input argument to umap() is treated as a data matrix or as a distance matrix
	#init: character or matrix. The default string "spectral" computes an initial embedding using eigenvectors of the connectivity graph matrix. An alternative is the string "random", which creates an initial layout based on random coordinates. This setting.can also be set to a matrix, in which case layout optimization begins from the provided coordinates.
	#min_dist: numeric; determines how close points appear in the final layout
	#set_op_ratio_mix_ratio: numeric in range [0,1]; determines who the knn-graph is used to create a fuzzy simplicial graph
	#local_connectivity: numeric; used during construction of fuzzy simplicial set
	#bandwidth: numeric; used during construction of fuzzy simplicial set
	#alpha: numeric; initial value of "learning rate" of layout optimization
	#beta: numeric; determines, together with alpha, the learning rate of layout optimization
	#negative_sample_rate: integer; determines how many non-neighbor points are used per point and per iteration during layout optimization
	#a: numeric; contributes to gradient calculations during layout optimization. When left at NA, a suitable value will be estimated automatically.
	#b: numeric; contributes to gradient calculationss during layout optimization.
	#spread: numeric; used during automatic estimation of a/b parameters.
	#random_state: integer; seed for random number generation used during umap()
	#transform_state: integer; seed for random number generation used during predict()
	#knn.repeat: number of times to restart knn search
	#verbose: logical or integer; determines whether to show progress messages
	#umap_learn_args: vector of arguments to python package umap-learn 
	require(umap)
	umap.config<-umap.defaults
	umap.config$init<-init
	umap.config$n_neighbors<-n_neighbors
	umap.config$n_components<-n_components
	umap.config$n_epochs<-n_epochs
	umap.config$min_dist<-min_dist
	umap.config$set_op_mix_ratio<-set_op_mix_ratio
	umap.config$local_connectivity<-local_connectivity
	umap.config$bandwidth<-bandwidth
	umap.config$alpha<-alpha
	umap.config$gamma<-gamma
	umap.config$negative_sample_rate<-negative_sample_rate
	umap.config$spread<-spread
	umap.config$random_state<-random_state
	umap.config$transform_state<-transform_state
	umap.config$knn_repeats<-knn_repeats
	umap.config$verbose<-verbose
	if(transpose) data<-t(data)
	if(PCAfirst){
		data<-ACP(data,transpose = FALSE,scale = FALSE)$x
		if(!is.null(nDimPCA)){
			data<-data[,1:nDimPCA]
		}
	}
	umap(as.matrix(data),config = umap.config)
}

#Return gene id correspondance, GO species code and KEGG species code
getSpeciesData<-function(sample.species="Human",genes,updateSpeciesPackage=FALSE){
	require(gage)
	require(AnnotationDbi)
	data(bods)
	
	species<-list()
	species.data<-data.frame(bods)
	species.index<-which(species.data$species==sample.species)
	if(len(species.index)!=1) stop("Wrong species name, type \ndata(bods)\nbods[,\"species\"]\nto see available species")
	species$package<-as.character(species.data[species.index,"package"])
	species$kegg<-as.character(species.data[species.index,"kegg.code"])
	species$go<-strsplit(as.character(species$package),split = ".",fixed = TRUE)[[1]][2]
	if(updateSpeciesPackage | !(require(species$package,character.only = TRUE))){
		require("BiocManager")
		print(paste0("Downloading species package: ",species.data$package))
		install(species$package, update=FALSE)
	}
	require(species$package,character.only = TRUE)
	
	print(genes)
	species$GeneIdTable<-AnnotationDbi::select(get(species$package),genes, "ENTREZID","SYMBOL")
	species$species<-sample.species
	return(species)
}

#use reshape2::melt instead
meltMatrix<-function(x, metricName="expression", colTitle="sample", rowTitle="gene"){
	tabGraph<-data.frame(unlist(data.frame(x)))
	colnames(tabGraph)<-metricName
	tabGraph[,rowTitle]<-rep(rownames(x),ncol(x))
	tabGraph[,colTitle]<-rep(colnames(x),each=nrow(x))
	tabGraph
}

#apply a function that take 2 argument on a vactor of n element, recursively, ! Deprecated, use Reduce function instead !
recursiveFun<-function(argList,fun){
	if(len(argList)==2) return(fun(argList[1][[1]],argList[2][[1]]))
	return(fun(argList[1][[1]],recursiveFun(argList[-1],fun)))
}



#' Easy and	quick view of expression
#' @param expr dataframe or matrix. Expression.
#' @param conditions dataframe of factors.	Sample table with sample in row and annotations in column.
#' @param legendName character. Custom legend name.
#' @param errorBar character. What represent the bars around each point. "se" : standard error mean, "ci": confidance interval (distribution must follow a normal law) "na":none
#' @param ciRate numeric. Confidance interval rate if errorBar = ci, 0.95 = CI 95%
#' @param geom character. GGplot function name for representation, Ex : geom="point" if geom_point is wanted.
#' @param addLine logical. Add extra line to the plot ?
#' @param xaxis character. 'gene' for representing genes and 'annot' for annotation on x axis (remaining parameter will be represented with color)
#' @param negValue logical. error bars sub zero ?
#' @param scale character. 'identity' or 'log' scaled ?
#' @param breaks numeric. Position of breaks on y axis.
#' @param xRotationLab numeric. Rotation of x axis labels
#' @param hjust numeric. Horizontal justification
#' @param main character. Main title
#' @param xLabelSize numeric. Size of x-axis labels
#' @param colorType character. How ggplot color graph, 'fill', 'contour' or 'both'
#' @param returnTab logical. Return processed data ready to plot by GGplot ?
#' @param axis.names character. A two element vector containing axis names
#' @param colorScale list or vector. If condition has one column, a character specifying color, else a list of character. 
#' @param ... list. Additional parameter passed to geom function
#' @return Nothing
#' @examples
#' genes<-c("geneA","geneB")
#' samples<-paste0("sample",as.character(1:10))
#' #Generation of random expression
#' expr<-rbind(abs(c(rnorm(5,6,4),rnorm(5,8,2))),abs(c(rnorm(5,100,15),rnorm(5,50,20))))
#' rownames(expr)<-genes; colnames(expr)<-samples
#'
#' #Generation of sample table
#' sampleAnnot<-data.frame(group=as.factor(c(rep("g1",5),rep("g2",5))),row.names=samples)
#'
#' plotExpr(expr=expr,conditions=sampleAnnot,scale="log")
plotExpr<-function(expr,conditions=NULL,legendName=NULL,errorBar="se", ciRate=0.95,	geom="point", addLine=FALSE,	xaxis="gene", negValue=FALSE, 
									 scale="identity",breaks = waiver(),xRotationLab=0, hjust=0.5,main=NULL,xLabelSize=10, colorType="contour",
									 returnTab=FALSE,returnGraphList=FALSE,printGraph=TRUE,colorScale=NULL,axis.names=NULL,...){
	
	require(grid)
	require(ggplot2)
	
	if(!errorBar%in%c("se","ci","na")) stop("errorBar must be 'na', 'ci' or 'se'")
	if(!xaxis%in%c("gene","annot")) stop("xaxis must be 'gene' or 'annot'")
	if(!scale%in%c("identity","log")) stop("scale must be 'identity' or 'log'")
	if(!colorType%in%c("contour","fill","both")) stop("scale must be 'identity' or 'log'")
	
	colorContour<-FALSE
	colorFill<-FALSE
	nullConditions<-is.null(conditions)
	nullColorScale<-is.null(colorScale)
	nullLegendName<-is.null(legendName)
	colorScaleList<-NULL
	if(colorType=="contour" | colorType=="both") colorContour=TRUE
	if(colorType=="fill" | colorType=="both") colorFill=TRUE
	
	if(is.vector(expr))	expr<-t(data.frame(expr=expr,row.names = names(expr)))
	if(nullConditions){
		conditions<-as.factor(colnames(expr))
		errorBar<-"na"
	}
	if(!(is.data.frame(conditions)|is.matrix(conditions)))	conditions<-data.frame(cond=conditions,row.names = colnames(expr))
	numPlots = ncol(conditions)
	if(numPlots>1 & is.list(colorScale)) colorScaleList<-colorScale
	cols<-floor(sqrt(numPlots))
	layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
									 ncol = cols, nrow = ceiling(numPlots/cols))
	
	if(printGraph){
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
	}
	
	ciFactor<-1
	if(errorBar=="ci") ciFactor<-qnorm(ciRate+(1-ciRate)/2)
	
	graphList<-list()
	
	for(condIndex in 1:ncol(conditions)){
		conditionName<-colnames(conditions)[condIndex]
		if(!is.null(colorScaleList)) colorScale<-colorScaleList[[conditionName]]
		cond<-unlist(conditions[,condIndex])
		names(cond)<-rownames(conditions)
		if(!is.factor(cond)) stop("You must give a factor dataframe/matrix") #Qualitatif
		tabGraph<-data.frame(matrix(nrow = length(levels(cond))*nrow(expr),ncol=5))
		colnames(tabGraph)<-c(conditionName,"means","errBarMin","errBarMax","gene")
		tabGraph[,conditionName]<-rep(levels(cond),nrow(expr))
		tabGraph$gene<-rep(rownames(expr),each=length(levels(cond)))
		for(lvl in levels(cond)){
			for(exprIndex in 1:nrow(expr)){
				values<-unlist(expr[exprIndex,which(colnames(expr)%in%names(cond[which(cond==lvl)]))])
				nameExpr<-rownames(expr)[exprIndex]
				tabGraph[which(tabGraph[,conditionName]==lvl & tabGraph$gene == nameExpr),"means"]<-mean(values)
				tabGraph[which(tabGraph[,conditionName]==lvl & tabGraph$gene == nameExpr),"errBarMin"]<-	se(values)*ciFactor
				tabGraph[which(tabGraph[,conditionName]==lvl & tabGraph$gene == nameExpr),"errBarMax"]<-	se(values)*ciFactor
			}
		}
		tabGraph[,conditionName]<-factor(tabGraph[,conditionName],levels=levels(cond))
		tabGraph$gene<-factor(tabGraph$gene,levels=rownames(expr))
		
		if(!negValue)tabGraph$errBarMin[which(tabGraph$means-tabGraph$errBarMin<0)]<-tabGraph$means[which(tabGraph$means-tabGraph$errBarMin<0)]
		
		if(xaxis=="gene"){
			g<-"gene"
			x<-conditionName
		}else{
			g<-conditionName
			x<-"gene"
		}
		paramAES<-list(x=x,group=g)
		if(!nullConditions){
			if(colorContour) paramAES$colour= g
			if(colorFill) paramAES$fill= g
		}
		addParam<-list(...)
		if(geom=="bar"){
			if(is.null(addParam$stat)) addParam$stat = "identity";
			if(is.null(addParam$width)) addParam$width=.4;
			if(is.null(addParam$position)) addParam$position = "dodge"
		}
		paramAES$y="means"
		graph<-ggplot(data=tabGraph,mapping=do.call("aes_string",paramAES))+do.call(paste0("geom_",geom),addParam)
		if(addLine) graph<-graph+geom_line()
		
		if(errorBar != "na"){
			if(geom=="bar"){
				graph<-graph+geom_errorbar(width=.1, aes(ymin=means-errBarMin, ymax=means+errBarMax),position = position_dodge(width=.4))
			}else{
				graph<-graph+geom_errorbar(width=.1, aes(ymin=means-errBarMin, ymax=means+errBarMax))
			}
		}
		
		if(scale=="log"){
			graph<-graph+scale_y_log10(breaks=breaks)
		}else{
			graph<-graph+scale_y_continuous(breaks=breaks)
		}
		
		if(nullLegendName) legendName=g
		if(nullColorScale) colorScale<-ggplotColours(nlevels(tabGraph[,g]))
		if(colorFill) graph<-graph+scale_fill_manual(name=legendName,values=colorScale)
		if(colorContour) graph<-graph+scale_color_manual(name=legendName,values=colorScale)
		
		if(! is.null(main)) graph<-graph+ggtitle(main)
		if(is.null(axis.names)) axis.names<-c("","Expression")
		graph<-graph+xlab(axis.names[1]) + ylab(axis.names[2])
		graph<-graph+theme(axis.text.x = element_text(angle = xRotationLab, hjust = hjust,vjust=.3,size=xLabelSize),
											 panel.background = element_rect(fill = NA,colour="black"),
											 panel.grid.major = element_line(colour = "grey50"),
											 panel.grid.minor = element_line(colour = "grey50"))
		graphList[[conditionName]]<-graph
		if(printGraph){
			matchidx <- as.data.frame(which(layout == condIndex, arr.ind = TRUE))
			print(graph, vp = viewport(layout.pos.row = matchidx$row,layout.pos.col = matchidx$col))
		}
	}
	if(returnTab & returnGraphList) warnings("returnTab and returnGraphList is set to TRUE, only graph list will be returned")
	if(returnGraphList) return(graphList)
	if(returnTab) return(tabGraph)
}


# 3rd generation enrichment (fgsea algorithm)
#' @param x : vector or dataframe/matrix of one column. Values are used to classify the genes, example, it can be Log2(Fold-Change0). Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-rnorm(n = 4); names(x)<-c("GATA2","SOX17","KLF4","POU5F1")
#' @param corrIdGenes : dataframe of genes id (used to convert genes), automatically computed if not provided
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param minSize : mininmum number of gene in each term
#' @param maxSize : maximum number of gene in each term
#' @param nperm : number of permutation in the GSEA algorithm
#' @param customAnnot : custom annotation database, as a list af gene symbols, named by annotations name
#' @param returnLeadingEdge : return genes that were the most important for the enrichment of term
#' @param keggDisease : retain kegg disease term ?
#' @param species : species, example: "Rat", "Mouse", "Human"
#' @param db_terms : precomputed list of term database, automatically computed if not provided
#' @param ... : additionnal parameters that are passed to fgsea
#' @param speciesData : result of getSpeciesData function, automatically gathered if not provided
Enrich<-function(x, corrIdGenes=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
								 maxSize=500,minSize=2,nperm=1000,customAnnot=NULL,returnLeadingEdge=FALSE,
								 keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL,...){
	require(fgsea)
	require(AnnotationDbi)
	select<-AnnotationDbi::select
	
	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}
	
	if(class(x)!="numeric") stop("Values must be numerical")
	if(class(names(x))!="character") stop("Values must be named with genes symbol")
	if((!is.null(speciesData)) & is.null(corrIdGenes)) corrIdGenes<-speciesData$GeneIdTable
	if(is.null(corrIdGenes)) corrIdGenes<-getSpeciesData(species,names(x))$GeneIdTable
	geneSym<-names(x)
	geneEntrez<-ConvertKey(names(x),tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
	new_x<-x[!is.na(geneEntrez)]
	names(new_x)<-geneEntrez[!is.na(geneEntrez)]
	if(is.null(db_terms)) db_terms<-getDBterms(geneSym=geneSym,geneEntrez=geneEntrez, corrIdGenes=corrIdGenes,database=database,customAnnot=customAnnot,keggDisease=keggDisease,species=species)
	if(length(db_terms)==0) stop("Error, no term in any database was found")
	res<-list()
	for(db in names(db_terms)){
		res[[db]]<-fgsea(db_terms[[db]], new_x,nperm=nperm,minSize=minSize,maxSize=maxSize,...)
		res[[db]]<-res[[db]][order(res[[db]]$NES),]
		res[[db]]$database<-db
		res[[db]]$leadingEdge<- if(returnLeadingEdge) sapply(res[[db]]$leadingEdge,ConvertKey,tabKey=corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL") else NULL
	}
	res<-do.call("rbind", res)
	return(res)
}


# 2nd generation enrichment (fisher algorithm)
#' @param x : vector or dataframe/matrix of one column. Values are booleans and say if gene is from the list of interest or not. Genes are contained in the names/rownames of the vector/dataframe. Example of valid x: x<-c(TRUE,TRUE,FALSE,FALSE); names(x)<-c("GATA2","SOX17","KLF4","POU5F1"). In this case, GATA2, SOX17, KLF4, POU5F1 are the universe of gene and GATA2 and SOX17 are the genes of interest
#' @param corrIdGenes : dataframe of genes id (used to convert genes), automatically computed if not provided
#' @param database Which annotation database ? valid: database: kegg reactom goBP goCC goMF custom
#' @param minSize : mininmum number of gene in each term
#' @param maxSize : maximum number of gene in each term
#' @param customAnnot : custom annotation database, as a list af gene symbols, named by annotations name
#' @param returnGenes : return genes of interest that are in the term 
#' @param keggDisease : retain kegg disease term ?
#' @param db_terms : precomputed list of term database, automatically computed if not provided
#' @param species : species, example: "Rat", "Mouse", "Human"
#' @param speciesData : result of getSpeciesData2 function, automatically gathered if not provided
EnrichFisher<-function(x, corrIdGenes=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
											 minSize=2,maxSize=500,returnGenes=FALSE, keggDisease=FALSE,species="Human",
											 customAnnot=NULL,db_terms=NULL,speciesData=NULL){
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	if(sum(database%in%validDBs)==0) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
	if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")
	
	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}
	
	if(class(x)!="logical") stop("Values must be logical (TRUE or FALSE)")
	if(class(names(x))!="character") stop("Values must be named with genes symbol")
	
	if(is.null(db_terms)) db_terms<-getDBterms(geneSym=geneSym,geneEntrez=geneEntrez, corrIdGenes=corrIdGenes,database=database,customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
	
	nInterest<-length(which(x))
	nNotInterest<-length(which(!x))
	
	results<-list()
	for(db in names(db_terms)){
		len_term<-sapply(db_terms[[db]],length)
		db_terms[[db]]<-db_terms[[db]][len_term>=minSize & len_term<=maxSize]
		
		nGeneByterm<-sapply(db_terms[[db]],length)
		nGeneOfInterestByterm<-sapply( db_terms[[db]],function(term){
			return(length(which(x[term])))
		})
		results[[db]]<-data.frame(row.names = names(db_terms[[db]]))
		results[[db]]$pathway<- names(db_terms[[db]])
		results[[db]]$pval<-phyper(q = nGeneOfInterestByterm-0.5, m = nInterest,n = nNotInterest, k = nGeneByterm, lower.tail=FALSE)
		results[[db]]$padj<-p.adjust(results[[db]]$pval,method = "BH")
		results[[db]]$nGeneOfInterest<-nGeneOfInterestByterm
		results[[db]]$nGene<-nGeneByterm
		results[[db]]$database<-db
		if(returnGenes){
			results[[db]]$Genes<- sapply(db_terms[[db]],function(term){
				term[term%in%geneEntrez[x]]
			});
		}
	}
	do.call("rbind", results)
}

exportEnrich<-function(enrichResults,file,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE,geneCol="leadingEdge",...){
	enrichResults$Gene<-sapply(enrichResults[[geneCol]],function(x){ return(paste0(x,collapse=sep))})
	enrichResults[[geneCol]]<-NULL
	write.table(enrichResults,file,quote = quote,sep = sep,col.names = col.names,row.names = row.names,...)
}



getDBterms<-function(geneSym,geneEntrez=NULL, corrIdGenes=NULL, speciesData=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),customAnnot=NULL,
										 keggDisease=FALSE,species="Human",returnGenesSymbol=FALSE){
	require(AnnotationDbi)
	select<-AnnotationDbi::select
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	if(!(combineLogical(database%in%validDBs))) stop(paste0("Error, valid values for database are: ",paste0(validDBs,collapse=", ")))
	if(is.null(customAnnot) & "custom"%in%database) stop("You must give a value a list in customAnnot if database=custom")
	if(is.null(speciesData)){
		speciesData<-getSpeciesData2(species)
	}else{
		species<-speciesData$species
	}
	if(is.null(corrIdGenes)) corrIdGenes<-speciesData$GeneIdTable
	options(warn=-1)
	if(is.null(geneEntrez)){
		geneEntrez<-ConvertKey(geneSym,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
		geneEntrez<-geneEntrez[!is.na(geneEntrez)]
	}
	db_terms<-list()
	if(is.list(customAnnot)){
		db_terms$custom<-lapply(customAnnot,function(x){
			new_x<-ConvertKey(x,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID")
			new_x[!is.na(new_x)]
		})
	}
	if(!(length(database)<=1 & database[1]=="custom")){
		if("reactom"%in%database){ 
			require("fgsea")
			require("reactome.db")
			db_terms$reactom<- reactomePathways(geneEntrez)
			db_terms$reactom<-db_terms$reactom[unique(names(db_terms$reactom))]
		}
		if("kegg"%in%database){
			require("gage")
			kg.species <- kegg.gsets(speciesData$kegg, id.type="entrez")
			db_terms$kegg<- if(keggDisease) kg.species$kg.sets else kg.species$kg.sets[kg.species$sigmet.idx]
		}
		if("go"%in%substr(database,1,2)){
			require("gage")
			go.species <- go.gsets(tolower(species), id.type="entrez")
			if("goBP"%in%database) db_terms$goBP<-go.species$go.sets[go.species$go.subs$BP]
			if("goMF"%in%database) db_terms$goMF<-go.species$go.sets[go.species$go.subs$MF]
			if("goCC"%in%database) db_terms$goCC<-go.species$go.sets[go.species$go.subs$CC]
		}
	}
	
	for(db in names(db_terms)){
		db_terms[[db]]<-sapply(db_terms[[db]],function(term){
			term[term %in% geneEntrez]
		})
	}
	options(warn=0)
	
	if(returnGenesSymbol){
		lapply(db_terms,function(db) lapply(db,ConvertKey,tabKey=corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL"))
	}else{
		db_terms
	}
}

acp2d<-function(pca, comp=1:2,group=NULL, plotVars = FALSE, pointSize=2, plotText=FALSE,fixedCoord=FALSE,main=NULL,ellipse=FALSE,color=NULL,returnGraph=FALSE){
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(comp)!=2) stop("You must give a vector of 2 integer for comp parameter");
	percentVar <- pca$percentVar
	functPlot<-ifelse(plotText,geom_text,geom_point)
	coor=ifelse(plotVars,"rotation","x")
	
	if(is.null(group)){
		d <- data.frame(PC1=pca[[coor]][,comp[1]], PC2=pca[[coor]][,comp[2]]);
		graph<-ggplot(data=d, mapping = aes(x=PC1, y=PC2, label = rownames(d)))
	}else{
		d <- data.frame(PC1=pca[[coor]][,comp[1]], PC2=pca[[coor]][,comp[2]], group=group);
		graph<-ggplot(data=d, mapping = aes(x=PC1, y=PC2,colour=group, label = rownames(d)))
	}
	graph<-graph+functPlot(size=pointSize)+
		xlab(paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance")) +
		ylab(paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance")) 
	if(fixedCoord)graph <- graph + coord_fixed(ratio=percentVar[comp[2]]/percentVar[comp[1]])
	if(ellipse) graph<-graph+stat_ellipse()
	if(!is.null(main)) graph <- graph + ggtitle(main)
	if(!is.null(color)) graph <- graph + scale_color_manual(values=color)
	if(returnGraph){
		return(graph)
	}else{
		print(graph)
	}
}
