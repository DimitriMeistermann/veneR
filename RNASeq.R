#!/usr/bin/R

#need general.R

filterMostExpressedGenesBySample <-function(data, numberOfGenes=nrow(data),minFreqOfGene=1,maxFreqOfGene=ncol(data),threshold=min(data)){
	#Fonction de filtrage des counts, on prend les x gènes les plus exprimés pour chache échantillon puis on effectue une jointure complète
	#puis on filtre selon la fréquence du gène
	#data: dataframe des counts
	#numberOfGenes: cutoff sur les rank d'expression des gènes, plus il est élevé moins le filtre est stringeant
	#minFreqOfGene : nombre de fois minimum où le gène revient dans la liste des x gènes les plus exprimés (où le gène > threshold)
	#maxFreqOfGene : nombre de fois maximale où le gène est exprimé
	#threshold: l'expression du gène doit être plus grande que ce paramètre pour que la fréquence soit comptée comme 1
	mostExpressedGenes<-list()
	for(i in 1:ncol(data)){
		col<- data[,i]
		names(col)<-rownames(data)
		col<-col[which(col>threshold)];
		mostExpressedGenes[[i]]<-names(sort(col,decreasing = T))[1:min(numberOfGenes,length(col))]
	}
	rm(col)

	freqTable<-summary(as.factor(unlist(mostExpressedGenes)),maxsum=nrow(data))
	mostExpressedgenesVector<-names(freqTable[which(freqTable>=minFreqOfGene & freqTable<=maxFreqOfGene)])
	return(data[mostExpressedgenesVector,])
}

vectSup0<-function(x){
	a<-x[which(x>0)]
	return(length(a))
}

nbGeneSup0<-function(data){
	nS<-ncol(data)
	i<-apply(data,1,vectSup0)
	return(length(i[which(i>0)]))
}

filterMostExpressedGenesBySampleThres <-function(data, numberOfGenes=3000,maxGenes=nrow(data)/2){
	#Fonction de filtrage des counts, on prend les x gènes les plus exprimés pour chache échantillon puis on effectue une jointure complète
	#puis on filtre selon la fréquence du gène
	#data: dataframe des counts
	#numberOfGenes: cutoff des gènes, plus il est élevé moins le filtre est stringeant
	#maxGenes : nombre maximal de gène que l'on veut en sortie
	mostExpressedGenes<-data.frame(matrix(ncol = ncol(data), nrow=numberOfGenes)) #création d'un dataframe vide
	for(i in 1:ncol(data)){
	col<- data[,i]
	names(col)<-rownames(data)
	mostExpressedGenes[,i]<-names(sort(col,decreasing = T,method="shell"))[1:numberOfGenes]
	}
	rm(col)
	colnames(mostExpressedGenes)<-colnames(data)
	freqTable<-sort(summary(as.factor(unlist(mostExpressedGenes)),maxsum=nrow(data)),decreasing=TRUE,method="shell")
	maxGenes<-min(maxGenes,length(freqTable))
	mostExpressedgenesVector<-names(freqTable[1:maxGenes])
	return(data[mostExpressedgenesVector,])
}

interceptDistThres<-function(data,threshold=nrow(data)/3){
	print("construction of database");
	nrows<-nrow(data);
	ncols<-ncol(data);
	listOfgenesBySample=matrix("",threshold,ncols);
	for(i in 1:ncols){
		col<-data[,i];
		names(col)<-1:nrows;
		#tri alp ?
		listOfgenesBySample[,i]<-names(sort(col,decreasing = T,method="shell"))[1:threshold];
	}
	print("construction of dist matrix");
	distMat<-matrix(0,ncols,ncols)
	rownames(distMat)<-colnames(data);
	colnames(distMat)<-rownames(distMat);
	c<-1;
	for(i in 2:ncols){
		for(j in 1:c){
			distMat[i,j]<-length(setdiff(listOfgenesBySample[,i],listOfgenesBySample[,j]));
		}
		c<-c+1;
	}
	distMat<-as.dist(distMat);
	return(distMat);
}

interceptDistanceLevenstein<-function(data){
	print("construction of database");
	nrows<-nrow(data);
	ncols<-ncol(data);
	listOfgenesBySample=matrix("",nrows,ncols);
	for(i in 1:ncols){
		col<-data[,i];
		names(col)<-1:nrows;
		listOfgenesBySample[,i]<-names(sort(col,decreasing = T,method="shell"));
	}
	print("construction of dist matrix");
	distMat<-matrix(0,ncols,ncols)
	rownames(distMat)<-colnames(data);
	colnames(distMat)<-rownames(distMat);
	c<-1;
	for(i in 2:ncols){
		for(j in 1:c){
			distMat[i,j]<-LevenshteinDist(listOfgenesBySample[,i],listOfgenesBySample[,j]);
		}
		c<-c+1;
	}
	distMat<-as.dist(distMat);
	return(distMat);
}

compareVectReplace<-function(data,verbose=FALSE){
	#return a dist matrix
	#@param data : dataframe or matrx (rows = genes, cols = samples)
	# tri alph
	print("construction of database");
	nrows<-nrow(data);
	ncols<-ncol(data);
	listOfgenesBySample=matrix(0,nrows,ncols);
	for(i in 1:ncols){
		col<-data[,i];
		names(col)<-1:nrows;
		listOfgenesBySample[,i]<-as.integer(names(sort(col,decreasing = T,method="shell")));
	}
	print("construction of dist matrix");
	distMat<-matrix(0,ncols,ncols)
	rownames(distMat)<-colnames(data);
	colnames(distMat)<-rownames(distMat);
	tot<-length(as.dist(distMat))
	c<-1;
	out<-0;
	comptT<-0
	for(i in 2:ncols){
		for(j in 1:c){
			distMat[i,j]<-.Call("orderDist",as.integer(listOfgenesBySample[,i]),as.integer(listOfgenesBySample[,j]))
			comptT<-comptT+1;
		}
		if(verbose) print(paste0(comptT,"/",tot));
		c<-c+1;
	}
	distMat<-as.dist(distMat);
	return(distMat);
}

genesRankVSExpr<-function(data, numberOfGenes=nrow(data)){
	#Stats par gènes
	nsamples<-ncol(data)
	mostExpressedGenes<-data.frame(matrix(ncol = nsamples, nrow=numberOfGenes));
	for(i in 1:nsamples){
		col<- data[,i];
		names(col)<-rownames(data);
		mostExpressedGenes[,i]<-names(sort(col,decreasing = T))[1:numberOfGenes];
	}
	colnames(mostExpressedGenes)<-colnames(data);
	ranks<-data.frame(matrix(ncol = nsamples, nrow=numberOfGenes));
	rownames(ranks)<-mostExpressedGenes[,1]
	colnames(ranks)<-colnames(data);
	for(i in 1:nsamples){
		col<-1:numberOfGenes;
		names(col)<-mostExpressedGenes[,i];
		ranks[,i]<-col[rownames(ranks)]
	}
	meanRank<-apply(ranks,1,mean)
	sdRank<-apply(ranks,1,sd)
	minRank<-apply(ranks,1,min)
	maxRank<-apply(ranks,1,max)
	data<-data[names(meanRank),]
	sdExpr<-apply(data,1,mean)
	minExpr<-apply(data,1,sd)
	maxExpr<-apply(data,1,max)
	meanExpr<-apply(data,1,min)
	res<-data.frame(meanRank = meanRank,SDRank=sdRank, minRank=minRank, maxRank=maxRank, 
		meanExpr = meanExpr,sdExpr=sdExpr,maxExpr=maxExpr,minExpr=minExpr) ;
	res<-res[order(res$meanRank),]
	return(res);
}

examineRNAseqSamples<-function(x, uncenter= FALSE){
	#Donne différentes stat par échantillon 
	if(uncenter){
		x<-x-min(x)
		zero<-0
	}else{
		zero<-min(x)
	}
	mean<-apply(x,2,mean)
	sd<-apply(x,2,sd)
	count<-colSums(x)
	CV<-apply(x,2,cv)
	noGenEx<-rep(0,ncol(x))
	for(i in 1:ncol(x)) noGenEx[i]<-length(which(x[,i]>zero))

	return(data.frame(mean=mean, sd=sd,CV=CV,TotalGenEx=noGenEx,TotalCount=count))
}

retrieveSexHumanEmbryoKmeans<-function(d,group){
	group<-droplevels(as.factor(group))
	#return a matrix of count of "male" and "female" predicted cells in each embryo (Kmeans method)
	maleGene<-c("DDX3Y","EIF1AY","TTTY15","RPS4Y1")
	k<-kmeans(t(d[maleGene,]),2)
	mORf<-rowSums(k$centers)
	if(mORf[1]<mORf[2]){
		mf<-c("F","M")
	}else{
		mf<-c("M","F")
	}
	count<-as.factor(mf[k$cluster])
	embryos<-as.factor(levels(group))
	names(count)<-group
	res<-list()
	res$count<-data.frame(matrix(ncol=2,nrow=length(embryos)))
	colnames(res$count)<-c("Male","Female")
	rownames(res$count)<-embryos
	
	for(embryo in embryos){
		res$count[embryo,1]<-length(which(count[which(names(count)==embryo)]=="M"))
		res$count[embryo,2]<-length(which(count[which(names(count)==embryo)]=="F"))
	}
	res$freq<-res$count/rowSums(res$count)
	return(res)
}

retrieveSexHumanEmbryoACP<-function(d,patternLen=3){
	#return a freq matrix of "male" and "female" predicted cells in each embryo (ACP method)
	maleGene<-c("DDX3Y","EIF1AY","TTTY15","RPS4Y1")
	z<-colSums(d[maleGene,])
	#voir ici pour renvoyer tout mâle ou tout femelle
	#normer acp sans réduire ? et prendre seuil
	
	a<-ACP(t(d[maleGene,]))
	M<-as.factor(substr(names(a$x[which(a$x[,1]>0),1]),1,patternLen)) #sélection du nom de l'embryon seul, pas des cellules (substr)
	F<-as.factor(substr(names(a$x[which(a$x[,1]<0),1]),1,patternLen))
	embryos<-as.factor(unique(substr(colnames(d),1,patternLen)))
	res<-list()
	res$count<-data.frame(matrix(ncol=2,nrow=length(embryos)))
	colnames(res$count)<-c("Male","Female")
	rownames(res$count)<-embryos
	
	for(embryo in embryos){
		res$count[embryo,1]<-length(which(M==embryo))
		res$count[embryo,2]<-length(which(F==embryo))
	}
	res$freq<-res$count/rowSums(res$count)
	return(res)
}


plotExpression<-function(expr,group=NULL,log10Plus1yScale=NULL,violin=TRUE,boxplot=TRUE,dotplot=FALSE,
				 violinArgs=list(),boxplotArgs=list(),dotplotArgs=list(),colorScale=mostDistantColor2,
				 defaultGroupName="group",dodge.width=.9,printGraph=FALSE){
	require(ggplot2)
	require(reshape2)
	require(ggbeeswarm)
	barplotGraph=greyGraph=coloredGraph=FALSE
	
	if(is.vector(expr))expr<-t(as.matrix(expr))
	if(is.vector(group) | is.factor(group)){
		group = data.frame(group=group,stringsAsFactors = TRUE)
		colnames(group)<-defaultGroupName
		rownames(group)<-colnames(expr)
	}
		
	if(!is.matrix(expr)) expr<-as.matrix(expr)
	if(is.null(log10Plus1yScale)) log10Plus1yScale<-nrow(expr)>1 #if more than one gene log10Plus1yScale is turned on

	if(is.null(group)){
		ggData<-reshape2::melt(expr,value.name="expression",varnames=c("gene","sample"))
		ggData$sample<-factor(ggData$sample,levels = levels(ggData$sample)[order(t(expr),decreasing = TRUE)])
		if(nrow(expr)==1){
			g<-ggplot(ggData,mapping = aes(x=sample,y=expression))+geom_bar(stat = "identity")+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3))
			barplotGraph<-TRUE
		}else{
			ggData$gene<-factor(ggData$gene,levels = rownames(expr))
			g<-ggplot(ggData,mapping = aes(x=gene,y=expression))+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold.italic"))
			greyGraph<-TRUE
		}
	}else{ #group on
		if(ncol(expr)!=nrow(group)) stop("expr should have the same number of sample than group")
		if(ncol(group)>1) stop("Multiple group are not allowed in the same time")
		groupName <- colnames(group)
		if(nrow(expr)==1){
			ggData<-data.frame(expression=expr[1,],group)
			g<-ggplot(ggData,mapping = aes_string(x=groupName,y="expression"))+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold"))
			greyGraph<-TRUE
		}else{
			coloredGraph=TRUE
			ggData<-reshape2::melt(data.frame(t(expr),group),value.name="expression",variable.name="gene",id.vars = groupName)
			if(violin){
				factorSampling<-table(group[,1])
				factor2drop<-names(factorSampling)[factorSampling<3]
				if(length(factor2drop)>1) warning(paste0(factor2drop,collapse = " "),
					" were dropped (n<3 is not compatible with violin plot). You can deactivate violin layer by setting violin argument to FALSE")
				ggData<-ggData[!ggData[,groupName]%in%factor2drop,] #drop levels where n < 3
				
			}
			colors <- if(is.function(colorScale)) colorScale(nlevels(group[,1])) else colorScale
			g<-ggplot(ggData,mapping = aes_string(x="gene",y="expression",fill=groupName))+
				theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3,face = "bold.italic"))+
				scale_fill_manual(values = colors[!levels(group[,1])%in%factor2drop])
		}
	}
	if(greyGraph){
			if(is.null(violinArgs$fill)) violinArgs$fill<-"grey50"
			if(is.null(violinArgs$scale)) violinArgs$scale<-"width"
			if(is.null(boxplotArgs$width)) boxplotArgs$width<-.2
	}
	if(coloredGraph){
			if(is.null(violinArgs$scale)) violinArgs$scale<-"width"
			if(is.null(boxplotArgs$width)) boxplotArgs$width<-.2
			if(is.null(violinArgs$position)) violinArgs$position<-position_dodge(preserve = "total",width = dodge.width)
			if(is.null(boxplotArgs$position)) boxplotArgs$position<-position_dodge(preserve = "total",width = dodge.width)
			if(is.null(dotplotArgs$dodge.width)) dotplotArgs$dodge.width<-dodge.width
	}
	if(!barplotGraph){
		g<-ggBorderedFactors(g,borderColor="black",borderSize = .5)
		if(violin) g<-g+do.call("geom_violin",violinArgs)
		if(boxplot) g<-g+do.call("geom_boxplot",boxplotArgs)
		if(dotplot)	g<-g+do.call("geom_beeswarm",dotplotArgs)
	}
	if(log10Plus1yScale) g<-g+scale_y_continuous(trans=log10plus1())
	if(printGraph){
		print(g)
	}else{
		return(g)
	}
}

#Utile pour les packages qui demandes des fonctions comme argument
colMap <- function(x) {
	.col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
	return(.col[match(1:length(x), order(x))])
}
selectGene<- function(x){
	return(x==1)
}

###

CPM<-function(data){ #Normalisation CPM
	data.CPM <- sweep(data, 2, colSums(data),`/`)
	data.CPM <-data.CPM * 1000000
	return(data.CPM)
}

UMI2UPM<-CPM

TPMfullLength<-function(data, gene.length){
	gene.length.kb <- gene.length[rn(data)]/1000
	data<-sweep(data, 1, gene.length.kb,`/`)
	return(CPM(data))
}

RPKM<-function(data, gene.length){
	gene.length.kb <- gene.length[rn(data)]/1000
	data<-CPM(data)
	sweep(data, 1, gene.length.kb,`/`)
}

plotDistrib<-function(data,type="boxplot",conditions=NULL,main=NULL,conditionName="Batch"){
	require(grid)
	require(ggplot2)
	#data : matrix of expression data
	#type: 'violin' or 'boxplot' ?
	#conditions : vector of factors to color plot
	if(!type%in%c("violin","boxplot")) stop("type must be 'violin', 'hist' or 'boxplot'")
	
	vectorDat<-as.vector(as.matrix(data))
	tabGraph<-data.frame(val=vectorDat,sample=rep(cn(data),each=nrow(data)))
	if(!is.null(conditions)){
		tabGraph[,conditionName]=rep(conditions,each=nrow(data))
	}
	
	tabGraph$sample<-factor(tabGraph$sample,levels=cn(data))
	
	if(is.null(conditions)){
		graph<-ggplot(data = tabGraph,mapping = aes(sample,val))
	}else{
		graph<-ggplot(data = tabGraph,mapping = aes_string("sample","val",color=conditionName))
	}
	
		graph<-graph+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3))
	if(type=="violin") graph<-graph+geom_violin()
	if(type=="boxplot") graph<-graph+geom_boxplot()
	if(! is.null(main)) graph<-graph+ggtitle(main)
	
	print(graph)
}



#Return gene id correspondance, GO species code and KEGG species code
getSpeciesData2<-function(sample.species="Human",updateSpeciesPackage=FALSE){
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
	suppressMessages(species$GeneIdTable<-AnnotationDbi::select(get(species$package),keys = AnnotationDbi::keys(get(species$package),"SYMBOL") , "ENTREZID","SYMBOL"))
	species$species<-sample.species
	return(species)
}


#Take a gene list, add transcription factor data and long Gene name
detailOnGenes<-function(x,tfDat,speciesDat){
	require(AnnotationDbi)
	if(is.data.frame(x) | is.matrix(x)){
		geneSym<-rn(x)
		res<-data.frame(x)
	}else{
		if(is.null(names(x))){
			res<-data.frame(row.names=x)
			geneSym<-x
		}else{
			geneSym<-names(x)
			res<-data.frame(val=x,row.names=names(x))
		}
	}
	geneNames<-select(get(speciesDat$package),geneSym,"GENENAME","SYMBOL")
	res$LongName<-ConvertKey(geneSym,tabKey = geneNames,colOldKey = "SYMBOL",colNewKey = "GENENAME")
	genesTF<-intersect(rn(tfDat),geneSym)
	res$TFdegree<-0
	res[genesTF,"TFdegree"]<-tfDat[genesTF,"tf_degree"]
	return(res)
}


#Optimal cutree on hclust
best.cutree <- function(hc, min=3, max=20, loss=FALSE, graph=FALSE, ...){
	require(ggplot2)
	if (class(hc)!="hclust") hc <- as.hclust(hc)
	max <- min(max, length(hc$height))
	inert.gain <- rev(hc$height)
	intra <- rev(cumsum(rev(inert.gain)))
	relative.loss = intra[min:(max)]/intra[(min - 1):(max - 1)]
	best = which.min(relative.loss)
	names(relative.loss) <- min:max
	if (graph) {
	print(
		ggplot(data.frame(partition=min:max,relative.loss=relative.loss),aes(x=partition,y=relative.loss))+
		geom_point()+
		scale_x_continuous(breaks=min:max)
	)
	} else {
	if (loss)
		relative.loss
	else
		best + min - 1
	}
}

best.cutree2 <- function(hc, min=2, max=20, loss=FALSE, graph=FALSE, ...){
	require(ggplot2)
	if (class(hc)!="hclust") hc <- as.hclust(hc)
	max <- min(max, length(hc$height)-1)
	inert.gain <- rev(hc$height)
	intra <- rev(cumsum(rev(inert.gain)))
	relative.loss = intra[min:(max+1)]/intra[(min - 1):(max)]
	derivative.loss = relative.loss[2:length(relative.loss)]-relative.loss[1:(length(relative.loss)-1)]
	names(derivative.loss) <- min:max
	if (graph) {
	print(
		ggplot(data.frame(partition=min:max,derivative.loss=derivative.loss),aes(x=partition,y=derivative.loss))+
		geom_point()+
		scale_x_continuous(breaks=min:max)
	)
	} else {
	if (loss)
		derivative.loss
	else
		as.numeric(names(which.max(derivative.loss)))
	}
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
#' @param speciesData : result of getSpeciesData2 function, automatically gathered if not provided
Enrich.gsea<-function(x, corrIdGenes=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
		maxSize=500,minSize=2,customAnnot=NULL,returnGenes=FALSE,
		keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL,...){
	require(fgsea)
	
	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}
	
	if(class(x)!="numeric") stop("Values must be numerical")
	if(is.null(db_terms)) db_terms<-getDBterms(geneSym=geneSym,geneEntrez=geneEntrez, corrIdGenes=corrIdGenes,database=database,
																						 customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
	if(length(db_terms)==0) stop("Error, no term in any database was found")
	res<-list()
	for(db in names(db_terms)){
		res[[db]]<-fgseaMultilevel (db_terms[[db]], x ,minSize=minSize,maxSize=maxSize,eps = 0,...)
		res[[db]]<-res[[db]][order(res[[db]]$padj),]
		res[[db]]$database<-db
		res[[db]]$leadingEdge<-NULL
		if(returnGenes) res[[db]]$genes <- db_terms[[db]][res[[db]]$pathway]
	}
	res<-do.call("rbind", res)
	res$padj<-p.adjust(res$pval,method = "BH")
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
Enrich.simple<-function(x, corrIdGenes=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),
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
		results[[db]]$term <- names(db_terms[[db]])
		results[[db]]$pval<-phyper(q = nGeneOfInterestByterm-0.5, m = nInterest,n = nNotInterest, k = nGeneByterm, lower.tail=FALSE)
		results[[db]]$nGeneOfInterest<-nGeneOfInterestByterm
		results[[db]]$nGene<-nGeneByterm
		results[[db]]$database<-db
		if(returnGenes){
			results[[db]]$genes<- db_terms[[db]]
		}
	}
	results<-do.call("rbind", results)
	results$padj<-p.adjust(results$pval,method = "BH")
	return(results)
}


computeActivationScore<-function(expressionMatrix,corrIdGenes=NULL,
			database=c("kegg","reactom","goBP","goCC","goMF"),
			maxSize=500,minSize=2,nperm=1000,customAnnot=NULL,
			keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL){
	
	if(!class(expressionMatrix)[1]%in%c("data.frame","matrix")) stop("expressionMatrix should be a matrix or a dataframe")
	if(class(rownames(expressionMatrix))!="character") stop("rows of expression matrix should be named with genes symbol")
	if(is.null(db_terms)){
		db_terms<-getDBterms(geneSym=rownames(expressionMatrix), corrIdGenes=corrIdGenes,database=database,
												 customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
	}
	
	if(length(db_terms)==0) stop("Error, no term in any database was found")
	
	lapply(db_terms,function(database){
		database<-lapply(database,function(genesOfTerm) intersect(genesOfTerm,rownames(expressionMatrix)))
		nGenePerTerm<-sapply(database,length)
		database<-database[nGenePerTerm>minSize & nGenePerTerm<maxSize]
		resPerPathway<-lapply(database,function(genesOfTerm){
			eigengenes(expressionMatrix,genesOfTerm,returnContribution = TRUE)
		})
		list(
			eigen=t(sapply(resPerPathway,function(term) term$eigen)),
			contribution = lapply(resPerPathway,function(term) term$contribution)	 
		)
	})
}


GSDA<-function(geneSetEigens=NULL,expressionMatrix=NULL,colData,contrast, corrIdGenes=NULL,
		database=c("kegg","reactom","goBP","goCC","goMF"),
		maxSize=500,minSize=2,customAnnot=NULL,keggDisease=FALSE,species="Human",db_terms=NULL,speciesData=NULL){
	
	if(is.null(geneSetEigens) & is.null(geneSetEigens)) stop("At least expressionMatrix or geneSetEigens miiust be given")
	
	if(is.null(db_terms)){
		db_terms<-getDBterms(geneSym=rownames(expressionMatrix), corrIdGenes=corrIdGenes,database=database,
												 customAnnot=customAnnot,keggDisease=keggDisease,species=species,returnGenesSymbol = TRUE)
	}
	
	if(is.null(geneSetEigens)){
		geneSetEigens<-computeActivationScore(expressionMatrix=expressionMatrix,db_terms=db_terms)
	}
	
	res<-list()

	for(db in names(db_terms)){
		if(is.list(geneSetEigens[[db]])){
			eigenPerPathway<-geneSetEigens[[db]]$eigen
		}else{
			eigenPerPathway<-geneSetEigens[[db]]
		}
		
		db_terms[[db]]<-db_terms[[db]][rownames(eigenPerPathway)]
		res[[db]]<-dfres<-data.frame(term=names(db_terms[[db]]),testLinearModel(eigenPerPathway,colData,contrast),database=db,
			size=sapply(db_terms[[db]],length),sd=apply(eigenPerPathway,1,sd),row.names = NULL)
	}
	do.call("rbind", res)

}

GSDA.HeatmapAnnot<-function(contributions,maxGeneContribAtOneSide=3,width=unit(3,"cm"),fontsizeFactor=400){
	AnnotationFunction(fun = function(index, k, n) {
		pushViewport(viewport(xscale = c(0,10), yscale = c(0.5, length(index) + 0.5)))
		grid.lines(.5,c(0,1))
		i<-length(index)
		for(sel in index){
			contrib<-sort(contribPerPathway[[sel]])
			
			left<-names(contrib)[contrib<0]
			right<-names(contrib)[contrib>0]
			left<-left[1:min(maxGeneContribAtOneSide,length(left))]
			right<-right[max(1,(length(right)-maxGeneContribAtOneSide)+1):length(right)]
			if(!is.na(left[1])){
				grid.text(paste0(left,collapse = "  "),just = "right",
									x=4.5,y=i,default.units = "native",gp=gpar(fontsize=1/length(index)*fontsizeFactor,fontface="italic"))
			}
			if(!is.na(right[1])){
				grid.text(paste0(right,collapse = "  "),just = "left",
									x=5.5,y=i,default.units = "native",gp=gpar(fontsize=1/length(index)*fontsizeFactor,fontface="italic"))
			}
			i<-i-1
		}
		popViewport()
	},
	var_import = list(contribPerPathway = contributions, maxGeneContribAtOneSide=maxGeneContribAtOneSide,width=width,
										fontsizeFactor=fontsizeFactor),
	subsetable = FALSE,
	width = width,which = "row"
	)
}

getDBterms2<-function(geneSym,geneEntrez=NULL, corrIdGenes=NULL, speciesData=NULL,database=c("kegg","reactom","goBP","goCC","goMF"),customAnnot=NULL,
	keggDisease=FALSE,species="Human",returnGenesSymbol=TRUE){
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
	options(warn=0)
	
	if(returnGenesSymbol){
		lapply(db_terms,function(db) lapply(db,ConvertKey,tabKey=corrIdGenes,colOldKey = "ENTREZID",colNewKey = "SYMBOL"))
	}else{
		db_terms
	}
}

exportEnrich.2<-function(enrichResults,file,quote = FALSE,sep = "\t",col.names = TRUE,row.names = FALSE,geneCol="genes",...){
	enrichResults$Gene<-sapply(enrichResults[[geneCol]],function(x){ return(paste0(x,collapse=sep))})
	enrichResults[[geneCol]]<-NULL
	write.table(enrichResults,file,quote = quote,sep = sep,col.names = col.names,row.names = row.names,...)
}

calConsensusRankingOld<-function(genes,pvalues,logFoldChanges){
	pvalues <- 1-pvalues; abslogFoldChanges<-abs(logFoldChanges)
	dataRank<-data.frame(pval=pvalues,LFC=abslogFoldChanges,row.names=genes)
	rankRes<-apply(dataRank,1,gmean)
	rankRes[logFoldChanges<0]<- -rankRes[logFoldChanges<0]
	rankRes
}

calConsensusRanking<-function(genes,pvalues,logFoldChanges){
	InvPvalues <- 1-pvalues;
	InvPvalues[logFoldChanges<0]<- -InvPvalues[logFoldChanges<0]
	names(InvPvalues)<-genes
	InvPvalues
}

###ViewKEGG####
#x : 
viewKEGG<-function(x,pathway,corrIdGenes=NULL,species="Human",speciesData=NULL,directory=getwd(),...){	
	require(pathview)
	
	blacklist<-c("hsa04215 Apoptosis - multiple species")
	if(pathway%in%blacklist){
		warning(pathway," is blacklisted as it contains issues in vizualisation, it will not be rendered.")
		return(NULL)
	}
	
	if(is.data.frame(x) | is.matrix(x)){
		tempx<-x
		x<-tempx[,1]
		names(x)<-rownames(tempx)
	}
	if(is.null(speciesData)) speciesData<-getSpeciesData2(species)
	if(is.null(corrIdGenes)) corrIdGenes<-speciesData$GeneIdTable
	entrezId<-ConvertKey(keyList = names(x),tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID" );
	notNA<-which(!is.na(entrezId))
	if(length(notNA)>0){
		dat<-x[notNA];
		names(dat)<-entrezId[notNA]
		dat<-dat[takefirst(names(dat),returnIndex = T)]
		
		
		pathview(gene.data = dat, pathway.id = pathway, species = speciesData$kegg,kegg.native=TRUE,
			low="#4B9AD5",mid="white",high="#FAB517",na.col="grey75",kegg.dir=directory,...)
	}else{
		warning("no entrez id were found")
	}
}

####Convert 2 Entrez
Sym2Entrez<-function(x, corrIdGenes=NULL,species="Human"){
	require(AnnotationDbi)
	select<-AnnotationDbi::select
	validDBs<-c("kegg","reactom","goBP","goCC","goMF","custom")
	options(warn=-1)
	if(is.null(corrIdGenes)) corrIdGenes<-getSpeciesData2(species)$GeneIdTable
	return(ConvertKey(x,tabKey = corrIdGenes,colOldKey = "SYMBOL",colNewKey = "ENTREZID"))
}
	
distPathway<-function(db_pathway){
	genesInPathway<-sort(unique(unlist(db_pathway,recursive = TRUE)))
	pathDF<-data.frame(matrix(FALSE,ncol = length(db_pathway),nrow = length(genesInPathway),dimnames = list(genesInPathway,names(db_pathway))))
	for(i in 1:length(db_pathway)){
		pathDF[db_pathway[[i]],i]<-TRUE
	}

	combPath<-combn(1:ncol(pathDF),2)
	distMat<-matrix(data=0,nrow=ncol(pathDF),ncol=ncol(pathDF),dimnames=list(cn(pathDF),cn(pathDF)))
	for(comp in 1:ncol(combPath)){
		i<-combPath[1,comp]
		j<-combPath[2,comp]
		lenInter<-length(which(pathDF[,i] & pathDF[,j]))
		leni<-length(which(pathDF[,i]))
		lenj<-length(which(pathDF[,j]))
		distMat[i,j]<-1-(lenInter/min(leni,lenj))
		distMat[j,i]<-distMat[i,j]
	}
	as.dist(distMat,upper=T)
}	
	
###Fonction pour retrouver les branch "Heatmap branch" de monocle
retrieveBranch<-function(cds,branch_point){
	require(monocle)
	require(igraph)
	pr_graph_cell_proj_mst <- minSpanningTree(cds)
	
	root_cell <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
	root_state <- pData(cds)[root_cell, ]$State
	pr_graph_root <- subset(pData(cds), State == root_state)
	
	closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
	root_cell_point_in_Y <- closest_vertex[row.names(pr_graph_root), ]
	
	root_cell <- names(which(igraph::degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, 
																	mode = "all") == 1, useNames = T))[1]
	paths_to_root <- list()
	
	pr_graph_cell_proj_mst <- minSpanningTree(cds)
	
	mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
	branch_cell <- mst_branch_nodes[branch_point] #Récupérer nom échantillon branch point
	mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
	path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, 
																		 branch_cell, root_cell)
	path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
	for (backbone_nei in V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name) {
		descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei], 
											 unreachable = FALSE)
		descendents <- descendents$order[!is.na(descendents$order)]
		descendents <- V(mst_no_branch_point)[descendents]$name
		if (root_cell %in% descendents == FALSE) {
			path_to_root <- unique(c(path_to_ancestor, branch_cell, 
															 descendents))
			closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
			path_to_root <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_root)]
			closest_vertex <- cds@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
			path_to_root <- intersect(path_to_root, colnames(cds))
			paths_to_root[[backbone_nei]] <- path_to_root
		}
	}
	return(paths_to_root)
}



# x : matrix
hierarchicalClustering<-function(x,transpose=TRUE,method.dist="euclidean",method.hclust="ward.D2",
																 bootstrap=FALSE,nboot=10,PCAfirst=FALSE,nDimPCA=NULL){
	if(transpose) x<-t(x)
	if(PCAfirst){
		x<-ACP(x,transpose = FALSE,scale = FALSE)$x
		if(!is.null(nDimPCA)){
			x<-x[,1:nDimPCA]
		}
	}
	if(bootstrap){
		require(pvclust)
		resClust<-pvclust(t(x),nboot=nboot,method.hclust = method.hclust,parallel = TRUE,method.dist = method.dist)$hclust
	}else{
		if(method.dist=="pearson"){
			resDist<-corrDist(x)
		}else if(method.dist=="bicor"){
			require("WGCNA")
			resDist<-as.dist((1 - suppressWarnings(bicor(Matrix::t(x))))/2)
		}else{
			resDist<-dist(x, method = method.dist)
		}
		resClust<-stats::hclust(resDist,method = method.hclust)
	}
	return(resClust)
}

#for compatibility
unsupervisedClustering<-hierarchicalClustering

genTopAnnot<-function(annot,colorScales=NULL,border = TRUE,...){
	require(ComplexHeatmap)
	require(circlize)
	
	if(is.factor(annot)|is.data.frame(annot)) annot<-droplevels(annot)
	if(is.vector(annot) | is.factor(annot)){
		annot=data.frame(Annotation=annot)
		if(is.list(colorScales)) colnames(annot)<-names(colorScales)[1]
		if( (!is.null(colorScales)) & !is.list(colorScales)) colorScales<-list("Annotation"=colorScales)
	}
	colorScales<-genColorsForAnnots(annots = annot,colorScales = colorScales,returnContinuousFun = TRUE,...)
	HeatmapAnnotation(df=annot,col = colorScales,border = border)
}


heatmap.DM3<-function(matrix,preSet="expression",clustering_distance_rows=NULL,clustering_distance_columns="euclidean",clustering_method_columns="ward.D2",
		clustering_method_rows="ward.D2",autoFontSizeRow=TRUE,autoFontSizeColumn=TRUE,scale=FALSE,center=TRUE,returnHeatmap=FALSE,name=NULL,
		additionnalRowNamesGpar=NULL,additionnalColNamesGpar=list(),border=TRUE,
		colorScale=NULL,colorScaleFun=NULL,midColorIs0=NULL,probs=NULL,useProb=TRUE,minProb=0.05,maxProb=0.95,
		cluster_rows=NULL,cluster_columns=NULL,sampleAnnot=NULL,colorAnnot=NULL,showGrid=NULL,gparGrid=gpar(col="black"),showValues=FALSE,Nsignif=3,
		column_dend_reorder = FALSE, row_dend_reorder=FALSE, ...){
	require("ComplexHeatmap")
	require("circlize")
	
	args<-list()
	allowedPreSet<-c("default","expression","correlation","distance")
	if(!preSet%in%allowedPreSet){
		stop(paste0("preSet must equal to one of this value: ",paste0(allowedPreSet,collapse = ", ")))
	}
	
	if(preSet=="expression"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows = corrDistBicor
		if(is.null(name)) name="centered log expression"
		if(is.null(midColorIs0)) midColorIs0=TRUE
		if(is.null(colorScale)) colorScale=	c("darkblue","white","red2")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list(fontface="italic")
	}else if(preSet=="correlation"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows ="euclidean"
		if(is.null(clustering_distance_columns)) clustering_distance_columns ="euclidean"
		if(is.null(name)) name="Pearson\ncorrelation"
		if(is.null(midColorIs0)) midColorIs0=TRUE
		if(is.null(colorScale)) colorScale=c("darkblue","white","#FFAA00")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list()
	}else if(preSet=="distance"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows ="euclidean"
		if(is.null(name)) name="Euclidean\ndistance"
		if(is.null(midColorIs0)) midColorIs0=FALSE
		if(is.null(colorScale)) colorScale=c("white","yellow","red","purple")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list()
	}else	if(preSet=="default"){
		if(is.null(clustering_distance_rows)) clustering_distance_rows = corrDistBicor
		if(is.null(name)) name="matrix"
		if(is.null(midColorIs0)) midColorIs0=TRUE
		if(is.null(colorScale)) colorScale=c("#2E3672","#4B9AD5","white","#FAB517","#E5261D")
		if(is.null(additionnalRowNamesGpar)) additionnalRowNamesGpar=list()
	}
	
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
	if(is.null(colorScaleFun)){
		colorScaleFun<-computeColorScaleFun(colors = colorScale,values = unlist(matrix),useProb = useProb,probs = probs,minProb = minProb,
																				maxProb = maxProb, midColorIs0 = midColorIs0,returnColorFun = TRUE)
	}
	args$col<-colorScaleFun
	if(is.null(showGrid)){
		if(nrow(matrix)*ncol(matrix)<500){
			showGrid = TRUE
		}else{
			showGrid = FALSE
		}
	}
	if(showGrid){
		args$rect_gp = gparGrid
	}
	if(showValues){
		args$cell_fun = function(j, i, x, y, w, h, col) {
			#dark or light background .
			if(colSums(col2rgb(col))<382.5) col="white" else col="black" 
			grid.text(as.character(signif(matrix[i,j],Nsignif)),x,y,gp=gpar(col=col))
		}
	}
	if(autoFontSizeRow) args$row_names_gp=do.call("autoGparFontSizeMatrix",c(list(nrow(matrix)),additionnalRowNamesGpar))
	if(autoFontSizeColumn) args$column_names_gp=do.call("autoGparFontSizeMatrix",c(list(ncol(matrix)),additionnalColNamesGpar))
	
	if(!is.null(sampleAnnot)){
		args$top_annotation<-genTopAnnot(sampleAnnot,colorAnnot)
	}
	
	args$column_dend_reorder<-column_dend_reorder;args$row_dend_reorder<-row_dend_reorder
	args$matrix<-matrix
	args$name=name
	args$border<-border
	args<-c(args,list(...))
	
	ht<-do.call("Heatmap",args)
	if(returnHeatmap){
		return(ht)
	}else{
		print(ht)
	}
}


normDeseq<-function(expr){
	# PE = pseudo reference sample
	PE<-apply(expr,1,gmean,keepZero=T)
	PE<-PE[PE>0]
	genes<-names(PE)
	ratioMat<-sweep(expr[genes,],1,PE,"/")
	normFactors<-apply(ratioMat,2,median)
	sweep(expr,2,normFactors,"/")
}

testLinearModel<-function(exprData,sampleData,contrast){
	samples<-rn(sampleData)[sampleData[,contrast[1]]%in%contrast[2:3]]
	data<-exprData[,samples]
	groups<-droplevels(sampleData[samples,contrast[1]])
	logicGroup<-rep(F,len(groups))
	logicGroup[groups==contrast[2]]<-T
	regTabList<-apply(data,1,function(x){
		data.frame(data=x,group=logicGroup)
	})
	resList<-lapply(regTabList,function(regTab){
		summary(lm(data ~ group,data=regTab))$coefficients[2,c(1,4)]
	})
	res<-data.frame(do.call("rbind",resList));colnames(res)<-c("log2FoldChange","pval")
	res<-cbind(data.frame(baseMean=apply(exprData[,samples],1,mean)),res)
	res$padj<-p.adjust(res$pval,method = "BH")
	return(res)
}

getMostVariableGenes<-function(counts,normalize=TRUE,lib.size=NULL,plot=FALSE){
	#from http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
	require(DESeq2)
	require(statmod)
	counts<-as.matrix(counts)
	counts<-counts[rowMeans(counts)>0,]
	
	if(normalize){
		if(is.null(lib.size)) lib.size <- estimateSizeFactorsForMatrix(counts)
		ed <- t(t(counts)/lib.size)
	}else{
		ed<-counts
	}

	means <- rowMeans(ed)
	vars <- apply(ed,1,var)
	cv2 <- vars/means^2
	par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9)
	
	fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means ),cv2 )
	a0 <- unname( fit$coefficients["a0"] )
	a1 <- unname( fit$coefficients["a1tilde"])

	afit <- a1/means+a0
	varFitRatio <- vars/(afit*means^2)
	varorder <- order(varFitRatio,decreasing=T)
	oed <- ed[varorder,]
	
	df <- ncol(ed) - 1
	pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
	adj.pval <- p.adjust(pval,"fdr")
	sigVariedGenes <- adj.pval<1e-3;
	
	xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
	vfit <- a1/xg + a0
	
	if(plot){
		par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2)); lines( log(xg), log(vfit), col="black", lwd=3 ); 
		lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black"); lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black");
		# add top 100 genes
		points(log(means[sigVariedGenes]),log(cv2[sigVariedGenes]),col=2)
	}
	
	data.frame(means=means,cv2=cv2,varFitRatio=varFitRatio,pval.adj=adj.pval,is.signif=sigVariedGenes)
}

getMostVariableGenes2<-function(counts,minCount=0,plot=TRUE){
	counts<-counts[rowMeans(counts)>minCount,]
	dispTable<-data.frame(mu=log(rowMeans(counts)),disp=log(apply(counts,1,var)),row.names =rownames(counts))
	fit<-lm(disp~mu,data = dispTable)
	dispTable$residuals<-fit$residuals
	dispTable$cookD<-cooks.distance(fit)
	if(plot){
		require(ggplot2)
		print(ggplot(dispTable,aes(x=mu,y=disp,label=rownames(dispTable),color=residuals))+
						geom_point()+geom_abline(slope=fit$coefficients[2],intercept = fit$coefficients[1]))
	}
	dispTable
}

getMostVariableGenes3<-function(counts,minCount=0,plot=TRUE,topPlot=100){
	counts<-counts[rowMeans(counts)>minCount,]
	dispTable<-data.frame(mu=log(rowMeans(counts)),disp=log(apply(counts,1,var)),row.names =rownames(counts))
	fit<-loess(disp~mu,data = dispTable)
	dispTable$residuals<-fit$residuals
	dispTable$residuals2<-dispTable$residuals^2
	dispTable$fitted<-fit$fitted
	if(plot){
		require(ggplot2)
		require(ggrepel)
		require(circlize)
		print(
			ggplot(dispTable,aes(x=mu,y=disp,label=rownames(dispTable),fill=residuals))+
				geom_point(stroke=1/8,colour = "black",shape=21)+geom_line(aes(y=fitted),color="red",size=1.5)
		)
	}
	dispTable
}

getMostVariableGenes4<-function(counts,minCount=0.01,plot=TRUE,returnPlot=FALSE){
	counts<-counts[rowMeans(counts)>minCount,]
	dispTable<-data.frame(mu=rowMeans(counts),var=apply(counts,1,var),row.names =rownames(counts))
	dispTable$cv2<- dispTable$var / dispTable$mu^2
	sumNullvariance<-sum(dispTable$cv2 <= 0)
	if(sumNullvariance>0){
		warning(paste0(sumNullvariance, " have null variance and will be removed"))
		dispTable<-dispTable[dispTable$cv2 > 0,]
	}
	fit<-loess(cv2~mu,data = log10(dispTable[,c("mu","cv2")]))
	dispTable$residuals<-fit$residuals
	dispTable$residuals2<-dispTable$residuals^2
	dispTable$fitted<-10^fit$fitted
	if(plot){
		require(ggplot2)
		require(ggrepel)
		require(circlize)
		g<-ggplot(dispTable,aes(x=mu,y=cv2,label=rownames(dispTable),fill=residuals))+
			geom_point(stroke=1/8,colour = "black",shape=21)+geom_line(aes(y=fitted),color="red",size=1.5)+
			scale_x_log10()+scale_y_log10()
		if(returnPlot){
			return(g)
		}else{
			print(g)
		}
	}
	dispTable
}



getMarkers<-function(expressionMatrix,annotation){
	require(SC3)
	markers<-SC3::get_marker_genes(expressionMatrix, annotation)
	rownames(markers)<-rn(expressionMatrix)
	markers$clusts<-factor(markers$clusts,levels=(as.character(1:nlevels(annotation))))
	levels(markers$clusts)<-levels(annotation)
	markers
}


# By Miron Kursa https://mbq.me
auroc <- function(score, bool) {
	n1 <- sum(!bool)
	n2 <- sum(bool)
	U	<- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
	return(1 - U / n1 / n2)
}

getMarkers2<-function(expressionMatrix,group){
	if(length(group)!=ncol(expressionMatrix)) stop("group should be a vector with same length as number of column in expressionMatrix")
	group<-as.factor(group)
	group<-droplevels(group)
	binaryGroup<-vapply(levels(group),function(x){
		x==group & !is.na(group)
	},FUN.VALUE = numeric(length(group)))
	apply(binaryGroup,2,function(labels){
		apply(expressionMatrix,1,function(x){
			auroc(x,labels)
		})
	})
}


getMarkers3<-function(expressionMatrix,group,BPPARAM=NULL){
  require(BiocParallel)
  if(is.null(BPPARAM )) BPPARAM=bpparam()
  if(!is.matrix(expressionMatrix)) expressionMatrix<-as.matrix(expressionMatrix)
  if(length(group)!=ncol(expressionMatrix)) stop("group should be a vector with same length as number of column in expressionMatrix")
  
  group<-as.factor(group)
  group<-droplevels(group)
  binaryGroup<-lapply(levels(group),function(x){
    x==group & !is.na(group)
  });names(binaryGroup)<-levels(group)
  
  res<-as.matrix(data.frame(lapply(binaryGroup,function(labels){
    unlist(bplapply(seq_len(nrow(expressionMatrix)),function(i,expressionMatrix,labels,auroc){
      auroc(expressionMatrix[i,],labels)
    },labels=labels,expressionMatrix=expressionMatrix,auroc=auroc,BPPARAM=BPPARAM),recursive = FALSE)
  })))
  rownames(res)<-rownames(expressionMatrix)
  res
}


corGeneToOthers<-function(gene,expression,corFun=cor,...){
	expression<-as.matrix(expression)
	t(corFun(expression[gene,],t(expression),...))[,1]
}

eigengenes<-function(exprMatrix,genes,scale=F,center=T,returnContribution=F){
	pca<-prcomp(x = t(exprMatrix[genes,]),retx = T,center = center,scale = scale)
	eigen<-pca$x[,1]
	contribution<-pca$rotation[,1]
	if(cor(colMeans(exprMatrix[genes,]),eigen)<0){
		eigen<- -eigen
		contribution<- -contribution
	}
	if(returnContribution){
		list(eigengenes=eigen,contribution=contribution)
	}else{
		eigen
	}
}

reScale <- function(nonCorrected, corrected) {
	apply(nonCorrected, 1, function(x) max(x) - min(x))/
		(apply(corrected, 1, function(x) max(x)-min(x)))*
		(corrected-apply(corrected, 1, max))+apply(nonCorrected, 1, max)
}

volcanoPlot.DESeq2<-function(DEresult,formula,downLevel,upLevel,condColumn,padjThreshold,LFCthreshold,topGene=30){
	require(ggrepel)
	require(grid)
	DEresult<-DEresult[!is.na(DEresult$padj),]
	gene2Plot<-order(DEresult$padj)
	gene2Plot<-gene2Plot[DEresult[gene2Plot,"isDE"]!="NONE"]
	gene2Plot<-gene2Plot[1:min(topGene,length(gene2Plot))]
	g<-ggplot(DEresult,aes(x=log2FoldChange,y=-log10(padj)+0.01,color=isDE))+
		geom_point(size=1)+theme_bw()+scale_color_manual(values=c("#3AAA35","grey75","#E40429"))+
		geom_text_repel(data = DEresult[gene2Plot,],aes(x=log2FoldChange,y=-log10(padj),label=gene),
										inherit.aes = FALSE,color="black",fontface = "bold.italic",size=3)+
		ylab("-log10(adjusted pvalue)")+xlab(NULL)+
		geom_vline(xintercept = c(-LFCthreshold,LFCthreshold))+
		geom_hline(yintercept = -log10(padjThreshold)) + guides(color = FALSE)+
		ggtitle("Volcano plot")
	
	grid.newpage();
	pushViewport(viewport(x = 0, y = 0,
												width = .8, height = .1,
												just = c("left", "bottom")))
	filledDoubleArrow(x=.3,y=1,width = .3,just=c("left","center"),gp = gpar(fill="black"))
	grid.text(label = downLevel,x = .28,y = 1.03,just=c("right","center"),gp = gpar(fontface="bold"))
	grid.text(label = upLevel,x = .62,y = 1.03,just=c("left","center"),gp = gpar(fontface="bold"))
	grid.text(label = "log2(Fold-Change)",x = .45,y = 1.5,just = c("center","center"))
	popViewport()
	pushViewport(viewport(x = .65, y = 1,
												width = .35, height = 1,
												just = c("left", "top"),default.units = "npc"))
	grid.text(label = paste0("Experimental design:\n",formula),x = 0,y = .95,just=c("left","center"))
	grid.text(label = paste0("Results for\n",condColumn,":\n",downLevel," vs ",upLevel),x = .0,y = .8,just=c("left","center"))
	grid.text(label = paste0(sum(DEresult$isDE=="DOWNREG")," downreg. genes"),x = 0,y = .65,just=c("left","center"),
						gp = gpar(col="#3AAA35",fontface="bold"))
	
	grid.text(label = paste0(sum(DEresult$isDE=="UPREG")," upreg. genes"),x =.0,y = .55,just=c("left","center"),
						gp = gpar(col="#E40429",fontface="bold"))
	grid.text(label = paste0("From ",nrow(DEresult),"\ntested genes"),x = 0,y = .45,just=c("left","center"))
	popViewport()
	main_vp <- viewport(x = 0, y = 1,
											width = .8, height = .9,
											just = c("left", "top"))
	pushViewport(main_vp);print(g,vp=main_vp);popViewport()
}



customUpsetPlot<-function(genePerGroupList,universe=NULL){
	require(ComplexHeatmap)
	if(is.null(universe)) universe<- unique(unlist(genePerGroupList))
	isInGroupMatrix<-list_to_matrix(genePerGroupList,universal_set = universe)
	upsetMatrix<-make_comb_mat(isInGroupMatrix,mode = "intersect")
	upsetMatrix<-upsetMatrix[comb_degree(upsetMatrix) > 1] # retain only intersections of sets
	
	combsize = comb_size(upsetMatrix)
	setsize = set_size(upsetMatrix)
	
	#Are the intersections sets (or venn diagramm region) enriched or not ?
	regionEnrich<-sapply(comb_name(upsetMatrix),function(region){
		colOfcomp=which(strsplit(region,split = "")[[1]]=="1")
		intersectionEnrichment(isInGroupMatrix[,colOfcomp])
	})
	
	enrich_ha = HeatmapAnnotation(
		"enrichment" = anno_barplot(
			regionEnrich, gp = gpar(fill = "black"), height = unit(3, "cm"),axis_param = list(side = "left"),
			ylim = c(0, max(regionEnrich)*1.1)
		),
		annotation_name_side = "left", annotation_name_rot = 0,annotation_name_gp = gpar(fontface="bold"),
		annotation_label = "Enrichment\n(real/expected size)"
	)
	intersect_ha = HeatmapAnnotation(
		"intersection_size" = anno_barplot(
			combsize, gp = gpar(fill = "black"), height = unit(3, "cm"),axis_param = list(side = "left"),
			ylim = c(0, max(combsize)*1.1)
		),
		annotation_name_side = "left", annotation_name_rot = 0,annotation_name_gp = gpar(fontface="bold"),
		annotation_label = "Intersection\nsize"
	)
	set_size_ha = rowAnnotation(
		"set_size" = anno_barplot(
			setsize,gp = gpar(fill = "black"),width = unit(2, "cm"),
			ylim = c(0, max(setsize)*1.3)
		), 
		annotation_name_side = "bottom", annotation_name_rot = 0,annotation_name_gp = gpar(fontface="bold"),
		annotation_label = "Set\nsize"
	)
	
	
	ht = draw(UpSet(upsetMatrix,top_annotation = intersect_ha,bottom_annotation = enrich_ha,right_annotation = set_size_ha,
		border=TRUE,column_split=comb_degree(upsetMatrix),
		row_names_gp=gpar(fontsize=min(1/max(nchar(rownames(upsetMatrix)))*260,20))#automatic fontsize to avoid out of bound text
	)) 
	
	
	#Offset to counterbalance column split space
	colPerSplit=sapply(column_order(ht),length)
	offsetPerSplit=seq(0,length(colPerSplit)-1)
	offsets<-unlist(lapply(seq_along(colPerSplit),function(i) rep(offsetPerSplit[i],colPerSplit[i])),use.names = FALSE)
	
	rowOrder = rev(row_order(ht))
	columnOrder = unlist(column_order(ht))
	
	decorate_annotation("intersection_size", {
		grid.text(combsize[columnOrder], x = unit(seq_along(combsize),"native")+unit(offsets,"mm"), 
							y = unit(combsize[columnOrder], "native") + unit(6, "pt"), 
							default.units = "native", just = "center", gp = gpar(fontsize = 8))
	})
	decorate_annotation("enrichment", {
		grid.text(round(regionEnrich[columnOrder],2), x = unit(seq_along(regionEnrich),"native")+unit(offsets,"mm"), 
							y = unit(regionEnrich[columnOrder], "native") + unit(6, "pt"), 
							default.units = "native", just = "center", gp = gpar(fontsize = 8))
	})
	decorate_annotation("set_size", {
		grid.text(round(setsize[rowOrder],2), y = seq_along(setsize), x = unit(setsize[rowOrder], "native") + unit(7, "pt"), 
							default.units = "native", just = "center", gp = gpar(fontsize = 10),rot=-90)
	})
}


quickSCnorm<-function(rawCounts,returnLog=TRUE,sizeFactors=NULL,...){
	require(scran)
	sce <- SingleCellExperiment(assays=list(counts=rawCounts))
	if(!is.null(sizeFactors)){
		sizeFactors(sce)<-sizeFactors
	}else{
		sce <- computeSumFactors(sce,...)	
	}
	scater::normalizeCounts(sce,log=returnLog,pseudo_count = 0,size_factors = sizeFactors(sce))
}

