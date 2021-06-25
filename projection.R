plotPalette<-function(colorScale,continuousStep=NULL){
	if(is.null(continuousStep)){
		image(1:length(colorScale), 1, as.matrix(1:length(colorScale)), 
			col=colorScale,
			xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
		if(!is.null(names(colorScale))) axis(1,1:length(colorScale),names(colorScale))
	}else{
		br = round(seq(1,continuousStep,length.out = length(colorScale)))
		cols = colorRamp2(breaks = br,colors = colorScale)(1:continuousStep)
		image(1:continuousStep, 1, as.matrix(1:continuousStep), 
			col=cols,
			xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
	}
}


computeColorScaleFun<-function(colors,values,useProb=FALSE,probs = NULL,minProb=0.05,maxProb=0.95,
 midColorIs0=FALSE,returnColorFun=TRUE,returnGGscale=FALSE,geomAes="fill",geomArgument=list()){
	require(circlize)
	if(!useProb){
		breaks = seq(min(values),max(values),length.out = length(colors))
	}else{
		if(is.null(probs)){
			probs=seq(minProb,maxProb,length.out = length(colors))
		}
		breaks<-quantile(values,probs=probs)
	}
	if(midColorIs0 & (length(colors) %% 2 == 1)){
		breaks[ceiling(length(breaks)/2)]<-0
	}
	colorFun<-colorRamp2(breaks=breaks,colors=colors)
	if(returnGGscale){
		scaledBreaks<-linearScale(values,c(0,1),returnFunction = TRUE)(breaks)
		if(scaledBreaks[1]>0){
			scaledBreaks<-c(0,scaledBreaks)
			colors<-c(colors[1],colors)
		}
		if(scaledBreaks[length(scaledBreaks)]<1){
			scaledBreaks<-c(scaledBreaks,1)
			colors<-c(colors,colors[length(colors)])
		}
		
		geomArgument$values <- scaledBreaks
		geomArgument$colors <- colors
		return(do.call(paste0("scale_",geomAes,"_gradientn"),geomArgument))
	}
	if(returnColorFun){
		return(colorFun)
	}else{
		return(colorFun(values))
	}
}



genColorsForAnnots<-function(annots, colorScales=NULL,discreteFuns = list(mostDistantColor2),returnContinuousFun=FALSE ,
 continuousPalettes=list(c("#440154","#6BA75B","#FDE725"),c("#2EB538","#1D1D1B","#DC0900"),c("#FFFFC6","#FF821B","#950961")),...){
	
	if(is.null(colnames(annots))) stop("annots must have colnames.")
	annotNames <- colorScalesToGen <- colnames(annots)
	newColorScales<-list()
	
	if(!is.null(colorScales)){
		for(colorScaleName in names(colorScales)){
			if(!colorScaleName%in%colorScalesToGen) stop("Condition '",colorScaleName,"' does not match with existing condition names")
			colorScale<-colorScales[[colorScaleName]]
			annotVect<-annots[,colorScaleName]
			if(!is.null(names(colorScale))){ #factors
				if(is.numeric(annotVect)){
					warning(colorScaleName," is numeric but encoded as factors (color vector has names). It will be converted to factors.")
					annots[,colorScaleName]<-as.factor(as.character(sampleAnnot[,colorScaleName]))
					annotVect<-annots[,colorScaleName]
				}else if(!is.factor(annotVect)){
					stop(colorScaleName," is not factors or numeric, please check the sample annotation table.")
				}
				if(sum(!levels(annotVect) %in% names(colorScale))>0) stop("Levels of ",colorScaleName," are existing in sample annotation table but not in ",colorScaleFile)
				newColorScales[[colorScaleName]]<-colorScale[levels(annotVect)]
			}else{ #numeric
				if(!is.numeric(annotVect)) stop(colorScaleName," is not numeric but encoded as numeric (color vector has no names)")
				if(is.function(colorScale) & !returnContinuousFun) stop("You must not provide function in colorScales if returnContinuousFun=FALSE")
				if(!is.function(colorScale) & returnContinuousFun){
					newColorScales[[colorScaleName]]<-computeColorScaleFun(colorScale,values = annotVect,...)
				}else{
					newColorScales[[colorScaleName]]<-colorScale
				}
				
			}
		}
		colorScalesToGen<-setdiff(colorScalesToGen,names(newColorScales))
	}
	cN<-1
	cF<-1
	for(colorScaleName in colorScalesToGen){
		annotVect<-annots[,colorScaleName]
		if(is.numeric(annotVect)){
			newColorScales[[colorScaleName]]<-continuousPalettes[[cN]]
			cN<-cN+1
			if(cN>length(continuousPalettes)) cN<-1
		}else{
			annots[,colorScaleName]<-as.factor(as.character(annots[,colorScaleName]))
			annotVect<-annots[,colorScaleName]
			newColorScales[[colorScaleName]]<-discreteFuns[[cF]](nlevels(annotVect))
			names(newColorScales[[colorScaleName]])<-levels(annotVect)
			cF<-cF+1
			if(cF>length(discreteFuns)) cF<-1
		}
	}
	newColorScales
}

ggBorderedFactors<-function(gg,borderSize=.75,borderColor="grey75"){
	nX<-nlevels(as.factor(gg$data[,quo_name(gg$mapping$x)]))
	gg+geom_vline(xintercept = seq(1.5,nX -0.5, 1),size=borderSize,color=borderColor)+
		scale_x_discrete(expand = c(0,0.5, 0, 0.5))+
		theme(
			panel.grid.major.x = element_line(colour = NA),
			panel.grid.minor.x = element_line(colour = NA),
		)
}


PCA<-function(d,transpose=T,scale=F,center=T) {
	if(transpose) d<-t(d);
	means<-0;sdeviations<-1
	if(center){
		means<-apply(d,2,mean)
		d<-sweep(d,2,means,"-")
	}
	if(scale){
		sdeviations<-apply(d,2,sd)
		d<-sweep(d,2,sdeviations,"/")
	}
	resacp <-prcomp(x = d,retx = T,center = FALSE,scale = FALSE);
	resacp$n.obs<-dim(d)[1];
	resacp$percentVar<- resacp$sdev^2 / sum( resacp$sdev^2 )
	resacp$scale<-scale
	resacp$center<-center
	resacp$transform<-list(sdeviations=sdeviations,means=means)
	return(resacp);
}

ACP<-PCA

fastPCA <- function(d,transpose=T,scale=F,center=T,nPC=NULL,
		weight.by.var = TRUE, ...) {
	require(irlba)
	if(transpose) d<-t(d);
	means<-0;sdeviations<-1
	if(center){
		means<-apply(d,2,mean)
		d<-sweep(d,2,means,"-")
	}
	if(scale){
		sdeviations<-apply(d,2,sd)
		d<-sweep(d,2,sdeviations,"/")
	}
	if(is.null(nPC | nPC>=min(nrow(d),ncol(d)))) {
		nPC <- min(ncol(d) - 1, nrow(d) -1)
	}
	
	resacp<-list()
	resacp$n.obs<-dim(d)[1];
	resacp$scale<-scale
	resacp$center<-center
	resacp$transform<-list(sdeviations=sdeviations,means=means)
	
	irlbaResults <- irlba(A = d, nv = nPC, ...)
	rotation <- irlbaResults$v
	resacp$sdev <- irlbaResults$d/sqrt(max(1, nrow(d) - 1))
	if (weight.by.var) {
		reducedSpace <- irlbaResults$u %*% diag(irlbaResults$d)
	} else {
		reducedSpace <- irlbaResults$u
	}
	rownames(rotation) <- colnames(d)
	colnames(rotation) <- paste0("PC", 1:nPC)
	rownames(reducedSpace) <- rownames(d)
	colnames(reducedSpace) <- colnames(rotation)
	resacp$x<-reducedSpace
	resacp$rotation<-rotation
	resacp$percentVar<- resacp$sdev^2 / sum( resacp$sdev^2 )
	resacp
}

pcaAddSamples<-function(pca,newSamplesMatrix,transpose=TRUE,combineMat=TRUE){
	if(transpose) newSamplesMatrix<-t(newSamplesMatrix)
	newSamplesMatrix<-newSamplesMatrix[,rownames(pca$rotation),drop=FALSE]
	
	if(pca$center) newSamplesMatrix<-sweep(newSamplesMatrix,2,pca$transform$means,"-")
	if(pca$scale) newSamplesMatrix<-sweep(newSamplesMatrix,2,pca$transform$sdeviations,"/")
	
	newSamplesCoord<-t(apply(newSamplesMatrix,1,function(x){
		colSums(x*pca$rotation)
	}))
	if(combineMat){
		pca$x<-rbind(pca$x,newSamplesCoord)
		return(pca)
	}else{
		return(newSamplesCoord)
	}
}


PCR<-function(pca,annotationDF,nComponent=10){
	require(reshape2)
	annots<-colnames(annotationDF)
	rSquaredMat<-matrixFromDimnames(cn(pca$x),annots,value = NA)
	for(x in cn(pca$x)){
		for(annot in annots){
			rSquaredMat[x,annot]<-summary(lm(formula(paste0(x,"~",annot)),data = data.frame(annotationDF[,annot,drop=F],pca$x[,x,drop=F])))$r.squared
		}
	}
	retDt<-melt(rSquaredMat[1:nComponent,], value.name = "Rsquared",varnames=c("PC","Annotation"))
	retDt$PC<-factor(retDt$PC,levels=cn(pca$x)[1:nComponent]) #so the levels of PCs are well ordered

	retDt
}
	
# Fonctions présentes dans R pour interpréter une ACP :
# soient "acp" un objet de type acp et T la table des données
#
# coefficients de corrélation entre les variables cor(T)
# Inertie des axes, axe i : acpn$sdev^2 , acpn$sdev[i]^2
# Nouvelles coordonnées des individus : acpn$scores
# De l'individu i sur CP j : acpn$scores[i,j]
# Graphique des inerties (valeurs propres) : plot(acp)
# Plan principal : biplot(acp)
# Plan i x j : biplot(acp,c(c1=i, c2=j))
#Tableau des inerties et des pourcentages cumulés
# Paramètre : résultat ACP
# Sortie : tableau des inerties, pourcentages et pourcentages cumulés
VP <- function(resacp) {
	N <- length(resacp$sdev);
	tab <- matrix(nrow=N,ncol=3);
	s <- sum(resacp$sdev^2);
	s1 <- 0;
	for (i in 1:N) {
		tab[i,1] <- resacp$sdev[i]^2;
		tab[i,2] <- resacp$sdev[i]^2/s;
		s1 <- s1+resacp$sdev[i]^2;
		tab[i,3] <- 100*s1/s;
	};
	return(tab)}
# Corrélations entre les axes et les variables initiales
# Paramètres :	table R des données
#		résultat ACP (produit par princomp)
# Sortie :	la matrice des corrélations
AXEVAR <- function(resacp) {return(resacp$rotation)}
# Corrélations entre les k premiers axes et les variables initiales
# Paramètres :	table R des données
#		résultat ACP (produit par princomp)
#		nombre axes
# Sortie :	la matrice des corrélations
AXEVARk <- function(d,resacp,k) {
	return(resacp$rotation[,1:k])
}
# Contribution de la ligne i à l'inertie de l'axe j
# Paramètres :	résultat ACP
#		numéro ligne
#		numéro axe
# Sortie :	pourcentage de la contribution
CTRij <- function(resacp,i,j) {
	x <- resacp$rotation[i,j]^2/(resacp$n.obs * resacp$sdev[j]^2);
        x <- 100*x;
	return(x)}
# Tableau des contribution des lignes aux axes
# Paramètres :	résultat ACP
#		nombre axes
# Sortie :	tableau des pourcentages des contributions
CTR <- function(resacp, nbax) {
      	matrice <- matrix(nrow=resacp$n.obs,ncol=nbax);
        row.names(matrice) <- row.names(resacp$x);
        for (j in 1:nbax) 
		for (i in 1:resacp$n.obs) matrice[i,j] <- CTRij(resacp,i,j);
     
       return(matrice)}
# Fonction utilitaire
SOMME2 <- function(resacp) {
	N <- resacp$n.obs ; M <- ncol(resacp$x);
	s2 <- vector (mode = "numeric", N);
	for (i in 1:N)
		for (j in 1:M) s2[i] <- s2[i] + resacp$x[i,j]^2;
	return(s2)
}
# Cosinus ** 2 des angles de projection
# Paramètres :	résultat ACP
#		nombre axes
# Sortie :	tableau des cos2 des angles de projection
COS2TETA <- function(resacp, nbax) {
	N <- resacp$n.obs ; 
	c2teta <- matrix(nrow=N,ncol=nbax);
	row.names(c2teta) <- row.names(resacp$x);
	s2 <- SOMME2(resacp);
	for (i in 1:N)
		for (j in 1:nbax) c2teta[i,j] <- resacp$x[i,j]^2 / s2[i];
	return(c2teta)
}

# visualisation 3d ggplot
acp3d<-function(pca, comp=1:3, group=rep(1,pca$n.obs), plotVars = FALSE, pointSize=2, plotText=FALSE){
	if(!require("rgl")) stop("You must install rgl");
	if(length(comp)!=3) stop("You must give a vector of 3 integer for comp parameter")
	if(!plotVars){
		x<-pca$x
	}else{
		x<-pca$rotation
	}
	if(is.null(levels(group))){ colors="black"}
	else{
		hashCol<-rainbow(nlevels(group))
		names(hashCol)<-levels(group)
		colors<-hashCol[group]
	}
	percentVar <- pca$percentVar
	plot3d(x[,comp[1]],x[,comp[2]],x[,comp[3]],
		xlab=paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance"), 
		ylab=paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance"), 
		zlab=paste0("PC",comp[3],": ",round(percentVar[comp[3]] * 100),"% variance"),
	col=colors,size=pointSize,type=ifelse(plotText,"n","p"))
	
	legend3d("topright", legend = names(hashCol), pch = 16, col = hashCol, cex=1, inset=c(0.02))
	
	if(plotText) text3d(x[,comp[1]],x[,comp[2]],x[,comp[3]],texts=rownames(x),cex=pointSize,col=colors)
	if(plotVars) spheres3d(x=0,y=0,z=0, radius = 1,alpha=0.5,color="white")
	spheres3d(x=0,y=0,z=0, radius = 0.005,alpha=1,color="red")
}
# visualisation 2d ggplot

pca2d<-function(pca, comp=1:2,group=NULL, plotVars = FALSE, pointSize=2, plotText=FALSE,fixedCoord=FALSE,main=NULL,
								ellipse=FALSE,colorScales=NULL,returnGraph=FALSE, outlierLabel=FALSE,outlierLabelThres=NULL,
								outlierLabelSize=3,...){
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(comp)!=2) stop("You must give a vector of 2 integer for comp parameter");
	percentVar <- pca$percentVar
	dt<-pca[[ifelse(plotVars, "rotation","x")]]
	graph<-proj2d(coord = dt,axis = comp,group = group,returnGraph = T,colorScale = colorScales,main = main,plotText = plotText,
								ellipse = ellipse,pointSize=pointSize,...)+
		xlab(paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance")) +
		ylab(paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance"))
	if(fixedCoord)graph <- graph + coord_fixed(ratio=percentVar[comp[2]]/percentVar[comp[1]])
	
	if(outlierLabel){
		require(ggrepel)
		dtOutlier<-data.frame(x=dt[,comp[1]],y=dt[,comp[2]],label=rn(dt))
		dens<-pointdensity.nrd(dtOutlier[,1:2]);
		if(is.null(outlierLabelThres)) outlierLabelThres=Mode(dens)/3
		dtOutlier<-dtOutlier[dens/max(dens)<outlierLabelThres,]
		graph<-graph+geom_text_repel(data = dtOutlier,mapping = aes(x=x,y=y,label=label),
																 size=outlierLabelSize,fontface="italic",inherit.aes = FALSE)
	}
	
	if(returnGraph){
		return(graph)
	}else{
		print(graph)
	}
}

acp2d.2<-pca2d



##########################
#Autre type de projection#
##########################
merge0dist<-function(disMat){
	mat<-as.matrix(disMat)
	merged<-list()
	found<-TRUE
	while(found==TRUE){
	  found<-FALSE
	  for(i in 2:nrow(mat)){
		for(j in 1:(i-1)){
		  if(mat[i,j]==0){
			newNames<-rownames(mat)
			newNames<-newNames[-i]
			newMat<-mat[-i,-i]
			colnames(newMat)<-rownames(newMat)<-newNames
			merged[[rownames(mat)[j]]]<-c(merged[[rownames(mat)[j]]],rownames(mat)[i])
			mat<-newMat
			found<-TRUE
			break
		  }
		}
		if(found) break
	  }
	}
	return(list(distMat=as.dist(mat),merged=merged))
}
NMDS<-function(data,transpose=TRUE,scale=FALSE,center=FALSE,metric=dist,ndim=2,maxit=100){
	merged<-FALSE
	require(MASS)
	if(transpose) data <- t(data)
	d <- metric(data)  # euclidean distances between the rows
	if(min(d,na.rm=TRUE)==0){
		merged<-TRUE
		md<-merge0dist(d)
		d<-md$distMat
		mergedSample<-md$merged	
	}
	fit <- isoMDS(d, k=ndim, maxit=maxit) # k is the number of dim
	fit$coord<-fit$points
	fit$points<-NULL
	if(merged){
		for(sple in names(mergedSample)){
			values<-matrix(rep(fit$coord[sple,],length(mergedSample[[sple]])),nrow=length(mergedSample[[sple]]),byrow = TRUE)
			rownames(values)<-mergedSample[[sple]]
			fit$coord<-rbind(fit$coord,values)
		}
	}
	return(fit)
}

proj2d<-function(coord,group=NULL,axis=1:2, pointSize=3, plotText=FALSE,main=NULL,alpha=9/10, 
		ellipse=FALSE,emph=NULL,colorScale=NULL,returnGraph=FALSE,legendTitle="Values",axis.names=NULL,
		na.color="grey50",na.bg=TRUE,obj=NULL,plotFactorsCentroids=FALSE,
		pointStroke=1/8,strokeColor="black",funAutoColorScale=ggplotColours,fixedCoord=FALSE,plotLabelRepel=FALSE,labelSize=3,
		sizeProp2Dens=FALSE,densEps=1,nnMatrix=NULL,nnSegmentParam=list(alpha=.75,size=.1)){
	coord<-data.frame(coord)
	if(is.null(coord)) coord<-obj$coord
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(axis)!=2) stop("You must give a vector of 2 integer for axis parameter");
	if((!is.null(axis.names)) & length(axis.names)!=2) stop("Error, if given axis.names parameter must contain 2 values.");
	d <- data.frame(lab=rownames(coord),Axis1=coord[,axis[1]], Axis2=coord[,axis[2]],sizeMultiplier=rep(pointSize,nrow(coord)));
	if(sizeProp2Dens){
		dens<-pointdensity.nrd(d[,c("Axis1","Axis2")],eps = densEps)
		d$sizeMultiplier<-d$sizeMultiplier*(1-dens/max(dens))*2
	}
	if(is.null(group)){
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab))
	}else{
		if(is.data.frame(group) | is.matrix(group)){
			if(!is.null(colnames(group))) legendTitle<- colnames(group)[1]
			group<-group[,1]
		}
		
		if(is.character(group)) group<-as.factor(group)
		d$group<-group;
		groupIsFactor<-is.factor(group)
		if(!is.null(emph)){
			if(!emph%in%levels(group)) stop("emph value not in levels of group")
			d$group<-as.character(d$group)
			d$group[which(d$group!=emph)]<-NA
			d$group<-as.factor(d$group)
			d$sizeMultiplier[which(d$group==emph)]<-d$sizeMultiplier[which(d$group==emph)]
		}
		if(na.bg){
			indexNA<-which(is.na(d$group))
			indexNotNA<-which(!is.na(d$group))
			if(length(indexNA)>0){
				tempd<-d
				tempd[1:length(indexNA),]<-d[indexNA,]
				tempd[(length(indexNA)+1):nrow(d),]<-d[indexNotNA,]
				d<-tempd
			}
		}
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab,color=group,fill=group))
		if(is.null(colorScale)){
			colorScale<-c("grey","red")
			if(groupIsFactor) colorScale<-funAutoColorScale(nlevels(group))
		}
		funColorScaleFill<-paste0("scale_fill_",ifelse(groupIsFactor,"manual","gradientn"))
		funColorScaleColor<-paste0("scale_color_",ifelse(groupIsFactor,"manual","gradientn"))
		paramColorScale<-list("name"=legendTitle,na.value=na.color)
		paramColorScaleType<-ifelse(groupIsFactor,"values","colors");paramColorScale[[paramColorScaleType]]<-colorScale
		graph<-graph+do.call(funColorScaleFill,paramColorScale)
		graph<-graph+do.call(funColorScaleColor,paramColorScale)
	}
	if(!is.null(nnMatrix)){
		if(!is.matrix(nnMatrix)) stop("If nnMatrix is not null, it should be a matrix !")
		retainCoord<-as.matrix(coord[,axis])
		segmentsCoord<-data.frame(do.call("rbind",lapply(1:nrow(retainCoord),function(i){
			t(vapply(nnMatrix[i,],FUN.VALUE = numeric(4),function(x){
				c(retainCoord[i,],retainCoord[x,])
			}))
		})));colnames(segmentsCoord)<-c("x","y","xend","yend")
		graph<-graph+do.call(
			"geom_segment",
			c(list(data = segmentsCoord,mapping = aes(x=x,y=y,xend=xend,yend=yend),inherit.aes = FALSE),nnSegmentParam)
		)
	}
	if(plotText){
		graph<-graph+geom_text(alpha=alpha,size=d$sizeMultiplier)
	}else{
		if(is.null(group)){
			graph<-graph+geom_point(stroke=pointStroke,colour = strokeColor,shape=21,alpha=alpha,fill="black",size=d$sizeMultiplier)
		}else{
			graph<-graph+geom_point(stroke=pointStroke,colour = strokeColor,shape=21,alpha=alpha,size=d$sizeMultiplier)
		}	
	}
	if(plotLabelRepel){
		require(ggrepel)
		graph<-graph+geom_text_repel(color="black",size=labelSize)
	}
	if(is.null(axis.names)) axis.names<-paste0("Axis",axis)
	graph<-graph+xlab(axis.names[1]) + ylab(axis.names[2])
	if(!is.null(main)) graph <- graph + ggtitle(main)
	if(ellipse) graph<-graph+stat_ellipse(size=.5)
	if(fixedCoord)graph <- graph + coord_fixed()
	if(plotFactorsCentroids){
		samplePerFactor<-lapply(levels(group),function(x) which(group==x))
		names(samplePerFactor)<-levels(group)
		centroids<-data.frame(t(sapply(samplePerFactor,function(x) colMeans(d[x,c("Axis1","Axis2")]))),groupName=names(samplePerFactor))
		graph<-graph+geom_label(data = centroids,mapping = aes(x=Axis1, y=Axis2, label = groupName),inherit.aes = FALSE)
	} 

	graph<-graph+theme(
		panel.background = element_rect(fill = NA,colour="black"),
		panel.grid.major = element_line(colour = NA),
		panel.grid.minor = element_line(colour = NA)
	) + guides(size=FALSE)
	if(returnGraph){
		return(graph)
	}else{
		print(graph)
	}
}

proj_densHex<-function(obj=NULL,coord=NULL, axis=1:2,group=NULL, main=NULL,bins=30,
 emph=NULL,colorScale=NULL,returnGraph=FALSE,legendTitle="Values",axis.names=NULL,na.color="grey50",na.bg=TRUE){
	if(!(is.factor(group)|is.null(group))) stop("Error, group must be factor or null.")
	if(is.null(coord)) coord<-obj$coord
	if(!require("ggplot2")) stop("You must install ggplot2");
	if(length(axis)!=2) stop("You must give a vector of 2 integer for axis parameter");
	if((!is.null(axis.names)) & length(axis.names)!=2) stop("Error, if given axis.names parameter must contain 2 values.");
	d <- data.frame(lab=rownames(coord),Axis1=coord[,axis[1]], Axis2=coord[,axis[2]]);
	if(is.null(group)){
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab))
	}else{
		d$group<-group;
		groupIsFactor<-is.factor(group)
		if(is.null(emph))  stop("if group provided, you must provide emph")
		if(!emph%in%levels(group)) stop("emph value not in levels of group")
		d$group<-as.character(d$group)
		d$group[which(d$group!=emph)]<-NA
		d$group<-as.factor(d$group)
		
		if(na.bg){
			indexNA<-which(is.na(d$group))
			indexNotNA<-which(!is.na(d$group))
			if(length(indexNA)>0){
				tempd<-d
				tempd[1:length(indexNA),]<-d[indexNA,]
				tempd[(length(indexNA)+1):nrow(d),]<-d[indexNotNA,]
				d<-tempd
			}
		}
		graph<-ggplot(data=d, mapping = aes(x=Axis1, y=Axis2, label = lab,color=group))
	}
	if(is.null(colorScale)) colorScale=c("grey","red")
	graph<-graph+scale_fill_gradientn(colors=colorScale)
	graph<-graph+stat_bin_hex(bins=bins)
	if(is.null(axis.names)) axis.names<-paste0("Axis",axis)
	graph<-graph+xlab(axis.names[1]) + ylab(axis.names[2])
	if(!is.null(main)) graph <- graph + ggtitle(main)
	graph<-graph+theme(
		panel.background = element_rect(fill = NA,colour="black"),
		panel.grid.major = element_line(colour = NA),
		panel.grid.minor = element_line(colour = NA)
	) 
	if(returnGraph) return(graph)
	print(graph)
}

proj2dWideColor<-function(coord, axis=1:2,group=NULL, pointSize=3, plotText=FALSE,main=NULL,alpha=9/10, 
ellipse=FALSE,emph=NULL,returnGraph=FALSE,legendTitle="Values",axis.names=NULL,na.color="grey50",na.bg=TRUE,obj=NULL,funAutoColorScale=ggplotColours,...){
	proj2d(coord, axis=axis,group=group, pointSize=pointSize, plotText=plotText,main=main,alpha=alpha, 
	ellipse=ellipse,emph=emph,colorScale=c("white","yellow","red","purple"),returnGraph=returnGraph,
	legendTitle=legendTitle,axis.names=axis.names,na.color=na.color,na.bg=na.bg,obj=obj,funAutoColorScale=ggplotColours,...)
}

proj3d<-function(coord, axis=1:3, group=NULL, pointSize=5, plotText=FALSE,colorScale=NULL,na.color="grey50",na.bg=TRUE,alpha=1,...){
	if(!(is.factor(group)|is.numeric(group)|is.null(group))) stop("Error, group must be numerical, factor or null.")
	if(!require("rgl")) stop("You must install rgl");
	if(!require("circlize")) stop("You must install circlize");
	if(length(axis)!=3) stop("You must give a vector of 3 integer for axis parameter")
	
	if(is.null(group)){ colors="black"}
	else{
		if(na.bg){
			indexNA<-which(is.na(group))
			indexNotNA<-which(!is.na(group))
			if(length(indexNA)>0){
				tempc<-coord;tempc[1:length(indexNA),]<-coord[indexNA,];tempc[(length(indexNA)+1):nrow(coord),]<-coord[indexNotNA,];coord<-tempc;
				tempg<-group;tempg[1:length(indexNA)]<-group[indexNA];tempg[(length(indexNA)+1):length(group)]<-group[indexNotNA];group<-tempg;
			}
		}
		if(is.factor(group)){
			if(is.null(colorScale)) colorScale<-rainbow(nlevels(group))
			hashCol<-colorScale
			names(hashCol)<-levels(group)
			colors<-hashCol[group]
		}else{
			if(is.null(colorScale)) colorScale<-c("grey","red")
			funCol<-colorRamp2(seq(min(group),max(group),length.out=length(colorScale)),colorScale)
			colors<-funCol(group)
		}
	}
	colors[is.na(colors)]<-na.color
	plot3d(coord[,axis[1]],coord[,axis[2]],coord[,axis[3]],
		xlab=paste0("Axis",axis[1]), 
		ylab=paste0("Axis",axis[2]), 
		zlab=paste0("Axis",axis[3]),
		type="n",...)
	if(plotText){
		text3d(coord[,axis[1]],coord[,axis[2]],coord[,axis[3]],texts=rownames(coord),cex=pointSize,col=colors,alpha=alpha)
	}else{
		spheres3d(coord[,axis[1]],coord[,axis[2]],coord[,axis[3]],col=colors,radius=pointSize,alpha=alpha)
	}
	if(!is.null(group)){
		if(is.factor(group)) legend3d("topright", legend = names(hashCol), pch = 16, col = hashCol, cex=1, inset=c(0.02))
		if(is.numeric(group)) legend3d("topright", legend = c(min(group),max(group)), pch = 16, col = colorScale, cex=1, inset=c(0.02))
	}
}

#' Add grids to a scatterplot3d
#' 
#' @description The goal of this function is to add grids on an existing
#'  plot created using the package scatterplot3d
#' @param x,y,z numeric vectors specifying the x, y, z coordinates of points.
#'  x can be a matrix or a data frame containing 3 columns corresponding to
#'  the x, y and z coordinates. In this case the arguments y and z are optional
#' @param grid specifies the facet(s) of the plot on which grids should be drawn.
#'  Possible values are the combination of "xy", "xz" or "yz".
#'  Example: grid = c("xy", "yz"). The default value is TRUE to add grids only on xy facet.
#' @param col.grid,lty.grid color and line type to be used for grids
#' @param lab a numerical vector of the form c(x, y, len).
#'  The values of x and y give the (approximate) number of tickmarks on the x and y axes.
#' @param lab.z the same as lab, but for z axis
#' @param scale.y of y axis related to x- and z axis
#' @param angle angle between x and y axis
#' @param "xlim, ylim, zlim" the x, y and z limits (min, max) of the plot.
#' 
#' @note
#' Users who want to extend an existing scatterplot3d graphic with the
#'  function addgrids3d, should consider to set the arguments scale.y, angle, ...,
#'  to the value used in scatterplot3d.
#' 
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' 
#' @example
#' library(scatterplot3d)
#' data(iris)
#' scatterplot3d(iris[, 1:3], pch = 16, grid=T, box=F)
#' addgrids3d(iris[, 1:3], grid = c("xy", "xz", "yz"))
addgrids3d <- function(x, y=NULL, z=NULL, grid = TRUE,
                    col.grid = "grey", lty.grid = par("lty"),
                    lab = par("lab"), lab.z = mean(lab[1:2]),
                    scale.y = 1, angle = 40,
                    xlim=NULL, ylim=NULL, zlim=NULL){
  
  
  if(inherits(x, c("matrix", "data.frame"))){
    x <- as.data.frame(x)
    y <- unlist(x[,2])
    z <- unlist(x[,3])
    x <- unlist(x[,1])
  }
  
  p.lab <- par("lab")
  
  angle <- (angle%%360)/90
  yz.f <- scale.y * abs(if (angle < 1) angle else if (angle >3) angle - 4 else 2 - angle)
  yx.f <- scale.y * (if (angle < 2) 1 - angle else angle - 3)
  
  
  # x axis range
  x.range <- range(x[is.finite(x)], xlim)
  x.prty <- pretty(x.range, n = lab[1], min.n = max(1, min(0.5 *lab[1], p.lab[1])))
  x.scal <- round(diff(x.prty[1:2]), digits = 12)
  x <- x/x.scal
  x.range <- range(x.prty)/x.scal
  x.max <- ceiling(x.range[2])
  x.min <- floor(x.range[1])
  if (!is.null(xlim)) {
    x.max <- max(x.max, ceiling(xlim[2]/x.scal))
    x.min <- min(x.min, floor(xlim[1]/x.scal))
  }
  x.range <- range(x.min, x.max)
  
  # y axis range
  y.range <- range(y[is.finite(y)], ylim)
  y.prty <- pretty(y.range, n = lab[2], min.n = max(1, min(0.5 *lab[2], p.lab[2])))
  y.scal <- round(diff(y.prty[1:2]), digits = 12)
  y.add <- min(y.prty)
  y <- (y - y.add)/y.scal
  y.max <- (max(y.prty) - y.add)/y.scal
  if (!is.null(ylim))
    y.max <- max(y.max, ceiling((ylim[2] - y.add)/y.scal))
  
  # Z axis range
  z.range <- range(z[is.finite(z)], zlim)
  z.prty <- pretty(z.range, n = lab.z, min.n = max(1, min(0.5 *lab.z, p.lab[2])))
  z.scal <- round(diff(z.prty[1:2]), digits = 12)
  z <- z/z.scal
  z.range <- range(z.prty)/z.scal
  z.max <- ceiling(z.range[2])
  z.min <- floor(z.range[1])
  if (!is.null(zlim)) {
    z.max <- max(z.max, ceiling(zlim[2]/z.scal))
    z.min <- min(z.min, floor(zlim[1]/z.scal))
  }
  z.range <- range(z.min, z.max)
  
  # Add grid
  if ("xy" %in% grid || grid == TRUE) {
    i <- x.min:x.max
    segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + 
               z.min, col = col.grid, lty = lty.grid)
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min, x.max + 
               (i * yx.f), i * yz.f + z.min, col = col.grid, lty = lty.grid)
  }
   
  if ("xz" %in% grid) {
    i <- x.min:x.max
    segments(i + (yx.f * y.max), yz.f * y.max + z.min, 
             i + (yx.f * y.max), yz.f * y.max + z.max, 
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.max + temp,temp1 + i , col = col.grid, lty = lty.grid)
    
  }
  
  if ("yz" %in% grid) {
    i <- 0:y.max
    segments(x.min + (i * yx.f), i * yz.f + z.min,  
             x.min + (i * yx.f) ,i * yz.f + z.max,  
             col = col.grid, lty = lty.grid)
    temp <- yx.f * y.max
    temp1 <- yz.f * y.max
    i <- z.min:z.max
    segments(x.min + temp,temp1 + i, 
             x.min, i , col = col.grid, lty = lty.grid)
    }
  
}


#Rotate point
rotate<-function(coord,angle,center=c(0,0)){
  x<-coord[,1]
  y<-coord[,2]
  xP<-x*cos(angle)+y*sin(angle)
  yP<- -x*sin(angle)+y*cos(angle)
  centerRot<-c(center[1]*cos(angle)+center[2]*sin(angle),-center[1]*sin(angle)+center[2]*cos(angle))
  transX<-center[1]-centerRot[1]
  transY<-center[2]-centerRot[2]
  return(cbind(xP+transX,yP+transY))
}

rotateSubspace<-function(coord,angle,center=c(0,0),which.sample=rownames(coord)){
	row_names<-rownames(coord)
	coord2Rotate<-coord[which.sample,]
	coordNotRotated<-coord[!row_names%in%which.sample ,]
	rotatedCoord<-rotate(coord2Rotate,angle,center=center)
	newCoord<-rbind(rotatedCoord,coordNotRotated)
	newCoord[row_names,]
}



rotate3dim<-function(coord,angle,center=c(0,0,0),axes=c(1,2)){
	if(length(axes)!=2) stop("you must give two axes for rotating the 3 space dimension")
	rotCoord<-rotate(coord[,axes],angle,center = center[axes])
	newCoord<-coord
	for(i in 1:2) newCoord[,axes[i]]<-rotCoord[,i]
	newCoord
}

rotate3dimSubSpace<-function(coord,angle,center=c(0,0,0),axes=c(1,2),which.sample=rownames(coord)){
	row_names<-rownames(coord)
	coord2Rotate<-coord[which.sample,]
	coordNotRotated<-coord[!row_names%in%which.sample ,]
	rotatedCoord<-rotate3dim(coord2Rotate,angle,center=center,axes=axes)
	newCoord<-rbind(rotatedCoord,coordNotRotated)
	newCoord[row_names,]
}


circlePlot = function(
  adjacency,
  labels,
  order,
  startNewPlot = TRUE,
  plotBox = c(-1, 1, -1, 1),
  center = c(0,0), 
  radii = c(0.8, 0.8),
  startAngle = 0,
  variable.cex.labels = TRUE,
  min.cex.labels = 1,
  max.cex.labels = 1.5,
  variable.cex.points = TRUE,
  min.cex.points = 1,
  max.cex.points = 3,
  variable.line.width = TRUE,
  min.line.width = 1,
  max.line.width = 5,
  lineColors = grey2red(50, 0.6, 1),
  pch = 21,
  labelColors = "black",
  pointColors = "black",
  pointBg = "black",
  xMargin = 1-radii[1],
  yMargin = 1-radii[2],
  xLabelOffset = 0.01,
  yLabelOffset = 0.01,
  variableLabelAngle = TRUE,
  
  ...)
{

  if (startNewPlot)
    plot(plotBox[1:2], plotBox[3:4], axes = FALSE, type = "n", xlab = "", ylab = "", ...);

  # plot(c(-1-xMargin,1+xMargin), c(-1-yMargin,1+yMargin), axes = FALSE, type = "n", xlab = "", ylab = "", ...) 
  checkAdjMat(adjacency, min = -1)
  n = length(labels);
  angles = seq(from = startAngle, to = startAngle + 2*pi * (1-1/n), length.out = n);
  x = center[1] + radii[1] * sin(angles);  # This is intentional; top should correspond to angle=0
  y = center[2] + radii[2] * cos(angles);

  adjx = adjacency
  adjx[is.na(adjx)] = 0;
  connectivity = apply(abs(adjx), 2, sum)-diag(adjx)
  minConn = min(connectivity, na.rm = TRUE);
  maxConn = max(connectivity, na.rm = TRUE);

  if (length(pch)==1) pch = rep(pch, n);
  if (length(labelColors)==1) labelColors = rep(labelColors, n);
  if (length(pointColors)==1) pointColors = rep(pointColors, n);
  if (length(pointBg)==1) pointBg = rep(pointBg, n);
  if (length(xLabelOffset)==1) xLabelOffset = rep(xLabelOffset, n);
  if (length(yLabelOffset)==1) yLabelOffset = rep(yLabelOffset, n);

  oLabs = labels[order]
  oLColors = labelColors[order];
  oPColors = pointColors[order];
  oPBg = pointBg[order];
  oConn = connectivity[order];
  oAdj = adjx[order, order];
  oPch = pch[order];

  actualCexPts = rep(0, n);
  for (node in 1:n)
  {
    cex = min.cex.points;
    if (variable.cex.points)
       cex = min.cex.points + (max.cex.points - min.cex.points) * 
                                  (oConn[node] - minConn)/(maxConn - minConn)
    actualCexPts[node] = cex
  }

  diag(oAdj) = 0;
  maxA = max(abs(oAdj));
  if (sum(oAdj < 0) > 0)
  {
     adjCol = numbers2colors(oAdj, signed = TRUE, lim = c(-maxA, maxA));
  } else {
     adjCol = numbers2colors(oAdj, signed = FALSE, lim = c(0, maxA));
  }


  ltA = oAdj;
  diag(ltA) = NA;
  ltA[upper.tri(ltA)] = NA;

  adjOrder = order(c(abs(ltA)))
  rows = row(oAdj)[adjOrder];
  cols = col(oAdj)[adjOrder];

  nLines = n*(n-1)/2;
  for (line in 1:nLines)
  {
    n1 = rows[line];
    n2 = cols[line];
    a = oAdj[n1, n2];
    normA = abs(a)/maxA;

    w = min.line.width;
    if (variable.line.width)
      w = min.line.width + (max.line.width - min.line.width) * normA;

    #pRadius1 = par("cxy") * actualCexPts[n1]/35;  # Emprical fudge factor..
    #pRadius2 = par("cxy") * actualCexPts[n2]/35;
    lineLen = sqrt( (x[n1] - x[n2])^2 + (y[n1] - y[n2])^2);
    x1 = x[n1] #+ pRadius1[1] * (x[n2] - x[n1]) / lineLen
    y1 = y[n1] #+ pRadius1[1] * (y[n2] - y[n1]) / lineLen
    x2 = x[n2] #+ pRadius2[1] * (x[n1] - x[n2]) / lineLen
    y2 = y[n2] #+ pRadius2[1] * (y[n1] - y[n2]) / lineLen

    lines(c(x1,x2),c(y1, y2), lwd = w, col = adjCol[n1, n2]);
  }

  for (node in 1:n)
    points(x[node], y[node], pch = oPch[node], cex = actualCexPts[node], bg = oPBg[node], col = oPColors[node]);

  for (node in 1:n)
  {
    cex = min.cex.labels;
    if (variable.cex.labels)
       cex = min.cex.labels + (max.cex.labels - min.cex.labels) *
                                 (oConn[node] - minConn)/(maxConn - minConn)
    textWidth = strwidth(oLabs[node], cex = cex);
    textHeight = strheight(oLabs[node], cex = cex);
    if (variableLabelAngle)
    {
      ang = angles[node]/pi * 180;
      if (ang < 180) 
      {
         dir = 1;
      } else {
         dir = -1;
         ang = ang - 180;
      }
      ang = (90 - ang)/2
      xDir = 1;
      yDir = 1;
      cosAng = cos(ang/180*pi);
      sinAng = sin(ang/180*pi);
    } else {
      ang = 0;
      xDir = x[node];
      yDir = y[node];
      cosAng = 1;
      sinAng = 1;
      dir = 1;
    }
    angRad = ang/180*pi;
    pRadius = par("cxy") * actualCexPts[node]/5  ;  # Emprical fudge factor..
    effPointRadius = sqrt(sum(c(cosAng^2, sinAng^2) * pRadius^2));
    rotMat = matrix( c(cosAng, sinAng, -sinAng, cosAng), 2, 2);
    labelShift = rotMat %*% as.matrix(c(textWidth, textHeight));
    text(x[node] + dir * xDir * (labelShift[1]/2 + cosAng * effPointRadius + xLabelOffset[node]), 
         y[node] + dir * yDir * (labelShift[2]/2 + sinAng * effPointRadius + yLabelOffset[node]), 
         labels = oLabs[node], adj = c(0.5, 0.5), 
         cex = cex, col = oLColors[node], srt = ang, xpd = TRUE);
  }
}

make.umap2<-function(data,nDimPCA=NULL,transpose=TRUE,n_neighbors=NULL, n_components = 2,min_dist=0.01,
										 init = "laplacian", metric = "euclidean",ret_model=FALSE,ret_nn=FALSE,...){
	require(uwot)
	if(transpose) data<-t(data)
	if(is.null(n_neighbors)) n_neighbors=nrow(data)
	if(!is.null(nDimPCA)){
		data<-fastPCA(data,transpose = FALSE,scale = FALSE,nPC = nDimPCA)$x
	}
	res<-uwot::umap(as.matrix(data),n_neighbors = n_neighbors, n_components = n_components,
									min_dist=min_dist, init = init, metric = metric,ret_model=ret_model,ret_nn = ret_nn,...)
	if(!ret_model & !ret_nn) rownames(res)<-rownames(data)
	res
}


serializeCurve=function(x){
	coordList<-sapply(x,function(el){
		return(paste0("\t[\n",paste0(apply(el,1,function(row){
			return(paste0("\t\t{\"x\":",row[1],",\"y\":",row[2],"}"))
		}),collapse=",\n"),"\n\t]"))
	})
	return(paste0("[\n",paste0(coordList,collapse = ",\n"),"\n]"))
}


qplotDensity<-function(x,printGraph=TRUE,...){
	require(ggplot2)
	if(class(x)=="data.frame" | class(x)=="list") x<-unlist(x)
	if(class(x)=="matrix") x<-as.vector(x)
	g<-qplot(x,geom="density",...)
	if(printGraph){
		print(g)
	}else{
		g
	}
}

qplotBarplot<-function(y,x,printGraph=TRUE,...){
	require(ggplot2)
	if(is.null(names(y))){
		aesX = 1:length(y)
	}else{
		aesX=names(y)
	}
	g<-ggplot(data.frame(x=aesX,y=y),aes(x=x,y=y))+
		geom_bar(stat="identity")
	if(printGraph){
		print(g)
	}else{
		g
	}
}

qplotAutoX<-function(x,printGraph=TRUE,geom="point",...){
	require(ggplot2)
	g<-qplot(x=1:length(x),y=x,geom=geom,...)
	if(printGraph){
		print(g)
	}else{
		g
	}
}

barplotPercentVar<-function(pca,printGraph=TRUE,...){
	g<-qplotBarplot(pca$percentVar*200,printGraph = FALSE)+ylab("% variance explained")+xlab("Principal component")+
		scale_x_continuous(breaks = seq(1,length(pca$percentVar),by = 2))+
		theme(
			panel.background = element_rect(fill = NA,colour="black"),
			panel.grid.major = element_line(colour = "grey"),
			panel.grid.minor = element_line(colour = NA)
		)
	if(printGraph){
		print(g)
	}else{
		g
	}
}

filledDoubleArrow<-function(x=0,y=0,width=1,height=1,just = c("left", "bottom"),gp=gpar(col="black"),...){
	pushViewport(viewport(x=x,y=y,width = width, height = height,just=just,...))
	grid.polygon(x = c(0,.2,.2,.8,.8,1,.8,.8,.2,.2,0),y=c(.5,.7,.6,.6,.7,.5,.3,.4,.4,.3,.5),gp = gp)
	popViewport()
}
