\name{Modules}
\docType{data}
\alias{Modules}
\title{Example Module Dataset}
\description{
  This example dataset and sample code illustrate the use of \code{NetWeaver} to visualize the complex features of gene coexpression network modules.
}
\usage{Modules}
\details{
This dataset contains a data.frame called \code{Modules}. Each row is a module, with the module id in the first column. The second column is ranking score. The next 4 columns are coefficients of module-trait correlations. The rest columns are P values of enrichment for various gene signatures. The sample R code shows a pipeline for visualizing the module feature data using a circos style plot.
}
\author{
Minghui Wang <m.h.wang@live.com>, Bin Zhang <bin.zhang@mssm.edu>
}
\keyword{datasets}
\seealso{\code{\link{rc.initialize}}}
\examples{
library(NetWeaver)
######################################
#set color functions
colfuncBlue=colorRampPalette(c("white", "blue"))
colfuncBrown=colorRampPalette(c("white", "brown"))
colfuncHeat=function(n) rev(heat.colors(n))
nCol=50 #number of colors in a gradient
######################################
#load module data
data(Modules)
#Create fake cyto information from modules.
#Here we treat each module as a hypothetical chromsome of equal size 100.
cyto=data.frame(Chr=Modules$id, Start=1, End=100, Data=1, BandColor='black',
		Score=Modules$Score,stringsAsFactors=FALSE)

###################################################################
# initialize circos
###################################################################
#set number of tracks
#should be larger than the actual number of rings required
#in order to leave sufficient sapce in the middle
num.tracks=36
rc.initialize(cyto,num.tracks=num.tracks,params=list(chr.padding=0.1))
params=rc.get.params()

######################################
#make plot area
######################################
rc.plot.area(0.9)
######################################
#plot ideogram
track.num=1:2
chrom.alias=1:nrow(cyto)
names(chrom.alias)=cyto$Chr
rc.plot.ideogram(track.num, plot.band=FALSE, plot.chromosome.id=TRUE,
	cex.text=1.0, chrom.alias=chrom.alias, track.border=NA)
######################################
#plot module ranking
Rank=cyto
track.num=4
data.col='Score'
params$color.hist <- '#A349A4' #'green' #'#A349A4'
rc.reset.params(params)
rc.plot.histogram(cyto, track.num, data.col, fixed.height=FALSE,
	custom.track.height=params$track.height*3, track.border=NA)
#############################
#plot module-trait correlation
eigenCorCols=colnames(Modules)[grep('Rho',colnames(Modules))]
#set maximum correlation coefficient
maxCorr=1
track.border='#999999'
track.color='white'
for(eigenCorCol in eigenCorCols){
	track.num <- track.num+1
	data.col <- 4
	color.col <- 5
	HistData=data.frame(Chr=Modules$id,Start=1,End=100,Data=Modules[,eigenCorCol],
		stringsAsFactors=FALSE)
	#plot positive correlation
	HistData1=HistData[HistData$Data > 0,]
	HistData1$col=colfuncBrown(nCol)[pmax(1,floor(HistData1[,data.col]*nCol/maxCorr))]
	rc.plot.histogram(HistData1, track.num, data.col, color.col=color.col,
		fixed.height=TRUE, track.color=track.color, track.border=track.border)
	#plot negative correlation
	HistData2=HistData[HistData$Data <= 0,]
	HistData2$Data=abs(HistData2$Data)
	HistData2$col=colfuncBlue(nCol)[pmax(1,floor(HistData2[,data.col]*nCol/maxCorr))]
	rc.plot.histogram(HistData2, track.num, data.col, color.col=color.col, fixed.height=TRUE,
		track.border=track.border)
}
#add color gradient legend at specified location
y.cor=rc.get.coordinates(1,1,1)$y[1]-1
x.cor=params$radius*0.8
bht=0.6
bwt=0.2
rc.plot.grColLegend(x.cor, y.cor, colfuncBrown(nCol), at=c(1,floor(nCol/2),nCol),
	legend=signif(c(0,maxCorr/2,maxCorr),2), title=expression(italic(r)),
	width=bwt,height=bht, cex.text=0.8)
rc.plot.grColLegend(x.cor, y.cor-bht, rev(colfuncBlue(nCol)), at=c(0,floor(nCol/2)),
	legend=signif(c(-maxCorr,-maxCorr/2),2), title='', width=bwt, height=bht, cex.text=0.8)
##########################
#plot P value significance of enrichment for gene signatures
#here we use rc.plot.heatmap function to plot multiple tracks of heat-map
track.num=track.num+1
heatmapData=t(Modules[,grep('Enrichment.*.Pvalue',colnames(Modules))])
colnames(heatmapData)=Modules$id
#convert P value to -log10 scale
heatmapData[,]=as.integer(-log(heatmapData+1.0e-320,base=10))
#cap the maximum -log10 P value so that there is good contrast
maxLogPval=25
heatmapData[heatmapData>maxLogPval]=maxLogPval
#
rc.plot.heatmap(heatmapData, track.num+1, color.gradient=colfuncHeat(nCol),
		track.color=track.color, track.border=track.border)
track.num=track.num+nrow(heatmapData)
#add color gradient legend
rc.plot.grColLegend(x.cor+0.8, y.cor-bht, colfuncHeat(nCol), at=c(1,floor(nCol/2),nCol),
	legend=signif(c(0,maxLogPval/2,maxLogPval),2), title=expression(paste(-log[10],'(P)')),
	width=bwt, height=2*bht, cex.text=0.8)

######################
#label track id for ease of reading
rc.plot.track.id(4, labels=1, col='black', custom.track.height=params$track.height*2, cex=0.8)
rc.plot.track.id(seq(7, track.num-1,by=3), labels=seq(7,track.num-1,by=3)-3, col='black', cex=0.8)

}
