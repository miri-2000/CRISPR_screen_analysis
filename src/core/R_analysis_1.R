# R_analysis_1: Script that normalizes read counts and creates quality control plots
# Last modified 19.11.2023
# ----------------------------------------------------------------------------------
# O) Create plots using the raw data (optional)
# 1) Normalize the data and store them in a new file
# 2) Create new plots using the normalized data
# --------------------------------------------------

# Clear all objects from work space,as well as the console
rm(list=ls())
cat("\14")

library(dplyr)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(DESeq2)
library(gplots)
library("RColorBrewer")
library(RUnit)
library(reshape2)
library(ggplot2)

analysis_tool_path <- file.path(dirname(whereami::thisfile()), "analysis_tools.R")

source(analysis_tool_path)

# Function storing default colors for plots for a maximum of 19 colors
getDefaultColors <- function(n,warning=T) {
  # the order in the palette determines which colors are given for a specific number of color
  
  if (warning & n > 19)
    stop("Maximum 19 colors")
  cn <- pal_npg()(9) 
  ca <- pal_aaas("default")(8) # leave out the colors grey and black
  
  # remove from AAAS the colors dark red and dark blue
  rm <- c("#BB0021FF","#3B4992FF","#EE0000FF")
  ca <- ca[!(ca %in% rm)]
  
  # some colorblind friendly colors
  cb <- c("#E69F00","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  cs <- c(cn,ca,cb)
  cs <- cs[c(1:n)]
  
  return(cs)
}  

# Function storing default colors for plots
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Function that creates a bar plot with the counts per condition
createCountBarplot <- function(x,metric,Conditions,colValues=NULL,title="Counts per sample",doColor=T) { 

  # Choose the x and y values for the chosen metric
  if (metric=="t") {
    y <- colSums(x,na.rm = T)
    ylab <- "Total counts" 
  }else if (metric=="ma") {
    x <- as.matrix(x) # other wise cannot calculate colMedians
    y <- colMedians(x,na.rm = T)
    ylab <- "Median all counts"
  }else if (metric=="mna") {
    x <- as.matrix(x)
    y <- colMeans(x,na.rm = T)
    ylab <- "Mean all counts"
  }else if (metric=="mn") {
    x <- x[anno$type=="n",]
    x <- as.matrix(x) # other wise cannot calculate colMedians
    y <- colMedians(x)
    ylab <- "Median neg. control counts"    
  }else if (metric=="mnt") {
    x <- x[anno$type=="o",]
    x <- as.matrix(x) # other wise cannot calculate colMedians
    y <- colMedians(x)
    ylab <- "Median NT control counts"   
  }
  
  # Create a new data frame with the x- and y-values
  df <- cbind.data.frame(x=colnames(x),y=y,cond=Conditions)
  df$x <- factor(df$x,levels=df$x)
  # Create the plot
  plot <- ggplot(df,aes(x=x,y=y,fill=Conditions))
  plot <- plot + scale_fill_manual(values = colValues,breaks=names(colValues))

  plot <- plot + 
    geom_bar(stat="identity") +
    xlab("Sample ID")+
    ylab(ylab)+
    theme_bw() +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text= element_text(family="Courier New",face="bold",hjust=0,vjust=0),
          axis.text.x =element_text(angle=-90))
  
  return(plot)
}

# Plot a scatter plot showing the distribution of the read counts
createRankPlot <- function(data,output,Conditions,colValues,showMedian=F,threshold=2,roundDigits=0,doColor=T,textSize=14,Nrrows=5,Nrcols=6,plotTitle="Counts per sgRNA",plotformats=c("jpg"),log10Transform=T) {
  if (log10Transform){
    data <- log10(data)
    dataformat <- "log10Count"
  } else {
    dataformat <- "raw_count"
  }
  # Set the upper limit of the y-axis (store the maximum value rounded up to an integer)
  max <- ceiling(max(data,na.rm=T))
  
  # Create a variable used to store all plots for each class
  plots <- list()
  
  # Iterate over each column/class
  for (i in c(1:ncol(data))) {
    plot_data=data.frame(dataformat=data[,i])
    # Store the current class for assigning the colors
    if (doColor) {
      condition1 <- Conditions[i]
      plot_data <- cbind(plot_data,condition1=condition1)
    }
    
    plot_data <- plot_data[order(plot_data$dataformat),,drop=F]
    # Add a column indicating the position in ordered list
    plot_data <- cbind(plot_data,rank=c(1:nrow(plot_data))/nrow(plot_data))
    
    # Define the plot title
    if (showMedian) {
      median1 <- round(median(plot_data$dataformat,na.rm=T),2)
      title <- paste(colnames(data)[i],", median:",median1,sep="")
    } else {
      # Calculate the percentage of rows that have a value above a set threshold (here = 2)
      perc <- round(nrow(plot_data[plot_data$dataformat > threshold,])/nrow(data) * 100,digits = roundDigits)
      title <- paste(colnames(data)[i]," (>",signif(threshold,3),": ",perc,"%)",sep="")
    }
    
    # Omit all rows that are infinite or NaN
    plot_data <- plot_data[is.finite(plot_data$dataformat) & !is.nan(plot_data$dataformat), ]
    
    # Define the x- and y-axis and visual parameters of the plot
    plot <- ggplot(data=plot_data,aes(x=rank,y=dataformat)) + ggtitle(title) + scale_y_continuous(limits=c(-1,max),breaks=c(-1:max)) + ylab(dataformat)
    plot <- plot + geom_point(aes(color=condition1)) +
      scale_color_manual(values = colValues,breaks=names(colValues))

    
    # Add a line for the median when yMetric="m"
    if (showMedian)
      plot <- plot + geom_hline(yintercept = median1, linetype="dashed")
    
    plot=plot+theme(axis.text=element_text(size=textSize,face="bold"), # ticks
              axis.title=element_text(size=textSize,face="bold",vjust=0.3), # the labels
              plot.title=element_text(face="bold", size=textSize + 2,vjust=1), # title
              legend.title=element_text(face="bold", size=textSize),
              legend.text=element_text(size=textSize),
              panel.grid.major = element_line(colour = "grey90",linewidth = 0.2), 
              panel.grid.minor = element_line(colour = "grey98", linewidth = 0.5),
              panel.background = element_rect(fill ="white", colour = "grey")) 
    
    # Store the current plot in a list
    plots[[i]] <- plot
  }  
  
  # Arrange all plots in one figure
  summaryplot <- ggarrange(plotlist=plots, nrow=Nrrows, ncol=Nrcols, common.legend = T, legend="right") + bgcolor("White") + labs(title = plotTitle)
  
  
  # Store the plot in the specified format
  for (plotformat in plotformats){
    if (plotformat=="pdf"){
      width <- Nrcols*30
      height <- Nrrows*30
      ggsave(filename = paste(output,".pdf",sep=""),plot=summaryplot,width=width,height=height,limitsize = F,,bg="white",device = cairo_pdf)
    }else { # (plotformat=="tiff")
      width <- Nrcols*15
      height <- Nrrows*15
      ggsave(filename = paste(output,".tiff",sep=""), plot=summaryplot,width=600, height=500, units="mm", dpi=300, compression = "lzw")
    }
  }
}

# Normalize read counts
normalizeData <- function(reads,Conditions,sizefactor_file,normalized_file,sizeFactorsBasedOnAll=T){
  
  data <- reads %>%
    select_if(is.numeric)
  annotations <- reads %>%
    select_if(~ !all(is.numeric(.)))
  rownames(data) <- reads$sgRNA

  # Store the replicate number
  replicate_nr <- sapply(strsplit(colnames(data),"_"),function(x)paste(x[length(x)],collapse="_"))
  
  # Create a simple data set object for normalization
  colData <- data.frame(Conditions=as.factor(Conditions), replicate_nr=as.factor(replicate_nr))
  
  # Create a DeSeq data set
  dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = colData,
                                design = ~ Conditions)
  
  # Calculate the size factor using DeSeq2
  deseq2SizeFactors <- sizeFactors(estimateSizeFactors(dds))
  deseq2Sf1 <- data.frame(sample=names(deseq2SizeFactors),deseq2=deseq2SizeFactors)

  # Calculate the size factor
  
  # Calculate the column sums
  total_reads <- apply(data,MARGIN=2,FUN=sum)
  # Transform to natural log values
  total_reads <- log(total_reads)
  # Calculate the mean
  mean_reads <- mean(total_reads)
  # Subtract the mean
  total_reads <- total_reads - mean_reads
  # Reverse the log
  totSizeFactors <- exp(total_reads)
  # Create a new data frame with the computed size factors
  totSf1 <- data.frame(sample=names(totSizeFactors),total=totSizeFactors)
  # Merge the data sets together
  sf1 <- merge(deseq2Sf1,totSf1,by="sample")
  # Calculate the ratio between the size factor and the new total
  sf1$ratioDT <- sf1$deseq2 / sf1$total
  sf1[,-1] <- signif(sf1[,-1],3)
  write.table(sf1,sizefactor_file,sep=";",row.names = FALSE)
  
  # Select the size factor method used based on the 
  # parameter sizeFactorsBasedOnAll
  if (sizeFactorsBasedOnAll){
    sizeFactors <- totSizeFactors
  }else{
    sizeFactors <- deseq2SizeFactors
  }
  
  # Calculate the normalized values
  if (!is.null(sizeFactors)) {
    norm <- t(t(data)/sizeFactors)
    norm <- data.frame(Tag=rownames(norm),norm,stringsAsFactors=F,check.names=F)
  }else
    stop("The size factor values could not be computed. Please verify your input file.")
  
  # Update the data variable and add a pseudo count to avoid log issues
  data <- round(norm[,-1]) + 1
  reads <- cbind(annotations,data,stringsAsFactors=F)
  # Store the normalized data in the specified file
  write.table(reads,file=normalized_file,row.names=F,sep=";")
  
  return(list(reads=reads,annotations=annotations,data=data))
}

# Create histograms with the raw and normalized data
# input_raw: list with raw input data
# input_norm: list with normalized input data
# writeFile: if True a  file is created, otherwise the two plots are returned
# filename: Name of the output file
# metrics: list of metrics for which a plot should be created
# doColor: boolean, defines if plot is in color or black and white (default = True)
# textSize: size of the axis labels
# colValues: can be used to define color values for the resulting plots (default = NULL)
# nextToEach: if True, raw and normalized plots (if inclNorm=T) are in a row, all metrics are in a column
# plotformats: List of file formats that should be used to create the output files
createRawNormBarplots <- function(input_raw,input_norm,Conditions,writeFile=T,filename="RawNormCount_Barplot",metrics=c("t"),doColor=T,textSize=14,colValues=NULL,nextToEach=T,plotformats=c("pdf")) {
  
  plots <- list()
  plotIdx <- 0 
  for (metric in metrics) { 
    # title per individual plot
    
    if (metric=="t"){
      metricTxt <- "Total"
    }else if (metric=="ma"){
      metricTxt <- "Median all"
    }else if (metric=="mna"){
      metricTxt <- "Mean all"
    }else if (metric=="mn"){
      metricTxt <- "Median neg. control"
    }else if (metric=="mnt"){
      metricTxt <- "Median NT control"
    }
    
    # In case only one metric is used, the metric is mentioned in the 
    # overall Title
    if (length(metrics)==1) {
      titleOverall <- paste(metricTxt,"counts per sample")
    } else {
      titleOverall <- "Counts per sample"
    }  
    
    for (Plot in c(1,2)) {
      if (Plot==1) {
        x <- input_raw
        titleInd <- "Raw"
      }else {
        x <- input_norm
        titleInd <- "Normalized"
      }
      
      # In case more than one metric is used, adjust the text of the title to
      # current current metric + raw/normalized data
      if (length(metrics) > 1)  # prefix with metric
        titleInd <- paste(metricTxt, ", ",titleInd,sep="")
      
      # Create the bar plot with the raw counts
      plot <- createCountBarplot(x,metric=metric,Conditions=Conditions,title=titleInd,doColor=doColor,colValues = colValues)
      
      # Modify the y-limits for the best fit
      if (Plot==1) {
        res <- ggplot_build(plot)
        yLimits1 <- res$layout$panel_scales_y[[1]]$get_limits()
        # Add a buffer in case the normalized values are a bit higher as the max value in the raw data
        yLimits1[2] <- yLimits1[2] + 2
        yBreaks1 <-res$layout$panel_scales_y[[1]]$get_breaks()
        # Remove potential NA values
        yBreaks1 <- yBreaks1[!is.na(yBreaks1)]
      }
      
      plot <- plot + scale_y_continuous(limits=yLimits1,breaks=yBreaks1)
      plot <- plot + theme(axis.text=element_text(size=textSize,face="bold",family="mono"),
                     axis.title=element_text(size=textSize,face="bold",vjust=0.3),
                     plot.title=element_text(face="bold", size=textSize + 2,vjust=1),
                     legend.title=element_text(face="bold", size=textSize),
                     legend.text=element_text(size=12),
                     panel.grid.major = element_line(colour = "grey90",linewidth = 0.2), 
                     panel.grid.minor = element_line(colour = "grey98", linewidth = 0.5),
                     panel.background = element_rect(fill ="white", colour = "grey"))  
      
      plotIdx <- plotIdx + 1
      plots[[plotIdx]] <- plot
    }
  } 
  
  if (writeFile) {
    # Define the layout of the overall figure
    if (nextToEach) {
      Nrrows <- length(metrics) #raw or norm
      Nrcols <- 2 #raw and norm
    }else {
      Nrrows <- 2
      Nrcols <- length(metrics)
    }
    
    plot <- ggarrange(plotlist=plots, nrow=Nrrows, ncol=Nrcols, common.legend = T, legend="bottom")
    plot <- annotate_figure(plot,top = text_grob(titleOverall, face = "bold", size = 16))
                            
    for (plotformat in plotformats) {
      if (plotformat=="pdf"){
        width <- Nrcols*9
        height <- Nrrows*6.5
        ggsave(filename = paste0(filename,".pdf"),plot=plot,width=width,height=height,limitsize = F,,bg="white")
      }else { # (plotformat=="jpg"){
        width <- Nrcols*600
        height <- Nrrows*400
        jpeg(paste(filename,".jpg"),width=width,height=height) #in px
        parameters <- c(plots,nrow=Nrrows,ncol=Nrcols)
        do.call(grid.arrange,parameters)
        dev.off()
        }
      }
    }else {
    return(plots)
  }  
}

#Create density plots
createCountDensity <- function(data,Conditions,dataset,plotformats=c("jpg")){
  plots  <- list()
  plotIdx <- 0
  Nrcols <- length(Conditions)/NrConditions
  Nrrows <- NrConditions
  
  for (i in c(1:length(data))) {
    data1 <- data[,i,drop=F]
    colnames(data1) <- "value"
    plot <- ggplot() + geom_density(data=data1, aes(x=value)) +
      ggtitle(colnames(data)[i])
    
    limits <- c(0,max(data1))
    plot <- plot + scale_x_continuous(limits=limits)
    
    # Add the median to the plot
    m <- median(data1$value,na.rm=T)
    plot <- plot + geom_vline(xintercept = m,linetype="dashed")
    
    plotIdx <- plotIdx + 1
    plots[[plotIdx]] <- plot
  }
  
  for (plotformat in plotformats){
    if (plotformat=="pdf"){
      pdf(paste(filename,".pdf",sep=""),width=Nrcols*3,height=Nrrows*3) #in inch
    }else if (plotformat=="jpg"){
      jpeg(paste(filename,".jpg",sep=""),width=Nrcols*300,height=Nrrows*300) #in px
      parameters <- c(plots,nrow=Nrrows,ncol=Nrcols)
      do.call(grid.arrange,parameters)
      dev.off()}
  }
}

# Plot normal approximation figures
plotLog2ratio <- function(reads,filename="norm-approx",conditions1,conditions2,Nrrows,Nrcols,textSize=12){ #normalized data
  
  # Specify the annotation columns to use for the log2 data- frame
  logdf <- reads[,c("sgRNA","Gene","type")]
  
  plots <- list()
  pIdx <- 0
  if (length(conditions1)!=length(conditions2))
    stop("The lists containing the conditions to compare must be of equal length!")
  for (i in c(1:length(conditions1))) {
    condition_name1 <- conditions1[i]
    condition_name2 <- conditions2[i]
    # Select all replicates with the same condition (conditions1 and conditions2)
    condition1 <- reads[,grep(condition_name1,colnames(reads))]
    condition2 <- reads[,grep(condition_name2,colnames(reads))]  
    # Calculate the mean for each row (conditions1 and conditions2)
    condition_means1 <- apply(condition1,MARGIN=1,FUN=mean)
    condition_means2 <- apply(condition2,MARGIN=1,FUN=mean)
    # Calculate the log2 between the ratio of the means of conditions1 and conditions2
    logdf$log2ratio <- log2(condition_means1/condition_means2)
    
    # Calculate the mean and standard deviation of the log2 ratio
    m <-  mean(logdf$log2ratio,na.rm=T)
    s <- sd(logdf$log2ratio,na.rm=T)
    
    title <- paste(condition_name1,"vs",condition_name2)
    pIdx <- pIdx + 1
    # Create a plot showing the empirical data and the normal approximation
    plots[[pIdx]] <- ggplot(data=logdf,aes(x=log2ratio)) + geom_density(aes(linetype="e")) + # scale_x_continuous(limits = c(0,1.4))
      stat_function(fun = dnorm, aes(linetype = "n"), args = list(mean = m, sd = s)) +
      theme(axis.text=element_text(size=textSize,face="bold"), # ticks
            axis.title=element_text(size=textSize,face="bold",vjust=0.3), # the labels
            plot.title=element_text(face="bold", size=textSize + 2,vjust=1), # title
            legend.title=element_text(face="bold", size=textSize),
            legend.text=element_text(size=textSize),
            panel.grid.major = element_line(colour = "grey90",linewidth = 0.2),
            panel.grid.minor = element_line(colour = "grey98", linewidth = 0.5),
            panel.background = element_rect(fill ="white", colour = "grey"))  +
      scale_linetype_manual("Type", values = c("e" ="solid","n" = "dashed"),labels=c("empirical","normal\napproximation")) +
      ylab("probability") + ggtitle(title) # +geom_vline(xintercept=mark5,linetype="dotted")
  }
  

  plot <- ggarrange(plotlist=plots, nrow=Nrrows, ncol=Nrcols, common.legend = T, legend="top") 
  plot <- annotate_figure(plot,
                          top = text_grob("Distribution log2 fold change values", face = "bold", size = 14))
                          
  ggsave(plot=plot,file=paste(filename,".pdf",sep=""),width=Nrcols*3,height=Nrrows*3)
}



args <- commandArgs(trailingOnly = TRUE)

reads <- read.csv(as.character(args[1]),stringsAsFactors=F,check.names=F,sep=";")

dataset=as.character(args[2])
output_file <- paste(dataset,"data_prep_norm")
Conditions <- strsplit(args[3], ',')[[1]]
sizefactor_file <- paste(dataset,"_sizefactors.csv",sep="")
normalized_file <- paste(dataset,"_norm.csv",sep="")
levels <- as.character(args[4])
combined_replicates_file <- NULL
# specify which conditions should be compared
conditions1 <- strsplit(args[5], ',')[[1]]
conditions2 <- strsplit(args[6], ',')[[1]]
distribution_condition1 <- as.character(args[7])
distribution_condition2 <- as.character(args[8])

data <- reads %>%
  select_if(is.numeric)
annotations <- reads %>%
  select_if(~ !all(is.numeric(.)))
unique_conditions <- unique(Conditions)
NrConditions <- length(unique_conditions)

# Default parameter settings
doColor=T 
colValues=NULL
metric="t"
plotRawRankPlot=F 
plotNormRankPlot=T
plotRawNormBarplots=T
plotCountsDensity <- F
createQualityPlots=T
plotDistr=T
corrPlotPerRep=F
doCorrPlotAll=F

# Set the color values for the data points
if (NrConditions <= 19) {
  colValues <- getDefaultColors(NrConditions)
}else {
  colValues <- gg_color_hue(NrConditions)
}
names(colValues) <- unique_conditions


# Create plots with the raw data
if (plotRawRankPlot)
  createRankPlot(data,Conditions,colValues,doColor=T,output=paste(dataset,"_RawRankPlot",sep=""))

# Normalize the data
rawdata=data
updates=normalizeData(reads,Conditions,sizefactor_file,normalized_file)
reads <- updates$reads
annotations <- updates$annotations
data <- updates$data

# Create plots corresponding to the type of replicates
for (level in levels) {
  if (level=="technical") {
    prefix <- "tech"
  }else { # level=="biological"
    prefix <- "bio"

    if (!is.null(combined_replicates_file)) {
      # Calculate the geometric mean out of the technical replicates
      data <- as.data.frame(sapply(unique_conditions,function(x)rowMeans(data[,which(Conditions==x),drop=F])))
      # Round the values(DeSeq2 requires integers)
      data <- round(data)
      # Update the reads
      reads <- cbind(annotations,data,stringsAsFactors=F)
      write.table(target_genes,file="_norm_bio.csv",row.names=F,sep=";")
      # Update the raw/un-normalized reads as well
      rawdata <- as.data.frame(sapply(unique_conditions,function(x)rowMeans(rawdata[,which(Conditions==x),drop=F])))
    }
  }

  # Create bar plots comparing the raw and normalized read counts
  if (plotRawNormBarplots){
    createRawNormBarplots(input_raw=rawdata,input_norm=data,
                          Conditions=Conditions,metric="t",
                          colValues=colValues,filename=paste(
                            dataset,prefix,"RawNormCount_Barplot",sep="_"))
  }
  # Create a new histogram with the normalized data
  if (plotNormRankPlot) {
    if (length(conditions1)>6){
      Nrrows=ceiling(length(conditions1)/6)
      Nrcols=6
    } else{
      Nrrows=1
      Nrcols=length(conditions1)
    }
    createRankPlot(data,output=paste(dataset,prefix,"NormRankPlot",sep="_"),Conditions,colValues)
  }

  # Create a density plot of the data per condition (eg t1_r1)
  if (plotCountsDensity) {
    createCountDensity(data,Conditions,filename=paste(dataset,prefix,"counts_sgrna_distr",sep="_"))
  }

  # Create quality control plots
  if (createQualityPlots) {
    if (length(Conditions)>6){
      Nrrows=ceiling(length(Conditions)/6)
      Nrcols=6
    } else{
      Nrrows=1
      Nrcols=length(Conditions)
    }

    qualityControlPlots(reads,length(annotations),Conditions,basename=paste(dataset,prefix,"norm",sep="_"),sizeFactors=sizeFactors,typeColoring=T,plotformat="pdf",
                        correlation=T,heatmap=T,paired=F,nrow=Nrrows,ncol=Nrcols,checkReplicates=F,ctrlDistrNRow=1,ctrlDistrNCol=1,ctrlDistrFrom=-6,ctrlDistrTo=3,
                        ctrlWidthSizeFactor=5)

  }

  # Create a distribution plot to compare different time points
  if (plotDistr) {
    read_copy=reads
    # Define the values for the distribution plot between 
    # positive and negative controls
    exps <- distribution_condition1
    refs <- distribution_condition2
    Nrrows <- length(exps)
    Nrcols <- length(refs)

    oFileName <- paste(dataset,prefix,"norm","distibution",sep="_")
    res <- controlDistrs(x=read_copy,exps=exps,refs=refs,
                         annotation_columns=length(annotations),paired=F,
                         other=NULL,yAxisMax=NULL,nCol=Nrcols,nRow=Nrrows,
                         oFileName=oFileName,plotformat="pdf",
                         types = c("n","p"),drawMedian=T,writeToFile = T,
                         widthSizeFactor = 5)
  }

  # Create a correlation plot per replicate over the condition
  if (corrPlotPerRep) {
    data_copy=data
    # Calculate the log10
    data_copy[,-c(1:annotation_columns)] <- log10(data_copy[,-c(1:annotation_columns)])
    # Store the time points from each column
    timepoints <- substring(colnames(data_copy)[-c(1:annotation_columns)],1,2)
    # Create the correlation plots using the column "type" along with all data columns
    plots <- correlationPlots(x=data_copy[,-c(1:2,4)],Conditions=timepoints,outputFile=paste(dataset,prefix,"normlog10_correlation_correlations_rep",sep="_"),nrow=NULL,ncol=NULL,typeColoring=T,plotformat="pdf",units="log10Counts")
  }

  # Create pairwise correlation plots irrespective of the time point they were made
  if (doCorrPlotAll)  {
    outputFile <- paste(dataset,prefix,"normlog10_correlation_plots_all",sep="_")
    corrPlotAll(input=target_genes,annotation_columns=annotation_columns,outputFile=outputFile)

  }
}

# Remove all rows that are of type "o"
reads <- reads[reads$type != "o",]
print(length( reads[reads$type != "o",]))

# Calculate log2 ratios
if (length(conditions1)>6){
  Nrrows=ceiling(length(conditions1)/6)
  Nrcols=6
} else{
  Nrrows=1
  Nrcols=length(conditions1)
}
plotLog2ratio(reads,filename=paste(dataset,"norm-approx",sep="_"),conditions1,
              conditions2,Nrrows=Nrrows,Nrcols=Nrcols)

# Remove unwanted annotation columns and store the final table in a new file
reads <- reads[, !(colnames(reads) %in% c("type", "sequence"))]
write.table(reads,paste(dataset,"drugz-input.txt",sep="_"),row.names=F,sep="\t")




























  

  


