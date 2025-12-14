# R_analysis_2: Script that creates gene- and guide-level significance plots
# Last modified 19.11.2023
# --------------------------------------------------------------------------

# Clear all objects from work space,as well as the console
rm(list=ls())
cat("\14")

library(ggplot2)
library(ggrepel)
library(gridExtra)
library(scales)
library(ggpubr)

# Create the reverse log for the y-axis scale
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# Add a specific theme to the plot
addTheme <- function(plot,base) {
  plot <- plot + theme(axis.text=element_text(size=base,face="bold"), # ticks
                           axis.title=element_text(size=base,face="bold",vjust=0.3), # the labels
                           plot.title=element_text(face="bold", size=base + 2,vjust=1,hjust=0.5), # title
                           legend.title=element_text(face="bold", size=base),
                           legend.text=element_text(size=base),
                           panel.grid.major = element_line(colour = "grey90",linewidth = 0.2), 
                           panel.grid.minor = element_line(colour = "grey98", linewidth = 0.5),
                           panel.background = element_rect(fill ="white", colour = "grey")) 
}

# Create guide-level significance plots
create_gRNA_plots <- function(significant_genes,gRNA_comparison,sel,title,thFdr=-log10(0.05)) {
  # create a directory for each comparison to store store the gRNA plots in
  output_dir=paste(".//drugz//",gRNA_comparison,sep="")
  if (!file.exists(output_dir)) {
    # If it doesn't exist, create it
    dir.create(output_dir)
  } 
  gRNA_file= read.csv(paste("drugz//gRNA_",gRNA_comparison,".txt",sep=""),stringsAsFactors = F,sep="\t")
  
  gRNA_normZ=gRNA_file[gRNA_file$GENE %in% significant_genes$Gene,]
  for (gene in unique(gRNA_normZ$GENE)){
    gRNAs=gRNA_normZ[gRNA_normZ$GENE==gene,]
    xLabel <- "z-score"
    show.legend <- F
    xlimit=ceiling(max(abs(gRNA_file$zscore)))
    xlimits=c(-(xlimit),xlimit)
    xbreaks=c(-(xlimit):xlimit)
    
    # Calculate the p-value
    gRNA_file$p_value <- 2 * (1 - pnorm(abs(gRNA_file$zscore)))
    
    # Create a negative log10 p-value column
    gRNA_file$log10_p_value <- -log10(gRNA_file$p_value)
    gRNAs=gRNA_file[gRNA_file$GENE==gene,]
    
    ylimits=c(0,ceiling(max(gRNA_file$log10_p_value)))
    print(xlimits)
    plot <- ggplot(data=gRNA_file,aes(x=zscore,y=log10_p_value)) +
      geom_point(size=2,color="grey90") +
      scale_y_continuous() + #, ,,labels=yBreaks ,limits=yLimits,breaks=yBreaks log10_reverse_trans()
      scale_x_continuous(limits=xlimits,breaks=xbreaks) +
      geom_point(data=gRNAs,aes(x=zscore,y=log10_p_value,color=GENE),size=2,show.legend = show.legend) + #,alpha=0.4 
      geom_text_repel(data=gRNAs,aes(x=zscore,y=log10_p_value,label=X,color=GENE),min.segment.length=0,show.legend = FALSE,nudge_y = 0.05,max.overlaps = 50) +
      ggtitle(paste(title,": drugZ results",sep="")) +
      #theme(plot.title = element_text(hjust = 0.5)) + # + 
      geom_hline(yintercept = thFdr,linetype="dashed") +
      geom_vline(xintercept = 0,linetype="solid") +
      xlab(xLabel) +
      ylab("-log10(p-value)")
    
    plot <- addTheme(plot,base=10)
    
    if (sel=="neg"){
      selection="negative_selection_"
    }else{
      selection="positive_selection_"
    }
    ggsave(plot=plot,filename=paste(output_dir,"//",gRNA_comparison,"_gRNA_",selection,gene,"_plot.","pdf",sep=""),limitsize = F)
    
  }
}

# Create gene-level significance plots
create_plots <- function(data,target_samples,reference_samples,x.axis="normZ",threshold_fdr=0.25,top=15,inclXyplot=F,selections=c("pos","neg"),plotformat="pdf"){
  for (i in c(1:length(target_samples))) {
    target_sample=target_samples[i]
    reference_sample=reference_samples[i]
    
    target_name <- sub("_[tr|ut]+$", "",target_sample)
    reference_name <- sub("_[tr|ut]+$", "",reference_sample)
    # store the column name of interest
    comparison <- paste(target_sample,reference_sample,sep= ".")
    plots <- list()
    pIdx <- 0
    
    columns <- c("Gene")
    
    inclXyplot=F
    if (target_name==reference_name) {
      inclXyplot <- T
      target_sample=paste(target_name,"_tr",sep="")
      reference_sample=paste(target_name,"_ut",sep="")
      if (target_name=="t2_2d"){
        xCol2= paste(reference_sample,".t1_2d.Gene.log2fc",sep="") #tr
        yCol2= paste(target_sample,".t1_2d.Gene.log2fc",sep="") #ut
      }else { # target_name=="t2_3d"
        xCol2= paste(reference_sample,".t1_3d.Gene.log2fc",sep="")
        yCol2= paste(target_sample,".t1_3d.Gene.log2fc",sep="")
      }
      columns=c(columns,xCol2,yCol2)
    }else{
      if (length(x.axis)>1){
        stop("Can only do normZ and log2FoldChange one on one if log2FoldChange is a value calculated by drugZ. Otherwise there can be inconsistent results")
      }
    }
    
    
    if ("log2 fold-change" %in% x.axis){
      comparison.log2fc <- paste(comparison,"Gene.log2fc",sep= ".")
      columns=c(columns,comparison.log2fc)
    }
    if ("normZ" %in% x.axis){
      normZCol <- paste(comparison,"normZ",sep= ".")
      columns=c(columns,normZCol)
    }
    
    for (selection in selections) {
      columns_per_selection=columns
      print(columns_per_selection)
      print(setdiff(columns_per_selection, colnames(data)))
      threshold_fdr <- 0.25
      
      if (selection=="pos") {
        thLfc <- 1
        thDiff <- 1
        
        fdrCol <- paste(comparison,"fdr_supp",sep= ".")
        rankCol <- paste(comparison,"rank_supp",sep= ".")
        selection_name="pos. sel."
      }else{
        thLfc <- -1
        thDiff <- -1
        
        fdrCol <- paste(comparison,"fdr_synth",sep= ".")
        rankCol <- paste(comparison,"rank_synth",sep= ".")
        selection_name="neg. sel."
      }
      columns_per_selection=c(columns_per_selection,fdrCol,rankCol)
      
      
      # Select the columns of interest
      data1 <- data[columns_per_selection]
      data1 <- data1[order(data1[,rankCol],decreasing=F),]
      
      # Select rows with a fdr beneath the fdr threshold
      significant_hits <- data1[data1[, fdrCol] <= threshold_fdr, ]
      if (nrow(significant_hits) > top) {
        takenTop <- T
        significant_hits <- significant_hits[c(1:top),]
      }else {
        takenTop <- F
      }
      
      significant_hits$Gene <- factor(significant_hits$Gene,levels=significant_hits$Gene)
      
      title <- paste(target_sample," vs ",reference_sample,", ",selection_name,sep="")
      if (takenTop)
        title <- paste(title,"Top",top)
      
      for (x.ax in x.axis) { 
        # VOLCANO PLOT DRUGZ
        # Set the plot details
        if (x.ax=="normZ") {
          xCol <- normZCol
          xLabel <- "normZ"
          show.legend <- F
        }else { # log2 fold change
          xCol <- comparison.log2fc
          xLabel <- xCol 
          show.legend <- T
        } 
        
        if (selection=="neg") {
          xMin <- floor(min(data[,xCol]))
          xLimits <-  c(xMin,0) 
          xBreaks <- c(xMin:0)
        }else { # pos. selection
          xMax <- ceiling(max(data[,xCol]))
          xLimits <-  c(0,xMax) 
          xBreaks <- c(0:xMax)
        }
        
        yLimits <- c(1,0.1)
        yBreaks <- c(10:1)/10
        
        pIdx <- pIdx + 1
        plots[[pIdx]] <- local({
          data1=data1
          significant_hits=significant_hits
          xCol=xCol
          fdrCol=fdrCol
          xLimits=xLimits
          xBreaks=xBreaks
          
          plot <- ggplot(data=data1,aes(x=data1[,xCol],y=data1[,fdrCol])) +
            geom_point(size=2,color="grey90") +
	    scale_y_continuous(
              trans = reverselog_trans(10),
              breaks = log_breaks(base = 10),  # Automatically generate breaks based on log scale
              #labels = scales::trans_format("log10", scales::math_format(10^.x))  # Optional: format labels
            ) +
            #scale_y_continuous(trans=reverselog_trans(10)) +
            scale_x_continuous(limits=xLimits,breaks=xBreaks) +
            geom_point(data=significant_hits,aes(x=significant_hits[,xCol],y=significant_hits[,fdrCol],color=Gene),size=2,show.legend = show.legend) + #,alpha=0.4 
            geom_text_repel(data=significant_hits,aes(x=significant_hits[,xCol],y=significant_hits[,fdrCol],label=Gene,color=Gene),min.segment.length=0,show.legend = FALSE,nudge_x = 0.1,nudge_y = 0.1, max.overlaps = 50,force = 2,box.padding = 0.5,point.padding = 0.5,) +
            ggtitle(paste(title,": drugZ results",sep="")) +
            geom_hline(yintercept = threshold_fdr,linetype="dashed") +
            geom_vline(xintercept = 0,linetype="solid") +
            xlab(xLabel) +
            ylab("FDR")
          plot <- addTheme(plot,base=10)
        })
        
        if (nrow(significant_hits)!=0){
          gRNA_comparison= paste(target_sample,reference_sample,sep="-")
          create_gRNA_plots(significant_hits,gRNA_comparison,selection,title)
        }
      }  
      
      if (inclXyplot) {
        #XY PLOT
        pIdx <- pIdx + 1
        plots[[pIdx]] <- local({
          data=data
          significant_hits=significant_hits
          
          plot <- ggplot()+
            ggtitle(paste(title,": log2fc vs t0",sep="")) +
            #geom_ribbon(aes(x, ymin=ymin, ymax=y), fill="#e6ffe6") +
            geom_point(data=data,aes(x=data[,xCol2],y=data[,yCol2]),size=0.3,color="grey90") +
            #geom_point(data=data1[data1$type=="neg",],aes(x=x,y=y),color="blue",size=0.3) +
            geom_point(data=significant_hits,aes(x=significant_hits[,xCol2],y=significant_hits[,yCol2],color=Gene),size=2) +
            geom_text_repel(data=significant_hits,mapping=aes(x=significant_hits[,xCol2],y=significant_hits[,yCol2],label=Gene,color=Gene),size=4,fontface="bold",show.legend = F,point.padding = 0.2,max.overlaps=50) +
            geom_vline(xintercept=0) +
            geom_hline(yintercept=0) +
            geom_abline(linetype="dotted")  +
            geom_hline(yintercept =thLfc,linetype="dashed") +
            geom_abline(intercept=thDiff,linetype="dashed") +
            xlab(xCol2) +
            ylab(yCol2)
          plot <- addTheme(plot,base=10)
        })
      }
    }
    
    
    ncol <- length(x.axis)
    if (inclXyplot==T){
      ncol=ncol+1
    }
    nrow <- length(selections)
    plot <- ggpubr::ggarrange(plotlist = plots, nrow =nrow, ncol = ncol,common.legend=F,legend="right") 
    widthBase <- 6
    heightBase <- 4
    
    for (type in plotformat){
      comparison_filename <- gsub("\\.","-",comparison)
      ggsave(plot=plot,filename =paste("drugz/drugz_",comparison_filename,"_plot.",type,sep=""),width=ncol*widthBase,height=nrow * heightBase,limitsize = F)
    }  
  }  
}



args <- commandArgs(trailingOnly = TRUE)

data <- read.csv(as.character(args[1]),stringsAsFactors = F,sep=";")
print(colnames(data))

create_plots(data=data,target_samples=strsplit(args[2], ',')[[1]],reference_samples=strsplit(args[3], ',')[[1]],x.axis=as.character(args[4]),threshold_fdr=as.numeric(args[5]),top=as.numeric(args[6]))
  

