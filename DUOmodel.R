### Adaptation of tag analysis package to work for duo libraries
### Hannah Dewey
### 22.11.2019

# Load in required packages

library(dplyr)
library(tibble)
library(DESeq2)
# library(preprocessCore)
library(reshape2)
library(GGally)


### Adds Enhancers to the attributes table
## INPUTS:
##    duoAttr:    Attributes table for duo data
##    enhList:    List of enhancers in project
## OUTPUTS:
##    attr_update: Updated attributes table with duo syntax
duoAttr <- function(duoAttr, enhList = c()){
  attr_update <- data.frame()
  for(enhancer in enhList){
    attr_temp <- duoAttr
    attr_temp$ID <- paste(enhancer,attr_temp$ID, sep = "^")
    attr_update <- rbind(attr_update, attr_temp)
  }
  return(attr_update)
}

### Transfoms barcode level output from count pipeline to oligo level counts and calculates the 
### plasmid mean and number of barcodes for each oligo.
## INPUTS:
##    dataCount:    barcode level counts from count pipeline
##    dataCond:     condition matrix for data, should have one column named condition with associated replicates as row names
## OUTPUTS:
##    counts_oligo: oligo level counts table with plasmid mean and number of barcodes for each oligo
duoStats <- function(dataCount,dataCond){
  plmd_counts <- dataCount[,c("Oligo",rownames(dataCond)[which(dataCond$condition=="DNA")])]
  message(class(plmd_counts[,colnames(plmd_counts)!="Oligo"]))
  BCSum <- as.data.frame(table(plmd_counts$Oligo[rowMeans(as.matrix(plmd_counts[,colnames(plmd_counts)!="Oligo"])) > 0]))
  message("set column names")
  colnames(BCSum) <- c("Oligo","PlasmidsBCsum")
  message("aggregating oligos")
  tag_counts <- aggregate(. ~Oligo+library, data=dataCount[,-1], FUN = sum)
  message("adding bc sums")
  tag_counts <- merge(tag_counts,BCSum, by="Oligo")
  counts_oligo <- tag_counts[,-1]
  rownames(counts_oligo) <- tag_counts[,1]
  message("calculating plasmid mean")
  message(paste0(colnames(counts_oligo[,rownames(dataCond)[dataCond$condition=="DNA"]])), collapse = "\t")
  plmean <- as.data.frame(rowMeans(counts_oligo[,rownames(dataCond)[dataCond$condition=="DNA"]]))
  colnames(plmean) <- "plmean"
  counts_oligo <- merge(counts_oligo,plmean,by="row.names")
  rownames(counts_oligo) <- counts_oligo$Row.names
  counts_oligo <- counts_oligo[,-1]
  counts_oligo$library <- as.factor(counts_oligo$library)
  y <- grepl(",",rownames(counts_oligo))
  counts_oligo <- counts_oligo[!y,]
  celltypes <- as.factor(dataCond$condition)
  for(celltype in levels(celltypes)){
    message(celltype)
    counts_oligo[,celltype] <- 0
    for(oligo in 1:nrow(counts_oligo)){
      for(rep in rownames(dataCond)[dataCond$condition==celltype]){
        if(as.integer(counts_oligo[oligo,rep]) == 0){
          counts_oligo[oligo,celltype] <- counts_oligo[oligo,celltype] + 1
        }
      }
    }
  }
  return(counts_oligo)
}

### Identifies the common names between each library for 2 runs, runs filtering on oligos based on
### plasmid mean and the number of barcodes in the plasmid
## INPUTS:
##    dataCount1: count data as a dataframe, for one run being analyzed
##    dataCount2: count data as a dataframe, for the second run being analyzed
##    libExcl1  : OPTIONAL. A list of libraries to be excluded from the first run
##    libExcl2  : OPTIONAL. A list of libraries to be excluded from the second run
## OUTPUTS:
##    names_list: List of oligos in both runs separated by library
duoNames <- function(dataCount1, dataCount2, libExcl1=c(), libExcl2=c(), duoOnly=F){
  for(lib1 in levels(dataCount1$library)){
    if(lib1 %in% libExcl1) next
    message(lib1)
    dataCount1 <- dataCount1[dataCount1$plmean >= 20, ]
    dataCount1 <- dataCount1[dataCount1$PlasmidsBCsum >= 10,]   
  }
  message("filtered 1")
  for(lib2 in levels(dataCount2$library)){
    if(lib2 %in% libExcl2) next
    message(lib2)
    dataCount2 <- dataCount2[dataCount2$plmean >= 20,] 
    dataCount2 <- dataCount2[dataCount2$PlasmidsBCsum >= 10,]
  }
  message("filtered 2")
  names_list <- list()
  for(lib in levels(dataCount1$library)){
    if(lib %in% libExcl1 | lib %in% libExcl2) next
    message(lib)
    names_list[[lib]] <- intersect(row.names(dataCount1)[dataCount1$library==lib], rownames(dataCount2)[dataCount2$library==lib])
  }
  if(duoOnly==F){
    message(length(names_list[[1]]))
    message(length(names_list[[2]]))
    duoname_a <- expand.grid(names_list[[1]], names_list[[2]])
    duoname_a <- mutate(duoname_a, duo1 = paste(!!!rlang::syms(c("Var2", "Var1")), sep="^"))
    duoname_a <- mutate(duoname_a, duo2 = paste(!!!rlang::syms(c("Var1", "Var2")), sep="^"))
    for(lib1 in libExcl1){
      for(lib2 in libExcl2){
      # if(duoname_a$duo1[1] %in% rownames(dataCount2[dataCount2$library==lib1,])){
        if(rownames(dataCount2[dataCount2$library==lib1,])[1] %in% duoname_a$duo1){
          message(paste0(lib, ": duo1"))
          assembled_a <- subset(duoname_a, duoname_a$duo1 %in% rownames(dataCount2)[dataCount2$library==lib1] & duoname_a$duo2 %in% rownames(dataCount1)[dataCount1$library==lib2])
        }
      # if(duoname_a$duo2[1] %in% rownames(dataCount2[dataCount2$library==lib1,])){
        if(rownames(dataCount2[dataCount2$library==lib1,])[1] %in% duoname_a$duo2){
          message(paste0(lib, ": duo2"))
          assembled_a <- subset(duoname_a, duoname_a$duo2 %in% rownames(dataCount2)[dataCount2$library==lib1] & duoname_a$duo1 %in% rownames(dataCount1)[dataCount1$library==lib2])
        }      
      }
    }
    if(length(libExcl1)==0&length(libExcl2)==0){
      assembled_a <- rbind(duoname_a$duo1,duoname_a$duo2)
    }
    message(paste0("total possible duos: ",2*nrow(duoname_a)))
    message(paste0("total duos found: ", nrow(assembled_a)))
    for(lib in libExcl1){
      if(assembled_a$duo1[1] %in% rownames(dataCount2[dataCount2$library==lib,])){
        message(paste0(lib, ": duo1"))
        names_list[[lib]] <- subset(rownames(dataCount2)[dataCount2$library==lib], rownames(dataCount2)[dataCount2$library==lib] %in% assembled_a$duo1)
      }
      else{
        message(paste0(lib, ": duo2"))
        names_list[[lib]] <- subset(rownames(dataCount2)[dataCount2$library==lib], rownames(dataCount2)[dataCount2$library==lib] %in% assembled_a$duo2)
      }
    }
    for(lib in libExcl2){
      if(assembled_a$duo1 [1] %in% rownames(dataCount1[dataCount1$library==lib,])){
        message(paste0(lib, ": duo1"))
        names_list[[lib]] <- subset(rownames(dataCount1)[dataCount1$library==lib], rownames(dataCount1)[dataCount1$library==lib] %in% assembled_a$duo1)
      }
      else{
        message(paste0(lib, ": duo2"))
        names_list[[lib]] <- subset(rownames(dataCount1)[dataCount1$library==lib], rownames(dataCount1)[dataCount1$library==lib] %in% assembled_a$duo2)
      }
    }
  }
  
  return(names_list)
}

### Separates libraries to allow for normalization based on library
## INPUTS:
##    dataCount : Count data as a data frame
##    dataCond  : Condition matrix for data - should have 2 columns labeled condition and batch
##    run       : Run identifier - can be int or character
##    libExcl   : OPTIONAL. A list of libraries to be excluded. 
## OUTPUTS:
##    lib_list  : List of separated libraries
duoPrep <- function(dataCount, dataCond, run, namesList, libExcl = c()){
  lib_list <- list()
  for(lib in levels(dataCount$library)){
    if(lib %in% libExcl) next
    message(lib)
    lib_list[[paste0(lib, "_", run)]] <- dataCount[namesList[[lib]], rownames(dataCond)]
    message("complete")
  }
  return(lib_list)
}

### Performs the DESeq analysis on the data writes out results file of initial DESeq results as well as 
### results of mean shifting
## INPUTS:
##    dataCount : Count data as a data frame
##    dataCond  : Condition matrix for data - should have one column labeled condition w/rownames as count column names
##    run       : Run identifier - can be int or character
##    filePrefix: Identifier for the specific dataset
##    namesList : Output from the duoNames function
##    libExcl   : OPTIONAL. A list of libraries to be excluded.
##    negListE  : negative controls for the enhancer library
##    negListS  : negative controls for the silencer library
## OUTPUTS:
##    dds_list  : DESeqDataSet for the given run with each library an element of the list
duoSeq <- function(duoAttr, dataCount, dataCond, run, filePrefix, namesList, libExcl=c(), negListS = c(), negListE = c(), filterPost = F,tTest=T){
  `%notin%` <- Negate(`%in%`)
  lib_list <- duoPrep(dataCount,dataCond, run, namesList, libExcl)
  dataCond$condition <- as.factor(dataCond$condition)
  silencer_enhancer_split <- colsplit(duoAttr$ID, "\\^", names = c("enhancer","silencer"))
  num_enh <- length(unique(silencer_enhancer_split$enhancer))
  message(num_enh)
  dds_list <- list()
  dds_rna <- list()
  control_list <- list()
  for(lib in names(lib_list)){
    message(lib)
    message(class(lib))
    dataCount_lib <- (lib_list[[lib]])
    ## Perform DESeq analysis across all cell types
    dds_temp <- DESeqDataSetFromMatrix(countData = dataCount_lib,colData = dataCond, design = ~ condition)
    dds_temp$condition <- relevel(dds_temp$condition, "DNA")
    dds_temp_res <- DESeq(dds_temp, fitType = 'local')
    dds_list[[lib]] <- dds_temp_res
    ## Initialize RNA Library List
    dds_rna[[lib]] <- list()
    ## Preform celltype specific analysis
    ## Overall levels: 1 - Library, 2 - Celltype
    for(celltype in levels(dataCond$condition)){
      if(celltype=="DNA") next
      message(celltype)
      rna_cols <- subset(dataCond, condition==celltype)
      rna_count <- dataCount_lib[,rownames(rna_cols)]
      dds_rna_temp <- DESeqDataSetFromMatrix(countData = rna_count, colData = rna_cols, design = ~1)
      sizeFactors(dds_rna_temp) <- sizeFactors(dds_temp_res)[rownames(rna_cols)]
      dds_rna_temp <- estimateDispersions(dds_rna_temp)
      message(paste0(dim(dds_rna_temp),collapse = "\t"))
      dds_rna[[lib]][[celltype]] <- dds_rna_temp
    }
  }
  dds_list_og <- dds_list
  res_list <- list()
  res_dds <- list()
  still_to_check <- list()
  for(dds in names(dds_list)){
    for(celltype in levels(dataCond$condition)){
      if(celltype=="DNA") next
      message(celltype)
      dds_temp_out <- dds_list[[dds]]
      res_dds[[dds]] <- results(dds_temp_out, contrast = c("condition", celltype, "DNA"))
      ## If the enhancer negative contol(s) is/are in the library
      if(length(intersect(negListE,rownames(dds_list[[dds]])))!=0){
        for(celltype in levels(dataCond$condition)){
          if(celltype=="DNA") next
          message(dds)
          neg_E_index <- c()
          for(neg in 1:length(negListE)){
            x <- grep(negListE[neg], rownames(dds_temp_out))
            neg_E_index[neg] <- x
          }
          neg_reference <- results(dds_temp_out)[neg_E_index,]
          control_list[[dds]] <- rownames(neg_reference)
          log_FC_E <- mean(neg_reference$log2FoldChange, na.rm=T)
          neg_adj <- 2^log_FC_E
          ## Summit Shift Normalization
          size_factor_E <- c(neg_adj, neg_adj, neg_adj, neg_adj,1,1,1,1)
          sizeFactors(dds_temp_out)[which(dataCond$condition == celltype)] <- sizeFactors(dds_temp_out)[which(dataCond$condition == celltype)]*neg_adj
          ## Replace dispersions with celltype specific dispersions and redo wald test
          dds_neg_res <- duoSig(dds_temp_out, dds_rna[[dds]], dataCond)
          dds_list[[dds]] <- dds_neg_res
          res_list[[dds]] <- results(dds_neg_res)
        }
      }
      ## If the silencer negative control(s) is/are in the library
      if(length(intersect(negListS,rownames(dds_list[[dds]])))!=0){
        for(celltype in levels(dataCond$condition)){
          if(celltype=="DNA") next
          message(dds)
          neg_S_index <- c()
          for(neg in 1:length(negListS)){
            x <- grep(negListS[neg], rownames(dds_temp_out))
            neg_S_index[neg] <- x
          }
          reference <- results(dds_temp_out)[neg_S_index,]
          list_reference <- row.names(reference)
          control_list[[dds]] <- list_reference
          log_FC_S <- mean(reference$log2FoldChange)
          rand_adj <- 2^log_FC_S
          ## Summit Shift Nomalization
          # size_factor_S <- c(rand_adj, rand_adj, rand_adj, rand_adj, 1,1,1,1)
          sizeFactors(dds_temp_out)[which(dataCond$condition == celltype)] <- sizeFactors(dds_temp_out)[which(dataCond$condition == celltype)]*rand_adj
          ## Replace dispersions with celltype specific dispersions and redo wald test
          dds_rand_res <- duoSig(dds_temp_out, dds_rna[[dds]], dataCond)
          dds_list[[dds]] <- dds_rand_res
          res_list[[dds]] <- results(dds_rand_res)
        }
      }
      ## If the first/only enhancer negative control appears in an element of the library, but doesn't appear by itself
      if(length(grep(negListE[1], rownames(dds_list[[dds]]))) > 0 & negListE[1] %notin% rownames(dds_list[[dds]])){
        message(paste0("Checking ", dds, " in a minute"))
        still_to_check[[dds]] <- dds_list[[dds]]
      }
    }
  }
  ## Get all combos of negative controls
  neg_randoms_A <- expand.grid(negListE, negListS)
  neg_randoms_A <- mutate(neg_randoms_A, duo_neg = paste(!!!rlang::syms(c("Var2","Var1")), sep = "^"))
  duo_neg_list_A <- neg_randoms_A$duo_neg
  duo_neg_list_A <- duo_neg_list_A[duo_neg_list_A %in% namesList[[1]]]
  message(duo_neg_list_A[1])
  neg_randoms_B <- expand.grid(negListE, negListS)
  neg_randoms_B <- mutate(neg_randoms_B, duo_neg = paste(!!!rlang::syms(c("Var1","Var2")), sep = "^"))
  duo_neg_list_B <- neg_randoms_B$duo_neg
  duo_neg_list_B <- duo_neg_list_B[duo_neg_list_B %in% namesList[[1]]]
  message(duo_neg_list_B[1])
  for(lib in names(still_to_check)){
    dds_temp_out <- dds_list[[lib]]
    for(celltype in levels(dataCond$condition)){
      if(celltype=="DNA") next
      duo_ctrl_index <- c()
      ## Check the list the library is in
      if(duo_neg_list_A[1] %in% rownames(results(dds_temp_out))){
        message("listA")
        duo_ctrl_index <- intersect(duo_neg_list_A, rownames(results(dds_temp_out)))
      }
      if(duo_neg_list_B[1] %in% rownames(results(dds_temp_out))){
        message("listB")
        duo_ctrl_index <- intersect(duo_neg_list_B, rownames(results(dds_temp_out)))
      }
      duo_ctrl <- results(dds_temp_out, contrast = c("condition",celltype,"DNA"))[duo_ctrl_index,]
      control_list[[lib]] <- rownames(duo_ctrl)
      log_FC_duo <- mean(duo_ctrl$log2FoldChange, na.rm=T)
      duo_adj <- 2^log_FC_duo
      ## Summit Shift Normalization
      # size_factor_duo <- c(duo_adj, duo_adj, duo_adj, duo_adj, 1,1,1,1)
      message(celltype)
      sizeFactors(dds_temp_out)[which(dataCond$condition == celltype)] <- sizeFactors(dds_temp_out)[which(dataCond$condition == celltype)]*duo_adj
    }
    ## Replace dispersions with celltype specific dispersions and redo wald test
    dds_duo_res <- duoSig(dds_temp_out,dds_rna[[lib]],dataCond)
    dds_list[[lib]] <- dds_duo_res
    res_list[[lib]] <- list()
    if(filterPost==T){
      plas_counts <- counts(dds_list[[lib]])[which(dataCond$condition=="DNA"),]
      plas_means <- rowMeans(plas_counts)
      dds_list[[lib]] <- dds_list[[lib]][plas_means>10]
    }
    for(celltype in levels(dataCond$condition)){
      if(celltype=="DNA") next
      res_list[[lib]][[celltype]] <- results(dds_duo_res, contrast = c("condition",celltype,"DNA"))
    }
  }
  results_out_list <- list()
  res_dds_list <- list()
  res_dds_rownames <- data.frame()
  for(res in names(res_list)){
    counts_norm <- counts(dds_list[[lib]], normalized = T)
    message("expanding counts")
    norm_exp <- expandDups(counts_norm)
    outA <- list()
    full_output <- list()
    for(celltype in names(res_list[[res]])){
      # if(celltype != "GM12878") next
      outputA <- res_list[[lib]][[celltype]]

      message(celltype)
      message("Results of dds_results recieved")

      ctrl_cols <- row.names(dataCond)[which(dataCond$condition=="DNA")]
      exp_cols <- row.names(dataCond)[which(dataCond$condition==celltype)]

      message("Control and Experiment Columns Set")

      ctrl_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% ctrl_cols], na.rm = T)
      exp_mean <- rowMeans(counts_norm[, colnames(counts_norm) %in% exp_cols], na.rm = T)
      output_2 <- cbind(ctrl_mean,exp_mean,outputA[,-1])

      message("Removing Duplicates")
      dups_output<-expandDups(output_2)

      message("Writing standard results files")
      results_out<-merge(duoAttr, as.matrix(dups_output), by.x="ID", by.y="row.names", all.x=TRUE, no.dups=F)
      full_output[[celltype]] <- results_out
      write.table(results_out, file=paste0("results/",filePrefix,"_",celltype,"_results.run",run,".txt"), row.names=F, col.names=T, sep='\t', quote=F)

      if(tTest==T){
        message("Writing T-Test Results Files")
        outA[[celltype]]<-cellSpecificTtest(duoAttr, norm_exp, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols, n=num_enh, altRef=T)
        write.table(outA[[celltype]],paste0("results/", filePrefix, "_", celltype, "_emVAR.out"), row.names=F, col.names=T, sep="\t", quote=F)
      }
      
      message("Writing bed File")
      full_bed_outputA<-merge(duoAttr, as.matrix(dups_output),by.x="ID",by.y="row.names",all.x=TRUE,no.dups=FALSE)
      #printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]    
      printbed<-full_bed_outputA[,c("chr","start","stop","ID","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","project")]        
      printbed$score<-"."
      #printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","cigar","md-tag","project")]      
      #colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","cigar","md-tag","project")
      printbed<-printbed[,c("chr","start","stop","ID","score","strand","log2FoldChange","Ctrl.Mean","Exp.Mean","pvalue","padj","lfcSE","project")]      
      colnames(printbed)<-c("chr","start","stop","id","score","strand","log2fc","input-count","output-count","log10pval","log10fdr","lfc-se","project")
      printbed$strand[printbed$strand=="fwd"]="+"
      printbed$strand[printbed$strand=="rev"]="-"
      printbed$log10pval=-log10(printbed$log10pval)
      printbed$log10fdr=-log10(printbed$log10fdr)
      
      write.table(printbed,paste0("results/",filePrefix,"_",celltype,".bed"),row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
    }
  }
  
  ## Plotting Baseline vs. normalized for overall and negative controls
  for(lib in names(dds_list)){
    for (celltype in levels(dataCond$condition)) {
      if(celltype == "DNA") next

      temp_outputB <- results(dds_list_og[[lib]], contrast=c("condition",celltype,"DNA"))

      outputA <- results(dds_list[[lib]], contrast=c("condition",celltype,"DNA"))

      message("Plotting Normalization Curves")
      pdf(paste0("plots/Normalized_FC_Density_",celltype,"_",filePrefix,".pdf"),width=10,height=10)
      plot(density(temp_outputB[control_list[[lib]],]$log2FoldChange,na.rm=TRUE),xlim=c(-5,5),ylim=c(0,1.5),col="grey",main=paste0("Normalization - ",celltype))
      lines(density(temp_outputB$log2FoldChange,na.rm=TRUE),xlim=c(-5,5),col="black")
      lines(density(outputA$log2FoldChange,na.rm=TRUE),xlim=c(-5,5),col="red")
      lines(density(outputA[control_list[[lib]],]$log2FoldChange,na.rm=TRUE),xlim=c(-5,5),col="salmon")
      text(1.5,0.4,adj=c(0,0),labels="All - baseline",col="black")
      text(1.5,0.35,adj=c(0,0),labels="All - corrected",col="red")
      text(1.5,0.3,adj=c(0,0),labels="Negatives - baseline",col="grey")
      text(1.5,0.25,adj=c(0,0),labels="Negatives - corrected",col="salmon")
      abline(v=0)
      dev.off()
    }
  }
  return(dds_list)
  # return(outA)
}

### Replaces celltype agnostic dispersions with celltype specific dispersions
## INPUTS:
##    dds_results : DESeq object with adjusted size factors
##    dds_rna     : Celltype specific DESeq object list for only that library (level 2)
##    cond_data   : Condition dataframe
## OUTPUTS:
##    dds_results : DESeq object with the dispersions replaced with celltype specific dispersions
duoSig <- function(dds_results, dds_rna, cond_data){
  for(celltype in levels(cond_data$condition)){
    if(celltype == "DNA") next
    message(celltype)
    dispersions(dds_rna[[celltype]])[which(is.na(dispersions(dds_rna[[celltype]])))] <- 10 #max(dispersions(dds_results))
    mcols(dds_results)$dispersion <- dispersions(dds_rna[[celltype]])
    message("dispersions replaced")
    dds_results <- nbinomWaldTest(dds_results)
  }
  message(paste(dim(dds_results), collapse = "\t"))
  return(dds_results)
}

### Function to perform TTest on individual cell types
## INPUTS:
##    attributesData  : table of full attributes, columns should include: ID, SNP, Project, Window, Strand, Allele,
##                      Haplotype, Bash
##    counts_norm     : passed from the dataOut function 
##    dups_output     : passed from the expandDups function
##    ctrl_mean       : passed from the dataOut function
##    exp_mean        : passed from the dataOut function
##    ctrl_cols       : passed from the dataOut function
##    exp_cols        : passed from the dataOut function
##    altRef          : Logical, default T indicating sorting by alt/ref, if sorting ref/alt set to F
cellSpecificTtest<-function(attributesData, counts_norm, dups_output, ctrl_mean, exp_mean, ctrl_cols, exp_cols,n, altRef = T, correction="BH", cutoff=0.01, prior=F){
  snp_data <- subset(attributesData,allele=="ref" | allele=="alt")
  snp_data <- snp_data[order(snp_data$SNP),]
  snp_data$comb <- paste(snp_data$SNP,"_",snp_data$window,"_",snp_data$strand,"_",snp_data$haplotype,sep="")
  tmp_ct <- as.data.frame(table(as.character(snp_data$comb)))
  
  snp_data_pairs <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq==n*2,]$Var1,]
  
  snp_data_rejected <- snp_data[snp_data$comb %in% tmp_ct[tmp_ct$Freq!=n*2,]$Var1,]
  
  snp_data_ctdata_pairs <- merge(snp_data_pairs, counts_norm, by.x="ID", by.y="row.names", all.x=T, no.dups=F)
  snp_data_ctdata_pairs <- snp_data_ctdata_pairs[order( snp_data_ctdata_pairs$SNP, snp_data_ctdata_pairs$window, snp_data_ctdata_pairs$allele),]
  snp_data_expdata_pairs <- merge(snp_data_pairs,dups_output, by.x="ID", by.y="row.names", all.x=T, no.dups=F)
  snp_data_expdata_pairs <- snp_data_expdata_pairs[order(snp_data_expdata_pairs$SNP, snp_data_expdata_pairs$window, snp_data_expdata_pairs$allele),]
  
  if(altRef==T){
    evens <- c()
    odds <- c()
    for(i in (n+1):(2*n)){
      evens <- c(evens, seq(i, nrow(snp_data_pairs),by=n*2))
    }
    for(i in 1:n){
      odds <- c(odds, seq(i, nrow(snp_data_pairs), by=n*2))
    }
  }else{
    evens <- c()
    odds <- c()
    for(i in 1:n){
      evens <- c(evens, seq(i, nrow(snp_data_pairs),by=n*2))
    }
    for(i in (n+1):(2*n)){
      odds <- c(odds, seq(i, nrow(snp_data_pairs), by=n*2))
    }
  }
  

  A_Ctrl_Mean <- snp_data_expdata_pairs[evens, "ctrl_mean"]
  A_Exp_Mean <- snp_data_expdata_pairs[evens, "exp_mean"]
  A_log2FC <- snp_data_expdata_pairs[evens, "log2FoldChange"]
  A_log2FC_SE <- -log10(snp_data_expdata_pairs[evens, "lfcSE"])
  A_logP <- -log10(snp_data_expdata_pairs[evens, "pvalue"])
  A_logP[is.na(A_logP)] <- 0
  A_logP[A_logP == Inf] <- max(A_logP[is.finite(A_logP)])
  A_logPadj_BH <- -log10(snp_data_expdata_pairs[evens, "padj"]) #BH Correction
  A_logPadj_BH[A_logPadj_BH < 0]<-0
  A_logPadj_BH[A_logPadj_BH == Inf] <- max(A_logPadj_BH[is.finite(A_logPadj_BH)])
  A_logPadj_BF <- -log10(snp_data_expdata_pairs[evens, "pvalue"]*(nrow(snp_data_expdata_pairs)/2)) #BF Correction
  A_logPadj_BF[A_logPadj_BF < 0]<-0
  A_logPadj_BF[A_logPadj_BF == Inf] <- max(A_logPadj_BF[is.finite(A_logPadj_BF)])
  B_Ctrl_Mean <- snp_data_expdata_pairs[odds, "ctrl_mean"]
  B_Exp_Mean <- snp_data_expdata_pairs[odds, "exp_mean"]
  B_log2FC <- snp_data_expdata_pairs[odds, "log2FoldChange"]
  B_log2FC_SE <- -log10(snp_data_expdata_pairs[odds, "lfcSE"])
  B_logP <- -log10(snp_data_expdata_pairs[odds, "pvalue"])
  B_logP[is.na(B_logP)] <- 0
  B_logP[B_logP == Inf] <- max(B_logP[is.finite(B_logP)])
  B_logPadj_BH <- -log10(snp_data_expdata_pairs[odds, "padj"]) #BH Correction
  B_logPadj_BH[B_logPadj_BH < 0]<-0
  B_logPadj_BH[B_logPadj_BH == Inf] <- max(B_logPadj_BH[is.finite(B_logPadj_BH)])
  B_logPadj_BF <- -log10(snp_data_expdata_pairs[odds, "pvalue"]*(nrow(snp_data_expdata_pairs)/2)) #BF Correction
  B_logPadj_BF[B_logPadj_BF < 0]<-0
  B_logPadj_BF[B_logPadj_BF == Inf] <- max(B_logPadj_BF[is.finite(B_logPadj_BF)])

  out2 <- cbind(snp_data_expdata_pairs[evens,c(1:7,9)],A_Ctrl_Mean, A_Exp_Mean, A_log2FC, A_log2FC_SE,A_logP,A_logPadj_BH,A_logPadj_BF,B_Ctrl_Mean,B_Exp_Mean,B_log2FC,B_log2FC_SE,B_logP,B_logPadj_BH,B_logPadj_BF)

  odds_ctrl_means <- rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols])
  odds_exp_means <- rowMeans(snp_data_ctdata_pairs[odds,exp_cols])
  even_ctrl_means <- rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols])
  even_exp_means <- rowMeans(snp_data_ctdata_pairs[evens,exp_cols])
  
  odds_ctrl_na <- is.na(odds_ctrl_means)
  odds_exp_na <- is.na(odds_exp_means)
  even_ctrl_na <- is.na(even_ctrl_means)
  even_exp_na <- is.na(even_exp_means)
  
  maxes <- c(max(odds_ctrl_means[!odds_ctrl_na]),max(odds_exp_means[!odds_exp_na]),max(even_ctrl_means[!even_ctrl_na]),max(even_exp_means[!even_exp_na]))

  # Don't try to do the t test for ones with all zeros.
  ignore_idx <- which(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[odds, exp_cols]) < 10 | 
                        is.na(rowMeans(snp_data_ctdata_pairs[odds,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[odds, exp_cols])) | 
                        rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]) < 10 | rowMeans(snp_data_ctdata_pairs[evens, exp_cols]) < 10 | 
                        is.na(rowMeans(snp_data_ctdata_pairs[evens,ctrl_cols]))  | is.na(rowMeans(snp_data_ctdata_pairs[evens, exp_cols])) )
  
  # For the numerator, set zero values to 1 so that the log-ratio is defined.
  counts1 <- snp_data_ctdata_pairs
  counts1[counts1 == 0] <- 1
  
  # t test
  ratios_A <- log((counts1[evens, exp_cols]) / rowMeans(snp_data_ctdata_pairs[evens, ctrl_cols]))
  ratios_B <- log((counts1[odds, exp_cols]) / rowMeans(snp_data_ctdata_pairs[odds, ctrl_cols]))
  
  ratios_list <- list(ratios_A, ratios_B)
  
  t_pvalue <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else if(any(duplicated(as.numeric(ratios_A[i,]))) | any(duplicated(as.numeric(ratios_B[i,])))){NA} else{
    t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$p.value})
  t_stat <- sapply(1:nrow(ratios_A), function(i) if(i %in% ignore_idx){NA} else if(any(duplicated(as.numeric(ratios_A[i,]))) | any(duplicated(as.numeric(ratios_B[i,])))){NA} else{
    t.test(as.numeric(ratios_A[i,]), as.numeric(ratios_B[i,]), var.equal=F, paired=T)$statistic})

  out2$LogSkew <- out2$B_log2FC - out2$A_log2FC
  out2$Skew_logP <- ifelse(is.na(t_pvalue), 0, -log10(t_pvalue))

  OE_threshold <- -log10(.01)
  is_OE <- out2$A_logPadj_BF >= OE_threshold | out2$B_logPadj_BF >= OE_threshold
  out2$Skew_logFDR <- rep(NA, dim(out2)[1])
  q_idx <- intersect(which(is_OE), which(!is.na(t_pvalue)))
  out2$Skew_logFDR[q_idx] <- -log10(p.adjust(t_pvalue[q_idx],method="BH"))
  
  return(out2)
  # return(ratios_list)
}

### Expand IDs that denote duplicate oligos in count/DESeq results
# Requires the output_2 of dataOut function. 
# i.e: dataOut(countsData, attributesData, conditionData)$output_2
## Returns: Final output table with no duplicates
expandDups <- function(output){
  output_orig <- output
  split_row_names <- colsplit(rownames(output), "\\^", names=c("Enhancer","Silencer"))
  message("adding split rownames")
  output_new <- cbind(split_row_names, output)
  # Identify duplicates, if any exist
  dups <- output_new[grep("\\(.*\\)$",output_new$Silencer),]
  dups$Silencer <- gsub("^\\((.*)\\)$","\\1",dups$Silencer)
  # Add everything but the duplicates to the final output
  output_final <- output_new[-(grep("\\(.*\\)$",output_new$Silencer)),] 
  output_final <- output_final[!(is.na(output_final$Silencer)),]
  # If there are 1 or more duplicates.... Not clear what this step is doing????
  if(nrow(dups) > 0) {
    for(i in 1:nrow(dups)){
      dup_id <- unlist(strsplit(dups$Silencer[i],";"))
      dup_exp <- dups[rep(i, length(dup_id)), ]
      dup_exp$Silencer <- dup_id
      # dup_exp$Enhancer <- dups$Enhancer[i]
      output_final <- rbind(dup_exp,output_final)
    }
    rownames(output_final) <- paste(output_final$Enhancer,output_final$Silencer, sep = "^")
    output_final <- as.data.frame(output_final[,-c(1,2)])
    message(class(output_final))
  }else {
    output_final <- output_orig
  }
  return(output_final)
}

### Performs a quantile normalization on the log2FoldChange and returns the DESeq Results with the normalized l2FC
## INPUTS:
##    ddsList1  : Result list from duoSeq function from one run
##    ddsList2  : Result list from duoSeq function from second run
##    dataCond  : Condition matrix for data - should have one column labled condition w/rownames as count column names
##    filePrefix: Identifier for the specific dataset
##    pairsList : List of pairs to normalize against each other (i.e. A_2,A_3)
## OUTPUTS:
##    dds_norm  : List of lists of DESeq results from each library, first element will be normalized ddsList1 second element will be ddsList2
# duoNorm <- function(ddsList1, ddsList2, dataCond, filePrefix, pairsList = c()){
#   `%notin%` <- Negate(`%in%`)
#   dds_norm1 <- list()
#   dds_norm2 <- list()
#   for(i in 0:((length(pairsList)/2)-1)){
#     pair = 2*i + 1
#     message(pairsList[pair])
#     pair1 <- pairsList[pair]
#     message(pairsList[(pair+1)])
#     pair2 <- pairsList[(pair+1)]
#     
#     dds_res1 <- as.data.frame(results(ddsList1[[pair1]]))
#     dds_res2 <- as.data.frame(results(ddsList2[[pair2]]))
#     
#     ## Normalize log2FoldChange between pair elements, then calculate the factors to adjust the count data by
#     list1_pair1_l2fc <- as.matrix(dds_res1[,2])
#     list2_pair2_l2fc <- as.matrix(dds_res2[,2])
# 
#     merged <- cbind(list1_pair1_l2fc, list2_pair2_l2fc)
#     rownames(merged) <- rownames(results(ddsList1[[pair1]]))
#     mergedn <- normalize.quantiles(as.matrix(merged), copy = F)
#     factor1 <- (2^mergedn[,1])/(2^merged[,1])
#     factor2 <- (2^mergedn[,2])/(2^merged[,2])
#     factor1[is.na(factor1)] <- 1
#     factor2[is.na(factor2)] <- 1
#     if(pairsList[pair] != "PinotA_2"){
#       rownames(mergedn) <- rownames(merged)
#     }
#     
#     #Expected log2Fold for testing purposes
#     l2fc1 <- as.data.frame(mergedn[,1])
#     l2fc2 <- as.data.frame(mergedn[,2])
#     
#     if(pairsList[pair] == "PintoA_2"){
#       rownames(l2fc1) <- rownames(results(ddsList1[[pair1]]))
#       rownames(l2fc2) <- rownames(results(ddsList2[[pair2]]))
#     }
#     
#     message("writing tables")
#     write.table(l2fc1, paste0(names(ddsList1[pair1]),"_", filePrefix, "_normalized_log2FoldChange.txt"), col.names = F)
#     write.table(l2fc2, paste0(names(ddsList2[pair2]),"_", filePrefix, "_normalized_log2FoldChange.txt"), col.names = F)
#     
#     dds_res1$log2FoldChange <- l2fc1
#     dds_res2$log2FoldChange <- l2fc2
#     
#     message("adding to list")
#     dds_norm1[[pair1]] <- dds_res1
#     dds_norm2[[pair2]] <- dds_res2
#     
#   }
#   return(list(dds_norm1, dds_norm2))
# }

### Plots the correlation of the overall dataset, as well as within individual celltypes
## INPUTS:
##    dataCount : Count data (either normalized from DESeq or oligo level raw)
##    namesList : output from the duoNames function
##    dataCond  : Condition matrix for data - should have one column labeled condition w/rownames as count column names
##    filePrefix: Identifier for the specific dataset
##    run       : run number to use
##    libExcl   : OPTIONAL. A list of libraries to be excluded. (RAW counts only)
##    libIncl   : A list of libraries to be included. (Neccesary for raw counts)
## OUTPUTS:
##    Plots correlation between all replicates for all celltypes as well as specific celltypes
duoCor <- function(dataCount, namesList, dataCond, filePrefix, run, libExcl = c(), libIncl = c()){
  celltypes <- as.factor(dataCond$condition)
  message(levels(celltypes))
  if("library" %in% colnames(dataCount)){
    message(paste0(dim(dataCount),collapse = "\t"))
    lib_list <- duoPrep(dataCount,dataCond, run, namesList, libExcl)
    message(names(lib_list))
    message(paste0(dim(lib_list[[1]]), collapse = "\t"))
    count_data <- as.data.frame(lib_list[[libIncl]])
    message(paste0(dim(count_data),collapse = "\t"))
    message("count_data assigned")
  }
  else {
    count_data <- dataCount
    colnames(count_data) <- c("P_1","P_2","P_3","P_4","G_1","G_2","G_3","G_4","K_1","K_2","K_3","K_4","H_1","H_2","H_3","H_4","S_1","S_2","S_3","S_4")
  }
  # return(count_data)
  # initial correlation matrix
  # p1 <- ggpairs(count_data, title = filePrefix, diag = list("naDiag"))
  # # Colored correlation plot
  # p2 <- ggcorr(count_data)
  # 
  # g2 <- ggplotGrob(p2)
  # col_vals <- g2$grobs[[6]]$children[[3]]$gp$fill
  # 
  # # Add colors as background of upper tiles
  # idx <- 1
  # for(k1 in 1:(ncol(count_data)-1)){
  #   for(k2 in (k1+1):ncol(count_data)){
  #     plt <- getPlot(p1,k1,k2) + theme(panel.background = element_rect(fill = col_vals[idx], color = "white"),
  #                                      panel.grid.major=element_line(color = col_vals[idx]))
  #     p1 <- putPlot(p1,plt,k1,k2)
  #     idx <- idx+1
  #   }
  # }
  # png(paste0("plots/correlation_plots_colored/",filePrefix,".png"), width = 1024, height = 1024)
  # print(p1)
  # dev.off()
  for(celltype in levels(celltypes)){
    message(celltype)
    cell_cols <- subset(dataCond, condition==celltype)
    message(paste0(rownames(cell_cols), collapse = "\t"))
    cell_count <- count_data[,rownames(cell_cols)]
    message(paste0(colnames(cell_count),collapse = "\t"))
    png(paste0("plots/corrleation_",filePrefix,"_",celltype,".png"))
    print(ggpairs(cell_count, title = paste0("Correlation ", filePrefix," ", celltype), diag = list("naDiag")))
    dev.off()
  }
}

### Plots the correlation between the log2FoldChange of all celltypes, with the option to subset from a list of oligos
## INPUTS:
##    dds_list  : List of DESeqDataSet objects, output of duoSeq function
##    dataCond  : Condition matrix for data - should have one column labeled condition w/rownames as count column names
##    filePrefix: Identifier for the specific dataset
##    subsetList: OPTIONAL. List of oligos to compare, if blank all oligos considered.
## OUTPUTS:
##    Plots correlation between the log2FoldChange of all celltypes
duoLogCor <- function(dds_list, dataCond, filePrefix, subsetList = c()){
  celltypes <- as.factor(dataCond$condition)
  for(lib in names(dds_list)){
    index <- 0
    for(celltype in levels(celltypes)){
      if(celltype=="DNA") next
      index <- index + 1
      message(celltype)
      if(index == 1){
        l2fc_df <- as.data.frame(results(dds_list[[lib]], contrast=c("condition",celltype,"DNA"))$log2FoldChange)
        rownames(l2fc_df) <- rownames(results(dds_list[[lib]]))
        colnames(l2fc_df) <- c(paste0(celltype,"_l2fc"))
      }
      if(index > 1){
        temp <- as.data.frame(results(dds_list[[lib]], contrast=c("condition",celltype,"DNA"))$log2FoldChange)
        rownames(temp) <- rownames(results(dds_list[[lib]]))
        colnames(temp) <- c(paste0(celltype,"_l2fc"))
        l2fc_df <- merge(l2fc_df, temp, by="row.names", all=T)
        row_names <- l2fc_df$Row.names
        l2fc_df <- l2fc_df[,-1]
        rownames(l2fc_df) <- row_names
      }
    }
  }
  if(length(subsetList) > 0){
    l2fc_df <- subset(l2fc_df, rownames(l2fc_df) %in% subsetList)
  }
  return(l2fc_df)
  # png(paste0("plots/l2fc_correlation_",filePrefix,".png"))
  # print(ggpairs(l2fc_df, title = paste0("Correlation of log2FoldChange across Cell Types ",filePrefix)), diag = list("naDiag"))
  # dev.off()
}

### Plots the log2FoldChange density for each enhancer provided in one plot (currently limited to 5 enhancers) with a subset of oligos 
### indicated in the oligo name (i.e. an element consistently in the oligo name)
## INPUTS:
##    dds_list  : List of DESeqDataSet objects, output of duoSeq function
##    dataCond  : Condition matrix for data - should have one column labeled condition w/rownames as count column names
##    enhList   : character vector of Enhancers to include
##    negIndc   : String indicating the subset of oligos to pull, string should be present in the oligo name.
## OUTPUTS:
##    Plots density of log2FoldChange for each enhancer within a celltype, with a subset of oligos based on the oligo name
duol2fcDensity <- function(duoAttr, dds_list, dataCond, filePrefix, enhList, negIndc, filterList){
  celltypes <- as.factor(dataCond$condition)
  for(lib in names(dds_list)){
    temp <- dds_list[[lib]]
    for(celltype in levels(celltypes)){
      if(celltype=="DNA") next
      message(celltype)
      res_ct <- results(dds_list[[lib]], contrast = c("condition", celltype, "DNA"))
      oligo_split <- colsplit(rownames(res_ct), pattern="\\^", names = c("col1","col2"))
      labs <- c()
      indic <- 0
      message(length(filterList[[celltype]]))
      attr_sub <- duoAttr[which(duoAttr$SNP %in% filterList[[celltype]]$V1),]
      message(paste0(dim(attr_sub), collapse = "\t"))
      sub_freq <- as.data.frame(table(attr_sub$SNP))
      message(paste0(colnames(sub_freq),collapse = "\t"))
      snp_to_keep <- sub_freq$Var1[which(sub_freq$Freq==10)]
      message(paste0(head(snp_to_keep), collapse = "\t"))
      pdf(paste0("plots/",filePrefix,"_",lib,"_l2fc_",celltype,"_with_",negIndc,".pdf"), height=10, width = 10)
      for(e in enhList){
        color_list <- c("red","salmon","blue","dodgerblue","green","limegreen","black","grey", "plum", "palevioletred")
        indic <- indic +1
        if(enhList[1] %in% oligo_split$col1){
          enh_l2fc <- res_ct[oligo_split$col1 == e,]
          enh_l2fc <- enh_l2fc[rownames(enh_l2fc) %in% attr_sub$ID[which(attr_sub$SNP %in% snp_to_keep)],]
          negs <- grepl(negIndc, rownames(enh_l2fc))
          enh_neg_l2fc <- enh_l2fc[negs,]
          enh_l2fc <- enh_l2fc[!negs,]
          neg_snps <- attr_sub$SNP[which(attr_sub$ID %in% rownames(enh_neg_l2fc))]
          ref_snps <- attr_sub$SNP[which(attr_sub$ID %in% rownames(enh_l2fc))]
          overlap_snps <- intersect(neg_snps,ref_snps)
          enh_neg_l2fc <- enh_neg_l2fc[rownames(enh_neg_l2fc) %in% attr_sub$ID[attr_sub$SNP %in% overlap_snps],]
          enh_l2fc <- enh_l2fc[rownames(enh_l2fc) %in% attr_sub$ID[attr_sub$SNP %in% overlap_snps],]
          labs <- c(labs, paste0(e," (", nrow(enh_l2fc), ")"), paste0(e, " - ",negIndc," (",nrow(enh_neg_l2fc),")"))
          if(indic == 1){
            plot(density(enh_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), ylim=c(0,1), col=color_list[indic], main=paste0("Density of enhancers over ", celltype))
            lines(density(enh_neg_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), col=color_list[2*indic])
          }
          if(indic > 1){
            lines(density(enh_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), col=color_list[2*indic-1])
            lines(density(enh_neg_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), col=color_list[2*indic])
          }
        }
        if(enhList[1] %in% oligo_split$col2){
          enh_l2fc <- res_ct[oligo_split$col2 == e & oligo_split$col1 %in% filterList[[celltype]]$V1,]
          enh_l2fc <- enh_l2fc[rownames(enh_l2fc) %in% attr_sub$ID[which(attr_sub$SNP %in% snp_to_keep)],]
          negs <- grepl(negIndc, rownames(enh_l2fc))
          enh_neg_l2fc <- enh_l2fc[negs,]
          enh_l2fc <- enh_l2fc[!negs,]
          neg_snps <- attr_sub$SNP[which(attr_sub$ID %in% rownames(enh_neg_l2fc))]
          ref_snps <- attr_sub$SNP[which(attr_sub$ID %in% rownames(enh_l2fc))]
          overlap_snps <- intersect(neg_snps,ref_snps)
          enh_neg_l2fc <- enh_neg_l2fc[rownames(enh_neg_l2fc) %in% attr_sub$ID[attr_sub$SNP %in% overlap_snps],]
          enh_l2fc <- enh_l2fc[rownames(enh_l2fc) %in% attr_sub$ID[attr_sub$SNP %in% overlap_snps],]
          labs <- c(labs, paste0(e," (", nrow(enh_l2fc), ")"), paste0(e, " - ",negIndc," (",nrow(enh_neg_l2fc),")"))
          if(indic == 1){
            plot(density(enh_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), ylim=c(0,1), col=color_list[indic], main=paste0("Density of enhancers over ", celltype))
            lines(density(enh_neg_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), col=color_list[2*indic])
          }
          if(indic > 1){
            lines(density(enh_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), col=color_list[2*indic-1])
            lines(density(enh_neg_l2fc$log2FoldChange, na.rm=T), xlim=c(-7,9), col=color_list[2*indic])
          }
        }
      }
      bottom <- 0.25
      for(lab in 1:length(labs)){
        yval <- ((length(labs)-lab)*0.05+.25)
        text(x=-5, y=yval, labels = labs[lab], col = color_list[lab])
      }
      dev.off()
    }
  }
}

