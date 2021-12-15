library(graphics)
library(stringr)
library(methods)
library(dplyr)
library(psych)
library(clusterProfiler)
library(stats)
library(magrittr)
library(utils)
library(circlize)
library(gridBase)
library(grid)
library(ComplexHeatmap)
library(pheatmap)
library(ggalluvial)
library(ggplot2)
library(RColorBrewer)
library(jsonlite)
library(networkD3)
library(reshape2)
library(tidyr)
library(enrichplot)
library(ggrepel)
library(ggridges)
library(DOSE)

#' transform count to RPKM or TPM within Hsapiens and Mmusculus
#' @param data a dataframe of count with row of gene and column of sample
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @param toType whether RPKM or TPM of result you wanted, eg "RPKM", "TPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @return the RPKM or TPM of \code{data}
counts2normalized_smartseq2 <- function(data, Org, toType, scale.factor=10^6){
  if(Org=="Homo sapiens"){
    f.tmp <- system.file("extdata", "homo/transcriptmaxLength.Rdata", package="cellcallEXT")
    load(f.tmp)
  }else if(Org=="Mus musculus"){
    f.tmp <- system.file("extdata", "mmus/transcriptmaxLength.Rdata", package="cellcallEXT")
    load(f.tmp)
  }
  # head(g_l)
  a <- data
  # a[1:4,1:4]
  ng=intersect(rownames(a),g_l$symbol) 
  
  exprSet=a[ng,]
  lengths=g_l[match(ng,g_l$symbol),2]
  head(lengths)
  head(rownames(exprSet))
  
  # exprSet[1:4,1:4]
  if(toType=="RPKM"){
    total_count<- colSums(exprSet)
    head(total_count)
    head(lengths)
    
    rpkm <- t(do.call( rbind,
                       lapply(1:length(total_count),
                              function(i){
                                10^9*exprSet[,i]/lengths/total_count[i]
                              }) ))
    return(rpkm)
  }else if(toType=="TPM"){
    rpk <- exprSet/lengths
    total_count<- colSums(rpk)
    tpm <- t(do.call( rbind,
                      lapply(1:length(total_count),
                             function(i){
                               scale.factor*rpk[,i]/total_count[i]
                             }) ))
    tpm <- data.frame(tpm, stringsAsFactors = FALSE)
    colnames(tpm) <- colnames(exprSet)
    rownames(tpm) <- rownames(exprSet)
    
    return(tpm)
  }
  
}

#' transform count to CPM
#' @param data a dataframe with row of gene and column of sample
#' @param toType whether RPKM or TPM of result you wanted, eg "CPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @return the CPM of \code{data}
counts2normalized_10X <- function(data, toType, scale.factor=10^6){
  exprSet <- data
  if(toType=="CPM"){
    total_count<- colSums(exprSet)
    
    cpm <- t(do.call( rbind,
                      lapply(1:length(total_count),
                             function(i){
                               scale.factor*exprSet[,i]/total_count[i]
                             }) ))
    cpm <- data.frame(cpm, stringsAsFactors = FALSE)
    colnames(cpm) <- colnames(exprSet)
    rownames(cpm) <- rownames(exprSet)
    return(cpm)
  }
}
#' Create Nichobject from Seurat object
#' @param Seurat.object  The Seurat object which stores the expression matrix and cell type information
#' @param slot The name of slot which contains expression matrix, default "counts".
#' @param cell_type The name of specific column which contains cell type information, default "orig.ident".
#' @param source the type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @param Org the species of this scRNA-seq
#' @param status the status of cell disease or control
#' @param expDirection the direction of expression change UP DOWN or BOTH
#' @importFrom Seurat GetAssayData
#' @export
#'
CreateObject_fromSeurat <- function(Seurat.object, slot="counts", status = NULL, expDirection = "BOTH", ##SG   add onre more parameters
                                    cell_type="orig.ident", data_source="UMI", scale.factor=10^6,
                                    Org="Homo sapiens"){
  Seurat.object <- Seurat.object
  slot <- slot # counts, data
  cell_type <- cell_type
  status <- status ##SG   add onre more parameters
  your_source <- data_source # "UMI", "fullLength", "TPM", or "CPM".
  
  is_Seurat <- is(Seurat.object, "Seurat")
  
  if(is_Seurat){
    if(require(Seurat)){
      
      #tryCatch({
      data <-  as.matrix(Seurat::GetAssayData(object = Seurat.object, slot = slot, assay = "RNA"))
      #},error=function(e){
      #  stop("there is no slot: ", slot," in your seurat object in RNA assay")
      #})
      
      #tryCatch({
      myCelltype <- as.character(Seurat.object@meta.data[,cell_type])
      #},error=function(e){
      #  stop("there is no columns: ", cell_type," in meta.data of Seurat.")
      #})
      
      data <- as.data.frame(data)
      colnames(data) <- paste(1:length(myCelltype), myCelltype, sep = "_")
      
      mt <- CreateNichConObject(data=data, min.feature = 0,
                                names.field = 2,
                                names.delim = "_",
                                source = your_source, # fullLength, UMI, TPM
                                scale.factor = scale.factor,
                                Org = Org,
                                status = status, ##SG   add onre more parameters
                                expDirection = expDirection, ##SG   add onre more parameters
                                project = "Microenvironment")
      
      return(mt)
      
    }else{
      stop("No package Seurat! Fail to extract data from Seurat object.")
    }
  }else{
    stop("It's not Seurat object!")
  }
  
}

#' melt enrich dataframe for bubble graph
#' @param pathway.hyper.list a list of pathway enrichment pvalue and NES
#' @param Org the species of this scRNA-seq
#' @return the result for pbubble graph
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join inner_join
#' @importFrom utils read.table
#' @export

getForBubble <- function(pathway.hyper.list = pathway.hyper.list, cella_cellb, Org="Homo sapiens"){
  
  pathway.hyper.df <- Reduce(function(x,y) full_join(x, y, by=c("pathway")), pathway.hyper.list, accumulate =FALSE)
  rownames(pathway.hyper.df) <- pathway.hyper.df$pathway
  pathway.hyper.df$pathway <- NULL
  
  pathway.hyper.df.pvalue <- pathway.hyper.df[, seq(from=1, to=ncol(pathway.hyper.df), by=3)]
  colnames(pathway.hyper.df.pvalue) <- cella_cellb
  rownames(pathway.hyper.df.pvalue) <- rownames(pathway.hyper.df)
  pathway.hyper.df.pvalue[is.na(pathway.hyper.df.pvalue)] <- 1
  
  pathway.hyper.df.NES <- pathway.hyper.df[, seq(from=3, to=ncol(pathway.hyper.df), by=3)]
  colnames(pathway.hyper.df.NES) <- cella_cellb
  rownames(pathway.hyper.df.NES) <- rownames(pathway.hyper.df)
  pathway.hyper.df.NES[is.na(pathway.hyper.df.NES)] <- 0
  
  pathway.hyper.df.NES <- pathway.hyper.df.NES[rowSums(pathway.hyper.df.NES)!=0,]
  pathway.hyper.df.pvalue <- pathway.hyper.df.pvalue[rowSums(pathway.hyper.df.NES)!=0,]
  
  if(Org=="Homo sapiens"){
    f.tmp <- system.file("extdata", "KEGGREACTOME_SYMBOL_ID.txt", package="cellcallEXT")
  }else{
    Org<-"Mus musculus"
    f.tmp <- system.file("extdata", "KEGGREACTOME_SYMBOL_ID_homology.txt", package="cellcallEXT")
  }
  pathway.info <- read.table(f.tmp, sep = '\t', quote = "", header = FALSE, stringsAsFactors = F)
  colnames(pathway.info) <- c('id', 'name', "main.object.name", "sub.object.name")
  rownames(pathway.info) <- pathway.info$id
  
  pathway.hyper.df.pvalue$pathway <- rownames(pathway.hyper.df.pvalue)
  pathway.hyper.df.pvalue$pathway <- pathway.info[pathway.hyper.df.pvalue$pathway, 'name']
  pvalue.melted <- melt(pathway.hyper.df.pvalue, id=c('pathway'),variable.name="cc",value.name="pvalue")
  
  pathway.hyper.df.NES$pathway <- rownames(pathway.hyper.df.NES)
  pathway.hyper.df.NES$pathway <- pathway.info[pathway.hyper.df.NES$pathway, 'name']
  NES.melted <- melt(pathway.hyper.df.NES, id=c('pathway'),variable.name="cc",value.name="NES")
  NES.melted$NES[NES.melted$NES > 2]=2
  NES.melted$NES[NES.melted$NES < 0]=0
  
  myPub.df <- inner_join(pvalue.melted, NES.melted, by=c('pathway', 'cc'))
  myPub.df$Org<-Org####Shouguo need to add one column for Org
  return(myPub.df)
}
#' create a Cellwave objects
#' @param object a Cellwave objects
#' @param probs Set the percentile of gene expression in one celltype to represent mean value, when use.type="median".
#' @param use.type With parameter "median", CellCall set the mean value of gene as zero, when the percentile of gene expression in one celltype below the parameter "probs". The other choice is "mean" and means that we not concern about the percentile of gene expression in one celltype but directly use the mean value.
#' @param pValueCor firlter target gene of TF with spearson, p > pValueCor, default is 0.05
#' @param CorValue firlter target gene of TF with spearson, value > CorValue, default is 0.1
#' @param topTargetCor use topTargetCor of candidate genes which has firlter by above parameters, default is 1, means 100%
#' @param p.adjust gsea pValue of regulons with BH adjusted threshold, default is 0.05
#' @param method "weighted", "max", "mean", of which "weighted" is default. choose the proper method to score downstream activation of ligand-receptor all regulons of given ligand-receptor relation
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @param IS_core logical variable ,whether use reference LR data or include extended datasets
#' @return the result dataframe of \code{cell communication}
#' @importFrom stringr str_split
#' @importFrom stats quantile median
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom utils read.table head
#' @export

ConnectProfile <- function(object, pValueCor=0.05, CorValue=0.1, topTargetCor=1, method="weighted", p.adjust=0.05, use.type="median", probs = 0.9, Org = 'Homo sapiens', IS_core = TRUE){
  options(stringsAsFactors = F) 
  
  if(Org == 'Homo sapiens'){
    f.tmp <- system.file("extdata", "rectome_ligand_receptor_TFs_withReactome.txt", package="cellcallEXT")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_extended.txt", package="cellcallEXT")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target.txt", package="cellcallEXT")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }else if(Org == 'Mus musculus'){
    f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology.txt", package="cellcallEXT")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology_extended.txt", package="cellcallEXT")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target_homology.txt", package="cellcallEXT")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }
  
  triple_relation$pathway_ID <- NULL
  print(triple_relation[1:4,])
  complex_tmp <- triple_relation$Receptor_Symbol[grep(",",triple_relation$Receptor_Symbol)] %>% unique()
  tmp_complex_symbol <- triple_relation$Receptor_Symbol[grep(",",triple_relation$Receptor_Symbol)] %>% unique() %>% str_split(",") %>% unlist %>% unique()
  all.gene.needed <- unique(as.character(c(triple_relation$Ligand_Symbol, triple_relation$Receptor_Symbol, triple_relation$TF_Symbol, target_relation$TF_Symbol, target_relation$Target_Symbol,tmp_complex_symbol)))
  # triple_relation[1:4,1:4]
  # target_relation[1:4,1:4]
  my_Expr <- object@data$withoutlog
  colnames(my_Expr) <- as.character(object@meta.data$celltype)
  my_Expr[1:4,1:4]
  detect_gene <- rownames(my_Expr)
  
  expr_set <- my_Expr[intersect(detect_gene, all.gene.needed),]  
  detect_gene <- rownames(expr_set)
  cell_type = unique(colnames(expr_set))
  expr.fc <- object@data$withoutlog[detect_gene,]
  colnames(expr.fc) <- colnames(expr_set)
  
  complex_matrix <- matrix(ncol = length(colnames(expr_set)))
  complex_matrix <- as.data.frame(complex_matrix)
  colnames(complex_matrix) <- colnames(expr_set)
  myrownames <- c()
  
  complex <- complex_tmp
  if(length(complex)>0){
    for(i in 1:length(complex)){
      i_tmp = strsplit(complex[i], ',')
      # print(i_tmp)
      if( sum(i_tmp[[1]] %in% detect_gene) == length(i_tmp[[1]]) ){
        tmp_df <- expr_set[i_tmp[[1]],]
        tmp_mean <- colMeans(tmp_df)
        tmp_index <- unique(unlist(apply(tmp_df, 1,function(x) {which(x==0)})))
        tmp_mean[tmp_index] <- 0
        
        # print(res_tmp)
        complex_matrix <- rbind(complex_matrix, tmp_mean)
        myrownames <- c(myrownames, complex[i])
      }
    }
    
    complex_matrix <- complex_matrix[-1,]
    
    ## 把complex的联合表达值加上
    if(nrow(complex_matrix) > 0){
      rownames(complex_matrix) <- myrownames
      expr_set <- rbind(expr_set, complex_matrix)
    }
  }
  
  expr_set <- expr_set[apply(expr_set, 1, function(x){sum(x!=0)})>0,]
  detect_gene <- rownames(expr_set)
  # expr_set[1:4,1:4]
  
  print("step1: compute means of gene")
  expr_mean <- matrix(nrow = nrow(expr_set), ncol = length(cell_type))
  myColnames <- c()
  for (i in 1:length(cell_type)) {
    myCell <- cell_type[i]
    myMatrix <- expr_set[,colnames(expr_set)==myCell,drop=F]
    if(use.type=="mean"){
      myMatrix_mean <- as.numeric(apply(myMatrix, 1, mean))
    }else if(use.type=="median"){
      quantil.tmp <- as.numeric(apply(myMatrix, 1, function(x){
        quantile(x, probs = probs,names=FALSE)
      }))
      mean.tmp <- rowMeans(myMatrix)
      mean.tmp[which(quantil.tmp==0)]<-0 
      myMatrix_mean <- mean.tmp
    }
    expr_mean[,i] <- myMatrix_mean
    myColnames <- c(myColnames, myCell)
    # print(myCell)
  }
  expr_mean <- data.frame(expr_mean)
  colnames(expr_mean) <- myColnames
  rownames(expr_mean) <- rownames(expr_set)
  
  expr_mean <- expr_mean[apply(expr_mean, 1, function(x){sum(x!=0)})>0,]
  detect_gene <- rownames(expr_mean)
  
  
  #if(use.type=="median"){
  #  fc.list <- mylog2foldChange.diy(inData = expr.fc, cell.type = cell_type, method="median", probs = probs)
  #}else{
  #  fc.list <- mylog2foldChange.diy(inData = expr.fc, cell.type = cell_type, method="mean", probs = probs)
  #}
  ####Shouguo
  if(use.type=="median"){
    fc.list <- mylog2foldChange.diy.casecontrol(inData = expr.fc, cell.type = cell_type, method="median", probs = probs, status = object@status)
  }else{
    fc.list <- mylog2foldChange.diy.casecontrol(inData = expr.fc, cell.type = cell_type, method="mean", probs = probs, status = object@status)
  }
  ####Shouguo
  print("step2: filrter tf-gene with correlation, then score regulons")
  tfs_set <- unique(triple_relation$TF_Symbol)
  regulons_matrix <- matrix(data = 0, nrow = length(tfs_set), ncol = length(cell_type))
  
  regulons_matrix <- as.data.frame(regulons_matrix)
  rownames(regulons_matrix) <- tfs_set
  colnames(regulons_matrix) <- cell_type
  my_minGSSize <- 5
  
  for (i in cell_type) {
    print(i)
    tf_val <- lapply(tfs_set, function(x) {
      if(x %in% detect_gene){
        targets <- target_relation[which(target_relation$TF_Symbol==x),2]
        targets <- targets[targets %in% detect_gene]
        if(length(targets)<=0){
          return(0)
        }
        corGene_tmp <- getCorrelatedGene(data = expr_set,cell_type = i,tf=x, target_list=targets, pValue=pValueCor, corValue=CorValue,topGene=topTargetCor)
        common_targets_tmp <- intersect(corGene_tmp, targets)
        if(length(common_targets_tmp)==0){
          return(0)
        }
        
        gene.name.tmp <- common_targets_tmp
        term_gene_list.tmp <- data.frame(term.name=rep(1, length(gene.name.tmp)), gene=gene.name.tmp)
        
        if(length(gene.name.tmp)<my_minGSSize){
          return(0)
        }
        
        tryCatch({
          nes.tmp <- getGSEA(term_gene_list = term_gene_list.tmp,
                             FC_OF_CELL = fc.list[[i]], minGSSize=my_minGSSize, maxGSSize=500)
          if(length(nes.tmp@result$NES)>0 & length(nes.tmp@result$p.adjust)>0 & expr_mean[x,i]>0){
            #if(nes.tmp@result$p.adjust<p.adjust & nes.tmp@result$NES>0){
            #  tf.val.enriched <- nes.tmp@result$NES
            #  return(tf.val.enriched)
            ##Shouguo changed
            if(nes.tmp@result$p.adjust<p.adjust){
              if(object@expDirection=="UP"){
                if(nes.tmp@result$NES>0){
                  cat("up FOUND\n")
                  tf.val.enriched <- nes.tmp@result$NES
                  return(tf.val.enriched)
                }else{
                  return(0)
                }
              }else if(object@expDirection=="DOWN"){
                if(nes.tmp@result$NES<0){
                  cat("down FOUND\n")
                  tf.val.enriched <- nes.tmp@result$NES * (-1)
                  return(tf.val.enriched)
                }else{
                  return(0)
                }
              }else{
                cat("both FOUND\n")
                tf.val.enriched <- abs(nes.tmp@result$NES)
                return(tf.val.enriched)
              }
              ##Shouguo changed
            }else{
              return(0)
            }
          }else{
            return(0)
          }
        },error=function(e){
          return(0)
        })
        
      }else{
        return(0)
      }
    })
    tf_val <- unlist(tf_val)
    regulons_matrix[,i] <- tf_val
    print(sum(tf_val>0))
  }
  
  gsea.list <- list()
  gsea.genes.list <- list()
  for (i in cell_type) {
    print(i)
    print(length(which(regulons_matrix[,i]!=0)))
    if(length(which(regulons_matrix[,i]!=0))!=0){
      tfs_set.tmp <- tfs_set[which(regulons_matrix[,i]!=0)]
      tf_val.df <- do.call(rbind, lapply(tfs_set.tmp, function(x) {
        targets <- target_relation[which(target_relation$TF_Symbol==x),2]
        targets <- targets[targets %in% detect_gene]
        if(length(targets)<=0){
          return(NULL)
        }
        
        corGene_tmp <- getCorrelatedGene(data = expr_set,cell_type = i,tf=x, target_list=targets, pValue=pValueCor, corValue=CorValue,topGene=topTargetCor)
        common_targets_tmp <- intersect(corGene_tmp, targets)
        
        if(length(common_targets_tmp)==0){
          return(NULL)
        }
        # print(length(common_targets_tmp))
        gene.name.tmp <- common_targets_tmp
        if(length(gene.name.tmp)<my_minGSSize){
          return(NULL)
        }
        
        term_gene_list.tmp <- data.frame(term.name=rep(x, length(gene.name.tmp)), gene=gene.name.tmp)
        return(term_gene_list.tmp)
      }))
      nes.tmp <- getGSEA(term_gene_list = tf_val.df, FC_OF_CELL = fc.list[[i]], minGSSize=my_minGSSize, maxGSSize=500)
      
      gsea.list <- c(gsea.list, list(nes.tmp))
    }else{
      gsea.list <- c(gsea.list, list(vector()))
    }
  }
  names(gsea.list) <- cell_type
  
  
  print("step3: get distance between receptor and tf in pathway")
  
  DistanceKEGG <- getDistanceKEGG(data = triple_relation,method = "mean")
  
  print("step4: score downstream activation of ligand-receptor all regulons of given ligand-receptor relation (weighted, max, or mean) ####")
  l_r_inter <- unique(triple_relation[,5:6])
  expr_r_regulons <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type)) ## A->A,A->B,A->C,,,,C->C
  expr_r_regulons <- as.data.frame(expr_r_regulons)
  rownames(expr_r_regulons) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  colnames(expr_r_regulons) <- cell_type
  
  for (n in 1:nrow(l_r_inter)) {
    sender_tmp <- l_r_inter[n,1]
    receiver_tmp <- l_r_inter[n,2]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")
    # print(n)
    
    val_tmp = 0
    if( sum(l_r_inter[n,] %in% detect_gene)==2 ){
      
      info_tmp <- dplyr::filter(triple_relation, Ligand_Symbol==sender_tmp & Receptor_Symbol==receiver_tmp)[,5:7]
      tfs_tmp <- info_tmp$TF_Symbol[info_tmp$TF_Symbol %in% detect_gene]
      if(length(tfs_tmp) > 0){
        regulon_tmp_df <- regulons_matrix[tfs_tmp,]
        if(method=='max'){
          expr_r_regulons[row_index,] = as.numeric(apply(regulon_tmp_df, 2, function(x){max(x)}))
        }else if(method=="weighted"){
          distance2w_tmp <- (1/DistanceKEGG[row_index,tfs_tmp])
          w_tmp<- distance2w_tmp/sum(distance2w_tmp)
          expr_r_regulons[row_index,] = as.numeric(apply(regulon_tmp_df, 2, function(x){
            sum(w_tmp*x)
          }))
        }else if(method=="mean"){
          expr_r_regulons[row_index,] = as.numeric(apply(regulon_tmp_df, 2, function(x){mean(x)}))
        }
      }
    }
  }
  
  print("step5: softmax for ligand")
  # softmax for ligand
  ligand_symbol <- unique(triple_relation$Ligand_Symbol)
  softmax_ligand <- expr_mean[intersect(ligand_symbol, detect_gene),]
  colnames(softmax_ligand) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_ligand)
  
  softmax_ligand <- do.call(rbind,lapply(1:nrow(softmax_ligand), function(i){
    softmax_ligand[i,]/rowCounts[i]
  }))
  
  # softmax for receptor
  receptor_symbol <- unique(triple_relation$Receptor_Symbol)
  softmax_receptor <- expr_mean[intersect(receptor_symbol, detect_gene),]
  colnames(softmax_receptor) <- colnames(expr_mean)
  rowCounts <- rowSums(softmax_receptor)
  
  softmax_receptor <- do.call(rbind,lapply(1:nrow(softmax_receptor), function(i){
    softmax_receptor[i,]/rowCounts[i]
  }))
  
  #  l-r in cell type level
  print("step6: score ligand-receptor relation (weighted, max, or mean) ####")
  
  l_r_inter <- unique(triple_relation[,5:6])
  expr_l_r <- matrix(data = 0,nrow = nrow(l_r_inter), ncol = length(cell_type)^2) ##A->A,A->B,A->C,,,,C->C
  expr_l_r <- as.data.frame(expr_l_r)
  rownames(expr_l_r) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  myColnames <- character()
  for (i in cell_type) {
    for (j in cell_type) {
      myColnames <- c(myColnames, paste(i,j,sep = "-"))
    }
  }
  colnames(expr_l_r) <- myColnames
  
  for (n in 1:nrow(l_r_inter)) {
    sender_tmp <- l_r_inter[n,1]
    receiver_tmp <- l_r_inter[n,2]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")
    # print(n)
    for (i in cell_type) {
      for (j in cell_type) {
        myColnames <- c(myColnames, paste(i,j,sep = "-"))
        val_tmp = 0
        if( sum(l_r_inter[n,] %in% detect_gene)==2 ){
          sender_val <- expr_mean[sender_tmp,i]
          receiver_val <- expr_mean[receiver_tmp,j]
          tf_val <- expr_r_regulons[row_index,j]
          
          if(tf_val > 0 & sender_val>0 & receiver_val >0){
            sender_val_weighted <- softmax_ligand[sender_tmp, i]
            receiver_val_weighted <- softmax_receptor[receiver_tmp, j]
            val_tmp <- 100*(sender_val_weighted^2 + receiver_val_weighted^2) * tf_val
          }else{
            val_tmp = 0
          }
        }else{
          val_tmp = 0
        }
        col_index_tmp <- paste(i,j,sep = "-")
        expr_l_r[n,col_index_tmp] <- val_tmp
      }
    }
  }
  
  expr_l_r <- expr_l_r[apply(expr_l_r, 1, function(x){sum(x!=0)})>0,]
  expr_l_r <- as.data.frame(expr_l_r)
  
  expr_l_r_log2 <- log2(expr_l_r+1)
  expr_l_r_log2_scale <- (expr_l_r_log2-min(expr_l_r_log2))/(max(expr_l_r_log2)-min(expr_l_r_log2))
  
  Result <- list(expr_mean = expr_mean,
                 regulons_matrix = regulons_matrix,
                 gsea.list = gsea.list,
                 fc.list = fc.list,
                 expr_r_regulons = expr_r_regulons,
                 softmax_ligand = softmax_ligand,
                 softmax_receptor = softmax_receptor,
                 expr_l_r =  expr_l_r,
                 expr_l_r_log2 = expr_l_r_log2,
                 expr_l_r_log2_scale = expr_l_r_log2_scale,
                 DistanceKEGG= DistanceKEGG,
                 object = object)  #SG add object
  
  return(Result)
}










#' get relation among specific TF
#' @param data a dataframe with row of gene and column of sample
#' @param cell_type specific cell type
#' @param tf choose one specific TF
#' @param target_list choose target genes corresponding TF
#' @param pValue set the pValue of spearman
#' @param corValue set the corValue of spearman
#' @param topGene use topTargetCor of candidate genes which has firlter by above parameters, default is 1, means 100%
#' @return result target genes in this regulon
#' @importFrom psych corr.test

getCorrelatedGene <- function(data, cell_type="", tf, target_list,
                              pValue=0.05, corValue=0, topGene=0.1
){
  my_Expr <- data
  expr_tmp <- my_Expr[,which(colnames(my_Expr)==cell_type)]
  if(tf %in% rownames(expr_tmp)){
    # library(psych)
    p_tmp <- psych::corr.test(t(expr_tmp[tf,]),t(expr_tmp[target_list,]),adjust = "none",use = "pairwise", method = "spearman")
    tmp_pValue = p_tmp$p
    tmp_corr = p_tmp$r
    tmp_pValue[is.na(tmp_pValue)]=1  
    tmp_corr[is.na(tmp_corr)]=0  
    corr_gene_tmp <- as.data.frame(t(tmp_corr[1,tmp_pValue<pValue,drop=F]))
    
    if(nrow(corr_gene_tmp)>1){
      corr_gene_tmp <- corr_gene_tmp[corr_gene_tmp>corValue,1,drop=F]
      corr_target_tmp <- rownames(corr_gene_tmp[order(corr_gene_tmp[,1],decreasing = T),1,drop=F])
      topTarget = floor(topGene * length(corr_target_tmp))
      
      if(length(corr_target_tmp)>(topTarget+1)){
        # print(topTarget)
        corr_target_tmp <- corr_target_tmp[2:(topTarget+1)]
      }
    }else{
      corr_target_tmp <- vector()
    }
  }else{
    corr_target_tmp <- vector()
  }
  return(corr_target_tmp)
}





#' abstract path from kegg pathway file
#' @param data a dataframe store L-R-TF relation in KEGG and the distance between LR and TF
#' @param method the method to calculate multiple path corresponding specific L-R-TF, default is "mean"
#' @return the distance between LR and TF
#' @importFrom stringr str_split
#' @importFrom stats median quantile
#' @importFrom magrittr %>%

getDistanceKEGG <- function(data, method){
  # library(stringr)
  triple_relation <- data
  l_r_inter <- unique(triple_relation[,5:6])
  dist_lr_tf <- matrix(data = 0, nrow = nrow(l_r_inter), ncol = length(unique(triple_relation$TF_Symbol)))
  dist_lr_tf <- as.data.frame(dist_lr_tf)
  rownames(dist_lr_tf) <- paste(l_r_inter$Ligand_Symbol, l_r_inter$Receptor_Symbol,sep = "-")
  colnames(dist_lr_tf) <- unique(triple_relation$TF_Symbol)
  
  for (i in 1:nrow(triple_relation)) {
    index_row_tmp <- paste(triple_relation[i,5],triple_relation[i,6],sep = "-")
    index_col_tmp <- triple_relation[i,7]
    tmp_mat <- str_split(triple_relation[i,4],",")[[1]] %>% str_split('_',simplify = T)
    if(method=="mean"){
      mean_dist <- mean(as.numeric(tmp_mat[,2]))
    }else if(method=="median"){
      mean_dist <- median(as.numeric(tmp_mat[,2]))
    }
    dist_lr_tf[index_row_tmp, index_col_tmp] <- mean_dist
  }
  return(dist_lr_tf)
}





#' plot pheatmap graph with communication profile
#' @param object a Cellwave objects
#' @param slot plot the graph with the data of specific slot
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param treeheight_row the height of a tree for rows, if these are clustered. Default value 0 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered. Default value 50 points.
#' @param cluster_rows boolean values determining if rows should be clustered or hclust object,
#' @param cluster_cols boolean values determining if columns should be clustered or hclust object.
#' @param fontsize base fontsize for the plot
#' @param angle_col angle of the column labels, right now one can choose only from few predefined options (0, 45, 90, 270 and 315)
#' @param color vector of colors used in heatmap.
#' @param main the title of the plot, default is "score".
#' @importFrom pheatmap pheatmap
#' @export

viewPheatmap <- function(object, slot="expr_l_r_log2_scale", show_rownames = T, show_colnames = T,
                         treeheight_row=0, treeheight_col=50, cluster_rows = T, cluster_cols = F,
                         fontsize = 12, angle_col = "45", color=NULL, main="score"){
  # library(pheatmap)
  dat <- object@data[[slot]]
  pheatmap::pheatmap(dat, show_rownames = show_rownames, show_colnames = show_colnames,
                     treeheight_row=treeheight_row, treeheight_col=treeheight_col,
                     cluster_rows = cluster_rows, cluster_cols = cluster_cols, fontsize = fontsize, angle_col = angle_col,
                     # color = colorRampPalette(colors = c('#FFFFFF','#ffffb2','#fd8d3c','#e31a1c'))(1000),
                     main=main)
}



#' enrich communication relation on the pathway
#' @param data a dataframe of communication score with row LR and column cellA-cellB
#' @param cella_cellb explore the LR between sender cellA and receiver cellB, eg: "A-B"
#' @param IS_core logical variable ,whether use reference LR data or include extended datasets
#' @param Org choose the species source of gene, only "Homo sapiens" in this version.
#' @return the dataframe with column: Pvalue, Jaccard, NES and pathway
#' @importFrom utils read.table
#' @importFrom stringr str_split
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_longer
#' @importFrom stats phyper sd
#' @export

getHyperPathway <- function(data, object, cella_cellb, IS_core=TRUE, Org="Homo sapiens"){
  n<- data
  Org <- Org
  mt <- object
  
  if(Org == 'Homo sapiens'){
    f.tmp <- system.file("extdata", "rectome_ligand_receptor_TFs_withReactome.txt", package="cellcallEXT")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_extended.txt", package="cellcallEXT")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target.txt", package="cellcallEXT")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }else if(Org == 'Mus musculus'){
    f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology.txt", package="cellcallEXT")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology_extended.txt", package="cellcallEXT")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target_homology.txt", package="cellcallEXT")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }
  
  tmp <- triple_relation
  tmp$triple <- paste(triple_relation$Ligand_Symbol, triple_relation$Receptor_Symbol, triple_relation$TF_Symbol, sep = "-")
  list.tmp <- apply(tmp, 1, function(x){
    str_split(x[5],",",simplify = F) %>% unlist()  -> path.tmp
    tmp.df.part <- data.frame(path=path.tmp, triple=rep(as.character(x[9]), length(path.tmp)))
  })
  path.triple.tmp.df <- do.call(rbind, list.tmp)
  
  path.list.tmp <- lapply(rownames(n), function(x){
    names.tmp <- colnames(n)[n[x,]>0]
    
    if(length(names.tmp)>0){
      # print(x)
      l.tmp <- str_split(x,"-")[[1]][1]
      r.tmp <- str_split(x,"-")[[1]][2]
      a.tmp <- triple_relation %>% dplyr::filter(Ligand_Symbol==l.tmp & Receptor_Symbol==r.tmp) %>% .[,c('pathway_ID','TF_Symbol')]
      recevier.tmp <- str_split(names.tmp, '-', simplify = T)[,2]
      df.tmp <- mt@data$regulons_matrix[a.tmp$TF_Symbol, recevier.tmp, drop=F]
      colnames(df.tmp) <- names.tmp
      df.tmp$path <- a.tmp$pathway_ID
      df.tmp$triple <- paste(x, a.tmp$TF_Symbol, sep = '-')
      df.longer.tmp <- df.tmp %>% tidyr::pivot_longer(cols = names.tmp, names_to = "cc", values_to = "count") %>% as.data.frame()
      df.longer.tmp <- df.longer.tmp[df.longer.tmp$count>0,]
      list.tmp <- apply(df.longer.tmp, 1, function(x){
        str_split(x[1],",",simplify = F) %>% unlist()  -> path.tmp
        tmp.df.part <- data.frame(triple=as.character(x[2]), path=path.tmp, cc=rep(as.character(x[3]), length(path.tmp)), count=rep(as.character(x[4]), length(path.tmp)))
      })
      list.tmp.df <- do.call(rbind, list.tmp)
      return(list.tmp.df)
    }
  })
  path.list.tmp.df <- do.call(rbind, path.list.tmp)
  
  apply(path.list.tmp.df, 1, function(x){
    row_index_tmp <- paste(str_split(x[1], '-')[[1]][1], str_split(x[1], '-')[[1]][2], sep = '-')
    col_index_tmp <- as.character(x[3])
    return(n[row_index_tmp, col_index_tmp])
  }) -> score
  path.list.tmp.df$lrscore <- score
  colnames(path.list.tmp.df)[4] <- "enrichment_score"
  
  cc.tmp <- cella_cellb
  cc.tmp.triple <- unique(dplyr::filter(path.list.tmp.df, cc==cc.tmp)[,1:2,drop=T] %>% as.data.frame())
  if(nrow(cc.tmp.triple)==0){
    res.df <- data.frame(Pvalue=1, Jaccard=0, NES=0, pathway="hsa04330")
    rownames(res.df) <- "hsa04330"
    return(res.df)
  }
  
  myHypergeometric <- function(intersect_mirna=intersect_mirna, mrna_has_mirna=mrna_has_mirna, all_mirna=all_mirna, lncrna_has_mirna=lncrna_has_mirna){
    q = intersect_mirna-1 
    m = mrna_has_mirna  
    n = all_mirna-mrna_has_mirna 
    k = lncrna_has_mirna 
    stats::phyper(q=q, m=m, n=n, k=k, log = FALSE, lower.tail = FALSE)
  }
  
  cc.tmp.triple$triple <- as.character(cc.tmp.triple$triple)
  cc.tmp.triple$path <- as.character(cc.tmp.triple$path)
  
  cc.tmp.pathway <- unique(cc.tmp.triple$path)
  do.call(rbind, lapply(cc.tmp.pathway, function(x){
    unique(dplyr::filter(path.triple.tmp.df, path==x)[,2]) -> pathway.triple.tmp
    unique(dplyr::filter(cc.tmp.triple, path==x)[,1]) -> cc.triple.tmp
    intersect_mirna <- length(intersect(pathway.triple.tmp, cc.triple.tmp))
    mrna_has_mirna <- length(pathway.triple.tmp)
    all_mirna <- length(unique(path.triple.tmp.df$triple))
    lncrna_has_mirna <- length(unique(cc.tmp.triple$triple))
    
    Hyper.pValue.tmp <- myHypergeometric(intersect_mirna = intersect_mirna,
                                         mrna_has_mirna = mrna_has_mirna,
                                         all_mirna = all_mirna,
                                         lncrna_has_mirna = lncrna_has_mirna)
    Jaccard.tmp <- length(intersect(pathway.triple.tmp, cc.triple.tmp))/length(unique(c(pathway.triple.tmp, cc.triple.tmp)))
    return(c(Hyper.pValue.tmp, Jaccard.tmp))
  })) -> res.df
  res.df <- as.data.frame(res.df)
  colnames(res.df) <- c("Pvalue", "Jaccard")
  rownames(res.df) <- cc.tmp.pathway
  
  all.path.tmp <- unique(path.triple.tmp.df$path)
  all.path.Jaccard.tmp <- do.call(rbind, lapply(all.path.tmp, function(x){
    unique(dplyr::filter(path.triple.tmp.df, path==x)[,2]) -> pathway.triple.tmp
    unique(dplyr::filter(cc.tmp.triple, path==x)[,1]) -> cc.triple.tmp
    
    Jaccard.tmp <- length(intersect(pathway.triple.tmp, cc.triple.tmp))/length(unique(c(pathway.triple.tmp, cc.triple.tmp)))
    return(Jaccard.tmp)
  })) %>% unlist()
  
  mean.tmp <- base::mean(all.path.Jaccard.tmp) 
  sd.tmp <- stats::sd(all.path.Jaccard.tmp)
  res.df$NES <- (res.df$Jaccard-mean.tmp)/sd.tmp
  
  res.df$pathway <- rownames(res.df)
  
  return(res.df)
}
#' calculate foldchange for each celltype
#' @param inData a dataframe of gene expression
#' @param cell.type the cell type which you want to calculate the foldchange
#' @param method default is "median",eg "median" or "mean"
#' @param probs the quantile for median calculate
#' @return the foldchange value of each celltype in \code{cellwave object}
#' @importFrom stats median

mylog2foldChange.diy<-function(inData, cell.type, method="median", probs = 0.75) # median or mean
{
  re.list <- list()
  cell.type <- cell.type
  print(cell.type)
  for (c in cell.type) {
    print(c)
    cell_fc_labels<-colnames(inData)
    cell_fc_labels[cell_fc_labels==c]<-0
    cell_fc_labels[cell_fc_labels!='0']<-1
    
    classLabel <- cell_fc_labels
    sampleIdsCase<-which(classLabel==0);#0 tumer
    sampleIdsControl<-which(classLabel==1);#1 normal
    probeFC <- as.numeric(apply(inData, 1, function(x){
      # print(x)
      if(method=="mean"){
        (mean(as.numeric(x[sampleIdsCase]))+1)/(mean(as.numeric(x[sampleIdsControl]))+1);
      }else if(method=="median"){
        quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)
        
        a.tmp <- quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)+1
        b.tmp <- quantile(as.numeric(x[sampleIdsControl]), probs = probs, names=FALSE)+1
        a.tmp/b.tmp
      }
    }))
    
    probeFC<-log(probeFC,base=2);
    fc_res<-probeFC;
    
    fc_res[is.infinite(fc_res)]<-0
    fc_res[is.na(fc_res)]<-0
    res = data.frame(gene_id = as.vector(rownames(inData)), log2fc = fc_res)
    # print(dim(res))
    re.list <- c(re.list, list(res))
  }
  names(re.list) <- cell.type
  return(re.list)
}

#' calculate foldchange between case and control for each celltype
#' @param inData a dataframe of gene expression
#' @param cell.type the cell type which you want to calculate the foldchange
#' @param method default is "median",eg "median" or "mean"
#' @param probs the quantile for median calculate
#' @return the foldchange value of each celltype in \code{cellwave object}
#' @importFrom stats median
#####Shouguo Gao  Add one parameters status
mylog2foldChange.diy.casecontrol<-function(inData, cell.type, method="median", probs = 0.75, status = NULL) # median or mean
{
  re.list <- list()
  cell.type <- cell.type
  print(cell.type)
  for (c in cell.type) {
    print(c)
    cell_fc_labels<-colnames(inData)
    #cell_fc_labels[cell_fc_labels==c]<-0
    #cell_fc_labels[cell_fc_labels!='0']<-1
    cell_fc_labels[cell_fc_labels==c & status =="CASE"]<-0
    cell_fc_labels[cell_fc_labels==c & status =="CONTROL"]<-1
    
    classLabel <- cell_fc_labels
    sampleIdsCase<-which(classLabel==0);#0 tumer
    sampleIdsControl<-which(classLabel==1);#1 normal
    probeFC <- as.numeric(apply(inData, 1, function(x){
      # print(x)
      if(method=="mean"){
        (mean(as.numeric(x[sampleIdsCase]))+1)/(mean(as.numeric(x[sampleIdsControl]))+1);
      }else if(method=="median"){
        quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)
        
        a.tmp <- quantile(as.numeric(x[sampleIdsCase]), probs = probs, names=FALSE)+1
        b.tmp <- quantile(as.numeric(x[sampleIdsControl]), probs = probs, names=FALSE)+1
        a.tmp/b.tmp
      }
    }))
    
    probeFC<-log(probeFC,base=2);
    fc_res<-probeFC;
    
    fc_res[is.infinite(fc_res)]<-0
    fc_res[is.na(fc_res)]<-0
    res = data.frame(gene_id = as.vector(rownames(inData)), log2fc = fc_res)
    # print(dim(res))
    re.list <- c(re.list, list(res))
  }
  names(re.list) <- cell.type
  return(re.list)
}

#' calculate enrich for each celltype with all regulons
#' @param term_gene_list all regulons of genes denotes geneSet
#' @param FC_OF_CELL the cell type which you want to calculate the foldchange
#' @param minGSSize set the min size of each geneSet
#' @param maxGSSize set the max size of each geneSet
#' @return the enrich result of each celltype with all regulons
#' @importFrom clusterProfiler GSEA
#' @importFrom dplyr arrange desc

getGSEA<-function(term_gene_list, FC_OF_CELL, minGSSize=5, maxGSSize=500) {
  term_gene <- term_gene_list
  # library(dplyr)
  FC_OF_CELL <- FC_OF_CELL
  geneList.sort <- dplyr::arrange(FC_OF_CELL, desc(FC_OF_CELL$log2fc))
  
  rownames(geneList.sort) <- geneList.sort[,1]
  gene.name <- geneList.sort[,1]
  geneList.sort$log2fc <- as.numeric(geneList.sort$log2fc)
  geneList<-geneList.sort[,-1]
  names(geneList) <- gene.name
  
  # if most fc are not high enough will to error like: GSEA Error in if (abs(max.ES) > abs(min.ES))
  # you can set exponent = 0
  y<-GSEA(geneList = geneList, TERM2GENE = term_gene, minGSSize = minGSSize,exponent = 1,
          maxGSSize = maxGSSize, pvalueCutoff = 1, pAdjustMethod = "BH",
          verbose = FALSE)
  return(y)
}







#' get score of Triple relation corresponding specific cellA-cellB
#' @param object a Cellwave objects
#' @param sender_cell the cell type of sender cell
#' @param recevier_cell the cell type of recevier cell
#' @param slot plot the graph with the data of specific slot
#' @param org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @param IS_core logical variable ,whether use reference LR data or include extended datasets
#' @importFrom utils read.table head
#' @importFrom stringr str_split
#' @importFrom dplyr filter
#' @export

LR2TF <- function(object, sender_cell, recevier_cell, slot="expr_l_r_log2_scale", org="Homo sapiens", IS_core=TRUE){
  options(stringsAsFactors = F) 
  
  myData <- object@data[[slot]]
  
  detect_gene <- rownames(object@data$expr_mean)
  regulons_matrix <- object@data$regulons_matrix
  DistanceKEGG <- object@data$DistanceKEGG
  expr_mean <- object@data$expr_mean
  sender_cell <- sender_cell
  recevier_cell <- recevier_cell
  # top_n <- top_n
  
  if(org == 'Homo sapiens'){
    f.tmp <- system.file("extdata", "rectome_ligand_receptor_TFs_withReactome.txt", package="cellcallEXT")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_extended.txt", package="cellcallEXT")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target.txt", package="cellcallEXT")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }else if(org == 'Mus musculus'){
    f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology.txt", package="cellcallEXT")
    triple_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
    
    if(IS_core){
    }else{
      f.tmp <- system.file("extdata", "new_ligand_receptor_TFs_homology_extended.txt", package="cellcallEXT")
      triple_relation_extended <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
      triple_relation <- rbind(triple_relation, triple_relation_extended)
    }
    
    f.tmp <- system.file("extdata", "tf_target_homology.txt", package="cellcallEXT")
    target_relation <- read.table(f.tmp, header = TRUE, quote = "", sep = '\t', stringsAsFactors=FALSE)
  }
  
  interest_group <- myData[,paste(sender_cell, recevier_cell, sep = "-"), drop=F]
  interest_df <- interest_group[order(interest_group,decreasing = T),,drop=F]
  interest_df <- interest_df[which(interest_df[,1]!=0),,drop=F]
  
  # library(stringr)
  my_ligand_receptor_df <- str_split(rownames(interest_df),"-",simplify = T)
  my_tfSet <- c()
  for (n in 1:nrow(interest_df)) {
    sender_tmp <- my_ligand_receptor_df[n,1]
    receiver_tmp <- my_ligand_receptor_df[n,2]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")
    
    info_tmp <- dplyr::filter(triple_relation, Ligand_Symbol==sender_tmp & Receptor_Symbol==receiver_tmp)[,c("Ligand_Symbol", "Receptor_Symbol", "TF_Symbol")]
    tfs_tmp <- info_tmp$TF_Symbol[info_tmp$TF_Symbol %in% detect_gene]
    tfs_tmp <- tfs_tmp[regulons_matrix[tfs_tmp, recevier_cell]>0]
    my_tfSet <- c(my_tfSet, tfs_tmp)
  }
  my_tfSet <- unique(my_tfSet)
  
  sankey_matrix <- matrix(0, ncol = 5)
  l_r_tf <- matrix(0, ncol = length(my_tfSet), nrow = nrow(interest_df))
  l_r_tf <- data.frame(l_r_tf)
  rownames(l_r_tf) <- rownames(interest_df)
  colnames(l_r_tf) <- my_tfSet
  for (n in 1:nrow(interest_df)) {
    print(n)
    sender_tmp <- my_ligand_receptor_df[n,1]
    receiver_tmp <- my_ligand_receptor_df[n,2]
    sender_val <- expr_mean[sender_tmp,sender_cell]
    receiver_val <- expr_mean[receiver_tmp,recevier_cell]
    row_index <- paste(sender_tmp, receiver_tmp,sep = "-")
    
    info_tmp <- dplyr::filter(triple_relation, Ligand_Symbol==sender_tmp & Receptor_Symbol==receiver_tmp)[,c("Ligand_Symbol", "Receptor_Symbol", "TF_Symbol")]
    my_tfSet_tmp <- info_tmp$TF_Symbol[info_tmp$TF_Symbol %in% detect_gene]
    my_tfSet_tmp <- my_tfSet_tmp[regulons_matrix[my_tfSet_tmp, recevier_cell]>0]
    
    tfVal_tmp <- regulons_matrix[my_tfSet_tmp,recevier_cell]
    distance2w_tmp <- (1/DistanceKEGG[row_index,my_tfSet_tmp])
    w_tmp<- distance2w_tmp/sum(distance2w_tmp)
    l_r_val_tmp <- myData[row_index,paste(sender_cell, recevier_cell, sep = "-")]
    
    # l_r_tf[row_index,my_tfSet_tmp] <- as.numeric(tfVal_tmp * w_tmp * l_r_val_tmp)
    sankey_tmp <- str_split( paste(sender_tmp, receiver_tmp, my_tfSet_tmp, l_r_val_tmp, tfVal_tmp * w_tmp )," ",simplify = TRUE)
    sankey_matrix <- rbind(sankey_matrix, sankey_tmp)
    
  }
  # l_r_tf <- l_r_tf[,apply(l_r_tf, 2, function(x){sum(x!=0)})>0] 
  sankey_matrix <- data.frame(sankey_matrix)
  colnames(sankey_matrix) <- c("Ligand", "Receptor",	"TF",	"weight1",	"weight2")
  sankey_matrix <- sankey_matrix[-1,]
  sankey_matrix$weight1 <- as.numeric(sankey_matrix$weight1)
  sankey_matrix$weight2 <- as.numeric(sankey_matrix$weight2)
  
  # object@reductions$LR_TF <- l_r_tf
  object@reductions$sankey <- sankey_matrix
  return(object)
}




#' plot enrichment graph
#' @param gsea.list a list of enrichment result from cellwave
#' @param myCelltype the cell type of receiver cell
#' @param fc.list foldchange list in the cellwave object
#' @param geneSetID the character of TF symbol, only significant activated can be inspected
#' @param selectedGeneID default is NULL, label the position of specific gene in FC flow.
#' @param mycol the color of each TF. the length is consistent with geneSetID
#' @importFrom enrichplot plot_grid
#' @importFrom ggplot2 ggplot geom_line scale_color_manual aes_ xlab ylab geom_hline theme_bw theme element_blank element_rect element_text margin geom_linerange scale_y_continuous geom_rect scale_fill_gradientn geom_segment geom_bar scale_fill_manual unit element_line
#' @importFrom ggrepel geom_text_repel
#' @export

getGSEAplot <- function(gsea.list, myCelltype, fc.list, geneSetID, selectedGeneID = NULL,
                        mycol = NULL)
{
  options(stringsAsFactors = FALSE)
  fc.df <- fc.list[[myCelltype]]
  
  fc.df.sorted <- fc.df[order(fc.df$log2fc, decreasing = T),]
  
  gsea.object <- gsea.list[[myCelltype]]
  
  if(is.null(mycol)){
    mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9",
               "#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D",
               "#7CC767")
  }
  x <- gsea.object
  geneList <- position <- NULL ## to satisfy codetool
  
  gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
  gsdata$gsym <- rep(fc.df$gene_id, length(geneSetID))
  
  p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
    scale_color_manual(values = mycol) +
    geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
    ylab("Enrichment\n Score") +
    theme_bw() +
    theme(panel.grid = element_blank())
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) + 
    theme(axis.text.y=element_text(size = 12, face = "bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b= 0, l=.2, unit="cm"))
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) +
    scale_color_manual(values = mycol) +
    
    theme_bw() +
    theme(panel.grid = element_blank()) +
    
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_y_continuous(expand=c(0,0))
  
  test <- gsdata[1:floor(nrow(gsdata)/length(geneSetID)),]
  test$xend = test$x+1
  
  p.grad <- ggplot(test)+ geom_rect(aes(xmin = x,xmax = xend , ymin = 0 , ymax = 1, fill=geneList))+
    scale_fill_gradientn(colours = c("#035BFD", "#397EFC", "#5B94FB","white", "#F77A7C", "#F45557", "#FB0407"), limits = c(-max(abs((test$geneList))),max(abs((test$geneList)))))+
    theme_bw() +
    theme(panel.grid = element_blank()) +
    
    theme(legend.position = "none",
          plot.margin = margin(t=0, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_y_continuous(expand=c(0,0))
  
  df2 <- p.res$data
  df2$y <- p.res$data$geneList[df2$x]
  df2$gsym <- p.res$data$gsym[df2$x]
  
  if(!is.null(selectedGeneID)){
    selectgenes <- data.frame(gsym = selectedGeneID)
    selectgenes <- merge(selectgenes, df2, by = "gsym")
    selectgenes <- unique(selectgenes[,c("x","y","gsym")])
    head(selectgenes)
    
    p.pos <- ggplot(selectgenes, aes(x, y, fill = "black", color = "black", label = gsym)) +
      geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                   color = "#80b1d3") +
      geom_bar(position = "dodge", stat = "identity",  width=0.5) +
      scale_fill_manual(values = "black", guide=FALSE) + 
      scale_color_manual(values = "black", guide=FALSE) +
      
      #scale_x_continuous(expand=c(0,0)) +
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
      ylab("Ranked list\n metric") +
      xlab("Rank in ordered dataset") +
      
      theme_bw() +
      theme(axis.text.y=element_text(size = 12, face = "bold"),
            panel.grid = element_blank()) +
      
      geom_text_repel(data = selectgenes,
                      show.legend = FALSE,
                      direction = "x", 
                      ylim = c(2, NA), 
                      angle = 90, 
                      size = 2.5, box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.3, "lines")) +
      theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  }else{
    p.pos <- ggplot() + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                                     color = "#80b1d3") +
      
      #scale_x_continuous(expand=c(0,0)) +
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
      ylab("Ranked list\n metric") +
      xlab("Rank in ordered dataset") +
      
      theme_bw() +
      theme(axis.text.y=element_text(size = 12, face = "bold"),
            panel.grid = element_blank())+  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  }
  
  rel_heights <- c(1.5, .5, .2, 1.5)
  plotlist <- list(p.res, p2, p.grad, p.pos)
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          axis.text.x = element_text(size = 12, face = "bold"))
  
  plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)
}




#' plot bubble graph
#' @param dat a dataframe of melted communication score with columns: pathway, cc, pvalue, NES, logPvalue
#' @importFrom ggplot2 ggplot aes geom_point scale_colour_gradientn scale_size labs theme element_text element_rect
#' @export

plotBubble <- function(dat){
  myPub.df <- dat
  myPub.df$logPvalue <- -log10(myPub.df$pvalue)
  myPub.df$logPvalue[ myPub.df$logPvalue>=2 ] = 2
  myPub.df$logPvalue[ myPub.df$logPvalue>=1.3 & myPub.df$logPvalue<2 ] = 1.3
  myPub.df$logPvalue[ myPub.df$logPvalue<1.3 ] = 1
  
  p = ggplot(myPub.df, aes(x = cc, y = pathway))
  
  pbubble = p + geom_point(aes(size=logPvalue,colour=NES))
  
  
  pr = pbubble+scale_colour_gradientn(colours=c('#1400FE', '#009FFE', '#D00B5D', '#FB6203','#F60404'),
                                      values = c(0,1/4,1/2,3/4,1)) +
    scale_size("Pvalue", breaks=c(1,1.3,2), range=c(1,5),
               labels = c("0.1", "0.05", "<0.01")) 
  
  pr = pr + labs(color="NES",
                 size="-log10 p-value",
                 x="CC",y="Pathway",title="")
  
  pr = pr + theme(axis.text.x = element_text(size = 8, color = "black", face = "plain", vjust = 1.0, hjust = 1, angle = 45),
                  panel.background = element_rect(fill = "white", colour = "black", size = 1))
  return(pr)
}

#' get score of Triple relation corresponding specific cellA-cellB
#' @param df a dataframe consist of five columns: "Ligand", "Receptor", "TF", "weight1", "weight2"
#' @return a dataframe consist of four columns: "Ligand", "Receptor", "TF", "weight"
#' @importFrom dplyr filter
#' @export

trans2tripleScore <- function(df){
  colnames(df) <- c("Ligand", "Receptor", "TF", "weight1", "weight2")
  df <- dplyr::filter(df, weight1 !=0 & weight2 !=0)
  df$Ligand <- paste("sender:", df$Ligand, sep = "")
  df$Receptor <- paste("receiver:", df$Receptor, sep = "")
  df$TF <- paste("TF:", df$TF, sep = "")
  a <- df
  sum_tf <- aggregate(a[,5],by=list(a$Ligand,a$Receptor),FUN=sum) 
  colnames(sum_tf) <- c("Ligand", "Receptor", 'Score')
  
  aa <- c()
  for (i in 1:nrow(a)) {
    tmp = dplyr::filter(sum_tf, Ligand == a[i,1] & Receptor == a[i,2])
    res_tmp = as.numeric(a[i,5])/as.numeric(tmp[,3][1])
    aa <- c(aa, res_tmp)
  }
  
  a$weight2 <- as.numeric(aa)
  a$value <- a$weight1*a$weight2
  a <- a[,c(1,2,3,6)]
  
  return(a)
}
#' The CellInter Class
#' The CellInter object is the center of each cell-cell communication analysis with the scRNA-seq.
#' @slot data  The normalized expression matrix, mean matrix, gsea result, lr score matrix, and lr score matrix(log-scale)
#' @slot meta.data cell type, barcode, nFeature and count of each cell
#' @slot reductions data.frame of L-R-TF with specific cellA-cellB, later for sankey plot
#' @slot project the name of this analysis
#' @slot Org the species of this scRNA-seq
#' @slot version the version of this R package
#' @slot status the disease status of different cells
#' @slot expDirection the direction of expression change UP DOWN or BOTH
#' @import methods
#' @name CellInter-class
#' @aliases CellInter-class
#' @exportClass CellInter
#' @export
#'
CellInter <- setClass("CellInter", slots = list(data = "list",
                                                meta.data = "data.frame",
                                                reductions = "list",
                                                project = "character",
                                                Org = "character",
                                                status = "character", ##SG   add onre more parameters
                                                expDirection = 'character', ##SG   add onre more parameters
                                                version = "character"),
                      
                      prototype = list(data = list(),
                                       meta.data = data.frame(),
                                       reductions = list(),
                                       project = "Microenvironment",
                                       Org = "Hsapiens",
                                       expDirection = "UP", ##SG   add onre more parameters
                                       status = c("CASE", "CONTROL"), ##SG   add onre more parameters
                                       version = "1.0")
)


#' create a Cellwave objects
#' @param data a dataframe with row of gene and column of sample
#' @param min.feature Include cells where at least this many features are detected
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the cell's column name. E.g. If your cells are named as BARCODE_CELLTYPE, set this to "_" to separate the cell name into its component parts for picking the relevant field.
#' @param names.field BARCODE_CELLTYPE in the input matrix, set names.field to 2 to set the initial identities to CELLTYPE.
#' @param project Project name for the this object
#' @param source the type of expression dataframe, eg "UMI", "fullLength", "TPM", or "CPM"
#' @param scale.factor set the scale factor, default "10^6"
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @param status the status of cell disease or control
#' @param expDirection the direction of expression change UP DOWN or BOTH
#' @return the value of \code{cellwave object}
#' @import graphics
#' @import methods
#' @importFrom stringr str_split
#' @export

CreateNichConObject <- function(data,
                                min.feature = 3,
                                names.field = 1,
                                names.delim = "_",
                                project = "Microenvironment",
                                source = "UMI", 
                                scale.factor = 10^6,
                                status = "status", ##SG   add onre more parameters
                                expDirection = "UP", ##SG   add onre more parameters
                                Org = "Homo sapiens"
)
{
  
  nichcon <- new("CellInter",project = "Microenvironment",Org = Org, status = status, expDirection =expDirection) ##SG   add onre more parameters
  
  if(!is.data.frame(data)){
    stop("data must be a dataframe.")
  }
  
  sampleID <- colnames(data)
  if(sum(duplicated(sampleID))>0){
    stop("The cell ID should be unique.")
  }
  
  cell_type <- stringr::str_split(colnames(data), names.delim, simplify = T)[,names.field]
  if(length(grep("-",unique(cell_type)))>0){
    stop("The cell ID can't contain '-'.")
  }
  if(length(grep("_",unique(cell_type)))>0){
    stop("The cell ID can't contain '_'.")
  }
  
  nFeature <- apply(data, 2, function(x) {sum(x>0)})
  nCounts <- colSums(data)
  
  init.meta.data <- data.frame(sampleID = sampleID, celltype = cell_type, nFeature = nFeature, nCounts = nCounts)
  
  if(source=="UMI"){
    data_cpm <- counts2normalized_10X(data, toType = "CPM", scale.factor=scale.factor)
  }else if(source=="fullLength"){
    data_cpm <- counts2normalized_smartseq2(data, Org, "TPM", scale.factor=scale.factor)
  }else if(source=="TPM"){
    data_cpm <- data
  }else if(source=="CPM"){
    data_cpm <- data
  }
  
  data.list = list()
  # data.list = c(data.list, list(count = data, data = data_cpm_log, withoutlog = data_cpm))
  data.list = c(data.list, list(count = data, withoutlog = data_cpm))
  
  nichcon@meta.data <- init.meta.data
  nichcon@data <- data.list
  
  return(nichcon)
}

#' get CommuProfile from a Cellwave objects
#' @param object a Cellwave objects
#' @param probs Percentile of gene expression in one cell type to represents this cell type
#' @param use.type the type of compute, default is "median"
#' @param pValueCor firlter target gene of TF with spearson, p > pValueCor, default is 0.05
#' @param CorValue firlter target gene of TF with spearson, value > CorValue, default is 0.1
#' @param topTargetCor use topTargetCor of candidate genes which has firlter by above parameters, default is 1, means 100%
#' @param p.adjust gsea pValue of regulons with BH adjusted threshold, default is 0.05
#' @param method "weighted", "max", "mean", of which "weighted" is default. choose the proper method to score downstream activation of ligand-receptor all regulons of given ligand-receptor relation
#' @param IS_core logical variable ,whether use reference LR data or include extended datasets
#' @param Org choose the species source of gene, eg "Homo sapiens", "Mus musculus"
#' @return the value of \code{cellwave object}
#' @import graphics
#' @import methods
#' @export

TransCommuProfile <- function(object,
                              pValueCor = 0.05,
                              CorValue = 0.1,
                              topTargetCor=1,
                              p.adjust=0.05,
                              use.type="median",
                              probs = 0.75,
                              method="weighted",  
                              Org = 'Homo sapiens',
                              IS_core = TRUE
)
{
  is_myObject <- is(object, "CellInter")
  if(is_myObject){
    profile <- ConnectProfile(object,
                              pValueCor = pValueCor,
                              CorValue = CorValue,
                              topTargetCor = topTargetCor,
                              p.adjust=p.adjust,
                              use.type=use.type,
                              probs = probs,
                              method = method,
                              IS_core = IS_core,
                              Org = Org)
    
    object@data = c(object@data, profile)
    return(object)
  }else{
    stop("object should be CellInter type.")
  }
}
#' plot circle graph with communication profile
#' @param object a Cellwave objects
#' @param font the size of font
#' @param cellColor a color dataframe, rownames is cell type, value is color
#' @param lrColor a color vector denotes the color of ligand and receptor, containing two elements, default is c('#D92E27', "#35C6F4")
#' @param order.vector default is null, a celltype vector with the order you want in the circle graph
#' @param trackhight1 Height of the outer track
#' @param trackhight2 Height of the inner track
#' @param linkcolor.from.sender logical value, whether the color of line correspond with color of sender cell
#' @param linkcolor one color you want link to be, only if parameter linkcolor.from.sender=FALSE
#' @param arr.type Type of the arrows, default value is big.arrow There is an additional option triangle
#' @param arr.length Length of the arrows, measured in 'cm'. If arr.type is set to big.arrow, the value is percent to the radius of the unit circle.
#' @param DIY logical value, if TRUE, the parameter object should be a dataframe, and set slot="expr_l_r_log2_scale". otherwise object should be a Cellwave objects.
#' @param slot plot the graph with the data of specific slot
#' @param gap.degree between two neighbour sectors. It can be a single value or a vector. If it is a vector, the first value corresponds to the gap after the first sector.
#' @param track.margin2 affect current track
#' @importFrom grid pushViewport unit upViewport viewport gpar
#' @importFrom graphics plot.new par
#' @importFrom gridBase gridOMI
#' @importFrom circlize circos.clear circos.par circos.initialize circos.trackPlotRegion circos.link get.cell.meta.data highlight.sector colorRamp2 rand_color
#' @importFrom stringr str_split
#' @importFrom magrittr %>% set_colnames
#' @importFrom dplyr filter
#' @importFrom ComplexHeatmap draw Legend packLegend
#' @export

ViewInterCircos <- function(object, font = 2, cellColor ,lrColor = NULL, order.vector = NULL,
                            trackhight1 = 0.05, linkcolor.from.sender = TRUE, linkcolor = NULL,
                            arr.type = "big.arrow",arr.length = 0.04, DIY = TRUE, gap.degree = NULL,
                            trackhight2 = 0.032, track.margin2 = c(0.01,0.12),slot="expr_l_r_log2_scale")
{
  
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size,
                        just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  
  
  circos.clear()  
  
  # cell.padding: the padding between sector and other sector
  
  if(is.null(gap.degree)){
    circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.01,0,0.01,0))
  }else{
    circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.01,0,0.01,0), gap.degree=gap.degree)
  }
  
  # library(stringr)
  if(DIY){
    a <- colSums(object)
  }else{
    a <- colSums(object@data[[slot]])
  }
  
  
  b <- stringr::str_split(names(a), "-", simplify = T)
  c <- data.frame(b, stringsAsFactors = FALSE)
  c$x3 <- as.numeric(a)
  test <- c
  colnames(test) <- c( "cell_from", "cell_to", "n")
  
  test$id1 <- paste("sender", test$cell_from, sep = "_", 1:nrow(test))
  test$id2 <- paste("recevier", test$cell_to, sep = "_", 1:nrow(test))
  
  a_tmp <- test[,c(1,4)]
  colnames(a_tmp) <- c('celltype', 'id')
  a_tmp$arrowType <-  "sender"
  
  b_tmp <- test[,c(2,5)]
  colnames(b_tmp) <- c('celltype', 'id')
  b_tmp$arrowType <-  "recevier"
  
  ab_tmp <- rbind(a_tmp,b_tmp)
  
  ab_tmp <- data.frame(ab_tmp, stringsAsFactors = FALSE)
  ab_tmp <- ab_tmp[order(ab_tmp$celltype,ab_tmp$id,decreasing = TRUE),]
  
  sector_id <- ab_tmp$id
  fa = ab_tmp$id
  if(is.null(order.vector)){
    fa = factor(fa,levels = fa)
  }else{
    fa.df <- str_split(fa, "_", simplify = T) %>% as.data.frame() %>% magrittr::set_colnames(c('sender_or_receiver', 'clltype', "index_number"))
    my.levels <- do.call(rbind,lapply(order.vector, function(x){
      fa.df[which(fa.df$clltype==x),]
    })) %>% apply(1, function(x){
      paste(x[1],x[2],x[3],sep='_')
    }) %>% unlist %>% as.character()
    
    fa = factor(fa,levels = my.levels)
  }
  circos.initialize(factors = fa, xlim = c(0,1)) # 
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = trackhight1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
      # print(sector.index)
      # print(xlim)
      # print(ylim)
    }
  )
  
  if(is.null(cellColor)){
    cell_color <- data.frame(color = rand_color(length(unique(ab_tmp$celltype)), luminosity = "light"), stringsAsFactors = FALSE)
    rownames(cell_color) <- unique(ab_tmp$celltype)
  }else{
    cell_color <- cellColor
  }
  cell_type <- unique(ab_tmp$celltype)
  for(i in 1:length(cell_type)){
    myCell_Type <-  cell_type[i]
    mySector_id <- as.character(ab_tmp[ab_tmp$celltype==myCell_Type,'id'])
    myColor <- as.character(cell_color[myCell_Type,1])
    highlight.sector(mySector_id, track.index = 1,
                     text = myCell_Type, text.vjust = -1,niceFacing = T, font = font, col = myColor)
    
  }
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = trackhight2,
    bg.border = NA,
    track.margin = track.margin2,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  
  if(is.null(lrColor)){
    ligand_receptor <- c('#D92E27', "#35C6F4")  #color of sender and receiver
  }else{
    ligand_receptor <- lrColor
  }
  
  for(i in 1:length(cell_type)){
    myCell_Type <-  cell_type[i]
    mySector_id_tmp <- ab_tmp[ab_tmp$celltype==myCell_Type,c('arrowType', 'id')]
    my_ligand <- dplyr::filter(mySector_id_tmp, arrowType == 'sender')
    my_receptor <- dplyr::filter(mySector_id_tmp, arrowType == 'recevier')
    
    if(length(as.character(my_ligand[,'id']))>0){
      highlight.sector(as.character(my_ligand[,'id']), track.index = 2,
                       text = '', niceFacing = F,col = ligand_receptor[1],text.col = 'white')
      
    }
    
    if(length(as.character(my_receptor[,'id']))>0){
      highlight.sector(as.character(my_receptor[,'id']), track.index = 2,
                       text = '', niceFacing = F,col = ligand_receptor[2],text.col = 'white')
      
    }
    
  }
  
  
  test$weighted_n <- (test$n-min(test$n))/(max(test$n)-min(test$n))
  test <- test[order(test$weighted_n, decreasing = F),]
  # print(test$weighted_n)
  # test$weighted_n <- as.numeric(test$n/sum(test$n))
  
  min_n <- min(test$weighted_n)
  max_n <- max(test$weighted_n)
  mean_n <- mean(min_n, max_n)
  
  if(linkcolor.from.sender){
    
    color.function <- apply(cell_color, 1, function(x){
      col_fun = colorRamp2(c(0, 1), c("#FFFFFF", x))
    })
    
    for(i in 1:nrow(test)){
      my_Line_color = as.character(cell_color[test[i,1],1])
      print(test$weighted_n[i])
      circos.link(sector.index1 = test[i,'id1'],
                  point1 = c(0,1),
                  sector.index2 = test[i,'id2'],
                  point2 = c(0,1),
                  directional = 1,
                  arr.type = arr.type,
                  # arr.width = 0,
                  arr.length = arr.length,
                  col=color.function[[test[i,'cell_from']]](test$weighted_n[i])
      )
      
    }
    
    upViewport()
    list.obj <- list()
    # discrete
    lgd_points = Legend(at = c('ligand', "receptor"), type = "points", gap = unit(2, "mm"),
                        legend_gp = gpar(col = lrColor), title_position = "topleft",
                        title = "")
    # discrete
    lgd_lines = Legend(at = rownames(cell_color), type = "lines", gap = unit(2, "mm"),
                       legend_gp = gpar(col = cell_color$color, lwd = 2), title_position = "topleft",
                       title = "Cell type")
    
    list.obj <- c(list.obj, list(lgd_lines))
    list.obj <- c(list.obj, list(lgd_points))
    # # continuous
    for (f in color.function) {
      lgd_links = Legend(at = c(0, 0.5, 1), col_fun = f, gap = unit(1, "mm"),direction="horizontal",
                         title_position = "topleft", title = "")
      list.obj <- c(list.obj, lgd_links)
    }
    
    lgd_list_vertical = packLegend(list = list.obj,
                                   gap = unit(2, "mm") # ligend distance
    )
    
    draw(lgd_list_vertical, x = circle_size, just = "left")
  }else{
    
    if(is.null(linkcolor)){
      col_fun = colorRamp2(c(0, 1), linkcolor)
    }else{
      col_fun = colorRamp2(c(0, 1), c("#FFFFFF", "#f349eb"))
    }
    
    for(i in 1:nrow(test)){
      my_Line_color = as.character(cell_color[test[i,1],1])
      print(test$weighted_n[i])
      circos.link(sector.index1 = test[i,'id1'],
                  point1 = c(0,1),
                  sector.index2 = test[i,'id2'],
                  point2 = c(0,1),
                  directional = 1,
                  arr.type = arr.type,
                  # arr.width = 0,
                  arr.length = arr.length,
                  col=col_fun(test$weighted_n[i])
      )
      
    }
    upViewport()
    
    # discrete
    lgd_points = Legend(at = c('ligand', "receptor"), type = "points", gap = unit(2, "mm"),
                        legend_gp = gpar(col = lrColor), title_position = "topleft",
                        title = "")
    # discrete
    lgd_lines = Legend(at = rownames(cell_color), type = "lines", gap = unit(2, "mm"),
                       legend_gp = gpar(col = cell_color$color, lwd = 2), title_position = "topleft",
                       background="#FFFFFF", title = "Cell type")
    # # continuous
    lgd_links = Legend(at = c(0, 0.5, 1), col_fun = col_fun, gap = unit(2, "mm"), direction="vertical",
                       background="#FFFFFF", title_position = "topleft", title = "Adjusted Score")
    
    lgd_list_vertical = packLegend(lgd_lines, lgd_points, lgd_links,
                                   gap = unit(4, "mm") # ligend distance
    )
    
    draw(lgd_list_vertical, x = circle_size, just = "left")
  }
  
  
}

####END####

