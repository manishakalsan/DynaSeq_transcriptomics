# this script will generate TSS regions for a list of DEGs (input as DESeq2 result file) 
# which are regulated by a TF enriched in their promoter regions and the DEGs which are 
# not regulated by the same TF
# This will also write whole enrichment results for TFs in the file enrichment_result
# .
# .
# the whole script is divided into three functions: enrichment_result, 
# tf_reg_unreg_genes, and extract_unique_tss
# .
# .
# USAGE: enrichment_result(data, pvalue_threshold_dge)
# parameter data the filename containing the DESeq2 result
# parameter pvalue_threshold_dge threshold value of pvalue to consider 
# the DGE result significant; is applied to adjusted p-value
# output writes the filtered DESeq2 result to file data_significant
# output writes file with enrichment result named enrichment_result
# return enrichment results
# .
# .
# USAGE: tf_reg_unreg_genes(enrichment_result,tf_row_id)
# parameter enrichment_result the df containing result from enrichr generated 
# with the function enrichment_result
# parameter tf_row_id numeric value of the row number corresponding to which 
# the list of genes regulated and remaining genes are to be calculated
# output writes two files named tf_regulated_genes and remaining_genes
# .
# .
# USAGE: extract_unique_tss(list_genes, mart, ucsc_genes, filename, chrom_sizes_file)
# parameter list_genes list of genes ()
# parameter mart the biomart mart to be used according to the organism
# parmeter ucsc_genes df corresponding to file containing all the metadata 
# of all genes for the organism derived from UCSC table browser
# parameter filename the name  to be used to write the file
# parameter chrom_sizes_file filename corresponding to text file containing list of 
# chromosomes and their sizes for the corresponding organism
# output writes a file with with value specified to filename containing unique TSS 
# for each gene in the list of genes
# .
# .
# Example
# library(biomaRt)
# library(enrichR)
# list_gse<-system("ls GSE* ",intern=T)
# enriched<-enrichment_result(list_gse[1],0.05)
# tf_reg_unreg_genes(enriched,0.05,1)
# ucsc_genes<-read.delim("genes.bed",header=T,sep="\t")
# ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# remaining_genes<-read.table("remaining_genes",header=F)
# tf_regulated_genes<-read.table("tf_regulated_genes",header=F)
# extract_unique_tss(remaining_genes, ensembl, ucsc_genes, filename="remaining_degs", chrom_sizes_file = "hg38.chrom.sizes")
# extract_unique_tss(tf_regulated_genes, ensembl, ucsc_genes, filename="tf_reg_degs", chrom_sizes_file = "hg38.chrom.sizes")


enrichment_result<-function(data, pvalue_threshold_dge){
   de_data<-read.table(data,header=F)
   pvalue_threshold_dge<-pvalue_threshold_dge
   data1<-de_data[complete.cases(de_data), ]
   data_significant<-data1[as.numeric(data1$V7) < pvalue_threshold_dge,]
   write.table(data_significant,"data_significant",sep="\t",row.names = F)
   enrichment_result <- as.data.frame(enrichr(data_significant[,1], "TRANSFAC_and_JASPAR_PWMs"))
   print("ENRICHMENT RESULTS written to enrichment_result")
   write.table(enrichment_result,"enrichment_result",row.names = F,sep="\t",quote=F)
   enrichment_result
}

tf_reg_unreg_genes<-function(enrichment_result,pvalue_threshold_enrich,tf_row_id){
   enrichment_result<-enrichment_result
   pvalue_threshold_enrich<-pvalue_threshold_enrich
   data_significant<-read.table("data_significant",header=T)
   tf_row_id<-tf_row_id
   padj_sig<-enrichment_result[enrichment_result$TRANSFAC_and_JASPAR_PWMs.Adjusted.P.value < pvalue_threshold_enrich,]
   hgnc_symbols <- unlist(strsplit(enrichment_result[as.numeric(tf_row_id),9],";"))
   tf_regulated_genes<-data_significant[grep(paste(hgnc_symbols,collapse = "|"),data_significant$V1),1]
   remaining_genes<-data_significant[-grep(paste(hgnc_symbols,collapse = "|"),data_significant$V1),1]
   write.table(tf_regulated_genes,"tf_regulated_genes",col.names = F,row.names = F,sep="\t")
   write.table(remaining_genes,"remaining_genes",col.names = F,row.names = F,sep="\t")
}

extract_unique_tss<-function(list_genes, mart, hg38_genes, filename, chrom_sizes_file){
   hgnc_symbols<-list_genes
   mart<-mart
   print(hgnc_symbols)
   #generating meta data for ID mapping from biomart
   id_info<-getBM(attributes = c("hgnc_symbol", "ensembl_gene_id_version",
                                 "ensembl_transcript_id_version",
                                 "entrezgene_id", "strand","transcription_start_site"),
                  filters = "hgnc_symbol",
                  values = hgnc_symbols,
                  mart = mart)
   print("metadata extracted ")
   print( head(id_info))
   #unique ensemble transcript IDs
   unique_entrez<-unique(id_info[,4])
   unique_entrez<-na.omit(unique_entrez)
   
   #generate one TSS per gene
   unique_gene_TSS<-lapply(1:length(unique_entrez), function(x){
      gene_info_mart<-id_info[grep(unique_entrez[x], id_info[,4]),]
      TSS<-strand<-c()
      if(gene_info_mart[1,5] == "-1"){
         TSS<-max(gene_info_mart$transcription_start_site)
      } else {
         TSS<-min(gene_info_mart$transcription_start_site)
      }
      c(gene_info_mart[1,c(1,3)], TSS)
   })
   unique_gene_TSS<-do.call(rbind,unique_gene_TSS)
   colnames(unique_gene_TSS)<-c("hgnc_symbol","ensembl_transcript_id_version","TSS")
   print(c("unique_gene_TSS", head(unique_gene_TSS)))
   unique_transcript_ids<-unlist(unique(unique_gene_TSS[,2]))
   
   #extract data for the selected genes
   unique_transcript_loc<-lapply(1:length(unique_transcript_ids), function(x){
      gene_info_mart<-id_info[id_info[,3] == unique_transcript_ids[x],]
      gene_info_ucsc<-hg38_genes[hg38_genes[,1]==unique_transcript_ids[x],]
      tss<-as.data.frame(unique_gene_TSS[unique_gene_TSS[,2] == unique_transcript_ids[x], ])
      one_gene<-c(gene_info_ucsc[1,c(2,8,3)],tss[1,3],gene_info_mart[1,1])
      names(one_gene)<-c("chr","exonCount","strand","tss","gene_name")
      one_gene
   })
   unique_transcript_loc<-do.call(rbind,unique_transcript_loc)
   
   #prepare bed file from info generated above
   TSS_coord_plus1<-as.numeric(unique_transcript_loc[,4]) +1
   genes_bed<-cbind(unique_transcript_loc[,c(1,4)],TSS_coord_plus1,unique_transcript_loc[,c(5,2,3)])

   #filter features by valid chromosome names
   chrom_sizes<-read.table(paste(chrom_sizes_file),header=F)
   selected_rows <- genes_bed[genes_bed[,1] %in% chrom_sizes[,1], ]
   write.table(selected_rows,file=paste(filename,"_filtered.bed",sep=""),col.names = F,row.names = F,sep="\t",quote=F)
   
   print(paste("Chromosome filtered coordinates written to ",paste(filename,"_filtered.bed",sep=""),sep=" "))
}



