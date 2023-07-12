# this script was written as a part of springer bookchapter titled 
# "Incorporating sequence-dependent DNA shape and dynamics into 
# transcriptome data analysis"
# Authors: Manisha Kalsan1, Almas Jabeen1, Shandar Ahmad1
# Affiliation: School of Computational and Integrative Sciences, 
# Jawaharlal Nehru University, Delhi, India
# Email: shandar@jnu.ac.in
# 
# This script takes input files available under input in this repository 
# to generate TSS bed files for a list of genes regulated by a TF enriched 
# in their promoter sites and the genes which are not regulated by it
# 
# Input files required are DESeq2 result file, UCSC meta data file, and 
# chromosome size files for hg38
# The input file used here as a result from DESeq2 have been downloaded 
# from GEO under the accession GSE138496 from a study by Weis-Banke et al 2020.
# 
# 
# 

library(enrichR)
list_gse<-system("ls GSE* ",intern=T)
de_data<-read.table(list_gse[1],header=F)

#remove rows with NA
data1<-de_data[complete.cases(de_data), ]

#take only results with p-adjusted less than 0.05
data_significant<-data1[as.numeric(data1$V7) < 0.05,]

#generate list of TFs regulating the genes
enrichment_result <- as.data.frame(enrichr(data_significant[,1], "TRANSFAC_and_JASPAR_PWMs"))
print(head(enrichment_result))

write.table(enrichment_result,"enrichment_result",row.names = F,sep="\t",quote=F)
print("ENRICHMENT RESULTS written to enrichment_result")

padj_sig<-enrichment_result[enrichment_result$TRANSFAC_and_JASPAR_PWMs.Adjusted.P.value <= 0.05,]

hgnc_symbols <- unlist(strsplit(enrichment_result[1,9],";"))
tf_regulated_genes<-data_significant[grep(paste(hgnc_symbols,collapse = "|"),data_significant$V1),1]
remaining_genes<-data_significant[-grep(paste(hgnc_symbols,collapse = "|"),data_significant$V1),1]

library(biomaRt)

# Connect to the appropriate database using the 'biomaRt' package
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

#load the data downloaded from UCSC table browser
hg38_genes<-read.delim("genes.bed",header=T,sep="\t")

extract_unique_tss<-function(list_genes, hg38_genes, filename, chrom_sizes_file){
   hgnc_symbols<-list_genes
   print(hgnc_symbols)
   #generating meta data for ID mapping from biomart
   id_info<-getBM(attributes = c("hgnc_symbol", "ensembl_gene_id_version",
                                 "ensembl_transcript_id_version",
                                 "entrezgene_id", "strand","transcription_start_site"),
                  filters = "hgnc_symbol",
                  values = hgnc_symbols,
                  mart = ensembl)
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

extract_unique_tss(remaining_genes, hg38_genes, filename="remaining_degs", chrom_sizes_file = "hg38.chrom.sizes")
extract_unique_tss(tf_regulated_genes, hg38_genes, filename="tf_reg_degs", chrom_sizes_file = "hg38.chrom.sizes")




