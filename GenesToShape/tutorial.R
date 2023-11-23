library(biomaRt)
library(enrichR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("biomaRt", "enrichR"))


tf_reg_genes<-function(enrichment_result, pvalue_threshold_enrich, tf_row_id){
  padj_sig<-enrichment_result[enrichment_result$TRANSFAC_and_JASPAR_PWMs.Adjusted.P.value < pvalue_threshold_enrich,]
  hgnc_symbols <- unlist(strsplit(enrichment_result[as.numeric(tf_row_id),9],";"))
  
  hgnc_symbols

}

# takes a list of genes, mart, genes data, chromosome size file and filename as input and 
# writes a file with filename_filtered.bed containing the coordinates of gene TSS
extract_unique_tss<-function(list_genes, mart, hg38_genes, filename, chrom_sizes_file){
  hgnc_symbols<-list_genes
  
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

generate_shape <- function(list_files, x){
  shape_res<-c()
  data<-read.table(list_files[x],header=F)
  r_seqid<-as.data.frame(data[seq(1,dim(data)[1],2),])
  #remove duplicate sequences if any
  uniq_seqid<-unique(r_seqid[,1]) #unique genomic indices
  r_seq<-lapply(1:length(uniq_seqid), function(id){
    genomic_seq<-data[(grep(uniq_seqid[id],data[,1]))+1,]
    if(length(genomic_seq) == 1){
      genomic_seq<-genomic_seq
    } else {
      genomic_seq<-genomic_seq[1]
    }
    genomic_seq
  })
  r_seq<-as.data.frame(unlist(r_seq))
  #remove sequences with N 
  if(length(grep("N", r_seq[,1]))>=1){
    n_seq<-grep("N", r_seq[,1]) # sequences with N
    r_seq<-as.data.frame(r_seq[-(n_seq),1]) 
    uniq_seqid<-uniq_seqid[-n_seq]
  }
  #generating pentamers of all sequences
  all_seq_sub_shape<-lapply(1:dim(r_seq)[1],function(y){
    aseq<-as.character(r_seq[y,1])
    aseq_sub<-lapply(1:(nchar(aseq)-(kmer_size-1)),function(l){
      one_sub<-toupper(substr(aseq,l,as.numeric(l+(kmer_size-1))))
      one_sub
    })
    aseq_sub<-unlist(aseq_sub)
    #generating DNA shape profiles for all sub sequences
    one_seq_shape<-lapply(1:length(aseq_sub),function(seq_id){
      ens_data<-norm_bin[grep(toupper(aseq_sub[seq_id]),rownames(norm_bin)),]
      avg_data<-avg_shape[grep(toupper(aseq_sub[seq_id]),rownames(avg_shape)),]
      seq_shape<-c(avg_data,ens_data)
      seq_shape
    })
    one_seq_shape<-do.call(rbind,one_seq_shape)
    #naming of rows
    shape_pos<-c(3:(length_seq-2))
    all_rows<-lapply(1:length(shape_pos), function(pos){
      pos_index<-paste(as.character(r_seqid[y,1]),"_pos",shape_pos[pos],"_",aseq_sub[pos],sep="")
      pos_index
    })
    rownames(one_seq_shape)<-all_rows
    one_seq_shape
  })
  shape_res<-do.call(rbind,all_seq_sub_shape)
  
  shape_res  
}

# generate average profiles of DNA shape generated 
mean_shapes <- function(file_name){
  shape_data<-as.matrix(read.table(file_name,header=T,stringsAsFactors = F))
  table<-c()
  num_pos<-(seq_length-4)
  #calculate mean per position
  col_mean_vals<-lapply(1:num_pos, function(i){
    posi<- seq(i, dim(shape_data)[1], num_pos)
    posi_data <- shape_data[posi,]
    mean_vals<-colMeans(posi_data)
    mean_vals
  })
  shape_means<-do.call(rbind,col_mean_vals)
  rownames(shape_means)<-paste("pos",3:(seq_length-2),sep="_")
  write.table(shape_means, file=paste("pos_mean_",dna_shape_files[x],sep=""),sep="\t")
  shape_means
}

# read the file containing list of genes
list_genes <- as.data.frame(read.table("list_genes", header = F))

# predict potential regulatory TFs for the genes
enrichment_result <- as.data.frame(enrichr(list_genes[,1], "TRANSFAC_and_JASPAR_PWMs"))
write.table(enrichment_result,"enrichment_result",row.names = F,sep="\t",quote=F)

# define which TF you wish to extract gene list for
tf_regulated_genes <- as.data.frame(tf_reg_genes(enrichment_result, 0.05, 1))
write.table(tf_regulated_genes,"tf_regulated_genes",col.names = F,row.names = F,sep="\t")

# read file containing gene data for reference genome 
ucsc_genes<-read.delim("genes.bed",header=T,sep="\t")

# load ensemble gene associated mart for reference genome
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# generate files for TSS coordinates of genes
extract_unique_tss(tf_regulated_genes, ensembl, ucsc_genes, 
                   filename="tf_reg_degs", 
                   chrom_sizes_file = "hg38.chrom.sizes")


# load DynaSeq libraries
norm_bin<-readRDS("DynaSeq65.rds")
avg_shape<-readRDS("DynaSeq13.rds")
kmer_size<-5
length_seq<-20

# list the fasta files
list_files<-system("ls *.fasta",intern=T)
print(list_files)

# generate the shape profiles
lapply(1:length(list_files),function(x){
  shape_profiles <- generate_shape(list_files, x)
  write.table(shape_profiles,
              paste("dna_shape_",strsplit(list_files[x],".fasta")[[1]][1],sep=""),quote=F)
  
})

#calculate averages per position
#for static shape data
params<-c("Shear", "Stretch", "Stagger", "Buckle", "PropTw", "Opening", 
          "Shift", "Slide", "Rise", "Tilt", "Roll", "Twist", "MGW")

seq_length<-20

dna_shape_files<-system("ls dna_shape_*",intern=T)

# generate mean profiles
mmat<-lapply(1:length(dna_shape_files),function(x){
  data<-mean_shapes(dna_shape_files[x])
})

list_pos_mean<-system("ls pos_mean*",intern=T)

data<-read.table(list_pos_mean[1], header=T)

pdf("plots_promoters_shape.pdf")
par(mfrow=c(2,2))
x_ticks<-seq(0,seq_length-4,2)[-1]
x_labels<-seq(-(seq_length-4),1,2)[-1]
width<-3
for(i in 1:13){
     plot(data[,i],xlab="position",main="DNA shape profiles", las=1,
          cex.main=1,type="l",ylab=params[i],col="red",lwd=3,cex.sub=0.8,
          sub="Promoter regions of DEGs regulated by TF", yaxt="n", xaxt="n")
     axis(1,at=x_ticks,labels=x_labels,cex.axis=0.8) #x axis labels
     axis(2,cex.axis=0.7,cex.lab=0.8,las=1) #y axis labels
}
dev.off()


#boxplots for ensemble data
dna_shape_files<-system("ls 5bin_dna_shape_*",intern=T)
num_bins<-as.numeric(5)
num_datasets<-as.numeric(1)
data_variable_names<-c(paste("data",1:num_datasets,sep="_"))

# # assign the shape profiles to variable names
# for (i in 1:length(data_variable_names)) {
#      assign(data_variable_names[i], read.table(dna_shape_files[i],header=T))
# }

pdf("ensemble_boxplots_5bin.pdf")
par(mfrow=c(2,2))
par_col<-seq(14,13+(13*num_bins),5)
for(i in 1:13){
     boxplot(list(data[,par_col[i]],data[,par_col[i]+1],
                  data[,par_col[i]+2],data[,par_col[i]+3],
                  data[,par_col[i]+4]),ylab=params[i], 
             xlab="ensemble bins", las=1)
     axis(1,at=seq(1,num_datasets*5,num_datasets),labels=c("bin 1", "bin 2", "bin 3", "bin 4", "bin 5"))
}
dev.off()


