# this script was written as a part of springer bookchapter titled 
# "Incorporating sequence-dependent DNA shape and dynamics into 
# transcriptome data analysis"
# Authors: Manisha Kalsan, Almas Jabeen, Shandar Ahmad
# Affiliation: School of Computational and Integrative Sciences, 
# Jawaharlal Nehru University, Delhi, India
# Email: shandar@jnu.ac.in
# 
# 
# this script uses the fasta files and generates their DNA shape profiles
# using the dictionaries from Dynaseq. It is a tool which offers dictionaries
# of 13 DNA confromational parameters for an exhaustive set of pentamers in 
# form of static shape which represents an average shape of the pentamer and
# shape ensemble which incorporates the dynamic aspect of structure to generate
# the DNA shape profiles. 
# 
# This tool is developed by SciWhyLab (http://www.sciwhylab.org/dynaseq/) 
# and the is available at the github repository DynaSeq at 
# https://github.com/manishakalsan/DynaSeq. To read more about it, refer to 
# 
# Andrabi M, Hutchins AP, Miranda-Saavedra D, Kono H, Nussinov R, 
# Mizuguchi K, Ahmad S. Predicting conformational ensembles and genome-wide 
# transcription factor binding sites from DNA sequences. Sci Rep. 2017 Jun 22;
# 7(1):4071. doi: 10.1038/s41598-017-03199-6. PMID: 28642456; PMCID: PMC5481346.
# 
# 
# 


list_files<-system("ls *.fasta",intern=T)
print(list_files)
norm_bin<-readRDS("dictionary_ensemble_5bin_5mer.rds")
avg_shape<-readRDS("dictionary_static_5mer.rds")
kmer_size<-5
length_seq<-20

all_shape<-lapply(1:length(list_files),function(x){
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
   write.table(shape_res,paste("5bin_dna_shape_",strsplit(list_files[x],".fasta")[[1]][1],sep=""),quote=F)
})

#calculate averages per position
#for static shape data
params<-c("Shear", "Stretch", "Stagger", "Buckle", "Prop-Tw", "Opening", "Shift", "Slide", "Rise", "Tilt", "Roll", "Twist","MGW")
seq_length<-20
dna_shape_files<-system("ls 5bin_dna_shape_*",intern=T)
mmat<-lapply(1:length(dna_shape_files),function(x){
   print(x)
   shape_data<-as.matrix(read.table(dna_shape_files[x],header=T,stringsAsFactors = F))
   table<-c()
   print(head(shape_data))
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
})

