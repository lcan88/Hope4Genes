options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
dataset_file<-args[1]
signature_file<-args[2]
class.number<-as.numeric(args[3])
if(class.number >= 3){
    gene.res<-as.numeric(as.numeric(args[4]))
    permutations.number<-as.numeric(args[5])
}else if(class.number < 2){
    print('error class number is minumum 2')
}else if(class.number == 2){
    permutations.number<-as.numeric(args[4])
}
source('hopfield_conv.R')
dataset<-as.matrix(read.delim(dataset_file,sep='\t',header=TRUE,row.names=1))
signature<-as.matrix(read.delim(signature_file,sep='\t',header=TRUE,row.names=1))
#dataset transformed in 1 and 0
for (i in 1:dim(dataset)[1]){
    for(j in 1:dim(dataset)[2]){
        if(dataset[i,j]>0){
            dataset[i,j]<- 1
        }else if(dataset[i,j] <= 0){
            dataset[i,j]<-(-1)
        }
    }
}
############### define memories
num.features<-length(signature[,1])
num.cls<-length(table(signature[,1]))
dimension<-as.matrix(table(signature))
mem_start<-matrix(nrow =num.features , ncol = num.cls)
row.names(mem_start)<-row.names(signature)
for(i in 1:num.features){
    for (j in 1:num.cls){
        if(signature[i,1] != j){
            mem_start[i,j] <- (-1)  
        }else if(signature[i,1] == j){
            mem_start[i,j] <- 1 
        }
    }
}  
final<-numeric(0)
E<-numeric(0)
for(permutation in 1:permutations.number){
    print(paste('permutation number=',permutation,sep=''))
   if (class.number >= 3){
       ####### memories resampling
       memories<-mem_start
       dimension2<-dimension
       for(cls in 1:num.cls){
           if(dimension[cls,1] > gene.res) {
               gene.del<-as.matrix((sum(dimension[1:cls-1,1])+1):sum(dimension[1:cls,1]))
               sample.genes<-sample(1:dimension[cls,1],gene.res)
               memories<-mem_start[-gene.del[-sample.genes,1],]
               dimension2[cls,1]<-gene.res
           }
       }
   }else{
       memories<-mem_start
   }
    #### match memories1 on the dataset
    p<-dim(memories)[2]
    m1<-match(rownames(memories),rownames(dataset))
    w1<-which(!is.na(m1))
    w1s<-which(is.na(m1))
    memories1<-as.matrix(memories[w1,])
    memories1_t<-t(memories1)
    #### construct the weight matrix
    W<-matrix(data = 0, nrow = dim(memories1)[1], ncol = dim(memories1)[1], byrow = FALSE,dimnames = NULL)
    sum<-0
    for(i in 1:dim(memories1)[1]){
        for(j in 1:dim(memories1)[1]){
            for(l in 1:dim(memories1)[2]){
                if(memories1[i,l]==memories1_t[l,j]){
                    sum<-sum + memories1[i,l]*memories1_t[l,j]
                }else{
                    sum<-sum + memories1[i,l]*memories1_t[l,j]*(p-1)
                }
            }
            W[i,j]<-(1/dim(memories1)[1])*sum
            sum<-0
        }
    }
    for(i in 1:dim(W)[2]){
        W[i,i]<-0
    }
    ##restriction dataset to only gene signatures(as memories1 are) and order dataset genes as memories1 nodes
    dataset_ord<-numeric(0)
    for(node in 1:dim(memories1)[1]){
        m<-match(rownames(dataset),rownames(memories1)[node])
        w<-which(!is.na(m))
        dataset_ord<-rbind(dataset_ord,dataset[w,])
    }
    ####classification
    sample<-sample(c(-1,1),length(as.matrix(dataset_ord[,1])),replace=T)
    E<- rbind(E,(-0.5*(t(sample) %*% W %*% sample)))
    final<-rbind(final,Hop_conv(sample, dataset_ord,memories1))
    print(final)
}
final<-as.matrix(final)
E<-as.matrix(E)
permutation.report<-cbind(final,E)
write.table(permutation.report,'permutation_report.txt',sep='\t',col.names=FALSE,row.names=FALSE)
