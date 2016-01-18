options(echo=TRUE) 
args <- commandArgs(trailingOnly = TRUE)
print(args)
dataset_file<-args[1]
signature_file<-args[2]
class.number<-as.numeric(args[3])
if(class.number >= 3){
    res.number<-as.numeric(args[4])
    gene.res<-as.numeric(args[5])
}else if(class.number < 2){
    print('error class number is minumum 2')
}
###
source('hopfield_conv.R')
########### maximum and corrispondent value function
smode<-function(x){
    xtab<-table(x)
    modes<-xtab[max(xtab)==xtab]
    mag<-as.numeric(modes[1]) #in case mult. modes, this is safer
    themodes<-as.numeric(names(modes))
    mout<-list(themodes=themodes,modeval=mag)
    return(mout)
}
##################
################################################################################################################################
########################### MAIN
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
if(class.number != num.cls){
    print('error class number')
}
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
summary<-numeric(0)
summary_E<-numeric(0)
summary_second_class<-numeric(0)
if(class.number >=3){
    for(rand in 1:res.number){
        print(paste('signature radomization number=',rand,sep=''))
        ####### memories resampling
        memories<-mem_start
        dimension2<-dimension
        #####
        for(cls in 1:num.cls){
            if(dimension[cls,1] > gene.res) {
                gene.del<-as.matrix((sum(dimension[1:cls-1,1])+1):sum(dimension[1:cls,1]))
                sample.genes<-sample(1:dimension[cls,1],gene.res)
                memories<-mem_start[-gene.del[-sample.genes,1],]
                dimension2[cls,1]<-gene.res
            }
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
        final<-numeric(0)
        E<-numeric(0)
        second_class<-numeric(0)
        for(col in 1:dim(dataset_ord)[2]){
            print(paste('column to be classified=',col,sep=''))
            sample<-as.matrix(dataset_ord[,col])
            E<- rbind(E,(-0.5*(t(sample) %*% W %*% sample)))
            final<-rbind(final,Hop_conv(sample, dataset_ord,memories1))
            ###################################################################################################################
            ################## secondary class
            if(num.cls > 2){
                delete<-as.numeric(final[dim(final)[1],1])
                ###delete memories of the convergence class
                memories2<-memories[,-delete]
                memories2<-memories2[-(sum(dimension2[1:delete-1,1])+1):-(sum(dimension2[1:delete,1])),]
                ###############
                m3<-match(rownames(memories2),rownames(dataset))
                w3<-which(!is.na(m3))
                w3s<-which(is.na(m3))
                memories2<-as.matrix(memories2[w3,])
                #### new weight matrix
                memories2_t<-t(memories2)
                W2<-matrix(data = 0, nrow = dim(memories2)[1], ncol = dim(memories2)[1], byrow = FALSE,dimnames = NULL)
                sum<-0
                for(i in 1:dim(memories2)[1]){
                    for(j in 1:dim(memories2)[1]){
                        for(l in 1:dim(memories2)[2]){
                            if(memories2[i,l]==memories2_t[l,j]){
                                sum<-sum + memories2[i,l]*memories2_t[l,j]
                            }else{
                                sum<-sum + memories2[i,l]*memories2_t[l,j]*(p-1)
                            }
                        }
                        W2[i,j]<-(1/dim(memories2)[1])*sum
                        sum<-0
                    }
                }
                for(i in 1:dim(W2)[2]){
                    W2[i,i]<-0
                }
                ##restriction dataset to only gene signatures(as memories2 are) and order dataset genes as memories2 nodes
                dataset_ord2<-numeric(0)
                for(node in 1:dim(memories2)[1]){
                    m<-match(rownames(dataset),rownames(memories2)[node])
                    w<-which(!is.na(m))
                    dataset_ord2<-rbind(dataset_ord2,dataset[w,])
                }
                l<-Hop_conv(sample, dataset_ord2,memories2)
                if(l >= delete){
                    l<-as.numeric(l)+1
                    print(l)
                }
                second_class<-rbind(second_class,l)
                #########################################################################################################################
            }
            final<-as.matrix(final)
            E<-as.matrix(E)
            second_class<-as.matrix(second_class)
            ####
            
        }
        summary<-cbind(summary,final)
        summary_E<-cbind(summary_E,E)
        summary_second_class<-cbind(summary_second_class,second_class)
    }  
    ######################### compute class as maximum on 100 resampling and energy as median on 100 resampling
    class<-vector(length=dim(summary)[1],mode="numeric")
    if(num.cls > 2){
        second_class2<-vector(length=dim(summary)[1],mode="numeric")
    }  
    energy<-vector(length=dim(summary)[1],mode="numeric")
    for (sample in 1:dim(summary)[1]){
        class_prov<-smode(as.matrix(summary[sample,]))$themodes
        if(length(class_prov)>=2){
            class[sample]<-smode(as.matrix(summary[sample,]))$themodes[1]
            second_class2[sample]<-smode(as.matrix(summary[sample,]))$themodes[1]
        }else{
            class[sample]<-smode(as.matrix(summary[sample,]))$themodes
            second_class2[sample]<-smode(as.matrix(summary_second_class[sample,]))$themodes
        }
        m<-match(as.matrix(summary[sample,]),class[sample])
        w<-which(!is.na(m))
        energy[sample]<-median(as.matrix(summary_E[sample,w]))
    }
    energy<-as.matrix(energy)
    class<-as.matrix(class)
    if(num.cls > 2){
        second_class2<-as.matrix(second_class2)
        report<-cbind(class,second_class2,energy)
    }
}else{
    memories<-mem_start
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
    final<-numeric(0)
    E<-numeric(0)
    second_class<-numeric(0)
    for(col in 1:dim(dataset_ord)[2]){
        print(paste('classified column=',col,sep=''))
        sample<-as.matrix(dataset_ord[,col])
        E<- rbind(E,(-0.5*(t(sample) %*% W %*% sample)))
        final<-rbind(final,Hop_conv(sample, dataset_ord,memories1))
    }
        final<-as.matrix(final)
        E<-as.matrix(E)
        report<-cbind(final,energy)
}    
###################################### 
write.table(report,'final_report_withoutFDR.txt',sep='\t',col.names=F,row.names=F)