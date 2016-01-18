#################### modified Hopefiled model implementation 
######### function for sample convergence
Hop_conv <- function(sample, TCGA_ord,memories){
    new_sample<-sample
    count<-TRUE
    indice<-0
    iteraz<-1
    contSI<-0
    contNO<-0
    count_countNO<-0
    while(count){ 
        node<-sample(1:dim(TCGA_ord)[1],1)
        z_1<-W[node,] %*% sample
        if(z_1>0){
            new_sample_node<-(1)
        }else if(z_1<0){
            new_sample_node<-0
        }else if(z_1==0){
            new_sample_node<-sample[node]
            #print('sn qui')
        } 
        
        if(new_sample[node]==new_sample_node){
            contSI<-contSI+1
        }else{
            contNO<-contNO+1
        }
        if(contNO==0){
            count_countNO<-count_countNO+1
        }else if(contNO!=0){
            count_countNO<-0
        }
        new_sample[node]<-new_sample_node
        riass_perc<-numeric(0)
        ###verifica classificazione
        stampa<-numeric(0)
        indice<-indice+1
        riga<-numeric(0)
        for(class in 1:dim(memories)[2]){
            hamm_dist<-0
            count_uni<-0
            for(i in 1:dim(memories)[1]){
                if(memories[i,class]==1){
                    count_uni<-count_uni+1
                    if(memories[i,class]==new_sample[i]){
                        hamm_dist<-hamm_dist+1
                    }  
                } 
            }
            if(hamm_dist>=(90*count_uni)/100){
                riga<-rbind(riga,cbind(class,hamm_dist))
                count<-FALSE
            }  
        }
        if(count==FALSE && dim(riga)[1]>1){
            c<-which(riga[,2] == max(riga[,2]))
            if(length(c)>1){
                finale<-riga[c[1],1]
            }else{
                finale<-riga[c,1]
            }
        }else if(count==FALSE && dim(riga)[1]==1) {finale<-riga[1,1]
                                                   
        }
        if(count_countNO >= dim(memories)[1]+(dim(memories)[1])/3 && count != FALSE){
            dist<-numeric(0)
            for(class in 1:dim(memories)[2]){
                hamm_dist<-0
                for(i in 1:dim(memories)[1]){
                    if(memories[i,class]==new_sample[i]){
                        hamm_dist<-hamm_dist+1
                    } 
                }
                dist<-rbind(dist,hamm_dist)  
            }  
            c<-which(dist == max(dist))
            if(length(c)>1){
                finale<-c[1]
            }else{
                finale<-c
            }
            count<-FALSE
        }  
        #######
        if(indice %% 100==0){
            contSI<-0
            contNO<-0
        }
        #####
        sample<-new_sample
        iteraz<-iteraz+1
    }
    return(finale)
}