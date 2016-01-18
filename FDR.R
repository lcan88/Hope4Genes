# options(echo=TRUE) 
# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# n.class<-as.numeric(args[1])
n.class<-4
####permutations results
permutation.report<-read.table('permutation_report.txt',sep='\t',header=F)
###classification results
report<-read.table('final_report_withoutFDR.txt',sep='\t',header=F)
####
report<-as.matrix(report)
permutation.report<-as.matrix(permutation.report)
### CALCOLO PVALUE
nominal.p_iteraz<-vector(length=dim(report)[1],mode="numeric")
for (i in 1:n.class){
    m1<-match(permutation.report[,1],i)
    w1<-which(!is.na(m1))
    energy.summary.perm<-as.matrix(permutation.report[w1,2]) #energia relativa a classe i
    for (sample in 1:dim(report)[1]){
        if(report[sample,1]== i) { 
            if(n.class >2){
                ranking_E<-rank(c(report[sample,3],energy.summary.perm),ties.method= "first")
            }  else{
                ranking_E<-rank(c(report[sample,2],energy.summary.perm),ties.method= "first")
            }
            nominal.p_iteraz[sample]<-ranking_E[1]/(length(ranking_E)*n.class)
        }
    }
}
#####clacolo FDR
BH.FDR_iteraz<-nominal.p_iteraz*dim(report)[1]/rank(nominal.p_iteraz,ties.method= "first")
BH.FDR_iteraz<-as.matrix(BH.FDR_iteraz)
FDR.report<-cbind(report,BH.FDR_iteraz)
write.table(FDR.report,'final_report_withFDR.txt',sep='\t',col.names=F,row.names=F)