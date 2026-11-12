library(CellChatESS)
f<-getFromNamespace('cellchat_prob_from_avg_cpp','CellChatESS')
ess<-getFromNamespace('cellchat_prob_ess_cpp','CellChatESS')
avg<-matrix(c(0.013,0.123,0.07,0,0.37,0.29,0.19,0),nrow=4)
li<-matrix(c(1L,0L,1L,3L,4L,0L),nrow=3,byrow=TRUE)
ri<-matrix(rep(c(2L,0L),3),nrow=3,byrow=TRUE)
empty<-matrix(0L,nrow=3,ncol=1)
flags<-rep(FALSE,3)
boot<-array(0,c(4,2,3));boot[,,1]<-avg;boot[,,2]<-avg*1.5;boot[,,3]<-avg*0.5
answer<-f(avg,boot,li,ri,empty,empty,empty,empty,flags,flags,0.5,1)
expected<-array(0,c(2,2,3))
for(lr in 1:2) for(s in 1:2) for(r in 1:2) {
 ligand<-if(lr==1) avg[1,s] else exp((log(avg[1,s])+log(avg[3,s]))/2)
 x<-ligand*avg[2,r];expected[s,r,lr]<-x/(0.5+x)
}
stopifnot(identical(as.numeric(answer$prob),as.numeric(expected)))
stopifnot(all(answer$pval[,,1:2]==1/3),all(answer$pval[,,3]==1))
group<-1:2;permutations<-matrix(rep(group,3),nrow=2)
call_ess<-function(threads) ess(avg,group,permutations,avg,li,ri,empty,empty,empty,empty,flags,flags,0.5,1,threads)
a<-call_ess(1L);b<-call_ess(4L)
stopifnot(identical(as.numeric(a$prob),as.numeric(expected)),all(a$pval[,,1:2]==0),
 all(a$pval[,,3]==1),identical(a$prob,b$prob),identical(a$pval,b$pval))
cat('Singleton direct value, multi-component geometric mean, strict ties, zero observed probability, dense and ESS kernels passed.\n')
cat('Synthetic requested=4 actual threads_used=',b$threads_used,'; representative server OpenMP validation remains required.\n',sep='')
