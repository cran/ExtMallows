
#HEMM model: hierarchical extened Mallows model
#the data has a subgroup of rankings that are over-correlated
#version 1.0
#last updated in 2018.6.28


#input arguments:
#1. rankings: a n*m data frame, with each column representing a ranking list, which ranks the items from the most preferred to the least preferred. For missing items, use 0 to denote them.
#2. num.kappa: number of over-correlated ranking groups
#3. is.kappa.ranker: a list of over-correlated ranking groups, with the k-th element denoting the column numbers of the rankings that belong to the k-th group
#4. initial.method: the method for initializing the value of pi0, with four options: mean, median, geometric and random (the mean of three randomly sampled ranking lists)


#return arguments:
#1. op.phi: optimal value of phi, the phi value in independent rankings
#2. op.phi1: optimal value of phi1, the phi value in over-correlated ranking groups
#3. op.omega: optimal value of omega
#4. op.alpha: optimal value of alpha
#5. op.pi0: optimal value of pi0, ranking the items from the most preferred to the least preferred
#6. op.kappa: optimal value of kappa, denoting the items from the most preferred to the least preferred
#7. logL.max: maximum value of log-likelihood


HEMM <- function(rankings, num.kappa, is.kappa.ranker, initial.method="mean", it.max=20){
    
    
    #######################################
    #
    #step 1. input the data
    #
    #######################################
    
    rank.list=list()
    num.ranker=ncol(rankings)
    num.top=rep(0,num.ranker)
    for(k in 1:num.ranker){
      temp=rankings[,k]
      rank.list[[k]]=temp[temp!=0]
      num.top[k]=length(rank.list[[k]])
      }
    
    
    item.name=as.character(sort(unique(unlist(rank.list))))
    num.item=length(item.name)
    item.name.ref=seq(1,num.item)
    names(item.name.ref)=item.name
    y=matrix(0,num.item,num.ranker)
    for(k in 1:num.ranker){
      y[1:num.top[k],k]=item.name.ref[as.character(rank.list[[k]])]
      }
    
    
    x=matrix(0,num.item,num.ranker)
    rownames(x)=item.name
    for(k in 1:num.ranker){
      x[y[1:num.top[k],k],k]=1:num.top[k]
      x[x[,k]==0,k]=num.top[k]+1
      }
    
    
    #######################################
    #
    #step 2. initialize the parameters
    #
    #######################################
    

    num.total.corr.ranker=0
    for(k in 1:num.kappa){num.total.corr.ranker=num.total.corr.ranker+length(is.kappa.ranker[[k]])}
    num.indep.ranker=num.ranker-num.total.corr.ranker+num.kappa

    
    op.phi=0
    op.phi.h1=rep(0, num.kappa)
    op.alpha=op.omega=rep(0,num.indep.ranker)
    op.kappa=op.kappa.rank=matrix(0, num.item, num.kappa)


    #2.1 initialize pi0
    
    score=rep(0,num.item)
    if(initial.method=="mean"){score=apply(x[, 1:12], 1, mean)}
    if(initial.method=="median"){
      for(i in 1:num.item){score[i]=quantile(x[i,],0.5)}
      }
    if(initial.method=="geometric"){
      for(i in 1:num.item){score[i]=(prod(x[i,]))^(1/num.ranker)}
      }
    if(initial.method=="random"){
      random.rankers=sample(1:num.ranker,3,replace=FALSE)
      score=apply(x[,random.rankers],1,mean)
      }
    
    temp=rank(score,ties.method="random")
    op.pi0=rep(0, num.item)
    op.pi0[temp]=1:num.item
    pi0.rank=seq(1,num.item)
    pi0.rank[op.pi0]=seq(1,num.item)
    pi0.prop=pi0.prop.rank=rep(0,num.item)
    
    
     #2.2 initialize kappa
     
    num.indep.top=c(rep(0, num.kappa), num.top[(num.total.corr.ranker+1):num.ranker])
    for(k in 1:num.kappa){
      index=as.vector(y[,is.kappa.ranker[[k]]])
      index1=sort(unique(index))
      if(index1[1]==0){index1=index1[-1]}
      num.indep.top[k]=length(index1)
      score=apply(x[index1,is.kappa.ranker[[k]]], 1, mean)
      temp=rank(score,ties.method="random")
      op.kappa[temp,k]=index1
      op.kappa.rank[index1,k]=temp
      op.kappa.rank[op.kappa.rank[,k]==0,k]=num.indep.top[k]+1
      }

    
    stage=rep(0, num.indep.ranker)
    for(k in 1:num.indep.ranker){stage[k]=min(num.item-1,num.indep.top[k])}
    
    v.indep=matrix(0,max(stage),num.indep.ranker)   #number of mismatch pairs
    tau=v.indep  #probability of latent variable=1
    
    
    #2.3 initialize phi, alpha, omega and phi.h1

    op.phi=0.5
    op.phi.h1=rep(0.3, num.kappa)
    op.alpha=rep(0.2,num.indep.ranker)
    op.omega=rep(0.8,num.indep.ranker)
    
    
    #2.4 initialize independent rankings
    
    y.indep=matrix(0,num.item,num.indep.ranker)
    x.indep=matrix(0,num.item,num.indep.ranker)
    
    
    for(k in 1:num.kappa){
      y.indep[,k]=op.kappa[,k]
      }
    
    
    y.indep[,(num.kappa+1):num.indep.ranker]=y[, (num.total.corr.ranker+1):num.ranker]
    
    
    for(k in 1:num.indep.ranker){
      x.indep[y.indep[1:num.indep.top[k],k],k]=1:num.indep.top[k]
      x.indep[x.indep[,k]==0,k]=num.indep.top[k]+1
     }
    

    #the derivative of logZ(theta), Z(theta)=1+theta+...+theta^n
    
    LogZDeriva <- function(theta,n){
      temp1=sum(theta^(1:n))+1
      temp2=sum(c(1, theta^(1:(n-1)))*(1:n))
      result=temp2/temp1
      return(result)
      }
    
    #Z''(theta)/Z(theta), Z(theta)=1+theta+...+theta^n
    
     ZDeriva2 <- function(theta,n){
       if(n>2){
        temp1=sum(theta^(1:n))+1
        temp2=(2:n)*(1:(n-1))
        temp3=sum(temp2*theta^(0:(n-2)))
        result=temp3/temp1
        }else if(n==2){
        result=2}else{result=0}

        return(result)
       }

    
    
    ############################################################
    #
    # Step 3. find the MLE of the parameters using grid method
    #
    #############################################################
    
    
    
    for(it in 1:it.max){
        
        #3.1 E-step for updating tau
        
        
        for(k in 1:num.indep.ranker){
          for(i in 1:stage[k]){
            index=y.indep[i,k]
            temp1=pi0.rank[index]
            temp2=pi0.rank[x.indep[,k]>i]
            v.indep[i,k]=sum(temp2<temp1)
            }
         }
        

        for(k in 1:num.indep.ranker){
          for(i in 1:stage[k]){
            temp1=op.phi*(1-op.alpha[k]^i)
            temp2=temp1^(1:(num.item-i))
            temp3=op.omega[k]*temp1^v.indep[i,k]/(sum(temp2)+1)
            temp4=(1-op.omega[k])/(num.item-i+1)
            tau[i,k]=temp3/(temp3+temp4)
            }
          }
        
        

        #3.2 update optimal phi, alpha, omega and pi0
        
        
        #3.2.1 update phi

        for(sub.it in 1:20){
          op.phi.init=op.phi
          temp1.sum=temp2.sum=0
          for(k in 1:num.indep.ranker){
            for(i in 1:stage[k]){
              temp1=1-op.alpha[k]^i
              temp2=LogZDeriva(op.phi*temp1,num.item-i)
              temp3=v.indep[i,k]/op.phi-temp2*temp1
              temp1.sum=temp1.sum+tau[i,k]*temp3
              
              temp4=-v.indep[i,k]/op.phi^2-temp1^2*(ZDeriva2(op.phi*temp1,num.item-i)-temp2^2)
              temp2.sum=temp2.sum+tau[i,k]*temp4
              }
           }
          
          delta.phi=-temp1.sum/temp2.sum    #use Newton-Rapshon method
          if(abs(delta.phi)<0.001){break}
          op.phi=op.phi+delta.phi
          if(op.phi<0){op.phi=abs(op.phi)}
          if(op.phi>1){
            if(temp1.sum>0){
                op.phi=1
              }else{op.phi=op.phi.init}
            }
         }
       
       
       
        # 3.2.2 update alpha
        
        for(k in 1:num.indep.ranker){

          for(subit in 1:20){
            op.alpha.init=op.alpha[k]
            temp1.sum=temp2.sum=0
            for(i in 1:stage[k]){
              temp1=1-op.alpha[k]^i
              temp2=LogZDeriva(op.phi*temp1,num.item-i)
              temp3=-(v.indep[i,k]/temp1-temp2*op.phi)*i*op.alpha[k]^(i-1)
              temp1.sum=temp1.sum+tau[i,k]*temp3
              
              if(i<2){
                  temp4=1/(1-op.alpha[k])^2
               }else{
                  temp4=(i-1)*op.alpha[k]^(i-2)/temp1+i*(op.alpha[k]^(i-1)/temp1)^2
               }
              temp5=ZDeriva2(op.phi*temp1,num.item-i)-temp2^2
              temp6=-i*v.indep[i,k]*temp4-(i*op.alpha[k]^(i-1))^2*temp5
              temp2.sum=temp2.sum+tau[i,k]*temp6
             }
            
            
            delta.alpha=-temp1.sum/temp2.sum    #use Newton-Rapshon method
            if(abs(delta.alpha)<0.001){break}
            op.alpha[k]=op.alpha[k]+delta.alpha
            if(op.alpha[k]<0){op.alpha[k]=abs(op.alpha[k])}
            if(op.alpha[k]>0.99){
              if(temp1.sum>0){
                  op.alpha[k]=0.99
              }else{op.alpha[k]=op.alpha.init}
              }
            
           }
         }
        
        
        
        # 3.2.3 update omega
        
        
        for(k in 1:num.indep.ranker){
          op.omega[k]=sum(tau[1:stage[k],k])/stage[k]
          }
        

        #3.2.4 iteratively update op.pi0 by swapping pairs
        
        pi0.curr=op.pi0
        for(test.i in 2:num.item){
            
          pi0.prop=op.pi0
          pi0.prop[c(test.i-1,test.i)]=op.pi0[c(test.i,test.i-1)]
          pi0.prop.rank[pi0.prop]=seq(1,num.item)
          index=op.pi0[c(test.i-1,test.i)]
          
          logL.swap1=0         #the change of log likelihood
          for(k in 1:num.indep.ranker){
              swap.index=x.indep[index,k]
              min.rank=min(swap.index)
              max.rank=max(swap.index)
              if(min.rank<max.rank){
                  phi=op.phi*(1-op.alpha[k]^min.rank)
                  diffv=sum(pi0.rank[x.indep[,k]>min.rank]<pi0.rank[y.indep[min.rank,k]])
                  temp=sum(phi^(1:(num.item-min.rank)))+1
                  bef.logL=log(op.omega[k]*phi^diffv/temp+(1-op.omega[k])/(num.item-min.rank+1))
                  if(swap.index[1]<swap.index[2]){
                      aft.logL=log(op.omega[k]*phi^(diffv+1)/temp+(1-op.omega[k])/(num.item-min.rank+1))
                  }
                  if(swap.index[2]<swap.index[1]){
                      aft.logL=log(op.omega[k]*phi^(diffv-1)/temp+(1-op.omega[k])/(num.item-min.rank+1))
                  }
                  logL.swap1=logL.swap1+aft.logL-bef.logL
              }
          }
            
        
        if(logL.swap1>0){
            
          op.pi0=pi0.prop
          pi0.rank=pi0.prop.rank
          
          #3.2.3 iteratively move the updated item forward until it stops
          
          for(test.j in max(test.i-1,2):2){
            
            pi0.prop=op.pi0
            pi0.prop[c(test.j-1,test.j)]=op.pi0[c(test.j,test.j-1)]
            pi0.prop.rank[pi0.prop]=seq(1,num.item)
            index=op.pi0[c(test.j-1,test.j)]
            
            logL.swap2=0         #the change of log likelihood
            for(k in 1:num.indep.ranker){
              swap.index=x.indep[index,k]
              min.rank=min(swap.index)
              max.rank=max(swap.index)
              if(min.rank<max.rank){
                phi=op.phi*(1-op.alpha[k]^min.rank)
                diffv=sum(pi0.rank[x.indep[,k]>min.rank]<pi0.rank[y.indep[min.rank,k]])
                temp=sum(phi^(1:(num.item-min.rank)))+1
                bef.logL=log(op.omega[k]*phi^diffv/temp+(1-op.omega[k])/(num.item-min.rank+1))
                if(swap.index[1]<swap.index[2]){
                  aft.logL=log(op.omega[k]*phi^(diffv+1)/temp+(1-op.omega[k])/(num.item-min.rank+1))
                  }
                if(swap.index[2]<swap.index[1]){
                  aft.logL=log(op.omega[k]*phi^(diffv-1)/temp+(1-op.omega[k])/(num.item-min.rank+1))
                  }
                logL.swap2=logL.swap2+aft.logL-bef.logL
                }
              }
            
            
            if(logL.swap2>0){
              op.pi0=pi0.prop
              pi0.rank=pi0.prop.rank
             }
            
            if(logL.swap2<0){break}
           }
            
          }
            
       }
        
        if(sum(pi0.curr==op.pi0)==num.item){break}
        
        ############################################################
        
        #estimate kappa
        
        
        for(k in 1:num.kappa){
            
          for(test.i in 1:(stage[k]-1)){

            kappa.prop=op.kappa[,k]
            kappa.prop[c(test.i+1,test.i)]=op.kappa[c(test.i,test.i+1),k]
            kappa.prop.rank=op.kappa.rank[,k]
            kappa.prop.rank[op.kappa[c(test.i,test.i+1),k]]=c(test.i+1,test.i)
            
            index=op.kappa[test.i,k]
            temp1=pi0.rank[index]
            temp2=pi0.rank[(1:num.item)[op.kappa.rank[,k]>=test.i+1]]
            err1=sum(temp2<temp1)
            
            index=op.kappa[test.i+1,k]
            temp1=pi0.rank[index]
            temp2=pi0.rank[(1:num.item)[op.kappa.rank[,k]>=min(test.i+2, num.item)]]
            err2=sum(temp2<temp1)
            
            index=kappa.prop[test.i]
            temp1=pi0.rank[index]
            temp2=pi0.rank[(1:num.item)[kappa.prop.rank>=test.i+1]]
            err3=sum(temp2<temp1)
            
            index=kappa.prop[test.i+1]
            temp1=pi0.rank[index]
            temp2=pi0.rank[(1:num.item)[kappa.prop.rank>=min(test.i+2, num.item)]]
            err4=sum(temp2<temp1)
            
            
            temp3=x[op.kappa[c(test.i,test.i+1),k], is.kappa.ranker[[k]]]
            err5=sum(temp3[2,]-temp3[1,]<0)
            err6=sum(temp3[2,]-temp3[1,]>0)
            
            temp.phi=op.phi*(1-op.alpha[k]^test.i)
            phi.sum=sum(temp.phi^(0:(num.item-test.i)))
            logL.swap1=log(op.omega[k]*temp.phi^err1/phi.sum+(1-op.omega[k])/(num.item-test.i+1))
            logL.swap2=log(op.omega[k]*temp.phi^err3/phi.sum+(1-op.omega[k])/(num.item-test.i+1))
            
            temp.phi=op.phi*(1-op.alpha[k]^(test.i+1))
            phi.sum=sum(temp.phi^(0:(num.item-test.i-1)))
            logL.swap1=logL.swap1+log(op.omega[k]*temp.phi^err2/phi.sum+(1-op.omega[k])/(num.item-test.i))
            logL.swap2=logL.swap2+log(op.omega[k]*temp.phi^err4/phi.sum+(1-op.omega[k])/(num.item-test.i))
            
            delta.logL=logL.swap2-logL.swap1+log(op.phi.h1[k])*(err6-err5)
            
            if(delta.logL>0){
              op.kappa[,k]=kappa.prop
              op.kappa.rank[,k]=kappa.prop.rank
              }
            
          }
            
        }
        
        
        #estimate phi1 for correlated rankings
        
        temp.err=rep(0,num.kappa)
        for(k in 1:num.kappa){
            for(w in is.kappa.ranker[[k]]){
              for(i in 1:num.top[w]){
                index=y[i,w]
                temp1=op.kappa.rank[index,k]
                temp2=op.kappa.rank[x[,w]>i,k]
                temp.err[k]=temp.err[k]+sum(temp2<temp1)
               }
             }
            
            
        for(sub.it in 1:20){
          op.phi.init=op.phi.h1[k]
          temp1.sum=temp2.sum=0
          temp.stage=num.top[is.kappa.ranker[[k]]]
          for(i in 1:max(temp.stage)){
            temp=LogZDeriva(op.phi.h1[k],num.item-i)
            temp1.sum=temp1.sum-temp*sum(temp.stage>=i)
            temp2.sum=temp2.sum-(ZDeriva2(op.phi.h1[k],num.item-i)-temp^2)*sum(temp.stage>=i)
            }

          temp1.sum=temp1.sum+temp.err[k]/op.phi.h1[k]
          temp2.sum=temp2.sum-temp.err[k]/(op.phi.h1[k]^2)
          
          delta.phi=-temp1.sum/temp2.sum    #use Newton-Rapshon method
          if(abs(delta.phi)<0.001){break}
          op.phi.h1[k]=op.phi.h1[k]+delta.phi
          if(op.phi.h1[k]<0){op.phi.h1[k]=abs(op.phi.h1[k])}
          if(op.phi.h1[k]>1){
              if(temp1.sum>0){
                  op.phi.h1[k]=1
                }else{op.phi.h1[k]=op.phi.init}
                
              }
            
          }
        
        }
        
        
        
        # print(it)
        
    }
    
    
    op.kappa.name=matrix(0,num.item,num.kappa)
    for(k in 1:num.kappa){
        op.kappa.name[1:num.indep.top[k],k]=as.character(item.name[op.kappa[1:num.indep.top[k],k]])
    }
    
    
    logL.max=0

    for(k in 1:num.indep.ranker){
      phi.vec=op.phi*(1-op.alpha[k]^(1:stage[k]))
      for(w in 1:stage[k]){
        temp=sum(phi.vec[w]^(0:(num.item-w)))
        logL.max=logL.max+log(op.omega[k]*phi.vec[w]^v.indep[w,k]/temp+(1-op.omega[k])/(num.item-w+1))
        }
      }
    
    
    for(k in 1:num.kappa){
        temp.stage=num.top[is.kappa.ranker[[k]]]
        for(i in 1:max(temp.stage)){
          logL.max=logL.max-log(sum(op.phi.h1[k]^(0:(num.item-i))))*sum(temp.stage>=i)
          }

       logL.max=logL.max+temp.err[k]*log(op.phi.h1[k])
      }

    
    return(list(op.phi=op.phi, op.phi1=op.phi.h1, op.omega=op.omega, op.alpha=op.alpha, op.pi0=item.name[op.pi0], op.kappa=op.kappa.name, logL.max=logL.max))
    
}


