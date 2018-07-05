
#MM: Mallows model
#the data has multiple independent rankings
#version 1.0
#last updated in 2018.6.28

#input arguments:
#1. rankings: a n*m data frame, with each column representing a ranking list, which ranks the items from the most preferred to the least preferred. For missing items, use 0 to denote them.
#2. initial.method: the method for initializing the value of pi0, with four options: mean, median, geometric and random (the mean of three randomly sampled ranking lists).


#return arguments:
#1. op.phi: optimal value of phi
#2. op.pi0: optimal value of pi0, ranking the items from the most preferred to the least preferred.
#3. logL.max: maximum value of log-likelihood


MM <- function(rankings, initial.method="mean", it.max=20){


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
    x[x[,k]==0,k]=num.top[k]+1       #the missing rank is set to be the number of top element +1
    }


  #######################################
  #
  #step 2. initialize the parameters
  #
  #######################################


  stage=rep(0,num.ranker)
  for(k in 1:num.ranker){stage[k]=min(num.item-1,num.top[k])}
  v=matrix(0,max(stage),num.ranker)   #number of mismatch pairs


  #2.1 initialize pi0

  score=rep(0,num.item)
  if(initial.method=="mean"){score=apply(x, 1, mean)}
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

  #2.2 initialize phi, alpha, omega and phi.h1, tau

  op.phi=0.2

  # the derivative of logZ(theta), Z(theta)=1+theta+...+theta^n

  LogZDeriva <- function(theta,n){
    temp1=sum(theta^(1:n))+1
    temp2=sum(c(1, theta^(1:(n-1)))*(1:n))
    result=temp2/temp1
    return(result)
    }


  # Z''(theta)/Z(theta), Z(theta)=1+theta+...+theta^n

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



  ########################################################################
  #
  # Step 3. find the MLE of the parameters using EM method
  #
  ########################################################################

  for(it in 1:it.max){

      #3.1 calculate v

      for(k in 1:num.ranker){
        for(i in 1:stage[k]){
          index=y[i,k]
          temp1=pi0.rank[index]
          temp2=pi0.rank[x[,k]>i]
          v[i,k]=sum(temp2<temp1)
          }
        }


    #3.2 update optimal phi and pi0


    #3.2.1 update phi


      v.total=sum(v)
      for(sub.it in 1:20){
        op.phi.init=op.phi
        temp1.sum=temp2.sum=0
        for(i in 1:stage[k]){
          temp=LogZDeriva(op.phi,num.item-i)
          temp1.sum=temp1.sum-temp*sum(stage>=i)
          temp2.sum=temp2.sum-(ZDeriva2(op.phi,num.item-i)-temp^2)*sum(stage>=i)
          }

        temp1.sum=temp1.sum+v.total/op.phi
        temp2.sum=temp2.sum-v.total/(op.phi^2)

        delta.phi=-temp1.sum/temp2.sum    #use Newton-Rapshon method
        if(abs(delta.phi)<0.001){break}
        op.phi=op.phi+delta.phi
        if(op.phi<0){
          op.phi=abs(op.phi)
          }
        if(op.phi>1){
          if(temp1.sum>0){
            op.phi=1
            }else{
            op.phi=op.phi.init
            }
          }

        }


    #3.2.2 iteratively update op.pi0 by swapping pairs

    pi0.curr=op.pi0
    for(test.i in 2:num.item){

      pi0.prop=op.pi0
      pi0.prop[c(test.i-1,test.i)]=op.pi0[c(test.i,test.i-1)]
      pi0.prop.rank[pi0.prop]=seq(1,num.item)
      index=op.pi0[c(test.i-1,test.i)]

      logL.swap1=0         #the change of log likelihood
      for(k in 1:num.ranker){
        swap.index=x[index,k]
        min.rank=min(swap.index)
        max.rank=max(swap.index)
        if(min.rank<max.rank){
          if(swap.index[1]<swap.index[2]){delta.logL=log(op.phi)}
          if(swap.index[2]<swap.index[1]){delta.logL=-log(op.phi)}
          logL.swap1=logL.swap1+delta.logL
        }
      }


      if(logL.swap1>0){

        op.pi0=pi0.prop
        pi0.rank=pi0.prop.rank

        #3.2.2.1 iteratively move the updated item forward until it stops

        for(test.j in max(test.i-1,2):2){

          pi0.prop=op.pi0
          pi0.prop[c(test.j-1,test.j)]=op.pi0[c(test.j,test.j-1)]
          pi0.prop.rank[pi0.prop]=seq(1,num.item)
          index=op.pi0[c(test.j-1,test.j)]

          logL.swap2=0         #the change of log likelihood
          for(k in 1:num.ranker){
            swap.index=x[index,k]
            min.rank=min(swap.index)
            max.rank=max(swap.index)
            if(min.rank<max.rank){
              if(swap.index[1]<swap.index[2]){delta.logL=log(op.phi)}
              if(swap.index[2]<swap.index[1]){delta.logL=-log(op.phi)}
              logL.swap2=logL.swap2+delta.logL
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
  }


  logL.max=0
  for(i in 1:stage[k]){
    temp=log(sum(op.phi^(0:(num.item-i))))
    logL.max=logL.max-temp*sum(stage>=i)
    }

  logL.max=logL.max+v.total*log(op.phi)

  return(list(op.phi=op.phi, op.pi0=item.name[op.pi0], logL.max=logL.max))

}


