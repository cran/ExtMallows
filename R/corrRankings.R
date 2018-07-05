
#function: corrRankings(), which caclulates the p values that measure the correlation of pairwise rankings.
#Note that the input rankings should have at least 8 rankings. When constructing the samples of rescaled V distance for a given rank position, the number of samples should at least be 28 and the number of rankings that have items up to this position should account for at least 2/3 of the total number of rankings, otherwise the p value calculation stops at this position.

#input argument:
#rankings: a n*m data frame, with each column representing a ranking list, which ranks the items from the most preferred to the least preferred. For missing items, use 0 to denote them.

#return argument:
#pair.pvalue: a symmetric matrix of p values, with the (i,j)-th element denoting the p value of the i,j-th rankings.


corrRankings <- function(rankings){
  
  rank.list=list() 
  num.ranker=ncol(rankings)
  if(num.ranker<8){
    stop("There are too few rankings.")
  }
  num.top=rep(0,num.ranker)
  for(k in 1:num.ranker){
    temp=rankings[,k]
    rank.list[[k]]=temp[temp!=0]
    num.top[k]=length(rank.list[[k]])
  }

  
  item.name=as.character(sort(unique(unlist(rank.list))))
  n=num.item=length(item.name)
  item.name.ref=seq(1,num.item)
  names(item.name.ref)=item.name
  d=max(num.top)
  y=matrix(0,d,num.ranker)
  for(k in 1:num.ranker){
    y[1:num.top[k],k]=item.name.ref[as.character(rank.list[[k]])]
  }
  
  
  x=matrix(0,num.item,num.ranker)
  rownames(x)=item.name
  for(k in 1:num.ranker){
    x[y[1:num.top[k],k],k]=1:num.top[k]
    x[x[,k]==0,k]=num.top[k]+1
  }
    
  
  pair.rescaled.dist.arr=array(0, dim=c(num.ranker,num.ranker,max(num.top)))
  for(i in 1:(num.ranker-1)){
    for(j in (i+1):num.ranker){
      d=min(num.top[i], num.top[j], n-1)
      rescaled.dist=rep(0,d)
      vec1=x[,i]
      vec2=x[,j]
      for(k in 1:d){
        temp1=seq(1,n)[vec1==k]
        temp2=seq(1,n)[vec1>k]
        temp3=vec2[temp1]
        temp4=vec2[temp2]
        temp5=temp4-temp3
        temp6=sum(temp5<0)  #incorrect
        temp7=sum(temp5>0)  #correct
        incorrect1=temp6/(temp6+temp7)

        temp1=seq(1,n)[vec2==k]
        temp2=seq(1,n)[vec2>k]
        temp3=vec1[temp1]
        temp4=vec1[temp2]
        temp5=temp4-temp3
        temp6=sum(temp5<0)  #incorrect
        temp7=sum(temp5>0)  #correct
        incorrect2=temp6/(temp6+temp7)

        rescaled.dist[k]=(incorrect1+incorrect2)/2
       }
      pair.rescaled.dist.arr[i,j,1:d]=rescaled.dist
    }
  }
  
  
  pooled.dist.sample=list()
  num.pair=num.ranker*(num.ranker-1)/2
  temp=round(2/3*num.ranker)
  min.comp.pair=temp*(temp-1)/2
  comp.d=min(max(num.top),n-1)
  for(w in 1:comp.d){
    count=1
    temp.count=rep(0,num.pair)
    for(i in 1:(num.ranker-1)){
      for(j in (i+1):num.ranker){
        d=min(num.top[i], num.top[j],n-1)
        if(d>=w){
          temp.count[count]=pair.rescaled.dist.arr[i,j,w]
          count=count+1
        }
      }
    }
    if(count>=28 & count>=min.comp.pair){pooled.dist.sample[[w]]=temp.count[1:(count-1)]}
    else{break}
  }
  
  
  sample.d=length(pooled.dist.sample)
  pair.pvalue.arr=array(0, dim=c(num.ranker,num.ranker,max(num.top)))
  for(i in 1:(num.ranker-1)){
    for(j in (i+1):num.ranker){
      d=min(num.top[i], num.top[j],sample.d)
      for(w in 1:d){
        temp1=pooled.dist.sample[[w]]
        temp2=sum(temp1<=pair.rescaled.dist.arr[i,j,w])
        pair.pvalue.arr[i,j,w]=temp2/length(temp1)
      }
    }
  } 
  
  pair.pvalue=matrix(0,num.ranker,num.ranker)
  for(i in 1:(num.ranker-1)){
    for(j in (i+1):num.ranker){
      temp1=pair.pvalue.arr[i,j,]
      temp2=temp1[temp1>0]
      temp3=length(temp2)
      temp4=sum(-log(temp2))
      pair.pvalue[i,j]=1-pgamma(temp4, shape=temp3, scale=1)
    }
  }  
  
  
  pair.pvalue=pair.pvalue+t(pair.pvalue)
  
  return(pair.pvalue)
  
}


