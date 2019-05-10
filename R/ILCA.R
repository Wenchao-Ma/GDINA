#' Iterative latent-class analysis
#'
#' This function implements an iterative latent class analysis (ILCA; Jiang, 2019) approach to estimating attributes for cognitive diagnosis.
#'
#' @param dat A required binary item response matrix.
#' @param Q A required binary item and attribute association matrix.
#' @param seed.num seed number
#' @return Estimated attribute profiles.
#'
#' @author {Zhehan Jiang, The University of Alabama}
#' @references
#' Jiang, Z. (2019). Using the iterative latent-class analysis approach to improve attribute accuracy in diagnostic classification models. \emph{Behavior research methods}, 1-10.
#' @export
#' @examples
#' ILCA(sim10GDINA$simdat, sim10GDINA$simQ)
#'
ILCA <- function(dat,Q,seed.num=5){
  ###Dependencies
  if (!requireNamespace(c("stringr", "poLCA"), quietly = TRUE)) {
    stop("stringr and poLCA needed for ILCA. Please install them.",
         call. = FALSE)
  }

  Na<-ncol(Q)
  Np<-nrow(dat)
  Ni<-ncol(dat)

  #Initial Values and Seeds
  mg.est <- GDINA(dat = dat, Q = Q, verbose = 0)
  posteriorAttr<-personparm(mg.est)
  class(posteriorAttr) <- "numeric"

  foursets_rate<-rep(0.5,(seed.num))
  foursets_seed<-vector(mode='list',length=(seed.num))
  ini_posteriorAttr<-posteriorAttr
  for ( i in 1:Ni){
    for(loop.a in 1:Na){
      sel.item<-unlist(lapply(1:Na,function(x){which(Q[,x]==1)})[[loop.a]])
      ini_posteriorAttr[,loop.a]  <- rbinom(Np,1,(apply(dat[,sel.item],1,mean)))
    }
    foursets_seed[[i]]<-ini_posteriorAttr
  }
  foursets_seed[[i+1]]<-posteriorAttr

  LCAgo<-function(posteriorAttr){
    Aname<-paste('A',1:Na,sep='')
    foursets_ini<-vector(mode='list',length=(Na))

    for(loop.a in 1:Na){
      sel.item<-unlist(lapply(1:Na,function(x){which(Q[,x]==1)})[[loop.a]])
      foursets_ini[[loop.a]]<-as.data.frame(dat[,sel.item])
    }
    IniAttr<-posteriorAttr
    NewAttr<-cbind(IniAttr[,-1],IniAttr[,1])
    iter_count<-1
    while(sum(IniAttr==NewAttr)!=Np*Na|iter_count==500){
      NewAttr<-IniAttr
      iter_count<-iter_count+1
      for(loop.a in 1:Na){
        #For an attribute loop.a
        group1id<-which(IniAttr[,loop.a]==1)
        group1<-cbind(IniAttr[IniAttr[,loop.a]==1,-loop.a],foursets_ini[[loop.a]][,
                                                                                  sample(1:ncol(foursets_ini[[loop.a]]),round(
                                                                                    ncol(foursets_ini[[loop.a]])/2
                                                                                  ))][IniAttr[,loop.a]==1,])
        colnames(group1)<- c(Aname[-loop.a],paste('x',1:(ncol(group1)-(Na-1)),sep=''))
        group1<-as.data.frame(group1)

        group0id<-which(IniAttr[,loop.a]==0)
        group0<-cbind(IniAttr[IniAttr[,loop.a]==0,-loop.a],foursets_ini[[loop.a]][,
                                                                                  sample(1:ncol(foursets_ini[[loop.a]]),round(
                                                                                    ncol(foursets_ini[[loop.a]])/2
                                                                                  ))][IniAttr[,loop.a]==0,])

        colnames(group0)<- c(Aname[-loop.a],paste('x',1:(ncol(group0)-(Na-1)),sep=''))
        group0<-as.data.frame(group0)

        group<-cbind(c(group1id,group0id),rbind(group1,group0))
        group<-group[order(group[,1]),-1]

        group<-apply( group,2,as.factor)

        result <- poLCA::poLCA(eval(parse(text=paste(
          'f<-cbind(',paste(colnames(group),collapse=','),')~1',collapse='')     )),
          as.data.frame(group),nclass=2,nrep=3,verbose =F)

        current.loop.a<-result$predclass-1
        IniAttr[,loop.a]<-current.loop.a
        assign('IniAttr',IniAttr)
      }
    }
    IniAttr
  }

  for(seed.loop in 1:seed.num){
    assign(paste('seed',seed.loop,sep=''),tryCatch(LCAgo(foursets_seed[[seed.loop]]),error = function(e){posteriorAttr}))
  }

  finalVote<-posteriorAttr
  for(i in 1:nrow(posteriorAttr)){
    for(j in 1:ncol(posteriorAttr)){
      for (seed.loop in 1:seed.num){
        finalVote[i,j]<-finalVote[i,j]+eval(parse(text=paste('seed',seed.loop,"[",i,",",j,"]",sep='')))
      }
    }
  }
  finalVote[finalVote<=(seed.num/2)]<-0;finalVote[finalVote>(seed.num/2)]<-1
  finalVote
}


