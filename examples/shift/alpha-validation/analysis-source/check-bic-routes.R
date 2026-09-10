library(kfl1ou)
tr <- ape::read.tree(text='((a:1,b:1):1,(c:1,d:1):1);')
y <- matrix(c(1,1.3,2.2,1.7),ncol=1,dimnames=list(tr$tip.label,'trait'))
maxgap <- 0
for(root in c('OUfixedRoot','OUrandomRoot')) for(se in c(0,.2)) for(floor in c(1e-7,.1)) {
 data <- adjust_data(tr,y,normalize=FALSE,quietly=TRUE,repair.tree=FALSE,drop.all.missing=FALSE,drop.invariant=FALSE)
 edge <- match(1L,data$tree$edge[,2])
 settings <- list(tree=data$tree,Y=data$Y,shift.configuration=edge,cr.regimes=list(0L,edge),root.model=root,
   alpha.lower=floor/2,alpha.upper=5,alpha.starting.value=.5,compute.hessian=FALSE,
   input_error=if(se==0) NULL else matrix(se^2,nrow=4,ncol=1,dimnames=dimnames(data$Y)))
 a<-suppressWarnings(do.call(fit_OU,c(settings,list(criterion='pBIC'))))
 b<-suppressWarnings(do.call(fit_OU,c(settings,list(criterion='BIC'))))
 gap<-abs(BIC(a)-b$score); maxgap<-max(maxgap,gap)
 if(gap>1e-6)stop('BIC route mismatch')
}
cat('8 BIC-route checks passed; max gap',format(maxgap,digits=17),'\n')
