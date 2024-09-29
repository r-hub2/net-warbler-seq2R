transform <-
function(x){

	if(is.list(x)) {x<-x[[1]]}
		n<-length(x)
		nseq<-1
		seqa=rep(NA,n)
		seqt=rep(NA,n)
		seqc=rep(NA,n)
		seqg=rep(NA,n)
		X<-c()

		res=list()
		#x[x!="a"&x!="t"&x!="c"&x!="g"]
		#table(nor[[1]])
		for(i in 1:n){
			if(x[i]=="a"){
				seqa[i]=1
				seqt[i]=0

				}else
			if(x[i]=="t"){
				seqa[i]=0
				seqt[i]=1
				}else
			if(x[i]=="c"){
				seqc[i]=1
				seqg[i]=0
				}else
			if(x[i]=="g"){
				seqc[i]=0
				seqg[i]=1
				}
		}

		AT=matrix(ncol=3,nrow=length(seqt))
		CG=matrix(ncol=3,nrow=length(seqg))


		Xlong_a=sum(na.omit(seqa[1:length(x)]=="0"|seqa[1:length(x)]=="1"))
		Xlong_g=sum(na.omit(seqg[1:length(x)]=="0"|seqg[1:length(x)]=="1"))
		Xlong_t=sum(na.omit(seqt[1:length(x)]=="0"|seqt[1:length(x)]=="1"))
		Xlong_c=sum(na.omit(seqc[1:length(x)]=="0"|seqc[1:length(x)]=="1"))
		X=c(1:length(seqt))
		AT=cbind(X,seqa,seqt)
		AT=na.omit(AT)

		colnames(AT)=c("X","A","T")
		AT=AT[1:Xlong_a,]
		X=c(1:length(seqg))
		CG=cbind(X,seqc,seqg)
		colnames(CG)=c("X","C","G")
		CG=na.omit(CG)
		CG=CG[1:Xlong_g,]
				res[[nseq]]<-list(AT,CG)
		return(list(AT=res[[nseq]][1],CG=res[[nseq]][2]))
}
