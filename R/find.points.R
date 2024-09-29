find.points <-
  function(x, kbin = 300, p = 3, bandwidth = -1, weights = 1, nboot = 100,
           kernel = "gaussian", n.bandwidths = 20, seed = NULL, ...){

	h <- bandwidth
	W <- weights
	nh <- n.bandwidths


  etiquetas<-unique(names(x))
	if(kernel=="gaussian")  kernel=2
	if(kernel=="epanech")      kernel=1
	if(kernel=="triang")     kernel=3
	if(length(h)==1) h=rep(h,2)
	t = matrix(ncol=2,nrow=kbin)
	m = matrix(ncol=2,nrow=kbin)
	lCI = matrix(ncol=2,nrow=kbin)
	uCI = matrix(ncol=2,nrow=kbin)
	m1 = matrix(ncol=2,nrow=kbin)
	lCI1 = matrix(ncol=2,nrow=kbin)
	uCI1 = matrix(ncol=2,nrow=kbin)
	critical = matrix(ncol=2,nrow=kbin)
	lcritical= matrix(ncol=2,nrow=kbin)
	ucritical= matrix(ncol=2,nrow=kbin)
	#mar = matrix(ncol=2,nrow=100)

	N=c()
	hhA=c()
	hhT=c()

	if (!is.null(seed)) {
	  set.seed(seed)
	}



	for(i in 1:2){

		A=x[[i]][[1]][,2] # sirve para Adeninas o para Guaninas
		T=x[[i]][[1]][,3] # sirve para Timinas o para Citosinas
		X=x[[i]][[1]][,1]
		n=length(X)
		umatrixA <- matrix(runif(n*nboot), ncol = nboot, nrow = n)
		umatrixT <- matrix(runif(n*nboot), ncol = nboot, nrow = n)
		#Y=rep(1,n)
		W=rep(1,n)
		#if(W==1){W=rep(1,n)}
	change.points<-.Fortran("change.points_",
			n = as.integer(n),
			hA= as.double(h[i]),
			hT=as.double(h[i]),
			p = as.integer(p),
			A = as.double(A),
			T = as.double(T),
			#Y = as.double(Y),
			X = as.double(X),
			W = as.double(W),
			Xb = as.double(rep(-1.0,kbin)),
			AT = as.double(rep(-1.0,kbin)),
			ATi = as.double(rep(-1.0,kbin)),
			ATs = as.double(rep(-1.0,kbin)),
			AT1 = as.double(rep(-1.0,kbin)),
			AT1i = as.double(rep(-1.0,kbin)),
			AT1s = as.double(rep(-1.0,kbin)),
			critical = as.double(rep(-1.0,kbin)),
			criticali = as.double(rep(-1.0,kbin)),
			criticals = as.double(rep(-1.0,kbin)),
			nboot=as.integer(nboot),
			kbin=as.integer(kbin),
			kernel=as.integer(kernel),
			nh=as.integer(nh),
			umatrixA <- array(umatrixA, c(n, nboot)),
			umatrixT <- array(umatrixT, c(n, nboot)),
			PACKAGE = "seq2R"
			#marta=as.double(rep(-1.0,100))
			#atboot=array(rep(-1.0),c(kbin,nboot)),
			)
	#mar[,i]<-change.points$marta
	t[,i]<-change.points$Xb
	m[,i]<-change.points$AT
	lCI[,i]<-change.points$ATi
	uCI[,i]<-change.points$ATs
	m1[,i]<-change.points$AT1
	lCI1[,i]<-change.points$AT1i
	uCI1[,i]<-change.points$AT1s
	critical[,i]<-change.points$critical
	lcritical[,i]<-change.points$criticali
	ucritical[,i]<-change.points$criticals
	N[i]<-change.points$n
	hhA[i]<-change.points$hA
	hhT[i]<-change.points$hT
	}

	hAT=hhA[1]
	hCG=hhA[2]
	#hCG=c(hhA[2],hhT[2])
	#hAT=c(hhA[1],hhT[1])
	res <- list(t_AT = t[,1],# t es el tiempo
			m_AT = m[,1],
			lCI_AT = lCI[,1],
			uCI_AT = uCI[,1],
			m1_AT = m1[,1],
			lCI1_AT = lCI1[,1],
			uCI1_AT = uCI1[,1],
			critical_AT = critical[,1][critical[,1]!=-1],
			lcritical_AT = lcritical[,1][lcritical[,1]!=-1],
			ucritical_AT = ucritical[,1][ucritical[,1]!=-1],
			t_CG = t[,2],
			m_CG = m[,2],
			lCI_CG = lCI[,2],
			uCI_CG = uCI[,2],
			m1_CG = m1[,2],
			lCI1_CG = lCI1[,2],
			uCI1_CG = uCI1[,2],
			critical_CG = critical[,2][critical[,2]!=-1],
			lcritical_CG = lcritical[,2][lcritical[,2]!=-1],
			ucritical_CG = ucritical[,2][ucritical[,2]!=-1],
			nboot=change.points$nboot,
			kbin=change.points$kbin,
			n_AT=N[1],
			n_CG=N[2],
			h=c(hAT,hCG),
			etiquetas=as.character(etiquetas),
			call=match.call())

		class(res) <- "change.points"
		return(res)

}
