critical <-
function(model,base.pairs=NULL){
	num_bp<- length(unique(c(model$n_AT,model$n_CG)))
	model$critical_AT[model$critical_AT==0]=NA
	model$lcritical_AT[model$lcritical_AT==0]=NA
	model$ucritical_AT[model$ucritical_AT==0]=NA
	model$critical_CG[model$critical_CG==0]=NA
	model$lcritical_CG[model$lcritical_CG==0]=NA
	model$ucritical_CG[model$ucritical_CG==0]=NA
	
	for(i in 1:length(model$critical_AT)){
		if(model$lcritical_AT[i]==model$ucritical_AT[i]) {
		model$lcritical_AT[i]=NA
		model$ucritical_AT[i]=NA
		model$critical_AT[i]=NA}}
		
		model$lcritical_AT=na.omit(model$lcritical_AT)
		model$ucritical_AT=na.omit(model$ucritical_AT)
		model$critical_AT=na.omit(model$critical_AT)
		
	for(i in 1:length(model$critical_CG)){
		if(model$lcritical_CG[i]==model$ucritical_CG[i]){
		model$lcritical_CG[i]=NA
		model$ucritical_CG[i]=NA
		model$critical_CG[i]=NA	}}		
		
		model$lcritical_CG=na.omit(model$lcritical_CG)
		model$ucritical_CG=na.omit(model$ucritical_CG)
		model$critical_CG=na.omit(model$critical_CG)
	
	
	
if(is.null(base.pairs)){
		res=matrix(NA,ncol=3,nrow=if(length(model$critical_AT)==0){1}else{length(model$critical_AT)})
			if(length(model$critical_AT)!=0){
			ii=8
			res[,1]=model[ii][[1]]
			res[,2]=model[ii+1][[1]]
			res[,3]=model[ii+2][[1]]
			}	
		colnames(res)=c("Critical","Low_CI","Up_CI")
 		
		res2=matrix(NA,ncol=3,nrow=if(length(model$critical_CG)==0){1}else{length(model$critical_CG)})
			if(length(model$critical_CG)!=0){
			iii=18
			res2[,1]=model[iii][[1]]
			res2[,2]=model[iii+1][[1]]
			res2[,3]=model[iii+2][[1]]
			}	
		colnames(res2)=c("Critical","Low_CI","Up_CI")
		return(list(AT=res,CG=res2))
	
}else if(base.pairs=="AT"){
		res=matrix(NA,ncol=3,nrow=if(length(model$critical_AT)==0){1}else{length(model$critical_AT)})
		if(length(model$critical_AT)!=0){
			ii=8
			res[,1]=model[ii][[1]]
			res[,2]=model[ii+1][[1]]
			res[,3]=model[ii+2][[1]]	
			}
		colnames(res)=c("Critical","Low_CI","Up_CI")
		return(res)
				
}else if(base.pairs=="CG"){
		res2=matrix(NA,ncol=3,nrow=if(length(model$critical_CG)==0){1}else{length(model$critical_CG)})
		if(length(model$critical_CG)!=0){
			iii=18
			res2[,1]=model[iii][[1]]
			res2[,2]=model[iii+1][[1]]
			res2[,3]=model[iii+2][[1]]
			}
		colnames(res2)=c("Critical","Low_CI","Up_CI")
		return(res2)
	}		
}
