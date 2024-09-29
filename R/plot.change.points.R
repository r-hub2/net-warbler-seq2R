plot.change.points <-
function(x=model,y=NULL,base.pairs=NULL,der=NULL,xlab="x",ylab="y",col="black",CIcol="black",main=NULL,type="l",CItype="l",critical=FALSE,CIcritical=FALSE,ylim=NULL,...)
{
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  model=x
	if(is.null(der)) der=5
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

	if (is.null(base.pairs) & der==5){



      par(mfrow=c(2,2))
			title=main
			if(is.null(main)) title=paste(model$etiquetas[1],"skew")
			plot(model$t_AT,model$m_AT,ylim=c(min(model$lCI_AT),max(model$uCI_AT)),type=type,xlab=xlab,ylab=ylab,col=col,main=title)
			lines(model$t_AT,model$lCI_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_AT,model$uCI_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
					if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
						graphics.off()
						stop("No exist any critical poits")
					}else{for(i in 1:length(model$critical_AT)){
					if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
					abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}}}
			if(CIcritical==TRUE){
					if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
						graphics.off()
						stop("No exist any critical poits")
					}else{for(i in 1:length(model$critical_AT)){
					if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
					lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],model$m_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],col=color,lwd=2.5)}}}
			title=main
			if(is.null(main)) title=paste(model$etiquetas[2],"skew")
			plot(model$t_CG,model$m_CG,ylim=c(min(model$lCI_CG),max(model												$uCI_CG)),type=type,xlab=xlab,ylab=ylab,col=col,main=title)
			lines(model$t_CG,model$lCI_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_CG,model$uCI_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)#abline(h=0,col="grey",lty=2)
			if(critical==TRUE){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model$m_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],col=color,lwd=2.5)}}


			title=main
			if(is.null(main)) title=paste(model$etiquetas[1],"skew")
			plot(model$t_AT,model$m1_AT,ylim=c(min(model$lCI1_AT),max(model$uCI1_AT)),type=type,xlab=xlab,ylab="First derivative",col=col,main=title)
			lines(model$t_AT,model$lCI1_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_AT,model$uCI1_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],model					$m1_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],col=color,lwd=2.5)}}
				abline(h=0,col="grey10",lwd=0.2)


			title=main
			if(is.null(main)) title=paste(model$etiquetas[2],"skew")
			plot(model$t_CG,model$m1_CG,ylim=c(min(model$lCI1_CG),max(model$uCI1_CG)),type=type,xlab=xlab,ylab="First derivative",col=col,main=title,...)
			lines(model$t_CG,model$lCI1_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_CG,model$uCI1_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$lcritical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model					$m1_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG										[i]],col=color,lwd=2.5)}}
			abline(h=0,col="grey10",lwd=0.2)

			}else if(is.null(base.pairs)&der == 0){


		par(mfrow=c(2,1))
			title=main
			if(is.null(main)) title=paste(model$etiquetas[1],"skew")
			#plot(model$t_AT,model$m_AT,ylim=c(min(model$m_AT),max(model$m_AT)),type=type,xlab=xlab,ylab=ylab,col=col,main=title,...)
			plot(model$t_AT,model$m_AT,ylim=c(min(model$lCI_AT),max(model$uCI_AT)),type=type,xlab=xlab,ylab=ylab,col=col,main=title,...)
			lines(model$t_AT,model$lCI_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_AT,model$uCI_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
					graphics.off()
					stop("No exist any critical poits")
				}else{for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}}}
			if(CIcritical==TRUE){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
				graphics.off()
				stop("No exist any critical poits")
				}else{for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],model					$m_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT										[i]],col=color,lwd=2.5)}}}

			title=main
			if(is.null(main)) title=paste(model$etiquetas[2],"skew")
			plot(model$t_CG,model$m_CG,ylim=c(min(model$lCI_CG),max(model$uCI_CG)),type=type,xlab=xlab,ylab=ylab,col=col,main=title,...)
			lines(model$t_CG,model$lCI_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_CG,model$uCI_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model					$m_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG										[i]],col=color,lwd=2.5)}}


	}else if(is.null(base.pairs) & der == 1){


		par(mfrow=c(2,1))

			title=main
			if(is.null(main)) title=paste(model$etiquetas[1],"skew")
			plot(model$t_AT,model$m1_AT,ylim=c(min(model$lCI1_AT),max(model$uCI1_AT)),type=type,xlab=xlab,ylab="First derivative",col=col,main=title,...)
			lines(model$t_AT,model$lCI1_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_AT,model$uCI1_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
					graphics.off()
					stop("No exist any critical poits")
				}else{for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}}}
			if(CIcritical==TRUE){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
					graphics.off()
					stop("No exist any critical poits")
				}else{for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],model					$m1_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],col=color,lwd=2.5)}}}
				abline(h=0,col="grey10",lwd=0.2)


			title=main
			if(is.null(main)) title=paste(model$etiquetas[2],"skew")
			plot(model$t_CG,model$m1_CG,ylim=c(min(model$lCI1_CG),max(model$uCI1_CG)),type=type,xlab=xlab,ylab="First derivative",col=col,main=title,...)
			lines(model$t_CG,model$lCI1_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_CG,model$uCI1_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$lcritical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model					$m1_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG										[i]],col=color,lwd=2.5)}}
			abline(h=0,col="grey10",lwd=0.2)




	}else if(der==5&base.pairs=="AT"){


		par(mfrow=c(2,1))
			title=main
			if(is.null(main)) title=paste(model$etiquetas[1],"skew")
			plot(model$t_AT,model$m_AT,ylim=c(min(model$lCI_AT),max(model$uCI_AT)),type=type,xlab=xlab,ylab=ylab,col=col,main=title,...)
			lines(model$t_AT,model$lCI_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_AT,model$uCI_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
				graphics.off()
				stop("No exist any critical poits")
				}else{for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$critical_AT)){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
					graphics.off()
					stop("No exist any critical poits")
				}else{if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],model					$m_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT										[i]],col=color,lwd=2.5)}}}

			title=main
			if(is.null(main)) title=paste(model$etiquetas[1],"skew")
			plot(model$t_AT,model$m1_AT,ylim=c(min(model$lCI1_AT),max(model$uCI1_AT)),type=type,xlab=xlab,ylab="First derivative",col=col,main=title,...)
			lines(model$t_AT,model$lCI1_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_AT,model$uCI1_AT,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],model					$m1_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],col=color,lwd=2.5)}}
				abline(h=0,col="grey10",lwd=0.2)



	}else if(der==5 &base.pairs=="CG"){

		par(mfrow=c(2,1))
			title=main
			if(is.null(main)) title=paste(model$etiquetas[2],"skew")
			plot(model$t_CG,model$m_CG,ylim=c(min(model$lCI_CG),max(model$uCI_CG)),type=type,xlab=xlab,ylab=ylab,col=col,main=title,...)
			lines(model$t_CG,model$lCI_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_CG,model$uCI_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
					graphics.off()
					stop("No exist any critical poits")
				}else{for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}}}
			if(CIcritical==TRUE){
				if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){
					graphics.off()
					stop("No exist any critical poits")
				}else{for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model					$m_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG										[i]],col=color,lwd=2.5)}}}

			title=main
			if(is.null(main)) title=paste(model$etiquetas[2],"skew")
			plot(model$t_CG,model$m1_CG,ylim=c(min(model$lCI1_CG),max(model$uCI1_CG)),type=type,xlab=xlab,ylab="First derivative",col=col,main=title,...)
			lines(model$t_CG,model$lCI1_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			lines(model$t_CG,model$uCI1_CG,lty=6,lwd=0.5,col=CIcol,type=CItype)
			if(critical==TRUE){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}}
			if(CIcritical==TRUE){
				for(i in 1:length(model$lcritical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model					$m1_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG										[i]],col=color,lwd=2.5)}}
			abline(h=0,col="grey10",lwd=0.2)

		}else{
if(base.pairs=="AT" & der==0){
	ii=c(2)
	iix=1
	ylab2=ylab
	if(is.null(ylim)) {ylim2=c(min(model[ii+1][[1]]),max(model[ii+2][[1]]))}else{ylim2=ylim}	}

if(base.pairs=="CG" & der==0){
	ii=c(12)
	iix=11
	ylab2=ylab
	if(is.null(ylim)) {ylim2=c(min(model[ii+1][[1]]),max(model[ii+2][[1]]))}else{ylim2=ylim}	}

if(base.pairs=="AT" & der==1){
	ii=c(5)
	iix=1
	ylab2="First derivative"
	if(is.null(ylim)) {ylim2=c(min(model[ii+1][[1]]),max(model[ii+2][[1]]))}else{ylim2=ylim}	}

if(base.pairs=="CG" & der==1){
	ii=c(15)
	iix=11
	ylab2="First derivative"
	if(is.null(ylim)) {ylim2=c(min(model[ii+1][[1]]),max(model[ii+2][[1]]))}else{ylim2=ylim}	}

		plot(model[iix][[1]],model[ii][[1]],ylim=ylim2,type=type,xlab=xlab,ylab=ylab2,col=col)
		lines(model[iix][[1]],model[ii+1][[1]],lty=6,lwd=0.5,col=CIcol,type=CItype)
		lines(model[iix][[1]],model[ii+2][[1]],lty=6,lwd=0.5,col=CIcol,type=CItype)
		if(critical==TRUE){
			if(base.pairs=="AT" & der==0 & length(model$critical_AT)==0){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="AT" & der==0 & length(model$critical_AT)>0){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}
			}else if(base.pairs=="CG" & der==0 & length(model$critical_CG)==0){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="CG" & der==0 & length(model$critical_CG)>0){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}
			}else if(base.pairs=="AT" & der==1 & length(model$critical_AT)==0){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="AT" & der==1 & length(model$critical_AT)>0){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				abline(v=model$t_AT[model$t_AT==model$critical_AT[i]],lty=2,col=color)}
			}else if(base.pairs=="CG" & der==1  & length(model$critical_CG)==0 ){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="CG" & der==1  & length(model$critical_CG)>0 ){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				abline(v=model$t_CG[model$t_CG==model$critical_CG[i]],lty=2,col=color)}}		}

		if(CIcritical==TRUE){
			if(base.pairs=="AT" & der==0 &  length(model $critical_AT)==0){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="AT" & der==0 &  length(model $critical_AT)>0){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT							[i]],model$m_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT							[i]],col=color,lwd=2.5)}
			}else if(base.pairs=="CG" & der==0 & length(model$critical_CG)==0){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="CG" & der==0 & length(model$critical_CG)>0){
				for(i in 1:length(model$critical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model					$m_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG										[i]],col=color,lwd=2.5)}
			}else if(base.pairs=="AT" & der==1 &  length(model $critical_AT)==0){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="AT" & der==1 &  length(model $critical_AT)>0){
				for(i in 1:length(model$critical_AT)){
				if(model$m1_AT[model$t_AT==model$critical_AT[i]]>0){color=4}else{color=2}
				lines(model$t_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT[i]],model					$m1_AT[model$t_AT>=model$lcritical_AT[i]&model$t_AT<=model$ucritical_AT										[i]],col=color,lwd=2.5)}
			}else if(base.pairs=="CG" & der==1 & length(model$critical_CG)==0){
				graphics.off()
				stop("No exist any critical poits")
			}else if(base.pairs=="CG" & der==1 & length(model$critical_CG)>0){
				for(i in 1:length(model$lcritical_CG)){
				if(model$m1_CG[model$t_CG==model$critical_CG[i]]>0){color=4}else{color=2}
				lines(model$t_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG[i]],model					$m1_CG[model$t_CG>=model$lcritical_CG[i]&model$t_CG<=model$ucritical_CG										[i]],col=color,lwd=2.5)}
				}}
}}
