print.change.points<-
function(x=model,...){
		model<-x
		cat("\nCall:\n")
		print(model $call)
		
		cat("\nNumber of A-T base pairs:")
		cat(format(model $n_AT))
		cat("\n")
		
		cat("\nNumber of C-G base pairs:")
		cat(format(model $n_CG))
		cat("\n")
	
		cat("\nNumber of binning nodes: ")
		cat(format(model $kbin))
		cat("\n")
		
		cat("\nNumber of bootstrap repeats: ")
		cat(format(model $nboot))
		cat("\n")
		
		cat("\nBandwidth: ")
		cat(format(model $h))
		cat("\n")
		
		cat("\nExists any critical point? ",if(length(model $critical_AT)==0&(length(model $critical_CG)==0)){"FALSE"}else{"TRUE"},"\n")
		cat("\n")
		
		}
