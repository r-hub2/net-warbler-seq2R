read.all <-
function(file=system.file(""),seqtype="DNA"){
	num=nchar(file)
	ii=substr(file,num,num)
	
	if(ii=="k"){
		
		input <- readLines(file)
    	head <- input[1]
   		head <- unlist(strsplit(head, split = " "))
    	head <- head[nchar(head) > 0]
    	seqname <- head[2]
    	seqsize <- as.integer(head[3])
    	outheader <- sprintf(">%s %d bp", seqname, seqsize)
    	sequences=list()
		sequences[[2]]=outheader 
	    debut <- which(substring(input, 1, 6) == "ORIGIN") + 1
   		 	if (length(debut) > 1)  stop("Multiple entries not yet implemented !")
    	fin <- which(substring(input, 1, 2) == "//") - 1
    		if (length(fin) > 1) stop("Multiple entries not yet implemented !")
    	input <- input[debut:fin]
		aux=unlist(strsplit(input, split = " "))
		aux=aux[nchar(aux) == 10]
		sequences[[1]]=paste(substr(aux,1,10),collapse="")
		sequences[[1]] <- lapply(sequences[[1]], s2c)
		seq=list()
		seq[[1]]=sequences[[1]][[1]]
		seq[[2]]=sequences[[2]]
		return(seq)
		}			
		else{
		lines<-readLines(file)
		ind<-which(substr(lines, 1L, 1L) == ">")
		nseq<-length(ind)
			if(nseq==0)  stop("no line starting with a > chatacter found")
		start<-ind + 1
		end<-ind - 1
		end <- c(end[-1], length(lines))
  		sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]],collapse = ""))
  		
  		sequences <- lapply(sequences, s2c)
  		nomseq <- lapply(seq_len(nseq), function(i) {
            firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
            substr(firstword, 2, nchar(firstword))
        })		
		names(sequences) <- nomseq
		return(sequences)
		}
}
