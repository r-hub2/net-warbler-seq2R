read.genbank <-
function (locus) {			
		X <- character()
   		URL <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",paste(locus, collapse = ","), "&rettype=gb&retmode=text",sep = "")
       X <- c(X, scan(file = URL, what = "", sep = "\n", quiet = TRUE))
   		prin <- grep("^ {0,}ORIGIN", X) + 1 
   		fin<- which(X == "//") - 1
    	res <- list()
     	code <- gsub("[[:digit:] ]", "", X[prin:fin])
      	res[[1]] <- unlist(strsplit(code, NULL))
   		res[[2]] <- locus
       code <- character()
       sp <- grep("ORGANISM", X)
     	code <- unlist(strsplit(X[sp], " +ORGANISM +"))[2]
       attr(res, "species") <- gsub(" ", "_", code)
  		return(res)
}
