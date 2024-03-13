#' @title dynGENIE3
#' 
#' @description \code{dynGENIE3} Infers a gene regulatory network (in the form of a weighted adjacency matrix) from time series of expression data, using ensembles of regression trees. The method can also be used to infer a network jointly from timeseries and steady-state expression data.
#'
#' @param TS.data List of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment. Each row of a matrix is a gene, each column is a time point. The number of time points does not need to be the same in all the experiments, but the genes must be the same.
#' @param time.points List of vectors, where the k-th vector contains the time points corresponding to the k-th time series experiment.
#' @param alpha Either "from.data" (default), or a vector containing the gene degradation rates, or a single number. When alpha is "from.data", the degradation rate of each gene is estimated from the data, by assuming an exponential decay between the highest and lowest observed expression values. When alpha is a single number, all the genes are assumed to have the same degradation rate alpha. 
#' @param SS.data Steady-state expression matrix (genes x samples). Every row is a gene, every column is a sample. The default value NULL means that no steady-state data are used for network inference.
#' @param tree.method Tree-based method used. Must be either "RF" for Random Forests (default) or "ET" for Extra-Trees.
#' @param K Number of candidate regulators randomly selected at each tree node (for the determination of the best split). Must be either "sqrt" for the square root of the total number of candidate regulators (default), "all" for the total number of candidate regulators, or a stricly positive integer.
#' @param ntrees Number of trees in an ensemble for each target gene. Default: 1000.
#' @param regulators Subset of genes used as candidate regulators. Must be either a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}. The default value NULL means that all the genes are used as candidate regulators.
#' @param ncores Number of cores to use for parallel computing. Default: 1.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: FALSE.
#' @param seed Random number generator seed for replication of analyses. The default value NULL means that the seed is not reset.
#'
#' @return 2 objects: (1) weight.matrix: the weighted adjacency matrix of the inferred network. Element w_ij (row i, column j) gives the importance of the link from regulatory gene i to target gene j. (2) alphas: the gene degradation rates
#' 
#' @examples
#' ## Generate fake data
#'
#' ## Time series data
#' data1 <- matrix(sample(1:10, 70, replace=TRUE), nrow=10)
#' rownames(data1) <- paste("Gene", 1:10, sep="")
#' colnames(data1) <- paste("Time point", 1:7, sep="")
#'
#' data2 <- matrix(sample(1:10, 50, replace=TRUE), nrow=10)
#' rownames(data2) <- paste("Gene", 1:10, sep="")
#' colnames(data2) <- paste("Time point", 1:5, sep="")
#'
#' TS.data <- list(data1,data2) 
#' 
#' ## Time points
#' tp1 <- seq(from=0, to=60, by=10)
#' tp2 <- seq(from=0, to=40, by=10) 
#' time.points <- list(tp1,tp2)
#'
#' ## Run dynGENIE3
#' res <- dynGENIE3(TS.data,time.points)
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(res$weight.matrix)
#' head(link.list)
#' @export
dynGENIE3 <- function(TS.data, time.points, alpha="from.data", SS.data=NULL, tree.method="RF", K="sqrt", ntrees=1000, regulators=NULL, ncores=1, verbose=FALSE, seed=NULL) {

	dyn.load("dynGENIE3.so")

	# check input arguments
	
	if (!is.list(TS.data) && !is.vector(TS.data)){
		stop("Parameter TS.data must be a list of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment.")
	}
	
	for (k in seq(from=1, to=length(TS.data))) {
		if (!is.matrix(TS.data[[k]]) && !is.array(TS.data[[k]])) {
			stop("Parameter TS.data must be a list of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment.")
		}
	
		if (length(dim(TS.data[[k]])) != 2) {
			stop("Parameter TS.data must be a list of expression matrices, where each matrix (genes x time points) corresponds to a time series experiment.")
		}
	
		if (is.null(rownames(TS.data[[k]]))) {
			stop("Each k-th matrix of TS.data must specify the names of the genes in rownames(TS.data[[k]]).")
		}
	}
	
	gene.names <- rownames(TS.data[[1]])
	num.genes <- length(gene.names)
	
	if (length(TS.data) > 1) {
		for (k in seq(from=2, to=length(TS.data))){
			if (length(union(gene.names, rownames(TS.data[[k]]))) != num.genes) {
				stop("The genes/rows must be the same in all the matrices of TS.data.")
			}
		}
	}
	
	
	if (!is.list(time.points) && !is.vector(time.points)){
		stop("Parameter time.points must be a list of vectors, where the k-th vector contains the time points related to the k-th time series experiment.")
	}
	
	if (length(time.points) != length(TS.data)) {
		stop("The number of vectors in time.points must be equal to the number of time series experiments in TS.data.")
	}
	
	for (k in seq(from=1, to=length(time.points))) {
		if (!is.vector(time.points[[k]])) {
			stop("Parameter time.points must be a list of vectors, where the k-th vector contains the time points related to the k-th time series experiment.")
		}
	}
	
	if (!is.numeric(alpha) && alpha != "from.data") {
		stop("Parameter alpha must be either 'from.data', a positive number or a vector of positive numbers.")	
	}
	
	if (is.numeric(alpha) && !is.vector(alpha)) {
		stop("Parameter alpha must be either 'from.data', a positive number or a vector of positive numbers.")
	}
	
	if (is.numeric(alpha)) {
		if (length(alpha) > 1) {
			if (length(alpha) != num.genes) {
				stop("When parameter alpha is a vector, this must be a vector of length p, where p is the number of genes.")
			}
			if (is.null(names(alpha))) {
				stop("When parameter alpha is a vector, the gene names must be specified.")
			}
		}
	
		for (i in seq(from=1, to=length(alpha))) {
			if (alpha[[i]] < 0) {
				stop("The gene degradation rates specified in parameter alpha must be positive.")
			}
		}	
	}
	
	if (!is.null(SS.data)) {
		
		if (!is.matrix(SS.data) && !is.array(SS.data)) {
			stop("Parameter SS.data must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
		}
	
		if (length(dim(SS.data)) != 2) {
			stop("Parameter SS.data must be a two-dimensional matrix where each row corresponds to a gene and each column corresponds to a condition/sample.")
		}
	
		if (is.null(rownames(SS.data))) {
			stop("SS.data must specify the names of the genes in rownames(SS.data).")
		}
		
		if (length(intersect(gene.names, rownames(SS.data))) != num.genes) {
			stop("The genes/rows of SS.data must be the same as the genes in TS.data.")
		}
		
	}
	
	if (tree.method != "RF" && tree.method != "ET") {
		stop("Parameter tree.method must be \"RF\" (Random Forests) or \"ET\" (Extra-Trees).")
	}
	
	if (K != "sqrt" && K != "all" && !is.numeric(K)) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}
	
	if (is.numeric(K) && K<1) {
		stop("Parameter K must be \"sqrt\", or \"all\", or a strictly positive integer.")
	}
	
	if (!is.numeric(ntrees) || ntrees<1) {
		stop("Parameter ntrees should be a stricly positive integer.")
	}
	
	if (!is.null(regulators)) {
		if (!is.vector(regulators)) {
			stop("Parameter regulators must be either a vector of indices or a vector of gene names.")
		}
		
		if (is.character(regulators) && length(intersect(regulators,gene.names)) == 0) {
			stop("The genes must contain at least one candidate regulator.")
		}
		
		if (is.numeric(regulators) && max(regulators) > num.genes) {
			stop("At least one index in regulators exceeds the number of genes.")
		}
	}
	
	if (!is.numeric(ncores) || ncores<1) {
		stop("Parameter ncores should be a stricly positive integer.")
	}
	
	
	
	# set random number generator seed if seed is given
    if (!is.null(seed)) {
        set.seed(seed)
    }
	
	num.exp <- length(TS.data)
    
    # transpose all expression matrices to (samples x genes)
    for (k in seq(from=1, to=num.exp)) {
    	TS.data[[k]] <- t(TS.data[[k]])
    }
	if (!is.null(SS.data)) {
		SS.data <- t(SS.data)
	}
	
	# re-order time points in increasing order
	for (k in seq(from=1, to=num.exp)) {
		tp <- time.points[[k]]
		tp_order <- order(tp)
		time.points[[k]] <- tp[tp_order]
		
		expr.data <- TS.data[[k]]
		TS.data[[k]] <- expr.data[tp_order,]
	}
	
	# degradation rates
	if (!is.numeric(alpha)) {
		alphas <- estimate.decay.rates(TS.data, time.points)
	} else if (length(alpha) == 1) {
		alphas <- rep(alpha,num.genes)
		alphas <- setNames(alphas,gene.names)
	} else {
		alphas <- alpha[gene.names]
	}
	
    # setup weight matrix
    weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
    rownames(weight.matrix) <- gene.names
    colnames(weight.matrix) <- gene.names
	
    # get names of input genes
    if (is.null(regulators)) {
        input.gene.names <- gene.names
    } else {
        # input gene indices given as integers
        if (is.numeric(regulators)) {
            input.gene.names <- gene.names[regulators]
        # input gene indices given as names
        } else {
            input.gene.names <- regulators
            # for security, abort if some input gene name is not in gene names
            missing.gene.names <- setdiff(input.gene.names, gene.names)
            if (length(missing.gene.names) != 0) {
                for (missing.gene.name in missing.gene.names) {
                    cat(paste("Gene ", missing.gene.name,
                              " was not in the expression matrix\n", sep=""))
                }
                stop("Aborting computation")
            }
        }
    }
	
	num.input.genes <- length(input.gene.names)
	
	# tree method
	if (tree.method == 'RF') {
		RF_randomisation <- 1
		ET_randomisation <- 0
		bootstrap_sampling <- 1
	} else {
		RF_randomisation <- 0
		ET_randomisation <- 1
		bootstrap_sampling <- 0
	} 
	
	if (verbose) {
        cat(paste("Tree method: ", tree.method, "\nK: ", K,
	              "\nNumber of trees: ", ntrees, "\nalpha min: ", min(alphas),
				  "\nalpha max: ", max(alphas), "\n\n",
                  sep=""))
        flush.console()
	}
	
	
	# generate time-lagged inputs and outputs
	h <- 1
	num.samples.time <- 0
	for (k in seq(from=1, to=num.exp)) {
		num.samples.time <- num.samples.time + dim(TS.data[[k]])[1]
	}
	
	input.matrix.time <- matrix(0.0, nrow=num.samples.time-h*num.exp, ncol=num.input.genes)
	output.matrix.time <- matrix(0.0, nrow=num.samples.time-h*num.exp, ncol=num.genes)
	
	num.samples.count <- 0
	
	for (k in seq(from=1, to=num.exp)){
		current.timeseries <- TS.data[[k]][,gene.names]
		current.time.points <- time.points[[k]]
		num.points <- length(current.time.points)
		time.diff.current <- current.time.points[(h+1):num.points] - current.time.points[1:(num.points-h)]
		current.timeseries.input <- current.timeseries[1:(num.points-h),input.gene.names]
		current.timeseries.output <- (current.timeseries[(h+1):num.points,] - current.timeseries[1:(num.points-h),]) / replicate(num.genes,time.diff.current) + t(replicate((num.points-h),alphas)) * current.timeseries[1:(num.points-h),]
		num.samples.current <- dim(current.timeseries.input)[1]
		input.matrix.time[(num.samples.count+1):(num.samples.count+num.samples.current),] <- current.timeseries.input
		output.matrix.time[(num.samples.count+1):(num.samples.count+num.samples.current),] <- current.timeseries.output
		num.samples.count <- num.samples.count + num.samples.current
	}
	
	# add steady-state data (if any)
	if (!is.null(SS.data)) {
		num.SS.exp <- dim(SS.data)[1]
		input.all <- rbind(SS.data[,input.gene.names],input.matrix.time)
		output.all <- rbind(SS.data[,gene.names]*t(replicate(num.SS.exp,alphas)),output.matrix.time)
	} else {
		input.all <- input.matrix.time
		output.all <- output.matrix.time
	}
	
	colnames(output.all) <- gene.names
	
	num.samples = dim(input.all)[1]
    
    # compute importances for every target gene
    
    x <- input.all
   
	if (ncores==1) {
		# serial computing
		if (verbose) {
		    cat("Using 1 core.\n\n")
		    flush.console()
		}
		
	    for (target.gene.idx in seq(from=1, to=num.genes)) {

            if (verbose) {	
                cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
                flush.console()
			 }

	        target.gene.name <- gene.names[target.gene.idx]
			y <- output.all[,target.gene.name]

		    # set mtry
		    if (class(K) == "numeric") {
		        mtry <- K
		    } else if (K == "sqrt") {
		        mtry <- round(sqrt(num.input.genes))
		    } else {
		        mtry <- num.input.genes
		    } 
		
			# some default parameters 
			nmin <- 1
			permutation_importance <- 0
		
	        im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
			          as.single(c(x)),as.single(c(y)),as.integer(nmin),
					  as.integer(ET_randomisation),as.integer(RF_randomisation),
					  as.integer(mtry),as.integer(ntrees),
					  as.integer(bootstrap_sampling),as.integer(permutation_importance),
					  as.double(vector("double",num.input.genes)))[[12]]
					  
		  	# some variable importances might be slighly negative due to some rounding error
		  	im[im<0] <- 0
					  
			im <- setNames(im,input.gene.names)
			im[target.gene.name] <- 0
			
			# normalize variable importances
			im <- im / sum(im)
			
	        weight.matrix[input.gene.names, target.gene.name] <- im[input.gene.names]
	    }
	} else {
		# parallel computing
	    library(doRNG); library(doParallel); registerDoParallel(); options(cores=ncores)
		
		if (verbose) {
		    message(paste("\nUsing", getDoParWorkers(), "cores."))
		}
		
	    weight.matrix.reg <- foreach(target.gene.name=gene.names, .combine=cbind) %dorng% 
	    {		
			y <- output.all[,target.gene.name]

		    # set mtry
		    if (class(K) == "numeric") {
		        mtry <- K
		    } else if (K == "sqrt") {
		        mtry <- round(sqrt(num.input.genes))
		    } else {
		        mtry <- num.input.genes
		    } 
			
			# some default parameters 
			nmin <- 1
			permutation_importance <- 0
		
	        im <- .C("BuildTreeEns",as.integer(num.samples),as.integer(num.input.genes),
			          as.single(c(x)),as.single(c(y)),as.integer(nmin),
					  as.integer(ET_randomisation),as.integer(RF_randomisation),
					  as.integer(mtry),as.integer(ntrees),
					  as.integer(bootstrap_sampling),as.integer(permutation_importance),
					  as.double(vector("double",num.input.genes)))[[12]]
					  
		  	# some variable importances might be slighly negative due to some rounding error
		  	im[im<0] <- 0
					  
		  	im <- setNames(im,input.gene.names)
		  	im[target.gene.name] <- 0
			
		  	# normalize variable importances
		  	im <- im / sum(im)
					  
			im[input.gene.names]
	    }
	    attr(weight.matrix.reg, "rng") <- NULL
	    weight.matrix[input.gene.names,] <- weight.matrix.reg
	}
	
	res.to.return <- list("weight.matrix" = weight.matrix, "alphas" = alphas)
	
    return(res.to.return)
}       
	



#' @title get.link.list
#' 
#' @description \code{get.link.list} Converts the weight matrix, as returned by \code{\link{dynGENIE3}}, to a sorted list of regulatory links (most likely links first).
#' 
#' @param weight.matrix Weighted adjacency matrix as returned by \code{\link{dynGENIE3}}.
#' @param report.max Maximum number of links to report. The default value NULL means that all the links are reported.
#' @param threshold Only links with a weight above the threshold are reported. Default: threshold = 0, i.e. all the links are reported.
#' 
#' @return List of regulatory links in a data frame. Each line of the data frame corresponds to a link. The first column is the regulatory gene, the second column is the target gene, and the third column is the weight of the link.
#'
#' @seealso \code{\link{dynGENIE3}}
#'
#' @examples
#' ## Generate fake data
#'
#' ## Time series data
#' data1 <- matrix(sample(1:10, 70, replace=TRUE), nrow=10)
#' rownames(data1) <- paste("Gene", 1:10, sep="")
#' colnames(data1) <- paste("Time point", 1:7, sep="")
#'
#' data2 <- matrix(sample(1:10, 50, replace=TRUE), nrow=10)
#' rownames(data2) <- paste("Gene", 1:10, sep="")
#' colnames(data2) <- paste("Time point", 1:5, sep="")
#'
#' TS.data <- list(data1,data2) 
#' 
#' ## Time points
#' tp1 <- seq(from=0, to=60, by=10)
#' tp2 <- seq(from=0, to=40, by=10) 
#' time.points <- list(tp1,tp2)
#'
#' ## Gene degradation rates
#' alpha <- rep(0.1,10) 
#' names(alpha) <- paste("Gene", 1:10, sep="")
#'
#' ## Run dynGENIE3
#' weight.matrix <- dynGENIE3(TS.data,time.points,alpha)
#' 
#' ## Get ranking of edges 
#' link.list <- get.link.list(weight.matrix)
#' head(link.list)
#' @export
get.link.list <- function(weight.matrix, report.max=NULL, threshold=0) {
    if(!is.numeric(threshold)) {
    	stop("threshold must be a number.")
    } 
	
	library(reshape2)
	
	# Only process weights that are not self-regulations
	input.genes <- rownames(weight.matrix)
	diag(weight.matrix[input.genes, input.genes]) <- NA
    link.list <- melt(weight.matrix, na.rm=TRUE)
    colnames(link.list) <- c("regulatory.gene", "target.gene", "weight")
    link.list <- link.list[link.list$weight>=threshold,]
    link.list <- link.list[order(link.list$weight, decreasing=TRUE),]
  
    if(!is.null(report.max)) {
    	link.list <- link.list[1:min(nrow(link.list), report.max),]
    } 
  
    rownames(link.list) <- NULL
  
    return(link.list)
}




read.expr.matrix <- function(filename, form="", sep="", default.gene.label="gene_", default.sample.label="sample_") {
    # Report when form is not correctly set
    if (form != "rows.are.genes" && form != "rows.are.samples") {
        stop("Parameter form must be set to \"rows.are.genes\" or \"rows.are.samples\"")
    }
    # read data
    m <- read.table(filename, sep=sep, as.is=TRUE)
    has.row.names <- FALSE
    has.col.names <- FALSE
    # have row and column names been recognized (cell (1,1) is empty in file) ?
    if (colnames(m)[1] != "V1" && rownames(m)[1] != "1") {
        has.row.names <- TRUE
        has.col.names <- TRUE
    }
    # is first column alphanumeric ?
    if (all(grepl("[[:alpha:]]", m[,1]))) {
        # check duplicated names
        if (any(duplicated(m[,1]))) {
            stop("Duplicated names in first column\n")
        }
        rownames(m) <- m[,1]
        m <- m[,-1]
        has.row.names <- TRUE
    }
    # is first row alphanumeric ?
    if (all(grepl("[[:alpha:]]", m[1,]))) {
        # check duplicated names
        if (any(duplicated(m[1,]))) {
            stop("Duplicated names in first row\n")
        }
        colnames(m) <- m[1,]
        m <- m[-1,]
        has.col.names <- TRUE
    }
    # coerce matrix data to numeric
    col.names <- colnames(m)
    row.names <- rownames(m)
    m <- as.matrix(m)
    m <- apply(m, 2, function(x) { as.numeric(x) })
    colnames(m) <- col.names
    rownames(m) <- row.names
    num.rows <- dim(m)[1]
    num.cols <- dim(m)[2]
    # fill in default gene names in rows if needed
    if (!has.row.names && form=="rows.are.genes") {
        rownames(m) <- paste(default.gene.label, seq(from=1, to=num.rows), sep="")
    }
    # fill in default sample names in rows if needed
    if (!has.row.names && form=="rows.are.samples") {
        rownames(m) <- paste(default.sample.label, seq(from=1, to=num.rows), sep="")
    }
    # fill in default sample names in columns if needed
    if (!has.col.names && form=="rows.are.genes") {
        colnames(m) <- paste(default.sample.label, seq(from=1, to=num.cols), sep="")
    }
    # fill in default gene names in columns if needed
    if (!has.col.names && form=="rows.are.samples") {
        colnames(m) <- paste(default.gene.label, seq(from=1, to=num.cols), sep="")
    }
    # transpose matrix to (genes x samples) if needed
    if (form == "rows.are.samples") m <- t(m)
    return(m)
}


estimate.decay.rates <- function(TS.data, time.points) {
	# For each gene, the decay rate is estimated by assuming that the gene expression x(t) follows:
    # x(t) =  A exp(-alpha * t) + C_min,
    # between the highest and lowest expression values.
    # C_min is set to the minimum expression value over all genes and all samples.
	
	num.exp <- length(TS.data)
	gene.names <- colnames(TS.data[[1]])
	num.genes <- length(gene.names)
	
	C_min <- min(TS.data[[1]])
	if (num.exp > 1) {
		for (k in seq(from=2, to=num.exp)) {
			C_min <- min(C_min,min(TS.data[[k]]))
		}
	}
	
	alphas <- matrix(0.0, nrow=num.exp, ncol=num.genes)
    colnames(alphas) <- gene.names
	
	for (k in seq(from=1, to=num.exp)) {
		current.timeseries <- TS.data[[k]]
		current.time.points <- time.points[[k]]
		
		for (target.gene.idx in seq(from=1, to=num.genes)) {
			target.gene.name <- gene.names[target.gene.idx]
			
			idx.min <- which.min(current.timeseries[,target.gene.name])
			idx.max <- which.max(current.timeseries[,target.gene.name])
			
			xmin <- current.timeseries[idx.min,target.gene.name]
			xmax <- current.timeseries[idx.max,target.gene.name]

			if (xmin != xmax) {
			
                tmin <- current.time.points[idx.min]
                tmax <- current.time.points[idx.max]

                xmin <- max(xmin-C_min,1e-6)
                xmax <- max(xmax-C_min,1e-6)

                xmin <- log(xmin)
                xmax <- log(xmax)

                alphas[k,target.gene.name] = (xmax - xmin) / abs(tmin - tmax)
                }
		}
	}	
	
	alphas <- apply(alphas,2,max)
	
	return(alphas)
}