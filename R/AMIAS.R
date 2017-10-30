AMIAS=function(y, D_type = "fused.1d", composite = FALSE, W_type = "identity", k1=NULL, k2=1, D=NULL, W=NULL, T1=NULL, T2=min(10,nmw[1]-1), rho1=n^2, rho2=n^2, h=5, tao=1, outer_itermax=20, select_max=min(20,nm[1]-1), eps=0.1, iter_max=10, smooth = TRUE, ...){

    # Checking
    if (storage.mode(y)=="integer") storage.mode(y) <- "double"
	if (any(is.na(y))) stop("Missing data (NA's) detected. Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing y to AMIAS")
	n = as.integer(length(y))
	nm = dim(D)
    if(!is.null(D)){
		if(!(nm[2]==n))stop("D and y do not have the same number of observations")
		if (storage.mode(D)=="integer") storage.mode(D) <- "double"
	}else{
		if(D_type == "user")stop("Choose a 'user' D_type without passing D to AMIAS")
	}
	nmw = dim(W)
	if(!is.null(W)){
		if(!(nmw[2]==n))stop("W and y do not have the same number of observations")
		if (storage.mode(W)=="integer") storage.mode(W) <- "double"
	}else{
		if(W_type == "user")stop("Choose a 'user' D_type without passing D to AMIAS")
	}
	if(D_type != "fused.1d" && D_type != "fused.tfk" && D_type != "user")stop("No defined method")
	if(W_type != "identity"&&W_type != "fused.1d" && W_type != "fused.tfk" && W_type != "user")stop("No defined method")

	# Smooth
	if(smooth)y <- my.rollmean(y, h = h)
	
	if(composite){
		if(!is.null(D))D_type = "user"
		if(!is.null(W))W_type = "user"
		if(is.null(T1)){
			if(W_type == "identity" && D_type == "fused.1d"){
				result = .Fortran("acomdifusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),rep(0,n-1),rep(0,n-1),as.integer(tao),as.integer(outer_itermax),
				as.integer(select_max),eps,rho1,rep(as.integer(1),outer_itermax),as.integer(iter_max),rep(0,n),rep(0,n),as.integer(T2),rho2,PACKAGE="AMIAS")
				result[[3]] = result[[3]][,1:result[[7]]]
				result[[11]] = result[[11]][1:result[[7]]]
			}else{
				if(D_type != "user"){
					D <- getM(D_type,k1,n)
					nm = dim(D)
				}
				if(W_type != "user"){
					W <- getM(W_type,k2,n)
					nmw = dim(W)
				}
				if(all(nm==nmw)){
					if(all(D==W))stop("D and W have the same value.")
				}
				DTD = D%*%t(D)
				WTW = W%*%t(W)
				result = .Fortran("acomggfusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),rep(0,nm[1]),rep(0,nm[1]),as.integer(tao),as.integer(outer_itermax),as.integer(select_max),eps,rho1,
				rep(as.integer(1),outer_itermax),as.integer(iter_max),rep(0,nmw[1]),rep(0,nmw[1]),DTD,D,as.integer(nm[1]),WTW,W,as.integer(nmw[1]),as.integer(T2),rho2,PACKAGE="AMIAS")
				result[[3]] = result[[3]][,1:result[[7]]]
				result[[11]] = result[[11]][1:result[[7]]]
			}
		}else{
			if(W_type == "identity" && D_type == "fused.1d"){
				result = .Fortran("comdifusedl0",n,y,rep(0,n),rep(0,n-1),rep(0,n-1),as.integer(T1),rho1,as.integer(1),
				as.integer(iter_max),rep(0,n),rep(0,n),as.integer(n),as.integer(T2),rho2,PACKAGE="AMIAS")
			}else{
				if(D_type != "user"){
					D <- getM(D_type,k1,n)
					nm = dim(D)
				}
				if(W_type != "user"){
					W <- getM(W_type,k2,n)
					nmw = dim(W)
				}
				if(all(nm==nmw)){
					if(all(D==W))stop("D and W have the same value.")
				}
				DTD = D%*%t(D)
				WTW = W%*%t(W)
				result = .Fortran("comggfusedl0",n,y,rep(0,n),rep(0,nm[1]),rep(0,nm[1]),as.integer(T1),rho1,as.integer(1),as.integer(iter_max),
				DTD,rep(0,nmw[1]),rep(0,nmw[1]),D,as.integer(nm[1]),WTW,W,as.integer(nmw[1]),as.integer(T2),rho2,PACKAGE="AMIAS")			
			}
		}
	}else{
		if(is.null(T1)){
			if(is.null(D)){
				if(D_type == "fused.1d"){
					result = .Fortran("afusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),y[2:n]-y[1:(n-1)],rep(0,n-1),as.integer(tao),
					as.integer(outer_itermax),as.integer(select_max),eps,rho1,rep(as.integer(1),outer_itermax),as.integer(iter_max),PACKAGE="AMIAS")
					result[[3]] = result[[3]][,1:result[[7]]]
					result[[11]] = result[[11]][1:result[[7]]]
				}else if(D_type == "fused.tfk"){
					if(is.null(k1))stop("Parameter k1 hasn't set to gennerate D (polynomial trend filtering of order k1)")
					k1 = as.integer(k1)
					stopifnot(k1>=1)
					vec = genDtf1d(k=k1)
					veclen = as.integer(length(vec))
					DTD = .Fortran("dtdmul",vec,n,veclen,matrix(0,nrow=2*veclen-1,ncol=n-veclen+1))[[4]]
					result = .Fortran("batfusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),rep(0,n-veclen+1),rep(0,n-veclen+1),as.integer(tao),as.integer(outer_itermax),
					as.integer(select_max),eps,rho1,rep(as.integer(1),outer_itermax),as.integer(iter_max),DTD,vec,veclen,PACKAGE="AMIAS")
					result[[3]] = result[[3]][,1:result[[7]]]
					result[[11]] = result[[11]][1:result[[7]]]
				}
			}else{
				D_type = "user"
				DTD <- D%*%t(D)
				result = .Fortran("agfusedl0",n,y,matrix(0,nrow=n,ncol=outer_itermax),rep(0,nm[1]),rep(0,nm[1]),as.integer(tao),as.integer(outer_itermax),
				as.integer(select_max),eps,rho1,rep(as.integer(1),outer_itermax),as.integer(iter_max),DTD,D,as.integer(nm[1]),PACKAGE="AMIAS")
				result[[3]] = result[[3]][,1:result[[7]]]
				result[[11]] = result[[11]][1:result[[7]]]
			}
		}else{
			if(is.null(D)){
				if(D_type == "fused.1d"){
					result = .Fortran("fusedl0",n,y,rep(0,n),y[2:n]-y[1:(n-1)],
					rep(0,n-1),as.integer(T1),rho1,as.integer(1),as.integer(iter_max),PACKAGE="AMIAS")
				}else if(D_type == "fused.tfk"){
					if(is.null(k1))stop("Parameter k1 hasn't set to gennerate D (polynomial trend filtering of order k1)")
					k1 = as.integer(k1)
					stopifnot(k1>=1)
					vec = genDtf1d(k=k1)
					veclen = as.integer(length(vec))
					DTD = .Fortran("dtdmul",vec,n,veclen,matrix(0,nrow=2*veclen-1,ncol=n-veclen+1),PACKAGE="AMIAS")[[4]]
					result = .Fortran("btfusedl0",n,y,rep(0,n),rep(0,n-veclen+1),rep(0,n-veclen+1),
					as.integer(T1),rho1,as.integer(1),as.integer(iter_max),DTD,vec,veclen,PACKAGE="AMIAS")
				}
			}else{
				D_type = "user"
				DTD = D%*%t(D)
				result = .Fortran("gfusedl0",n,y,rep(0,n),rep(0,nm[1]),rep(0,nm[1]),as.integer(T1),
				rho1,as.integer(1),as.integer(iter_max),DTD,D,as.integer(nm[1]),PACKAGE="AMIAS")
			}
		}
	}

	if (is.null(T1)){
			output = structure(list(call = match.call(),
				y = y,
                beta = result[[3]],
				composite = composite,
				k1 = k1,
				k2 = k2,
				alpha = result[[4]],
				u = result[[5]],
				gamma = if(!composite){NULL}else{result[[13]]},
				v = if(!composite){NULL}else{result[[14]]},
				df = seq(from=tao,by=tao,length.out=result[[7]]),
				T2 = T2,
                D_type = D_type,
				W_type = W_type,
				rho1 = rho1,
				rho2 = rho2,
				tao = tao,
				eps = eps,
				iter = result[[11]],
				iter_max = iter_max,
				smooth = smooth,
				outer_itermax = outer_itermax,
                select_max = select_max),
                class = 'AMIAS.1d')
	}else{
		output = structure(list(call = match.call(),
			y = y,
            beta = result[[3]],
			composite = composite,
			k1 = k1,
			k2 = k2,
			alpha = result[[4]],
			u = result[[5]],
			gamma = if(!composite){NULL}else{result[[11]]},
			v = if(!composite){NULL}else{result[[12]]},
			df = T1,
			T2 = T2,
            D_type = D_type,
			W_type = W_type,
			rho1 = rho1,
			rho2 = rho2,
			tao = tao,
			eps = eps,
			iter = result[[8]],
			iter_max = iter_max,
			smooth = smooth,
			outer_itermax = outer_itermax,
            select_max = select_max),
            class = 'AMIAS.1d')
	}
	output
}

print.AMIAS.1d = function(x, ...) {
	result = x
	cat("\nCall:\n", paste(deparse(result$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	selected = result$alpha
	if (sum(selected!=0)) {
		beta.selected <- result$alpha[result$alpha!=0]
		cat("Selected", length(beta.selected), "change points with nozero alpha(gap):\n")
		print.default(beta.selected, print.gap = 2L, quote = FALSE)
	} else {
		cat("No selected variables\n")
	}
	cat("\n")
	invisible(result)
}

plot.AMIAS.1d = function(x, s = length(x$df), betatype = "l", betacolor = "blue", betalwd = 2, ...) {
	s <- as.integer(s)
	stopifnot(s>=1&&s<=length(x$df))
	yname = "original-y"
    if(x$smooth)yname <- "smooth-y" 
	plot.default(x$y,ylab = yname,...)
	if(length(x$df) == 1){
		lines.default(x$beta, type = betatype, col = betacolor)
	}else{
		lines.default(x$beta[,s], type = betatype, col = betacolor, lwd = betalwd)
	}
}
