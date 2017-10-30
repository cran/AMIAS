AMIAS2d=function(Mat, D_type = "fused.2d", composite = FALSE, W_type = "identity", k1=NULL, k2=1, D=NULL, W=NULL, T1=NULL, T2=min(10,nmw[1]-1), rho1=(dim1*dim2)^6, rho2=(dim1*dim2)^6, h=5, tao=1, outer_itermax=20, select_max=min(20,nm[1]-1), eps=0.1, iter_max=10, ...){

	# Checking
	if (storage.mode(Mat)=="integer") storage.mode(Mat) <- "double"
	if (any(is.na(Mat))) stop("Missing data (NA's) detected. Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing Mat to AMIAS2d")
	dim1 = as.integer(nrow(Mat))
	dim2 = as.integer(ncol(Mat))
	nm = dim(D)
	nmw = dim(W)
    if(!is.null(D)){
		if(!(nm[2]==dim1*dim2))stop("D and Mat do not have the same number of observations")
	}else{
		if(D_type == "user")stop("Choose a 'user' D_type without passing D to AMIAS")
	}
	if(!is.null(W)){
		if(!(nmw[2]==dim1*dim2))stop("D and Mat do not have the same number of observations")
	}else{
		if(D_type == "user")stop("Choose a 'user' D_type without passing D to AMIAS")
	}
	if(D_type != "fused.2d" && D_type != "fused.tfk" && D_type != "user")stop("No defined method")
	if(W_type != "identity"&&W_type != "fused.2d" && W_type != "fused.tfk" && W_type != "user")stop("No defined method")

	if(composite){
		if(!is.null(D))D_type = "user"
		if(!is.null(W))W_type = "user"
		if(D_type != "user"){
			D <- getM2d(D_type,k1,dim1,dim2)
			nm = dim(D)
		}
		if(W_type != "user"){
			W <- getM2d(W_type,k2,dim1,dim2)
			nmw = dim(W)
		}
		if(all(nm==nmw)){
			if(all(D==W))stop("D and W have the same value.")
		}
		DTD = D%*%t(D)
		WTW = W%*%t(W)
		n = as.integer(dim1*dim2)
		if(is.null(T1)){
			result = .Fortran("acomggfusedl0",n,as.vector(Mat),matrix(0,nrow=n,ncol=outer_itermax),rep(0,nm[1]),rep(0,nm[1]),as.integer(tao),as.integer(outer_itermax),as.integer(select_max),eps,rho1,
			rep(as.integer(1),outer_itermax),as.integer(iter_max),DTD,D,as.integer(nm[1]),rep(0,nmw[1]),rep(0,nmw[1]),WTW,W,as.integer(nmw[1]),as.integer(T2),rho2,PACKAGE="AMIAS")
			result[[3]] = result[[3]][,1:result[[7]]]
			result[[11]] = result[[11]][1:result[[7]]]
			result = c(dim1,result)
		}else{
			result = .Fortran("comggfusedl0",n,as.vector(Mat),rep(0,n),rep(0,nm[1]),rep(0,nm[1]),as.integer(T1),rho1,as.integer(1),as.integer(iter_max),
			DTD,D,as.integer(nm[1]),rep(0,nmw[1]),rep(0,nmw[1]),WTW,W,as.integer(nmw[1]),as.integer(T2),rho2,PACKAGE="AMIAS")
			result = c(dim1,result)
		}
	}else{
		if(is.null(T1)){
			if(is.null(D)){
				if(D_type == "fused.2d"){
					vec = c(-1,1)
					veclen = as.integer(2)
					inv = as.matrix(inv2d(dim1,dim2))				
				}else if(D_type == "fused.tfk"){
					if(is.null(k1))stop("Parameter k1 hasn't set to gennerate D (polynomial trend filtering of order k1)")
					k1 = as.integer(k1)
					stopifnot(k1>=1)
					vec = genDtf2d(k=k1)
					veclen = as.integer(length(vec))
					inv = as.matrix(invtf2d(dim1,dim2,k1))
				}
				result = .Fortran("atfusedl02d",dim1,dim2,as.vector(Mat),matrix(0,nrow=dim1*dim2,ncol=outer_itermax),rep(0,nrow(inv)),rep(0,nrow(inv)),as.integer(tao),
				as.integer(outer_itermax),as.integer(select_max),eps,rho1,rep(as.integer(1),outer_itermax),as.integer(iter_max),inv,vec,veclen,PACKAGE="AMIAS")
				result[[4]] = result[[4]][,1:result[[8]]]
				result[[12]] = result[[12]][1:result[[8]]]
			}else{
				D_type = "user"
				DTD <- D%*%t(D)
				if(rcond(DTD)<1e-12){
					diag(DTD) = diag(DTD) + 1e-6
				}
				DTD = solve(DTD);gc()
				n  = as.integer(dim1*dim2)
				result = .Fortran("agfusedl0",n,as.vector(Mat),matrix(0,nrow=n,ncol=outer_itermax),rep(0,nm[1]),rep(0,nm[1]),as.integer(tao),as.integer(outer_itermax),
				as.integer(select_max),eps,rho1,rep(as.integer(1),outer_itermax),as.integer(iter_max),DTD,D,as.integer(nm[1]),PACKAGE="AMIAS")
				result[[3]] = result[[3]][,1:result[[7]]]
				result[[11]] = result[[11]][1:result[[7]]]
				result = c(dim1,result)
			}
		}else{
			if(is.null(D)){
				if(D_type == "fused.2d"){
					vec = c(-1,1)
					veclen = as.integer(2)
					inv = as.matrix(inv2d(dim1,dim2))
					result = .Fortran("tfusedl02d",dim1,dim2,as.vector(Mat),rep(0,dim1*dim2),rep(0,nrow(inv)),rep(0,nrow(inv)),
					as.integer(T1),rho1,as.integer(1),as.integer(iter_max),inv,vec,veclen,PACKAGE="AMIAS")
				}else if(D_type == "fused.tfk"){
					if(is.null(k1))stop("Parameter k1 hasn't set to gennerate D (polynomial trend filtering of order k1)")
					k1 = as.integer(k1)
					stopifnot(k1>=1)
					vec = genDtf2d(k=k1)
					veclen = as.integer(length(vec))
					inv = as.matrix(invtf2d(dim1,dim2,k1))
					result = .Fortran("tfusedl02d",dim1,dim2,as.vector(Mat),rep(0,dim1*dim2),rep(0,nrow(inv)),rep(0,nrow(inv)),
					as.integer(T1),rho1,as.integer(1),as.integer(iter_max),inv,vec,veclen,PACKAGE="AMIAS")
				}
			}else{
				D_type = "user"
				DTD = D%*%t(D)
				if(rcond(DTD)<1e-12){
					diag(DTD) = diag(DTD) + 1e-6
				}
				DTD = solve(DTD);gc()
				n = as.integer(dim1*dim2)
				result = .Fortran("gfusedl0",n,as.vector(Mat),rep(0,n),rep(0,nm[1]),rep(0,nm[1]),as.integer(T1),
				rho1,as.integer(1),as.integer(iter_max),DTD,D,as.integer(nm[1]),PACKAGE="AMIAS")
				result = c(dim1,result)
			}
		}
	}

	if (is.null(T1)){
			output = structure(list(call = match.call(),
				Mat = Mat,
                Mathat = if(length(result[[4]])>dim1*dim2){lapply(1:result[[8]],function(x)matrix(result[[4]][,x],nrow=dim1))}else{matrix(result[[4]],nrow=dim1)},
				composite = composite,
				k1 = k1,
				k2 = k2,
				alpha = result[[5]],
				u = result[[6]],
				gamma = if(!composite){NULL}else{result[[17]]},
				v = if(!composite){NULL}else{result[[18]]},
				df = seq(from=tao,by=tao,length.out=result[[8]]),
				T2 = T2,
                D_type = D_type,
				W_type = W_type,
				rho1 = rho1,
				rho2 = rho2,
				tao = tao,
				eps = eps,
				iter = result[[12]],
				iter_max = iter_max,
				outer_itermax = outer_itermax,
                select_max = select_max),
                class = 'AMIAS.2d')
	}else{
		output = structure(list(call = match.call(),
			Mat = Mat,
            Mathat = matrix(result[[4]],nrow=dim1),
			composite = composite,
			k1 = k1,
			k2 = k2,
			alpha = result[[5]],
			u = result[[6]],
			gamma = if(!composite){NULL}else{result[[17]]},
			v = if(!composite){NULL}else{result[[18]]},
			df = T1,
			T2 = T2,
            D_type = D_type,
			W_type = W_type,
			rho1 = rho1,
			rho2 = rho2,
			tao = tao,
			eps = eps,
			iter = result[[9]],
			iter_max = iter_max,
			outer_itermax = outer_itermax,
            select_max = select_max),
            class = 'AMIAS.2d')
	}
	output
}
