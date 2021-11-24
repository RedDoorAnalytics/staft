version 15.1

local RC real colvector
local RS real scalar

mata:

//calculate splines with provided knots
real matrix rcsgen_core(	real colvector variable,	///
                                real rowvector knots, 		///
                                real scalar deriv,|			///
                                real matrix rmatrix			///
                        )
{
	`RS' Nobs, Nknots, kmin, kmax, interior, Nparams
	real matrix splines, knots2

	//======================================================================================================================================//
	// Extract knot locations

		Nobs 	= rows(variable)
		Nknots 	= cols(knots)
		kmin 	= knots[1,1]
		kmax 	= knots[1,Nknots]
	
		if (Nknots==2) interior = 0
		else interior = Nknots - 2
		Nparams = interior + 1
		
		splines = J(Nobs,Nparams,.)

	//======================================================================================================================================//
	// Calculate splines

		if (Nparams>1) {
			lambda = J(Nobs,1,(kmax:-knots[,2..Nparams]):/(kmax:-kmin))
			knots2 = J(Nobs,1,knots[,2..Nparams])
		}

		if (deriv==0) {
			splines[,1] = variable
			if (Nparams>1) {
				splines[,2..Nparams] = (variable:-knots2):^3 :* (variable:>knots2) :- lambda:*((variable:-kmin):^3):*(variable:>kmin) :- (1:-lambda):*((variable:-kmax):^3):*(variable:>kmax) 
			}
		}
		else if (deriv==1) {
			splines[,1] = J(Nobs,1,1)
			if (Nparams>1) {
				splines[,2..Nparams] = 3:*(variable:-knots2):^2 :* (variable:>knots2) :- lambda:*(3:*(variable:-kmin):^2):*(variable:>kmin) :- (1:-lambda):*(3:*(variable:-kmax):^2):*(variable:>kmax) 	
			}
		}
		else if (deriv==2) {
			splines[,1] = J(Nobs,1,0)
			if (Nparams>1) {
				splines[,2..Nparams] = 6:*(variable:-knots2) :* (variable:>knots2) :- lambda:*(6:*(variable:-kmin)):*(variable:>kmin) :- (1:-lambda):*(6:*(variable:-kmax)):*(variable:>kmax) 	
			}
		}
		else if (deriv==3) {
			splines[,1] = J(Nobs,1,0)
			if (Nparams>1) {
				splines[,2..Nparams] = 6:*(variable:>knots2) :- lambda:*6:*(variable:>kmin) :- (1:-lambda):*6:*(variable:>kmax)
			}
		}
	
		//orthog
		if (args()==4) {
			real matrix rmat
			rmat = luinv(rmatrix)
			if (deriv==0) splines = (splines,J(Nobs,1,1)) * rmat[,1..Nparams]
			else splines = splines * rmat[1..Nparams,1..Nparams]
		}
		
		return(splines)
}

real rowvector GetNewKnots(	real colvector variable,	///
							real colvector d,			///
							real scalar df, 			///
							real colvector knotindex)
{
	`RC' knotsvar
	`RS' kmin, kmax, Nknots

	knotsvar = select(variable,d)
	kmin	= min(knotsvar)
	kmax	= max(knotsvar)		
	Nknots = df+1

	if (Nknots==2) return(kmin,kmax)
	else return(kmin,(sort(knotsvar,1)[knotindex])',kmax)
}

real matrix GetNewRmatrix(	real colvector variable,	///
							real rowvector knots)
{
	//======================================================================================================================================//
	// Extract knot locations

		Nobs 	= rows(variable)
		Nknots 	= cols(knots)
		kmin 	= knots[1,1]
		kmax 	= knots[1,Nknots]
	
		if (Nknots==2) interior = 0
		else interior = Nknots - 2
		Nparams = interior + 1
		
		x = J(Nobs,Nparams,.)

	//======================================================================================================================================//
	// Calculate splines

		if (Nparams>1) {
			lambda = J(Nobs,1,(kmax:-knots[,2..Nparams]):/(kmax:-kmin))
			knots2 = J(Nobs,1,knots[,2..Nparams])
		}

		x[,1] = variable
		if (Nparams>1) {
			x[,2..Nparams] = (variable:-knots2):^3 :* (variable:>knots2) :- lambda:*((variable:-kmin):^3):*(variable:>kmin) :- (1:-lambda):*((variable:-kmax):^3):*(variable:>kmax) 
		}
		
        meanx = mean(x)
        v = x :- meanx ,J(rows(x),1,1) 
        q = J(rows(v),0,.)
        R = J(cols(v),cols(v),0)
        R[cols(v),] = (meanx,1)
        for (i=1;i<=cols(x);i++){
			r = norm(v[,i])/sqrt(rows(v))
			q = q, (v[,i]:/ r)
			R[i,i] = r
			for (j = i + 1; j<=cols(x); j++){
				r = (q[,i]' * v[,j])/rows(v)
				v[,j] = v[,j] - r*q[,i]
				R[i,j] = r 
			}
        }
        return(R)		
}

end
