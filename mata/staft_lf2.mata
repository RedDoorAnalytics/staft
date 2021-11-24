*! version 1.0.0 MJC

version 15.1
mata:

void staft_lf2(	transmorphic scalar M, 		///
				real scalar todo, 			///
				real rowvector b, 			///
				real colvector lnfj, 		///
				real matrix S, 				///
				real matrix H)
{
	pointer(struct staft_main scalar) scalar pS
	struct staft_main scalar sS
	real colvector xb, rcscoeffs
	real matrix rcs, drcs, d2rcs, d3rcs
	real scalar ind, int0, i

	pS = &moptimize_util_userinfo(M,1)
	sS = *pS

	//linear predictors
	ind = 1
	xb = -moptimize_util_xb(M,b,ind++)				//negative imposed here
	if (sS.hastvcs) {	
		real colvector dxb
		dxb = -moptimize_util_xb(M,b,ind++)			//negative imposed here
		if (sS.delentry) {
			real colvector s0xb
			s0xb = -moptimize_util_xb(M,b,ind++)	//only need this when tvc's present as need tvc splines at _t0, o/w can use xb
		}
	}
	
	//spline terms
	rcscoeffs = J(sS.Nsplines,1,.)
	for (i=1;i<=sS.Nsplines;i++) rcscoeffs[i,1] = moptimize_util_xb(M,b,ind++)
	//intercept
	int0 = moptimize_util_xb(M,b,ind)
	pS->lntxb = sS.lntxb = log(sS.t:*exp(xb))

	if (sS.orthog) {
		rcsvars = rcsgen_core(sS.lntxb,sS.knots,0,sS.rmatrix)
		drcsvars = rcsgen_core(sS.lntxb,sS.knots,1,sS.rmatrix)
		d2rcsvars = rcsgen_core(sS.lntxb,sS.knots,2,sS.rmatrix)
		d3rcs = rcsgen_core(sS.lntxb,sS.knots,3,sS.rmatrix) * rcscoeffs
		rcs = rcsvars * rcscoeffs :+ int0
		drcs = drcsvars * rcscoeffs
		d2rcs = d2rcsvars * rcscoeffs
	}
	else {
		rcsvars = rcsgen_core(sS.lntxb,sS.knots,0)
		drcsvars = rcsgen_core(sS.lntxb,sS.knots,1)
		d2rcsvars = rcsgen_core(sS.lntxb,sS.knots,2)
		d3rcs = rcsgen_core(sS.lntxb,sS.knots,3) * rcscoeffs		
		rcs = rcsvars * rcscoeffs :+ int0
		drcs = drcsvars * rcscoeffs
		d2rcs = d2rcsvars * rcscoeffs
	}

	if (sS.delentry) {
		real colvector lnt0xb, s0rcs, ds0rcs, d2s0rcs
		real matrix s0rcsvars, ds0rcsvars
		if (sS.hastvcs) lnt0xb = log(sS.t0:*exp(s0xb[sS.t0index]))
		else lnt0xb = log(sS.t0:*exp(xb[sS.t0index]))
		if (sS.orthog) {
			s0rcsvars = rcsgen_core(lnt0xb,sS.knots,0,sS.rmatrix)
			s0rcs = s0rcsvars * rcscoeffs :+ int0
			ds0rcsvars = rcsgen_core(lnt0xb,sS.knots,1,sS.rmatrix)
			ds0rcs = ds0rcsvars * rcscoeffs
			d2s0rcs = rcsgen_core(lnt0xb,sS.knots,2,sS.rmatrix) * rcscoeffs
		}
		else {
			s0rcsvars = rcsgen_core(lnt0xb,sS.knots,0)
			s0rcs = s0rcsvars * rcscoeffs :+ int0
			ds0rcsvars = rcsgen_core(lnt0xb,sS.knots,1)
			ds0rcs = ds0rcsvars * rcscoeffs		
			d2s0rcs = rcsgen_core(lnt0xb,sS.knots,2) * rcscoeffs		
		}
	}

	if (sS.hastvcs) lnfj = sS.d:*(rcs + log(drcs) :+ log(dxb :+ 1)) :- exp(rcs)
	else lnfj = sS.d:*(rcs + log(drcs)) :- exp(rcs)

	if (sS.delentry) lnfj[sS.t0index] = lnfj[sS.t0index] :+ exp(s0rcs)
	if (todo==0) return

	//Gradient
	real scalar NS
	NS = sS.Nsplines + 2
	if (sS.hastvcs) {
		NS = NS + 1
		if (sS.delentry) NS = NS + 1
	}
	S = J(sS.Nobs,NS,0)

		ind = 1
		S[,ind] = -sS.d :* ( drcs :+ d2rcs:/drcs ) :+ exp(rcs) :* drcs																//d/dxb
		if (sS.delentry & !sS.hastvcs) S[sS.t0index,ind] = S[sS.t0index,ind] :- exp(s0rcs) :* ds0rcs
		ind++
		if (sS.hastvcs) {																											
			S[,ind++] = -sS.d :/ (dxb :+ 1)																							//d/ddxb
			if (sS.delentry) S[sS.t0index,ind++] = -exp(s0rcs) :* ds0rcs															//d/ds0xb
		}
		S[,ind..(ind-1+sS.Nsplines)] = sS.d:*(rcsvars :+ drcsvars:/drcs) :- exp(rcs) :* rcsvars										//d/dgammak
		if (sS.delentry) S[sS.t0index,ind..(ind-1+sS.Nsplines)] = S[sS.t0index,ind..(ind-1+sS.Nsplines)] :+ s0rcsvars :*exp(s0rcs)
		ind++
		S[,NS] = sS.d :- exp(rcs)																									//d/dintercept
		if (sS.delentry) S[sS.t0index,NS] = S[sS.t0index,NS] :+ exp(s0rcs)

	if (todo==1) return

	//Hessian
	real colvector Hxb, Hgamma0
	real matrix Hgammas, Hgammasxb
		
		//d2/dxb2
		Hxb = sS.d :* d2rcs :+ sS.d :* d3rcs :/ drcs :- sS.d :* d2rcs:^2 :/ drcs:^2 :- exp(rcs) :* (d2rcs :+ drcs:^2)
		if (sS.delentry & !sS.hastvcs) Hxb[sS.t0index,] = Hxb[sS.t0index,] :+ ds0rcs:^2 :* exp(s0rcs) :+ exp(s0rcs):*d2s0rcs
		//d2/ddxb2
		if (sS.hastvcs) {
			real matrix Hdxb
			Hdxb = -sS.d :* (dxb :+ 1):^(-2)
			//d2/ds0xb2
			if (sS.delentry) {
				real matrix Hs0xb
				Hs0xb = J(sS.Nobs,1,0)
				Hs0xb[sS.t0index,] = Hs0xb[sS.t0index,] :+ exp(s0rcs) :* d2s0rcs :+ ds0rcs:^2 :* exp(s0rcs)
			}
		}
		//d2/dgammak2
		Hgammas = -sS.d :* drcs:^(-2) :* drcsvars:^2 :- rcsvars:^2 :* exp(rcs)
		if (sS.delentry) Hgammas[sS.t0index,] = Hgammas[sS.t0index,] :+ s0rcsvars:^2 :* exp(s0rcs)
		//d2/dintercept2
		Hgamma0 = -exp(rcs)
		if (sS.delentry) Hgamma0[sS.t0index,] = Hgamma0[sS.t0index,] :+ exp(s0rcs)

		//d2/ddxb dxb
		if (sS.hastvcs) {
			real colvector Hdxbxb
			Hdxbxb = J(sS.Nobs,1,0)
			//d2/ds0xb dxb
			if (sS.delentry) {
				real matrix Hs0xbxb
				Hs0xbxb = J(sS.Nobs,1,0)
			}
		}		
		//d2/dgammakdxb
		Hgammasxb = -sS.d :* drcsvars :- sS.d :* (d2rcsvars :* drcs :- drcsvars :* d2rcs):/(drcs:^2) :+ exp(rcs) :* drcsvars :+ drcs :* exp(rcs) :* rcsvars
		if (sS.delentry & !sS.hastvcs) Hgammasxb[sS.t0index,] = Hgammasxb[sS.t0index,] :- exp(s0rcs) :* ds0rcsvars :- ds0rcs :* exp(s0rcs) :* s0rcsvars
		//d2/dgamma0dxb
		Hgamma0xb = exp(rcs) :* drcs
		if (sS.delentry & !sS.hastvcs) Hgamma0xb[sS.t0index,] = Hgamma0xb[sS.t0index,] :- ds0rcs :* exp(s0rcs)

		//d2/ds0xb ddxb
		//d2/dgammak ddxb
		//d2/dgamma0 ddxb
		if (sS.hastvcs) {			
			real colvector Hs0xbdxb, Hgamma0dxb
			real matrix Hgammasdxb
			Hs0xbdxb = Hgamma0dxb = J(sS.Nobs,1,0)
			Hgammasdxb = J(sS.Nobs,sS.Nsplines,0)
			
			//d2/dgammak ds0xb
			//d2/dgamma0 ds0xb
			if (sS.delentry) {
				real matrix Hgammass0xb
				real colvector Hgamma0s0xb
				Hgammass0xb = J(sS.Nobs,sS.Nsplines,0)
				Hgammass0xb[sS.t0index,] = Hgammass0xb[sS.t0index,] :- exp(s0rcs) :* ds0rcsvars :- ds0rcs :* exp(s0rcs) :* s0rcsvars
				Hgamma0s0xb = J(sS.Nobs,1,0)
				Hgamma0s0xb[sS.t0index,] = Hgamma0s0xb[sS.t0index,] :- exp(s0rcs) :* ds0rcs
			}
		}

		//d2/dgamma0 dgammak
		real matrix Hgamma0gammas
		Hgamma0gammas = -rcsvars :* exp(rcs)
		if (sS.delentry) Hgamma0gammas[sS.t0index,] = Hgamma0gammas[sS.t0index,] :+ s0rcsvars :* exp(s0rcs)
		
		// d2/ dgamma not k dgammak
		//!!

		//build H
			real matrix cov
			
			ind = 1
			H =	moptimize_util_matsum(M,ind,ind,Hxb,lnfj[1])
			ind++
			if (sS.hastvcs) {
				cov = moptimize_util_matsum(M,1,ind,Hdxbxb,lnfj[1])
				H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hdxb,lnfj[1])
				ind++
				if (sS.delentry) {	
					cov = moptimize_util_matsum(M,1,ind,Hs0xbxb,lnfj[1])\moptimize_util_matsum(M,2,ind,Hs0xbdxb,lnfj[1])
					H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hs0xb,lnfj[1])
					ind++
					
					//splines  	//note different to other spline sections due to delayed entry being handled with a separate equation
					for (i=1;i<=sS.Nsplines;i++) {
						cov = moptimize_util_matsum(M,1,ind,Hgammasxb[,i],lnfj[1])  		//d2/ dgammak dxb
						cov = cov\moptimize_util_matsum(M,2,ind,Hgammasdxb[,i],lnfj[1])  	//d2/ dgammak ddxb
						cov = cov\moptimize_util_matsum(M,3,ind,Hgammass0xb[,i],lnfj[1])  	//d2/ dgammak ds0xb
						if (i>1) {
							index = 4
							while (index<(i+3)) {		//i is k, (index-1) is not k
								index2 = index-3
								tempcov = -sS.d :* drcsvars[,i] :* drcs:^(-2) :* drcsvars[,index2] :- rcsvars[,i] :* exp(rcs) :* rcsvars[,index2] 
								tempcov[sS.t0index,] = tempcov[sS.t0index,] :+ s0rcsvars[,i] :* exp(s0rcs) :* s0rcsvars[,index2]
								cov = cov\moptimize_util_matsum(M,index,ind,tempcov,lnfj[1])
								index++
							}
						}
						H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hgammas[,i],lnfj[1])
						ind++
					}

					//gamma0
					cov = moptimize_util_matsum(M,1,ind,Hgamma0xb,lnfj[1])
					cov = cov\moptimize_util_matsum(M,2,ind,Hgamma0dxb,lnfj[1])
					cov = cov\moptimize_util_matsum(M,3,ind,Hgamma0s0xb,lnfj[1])
					for (i=1;i<=sS.Nsplines;i++) {
						cov = cov\moptimize_util_matsum(M,3+i,ind,Hgamma0gammas[,i],lnfj[1])
					}
					H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hgamma0,lnfj[1])
					
				}
				else {
					//splines
					for (i=1;i<=sS.Nsplines;i++) {
						cov = moptimize_util_matsum(M,1,ind,Hgammasxb[,i],lnfj[1])  		//d2/ dgammak dxb
						cov = cov\moptimize_util_matsum(M,2,ind,Hgammasdxb[,i],lnfj[1])  	//d2/ dgammak ddxb
						if (i>1) {
							index = 3
							while (index<(i+2)) {		//i is k, (index-1) is not k
								index2 = index-2
								cov = cov\moptimize_util_matsum(M,index,ind,-sS.d :* drcsvars[,i] :* drcs:^(-2) :* drcsvars[,index2] :- rcsvars[,i] :* exp(rcs) :* rcsvars[,index2],lnfj[1])
								index++
							}
						}
						H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hgammas[,i],lnfj[1])
						ind++
					}

					//gamma0
					cov = moptimize_util_matsum(M,1,ind,Hgamma0xb,lnfj[1])
					cov = cov\moptimize_util_matsum(M,2,ind,Hgamma0dxb,lnfj[1])
					for (i=1;i<=sS.Nsplines;i++) {
						cov = cov\moptimize_util_matsum(M,2+i,ind,Hgamma0gammas[,i],lnfj[1])
					}
					H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hgamma0,lnfj[1])
				}
			}
			else {
				//splines
				for (i=1;i<=sS.Nsplines;i++) {
					//d2/ dgammak dxb
					cov = moptimize_util_matsum(M,1,ind,Hgammasxb[,i],lnfj[1])  
					if (i>1) {
						index = 2
						while (index<(i+1)) {		//i is k, (index-1) is not k
							index2 = index-1
							if (sS.delentry) {
								tempcov = -sS.d :* drcsvars[,i] :* drcs:^(-2) :* drcsvars[,index2] :- rcsvars[,i] :* exp(rcs) :* rcsvars[,index2] 
								tempcov[sS.t0index,] = tempcov[sS.t0index,] :+ s0rcsvars[,i] :* exp(s0rcs) :* s0rcsvars[,index2]
								cov = cov\moptimize_util_matsum(M,index,ind,tempcov,lnfj[1])
							}
							else cov = cov\moptimize_util_matsum(M,index,ind,-sS.d :* drcsvars[,i] :* drcs:^(-2) :* drcsvars[,index2] :- rcsvars[,i] :* exp(rcs) :* rcsvars[,index2],lnfj[1])
							index++
						}
					}
					H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hgammas[,i],lnfj[1])
					ind++
				}
				
				//gamma0
				cov = moptimize_util_matsum(M,1,ind,Hgamma0xb,lnfj[1])
				for (i=1;i<=sS.Nsplines;i++) {
					cov = cov\moptimize_util_matsum(M,1+i,ind,Hgamma0gammas[,i],lnfj[1])
				}
				H = H,cov\cov',moptimize_util_matsum(M,ind,ind,Hgamma0,lnfj[1])
			}
}
end

