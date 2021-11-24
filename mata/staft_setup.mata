*! version 1.0.0 ?????2014 MJC

version 15.1

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix

mata:

struct staft_main {
	real colvector t, d, t0, t0index, knotindex, lntxb
	`SS' touse, t0ind
	`RS' Nsplines, delentry, hastvcs, orthog, Nobs, df, varyknots, index
	`NM' rmatrix
	real rowvector knots
}

void staft_setup()
{
	struct staft_main scalar S
	pointer scalar p

	//get tempname and create struct
	stata("tempname staft_struct")
	rmexternal(st_local("staft_struct"))
	p = crexternal(st_local("staft_struct"))

	S.touse = st_local("touse")	
	S.hastvcs = st_local("tvc")!=""
	
	S.t = st_data(.,"_t",S.touse)
	S.d = st_data(.,"_d",S.touse)
	S.Nobs = rows(S.t)
	
	S.orthog = st_local("noorthog")==""
	if (S.orthog) S.rmatrix = st_matrix(st_local("rmat"))

	S.varyknots = st_local("adapt")!=""
	S.lntxb = log(S.t)
	
	S.knots = strtoreal(tokens(st_local("ln_bknots")))
	S.index = 0
	if (S.varyknots) {
		//get indexes from df and Nobs
		S.df = cols(S.knots)-1
		if (S.df==0) {
			S.df=1
			S.knotindex = 0
		}
		
		if (S.df==2) S.knotindex = 0.50
		else if (S.df==3) S.knotindex = (0.33,0.67)'
		else if (S.df==4) S.knotindex = (0.25,0.50,0.75)'
		else if (S.df==5) S.knotindex = (0.20,0.40,0.60,0.80)'
		else if (S.df==6) S.knotindex = (0.17,0.33,0.50,0.67,0.83)'
		else if (S.df==7) S.knotindex = (0.14,0.29,0.43,0.57,0.71,0.86)'
		else if (S.df==8) S.knotindex = (0.125,0.25,0.375,0.50,0.625,0.75,0.875)'
		else if (S.df==9) S.knotindex = (0.111,0.222,0.333,0.444,0.556,0.667,0.778,0.889)'
		else if (S.df==10) S.knotindex = (0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90)'
			
		S.knotindex = round(S.knotindex * sum(S.d))	//restricted to events
	}
	else S.df = 0
	S.Nsplines = strtoreal(st_local("Nsplines"))	
	
	S.delentry = strtoreal(st_local("delentry"))
	if (S.delentry) {
		S.t0ind = st_local("t0ind")
		S.t0index = st_data(.,st_local("index"),S.t0ind)
		S.t0 = st_data(.,"_t0",S.t0ind)
	}

	//Done 	
	swap((*p), S)
}

void staft_prolog(real rowvector b,scalar M,real scalar lnl)
{
	pointer(struct staft_main scalar) scalar pS
	pS = &moptimize_util_userinfo(M,1)
	
	pS->knots = GetNewKnots((*pS).lntxb,(*pS).d,(*pS).df,(*pS).knotindex)
	if ((*pS).orthog) {
		pS->rmatrix = GetNewRmatrix((*pS).lntxb,(*pS).knots)
	}

}

end

exit
