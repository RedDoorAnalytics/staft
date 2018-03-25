//local drive Z:
local drive /Users/Michael/Documents


cd "`drive'/stdev/staft"
adopath ++ "./staft"
clear all

do ./build/buildmlib.do

webuse brcancer
stset rectime, f(censrec) scale(365.25)

/*
staft hormon, df(3) 

predict s1, survival
predict h1, hazard
predict ch1, cumhaz
predict af1, af

predict s12, survival ci
predict h12, hazard ci
predict ch12, cumhaz ci
predict af12, af ci

//expand 100
timer clear
timer on 1
staft hormon, df(3) //eform mlmethod(lf2debug)
timer off 1
timer list
*/

forvalues i=1/5 {
	qui staft hormon, df(3) tvc(hormon) dftvc(`i') 
	di e(AIC)
}
staft hormon, df(3)
staft hormon, df(3) tvc(hormon) dftvc(1) 

predict af1, af at(hormon 1) ci

scatter af1* _t 

twoway (rarea af1_lci af1_uci _t if _t<=7,sort)(line af1 _t if _t<=7,sort)	,	///
	xtitle("Follow-up time (years)") ytitle("Acceleration factor") yscale(log)	///
	ylabel(,angle(h) format(%2.1f))	xlabel(0(1)7)		///
	legend(order(2 "Acceleration factor" 1 "95% CI"))

//tr:staft hormon, df(9)  tvc(hormon) dftvc(3) debug
//staft hormon, df(9) 

stpm2 hormon, df(3) tvc(hormon) dftvc(1) scale(h)
