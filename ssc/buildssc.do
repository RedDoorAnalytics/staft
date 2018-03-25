//build new version for SSC and zip up
local drive Z:
cd `drive'/stdev/staft

local includemata = 0

local vsnumber 0_1_2										//remove _temp for release
local newssc ./ssc/version_`vsnumber'`temp'/
local newsscnoslash ./ssc/version_`vsnumber'

cap mkdir `newssc'

//description
copy .\ssc\staft_details.txt `newssc', replace

//data
//copy .\data\staft_example.dta `newssc', replace
	
//staft
copy .\staft\staft.ado `newssc', replace
copy .\staft\staft.sthlp `newssc', replace

if `includemata' {
	copy .\staft\staft_setup.mata `newssc', replace
	copy .\staft\staft_lf2.mata `newssc', replace
	copy .\staft\rcsgen_mata.mata `newssc', replace
}

//predictions
copy .\staft\staft_pred.ado `newssc', replace
copy .\staft\staft_postestimation.sthlp `newssc', replace
copy .\staft\rcsgen2.ado `newssc', replace

//mlib
copy .\lstaft.mlib `newssc', replace

//zip up
//zipfile `newssc', saving(`newsscnoslash',replace) 



