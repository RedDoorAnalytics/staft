//build new version for SSC and zip up
local drive /Users/Michael/Documents/reddooranalytics/products
cd `drive'/staft

local includemata = 0

local vsnumber 1_0_1										//remove _temp for release
local newssc ./ssc/version_`vsnumber'/

cap mkdir `newssc'

//description
copy ./ssc/staft_details.txt `newssc', replace

//data
//copy ./data/staft_example.dta `newssc', replace
	
//staft
copy ./staft.ado `newssc', replace
copy ./staft.sthlp `newssc', replace

if `includemata' {
	copy ./staft_setup.mata `newssc', replace
	copy ./staft_lf2.mata `newssc', replace
	copy ./rcsgen_mata.mata `newssc', replace
}

//predictions
copy ./staft_pred.ado `newssc', replace
copy ./staft_postestimation.sthlp `newssc', replace
copy ./rcsgen2.ado `newssc', replace

//mlib
copy ./lstaft.mlib `newssc', replace

