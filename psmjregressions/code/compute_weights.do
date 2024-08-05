/***
This do file computes weights and saves them as a tempfile `weights'.
Call this do file using do compute_weights "weightversion" where weightversion is
either OLS or Sample.

*/


preserve
if "`weighttype'" == "OLS" {
	eventstudyweights `rebatevars' [w=finlwt21] if insample`year' `rebatesample'  , cohort(intdate) rel_time(rel_interview)  absorb(`absorbvars') covariates(`controlvars' `extravar') 
	mat A = e(weights)

	clear
	svmat A, names(col)

	keep `rebatevars' rel_interview intdate
	rename rel_interview rint
	rename intdate interview
	
	local lagrbt = "lag1rbt" + "`year'" +"indicator"
	local rbt = "rbt" + "`year'" +"indicator"
	
	*Make sure weights sums to one
	cap egen sumweightlag = sum(`lagrbt') if rint == 3
	cap replace `lagrbt' = `lagrbt'/sumweightlag
	
	egen sumweight = sum(`rbt') if rint == 0
	replace `rbt' = `rbt'/sumweight 

	
	rename interview intdate

	cap replace `lagrbt' = . if rint != 3
	replace `rbt' = . if rint != 0
	collapse (sum)  `rebatevars', by(intdate )
	rename `rbt' wconOLS
	cap rename `lagrbt' wlagOLS
	tempfile weights
	save `weights', replace	
}

if "`weighttype'" == "Sample" {
	cd ../../did_imputation		

	
			
	*Weights in lag specification
	did_imputation `depvar' cuid intdate firstrbtintdate [w=finlwt21] if insample`year' `rebatesample', cluster(cuid) fe(`absorbvars') controls(`controlvars' `extravar') horizons(0 3) saveweights  autosample maxit(10000)
	
	gen treat_status = 0
	replace treat_status = 1 if firstrbtintdate == intdate
	replace treat_status = 2 if firstrbtintdate < intdate


	collapse (sum) __w_tau0 __w_tau3, by(intdate treat_status)
	cd ../psmjregressions/code
	
	rename __w_tau0 weightBJScon
	rename __w_tau3 weightBJSlag
		
	keep weight* intdate 
	collapse (max) weightBJScon weightBJSlag, by(intdate)
	egen totalweightcon = sum(weightBJScon)
	egen totalweightlag = sum(weightBJSlag)
	replace weightBJScon = weightBJScon/totalweightcon
	replace weightBJSlag = weightBJSlag/totalweightlag
	tempfile weights
	save `weights', replace		

}

restore





	
