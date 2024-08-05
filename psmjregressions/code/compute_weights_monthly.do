/***
This do file computes weights and saves them as a tempfile `weights'.
Call this do file using do compute_weights "weightversion" where weightversion is
either OLS or Sample.

*/


preserve
if "`weighttype'" == "OLS" {
	eventstudyweights `rebatevars' [w=finlwt21] if insample`year' `rebatesample'  , cohort(date) rel_time(rel_month)  absorb(`absorbvars') covariates(`controlvars' `extravar') 
	mat A = e(weights)

	clear
	svmat A, names(col)

	keep `rebatevars' rel_month date
	rename rel_month rint
	rename date month
	
	
	*Make sure weights sums to one
	local h = -1
	foreach rbt of varlist `rebatevars'{
		local hname "`h'"
		if "`rbt'" == "lead1rbtindicator"{
			local hname "lead1"
		}
		egen sumweight = sum(`rbt') if rint == `h'
		replace `rbt' = `rbt'/sumweight 
		replace `rbt' = . if rint != `h'
		
		rename `rbt' wOLS`hname'
		drop sumweight
		local h = `h' + 1
	}
	
	rename month date

	
	collapse (sum)  wOLS*, by(date)


	tempfile weights
	save `weights', replace	
}

restore





	
