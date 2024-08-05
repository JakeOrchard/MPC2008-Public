


/*******************************************************************************

Step 0: Load Data

*******************************************************************************/

clear all
set more off

use ../output/psmjsampleinterview_wlabels.dta, clear


/*******************************************************************************

Step 1: Aggregate to per household spending on vehicles

*******************************************************************************/

keep if insample

replace intdate = mofd(dofm(intdate))

local plotvars nipa_mv_parts

local mindate = ym(2007,12)
local treatdate = ym(2008,6)
local maxdate = ym(2009,3)

foreach var of varlist `plotvars' {

	reg `var' ibn.intdate#ibn.everrbtindicator [aw=finlwt21], nocons cluster(cuid)
	
	* tests for equality run separate regression for t-test
	local testexplhs "`treatdate'.intdate#0.everrbtindicator"
	local testexprhs "`treatdate'.intdate#1.everrbtindicator"
	test `treatdate'.intdate#0.everrbtindicator = `treatdate'.intdate#1.everrbtindicator
	forvalues tt = `=`treatdate'+1'(1)`=ym(2008,12)' {
		local testexplhs "`testexplhs' + `tt'.intdate#0.everrbtindicator"
		local testexprhs "`testexprhs' + `tt'.intdate#1.everrbtindicator"
		test `tt'.intdate#0.everrbtindicator = `tt'.intdate#1.everrbtindicator
	}
	disp "`testexplhs' = `testexprhs'"
	test `testexplhs' = `testexprhs'
	
	
	
	
	matrix coeff = J(`=`maxdate' - `mindate' + 1', 7, .)
	
	forvalues tt = `mindate'(1)`maxdate' {
		local hh = `tt' - `mindate' + 1
		
		forvalues aa = 0(1)1 {
			local select "`tt'.intdate#`aa'.everrbtindicator"
			
			matrix coeff[`hh', 1 + `aa'*3] = _b[`select']
			matrix coeff[`hh', 2 + `aa'*3] = _b[`select'] + 1.96 * _se[`select']
			matrix coeff[`hh', 3 + `aa'*3] = _b[`select'] - 1.96 * _se[`select']
		}
		matrix coeff[`hh', 7] = `tt'
	}
	
	matrix colnames coeff = control control_ciu control_cil treat treat_ciu treat_cil h 
	
	clear
	qui svmat coeff, names(col)
	
	tsset h, monthly
	
	graph twoway ///
	(tsline treat control , lwidth(medthick medthick) lcolor(blue red)) ///
	(rarea treat_ciu treat_cil h, fcolor(blue%30) lwidth(none)) ///
	(rarea control_ciu control_cil h, fcolor(red%30) lwidth(none)), ///
	ylabel(,nogrid tposition(outside) angle(horizontal))  ///
	legend(order(1 2) label(1 "Rebate Group") label(2 "Never Rebate Group")  cols(1) ring(0) position(2)) ///
	ytitle("Motor Vehicle Expenditure per Household")  ttitle("Interview Month") ///
	xline(`treatdate', lcolor(red) lpattern(dash)) tlabel(, format(%tmMonCCYY) ) ///
	graphregion(color(white))
	

	
	graph export ../output/motorvehicle_eventstudy.pdf, replace

}





// collapse (sum) `plotvars' finlwt21, by(intdate everrbtindicator)

// foreach var of varlist `plotvars' {
// 	replace `var' = `var' / finlwt21
// }

// tsset everrbtindicator intdate, monthly

