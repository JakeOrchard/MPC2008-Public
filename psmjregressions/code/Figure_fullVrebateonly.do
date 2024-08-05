clear all
set more off
set scheme s1color
cap ssc install eventstudyweights

cap cd ../psmjregressions/code





/*******************************************************************************

Step 1: 

*******************************************************************************/


local dependentvars "d_pce"

local baselinecontrols "age d_num_adults d_perslt18"
local extracontrollist `"" "lag1rbtindicator" "lag1rbtindicator _Iincdinte* lag_pce lag_mv"'
local extracontrollist `""'



foreach year in ""{
	
		local rebatevar = "rbt" + "`year'" + "indicator"
		local cohortvar = "firstrbt" + "`year'" + "intdate"
		
		*Rebate dates
		if "`year'" == ""{
			local beginrebatedate = 581
			local endrebatedate = 590
			local beginsampledate = 575
			local endsampledate = 590
			
		}
		else if "`year'" == "01"{
			local beginrebatedate = 499
			local endrebatedate = 503
			local beginsampledate = 494
			local endsampledate = 506
		}


foreach depvar in `dependentvars' {
	
	
	
	
	

	use ../output/psmjsampleinterview_wlabels.dta, clear
	xi gen i.incdinterview

	egen numrbts = sum(`rebatevar'), by(cuid)
	keep if numrbts <2
	
	tsset cuid intdate

	



	*For eventstudyinteract the treatment timing variable must be missing for non-treated
	gen `cohortvar'_esi = `cohortvar'
	gen control_cohort_esi = `cohortvar'_esi == .
	gen rel_interview = intdate-`cohortvar'
	replace rel_interview = -9 if rel_interview < -7 //Never treated

	
	keep if insample`year' 
	




	qui forval rint = -6(3)9{

		local hname = `rint'
		if `rint' == -9{
			local hname  p9
		}
		if `rint' == -6{
			local hname  p6
		}
		if `rint' == -3{
			local hname  p3
		}
		
		
		forval interview = `beginsampledate'(1)`endsampledate'{
			local absdate = `interview'-`rint'
			if `absdate' >= `beginrebatedate' & `absdate' <= `endrebatedate'{
				noi di `absdate'
				gen g_`hname'_`interview' = (intdate == `interview') & (rel_interview == `rint')
			}
		}
	}



	xtset cuid intdate



/*******************************************************************************

Step 2: OLS Weights

*******************************************************************************/
foreach rebatesample in "" "& ever`rebatevar'" {
	foreach extra in "`extracontrollist'"{

if "`rebatesample'" == ""{
	local rebatename "_full"
	
	}
	else{
	
	local rebatename = "_rbt" + "`year'" + "only"
	}
	if "`extra'" == ""{
	local extraname ""
	
	}
	if "`extra'" == "lag1rbtindicator"{
	
	local extraname "lag"
	}

	if "`extra'" == "lag1rbtindicator _Iincdinte* lag_pce lag_mv"{
	
	local extraname "extra"
	}
	*baseline
	reghdfe `depvar' `rebatevar' `baselinecontrols' `extra' [w=finlwt21] if insample`year' `rebatesample' ,  absorb(intdate) 

	if "`rebatesample'" == "" | "`year'" == "01"{
		
		local satvars "g_0* g_3* g_6* g_9* g_p3* g_p6*"
		if "`year'" == "01"{
			local satvars "g_0* g_3* g_6* g_p3* g_p6*"
		}
		
	}
	else{
		local satvars "g_0* g_3*  g_6_587 g_6_590 g_p3_578 g_p3_579 g_p3_580 g_p3_581 g_p3_582 g_p3_583 g_p3_584 "
		*g_p6 is the ommitted dummy. Starting in 585, g_p3 is the ommitted dummy. Starting in period 588, all household have already been treated 
		*so the past past treated (g_6_588, g_6_589) or the past past past treated (g_9_590) is the ommitted dummy.
			}
	*Saturated Regression 
	reghdfe `depvar'  `satvars'  `baselinecontrols'  `extra' [w=finlwt21] if insample`year' `rebatesample', nocons absorb(intdate)  //


	preserve

	gen coef = .
	gen interview = .
	gen rint = .
	local i = 1

	qui forval rint = -6(3)9{

		local hname = `rint'
		if `rint' == -9{
			local hname  p9
		}
		if `rint' == -6{
			local hname  p6
		}
		if `rint' == -3{
			local hname  p3
		}
		
		
		forval interview = `beginsampledate'(1)`endsampledate'{
			local absdate = `interview'-`rint'
			if `absdate' >= `beginrebatedate' & `absdate' < `endrebatedate'{
				replace interview = `interview' if _n == `i'
				replace rint = `rint' if _n == `i'
				cap replace coef = _b[g_`hname'_`interview']	if interview == `interview' & rint == `rint'
				local i = `i' + 1
			}
		}
	}
	keep if coef ~= .
	keep coef interview rint

	tempfile gamma
	save `gamma'


	restore



	*Derive OLS weights using eventstudyweights command (virtually equivalent to above)

	preserve	
		


	eventstudyweights `rebatevar' [w=finlwt21] if insample`year' `rebatesample'  , cohort(intdate) rel_time(rel_interview)  absorb(intdate) covariates(`baselinecontrols' `extra') 
	mat A = e(weights)

	clear

	svmat A, names(col)



	keep `rebatevar' rel_interview intdate
	rename `rebatevar' w

	rename rel_interview rint
	rename intdate interview

	merge 1:1 rint interview  using `gamma'



	gen w_coef= w*coef

	egen didcoef = sum(w_coef)
	sum didcoef

	*Test weight sums

	gen treat_status = 0
	replace treat_status = 1 if rint == 0
	replace treat_status = 2 if rint >0

	egen sumweight = sum(w), by(treat_status)
	bys treat_status: sum sumweight

	rename interview intdate


	collapse (sum) w_coef`rebatename'`extraname' = w_coef w`rebatename'`extraname'= w, by(intdate treat_status)

	gen coef`rebatename'`extraname' = w_coef`rebatename'`extraname'/w`rebatename'`extraname'
	
	
	tempfile weights`rebatename'`extraname'
	save `weights`rebatename'`extraname'', replace
	restore
}
}
				
	use `weights_full', clear
	merge 1:1 intdate treat_status using `weights_rbtonly', nogen

	

	/********************************************************************************
	Cumulative Weight graphs
	********************************************************************************/

	

	reshape long w_ w_coef_ coef_, i(intdate treat_status) j(method) string


	reshape wide w_ w_coef_ coef_ , i(intdate method) j(treat_status)

	
	lab val intdate time
	
	replace method = "1" if method == "full"
	replace method = "2" if method == "rbtonly"
	

	destring method, replace
	lab def met 1 "Full Sample" 2 "Rebate Only" //For weight Graph
	lab val method met
	
	

	*Weight each interview
	replace w_2 = 0 if w_2 == .
	gen tw = -(w_0+w_2)
	gen rw0 = w_0/tw
	gen rw1 = w_1
	gen rw2 = w_2/tw


	*Avg. Coefficient
	replace coef_2 = 0 if coef_2 == .
	gen pw_coef1 = coef_1 + rw0*coef_0
	gen pw_coef2 = rw2*coef_2
	
	replace pw_coef2 = 0 if w_coef_2 == .
	
	

	gen avg_coef= pw_coef1 + pw_coef2 
	
	*Contribution of treatment compared to not treated
	gen w_coef_diff = w_coef_1 + w_coef_0

	*Defines desired time interval
	if "`year'" == ""{
		local begincutoff = ym(2008,5)
		local endcutoff = ym(2008,11)
		lab def time 581 "Jun" 582 "Jul" 583 "Aug" 584 "Sep" 585 "Oct" 586 "Nov"

	}
	else if "`year'" == "01"{
		local begincutoff = ym(2001,7)
		local endcutoff = ym(2002,1)
		lab def time 499 "Aug" 500 "Sep" 501 "Oct" 502 "Nov" 503 "Dec" 

	}
		lab val intdate time


	*Weight

	graph bar  w_1 if intdate > `begincutoff' & intdate < `endcutoff', over(method) over(intdate)   ytitle("Interview Weight")  asyvars bar(1, color(black)) bar(2, color(gs8))  
	cap graph export ../output/Figure_fullVrebateonly_weight_`depvar'`year'..png, replace
	graph export ../output/Figure_fullVrebateonly_weight_`depvar'`year'.pdf, replace
	
	
	
	*Coefficient
	graph bar  avg_coef if intdate > `begincutoff' & intdate < `endcutoff', over(method) over(intdate)   ytitle("Treatment Effect, $ Expenditure")  asyvars bar(1, color(black)) bar(2, color(gs8))  
		cap graph export ../output/Figure_fullVrebateonly_coef_`depvar'`year'.png, replace
		graph export ../output/Figure_fullVrebateonly_coef_`depvar'`year'.pdf, replace

*Coefficeint
		graph bar  pw_coef1 pw_coef2 if intdate > `begincutoff' & intdate < `endcutoff', over(method,label(angle(90)))  over(intdate, label(angle(90)))   stack ytitle("Treatment Effect, $ Expenditure")  asyvars  bar(1, color(black))   bar(2, color(maroon))  legend(order(1 "Treated - Not Yet Treated" 2 "Past Treated")) 
		
		
		cap graph export ../output/Figure_fullVrebateonly_coefstack_`depvar'`year'.png, replace
		graph export ../output/Figure_fullVrebateonly_coefstack_`depvar'`year'.pdf, replace
	

}


}



