
set scheme s1color
set fredkey b2d9da12107485ec086c5691d41d626c


local yearlist `"2001" "2008"'

foreach year in "`yearlist'"{
	if "`year'" == "2008"{
		local spfdate = yq(2007,4)
	}
	else{
		local spfdate = yq(2001,2)
	}
	import fred PCE PCEPI,  clear 
	gen idate = date(datestr,"YMD")
	gen qdate = qofd(idate)
	gen pcec96 = 100*PCE/PCEPI //(REAL CONSUMPTION)
	collapse (mean) pcec96 , by(qdate)
	format qdate %tq
	tsset qdate
	gen _rpcebase = pcec96 if qdate == `spfdate'
	egen rpcebase = max(_rpcebase)
	gen delta_rpce = 100*((pcec96/rpcebase)-1)
	gen delta_yy = 100*((pcec96/l4.pcec96)-1)
	keep qdate delta* pcec96
	tempfile rpce
	save `rpce'



use ../output/spf_rcons.dta, clear



gen qdate = yq(YEAR,QUARTER)




/***************************************************************************
Generates Percent change for each individual forecaster
****************************************************************************/

forval i = 2/6{

gen rpc_CP`i' = 100*((RCONSUM`i'/RCONSUM2)-1)
gen rpc_PP`i' = 100*((RCONSUM`i'/RCONSUM1)-1)


}



keep if qdate == `spfdate'

reshape long rpc_CP rpc_PP RCONSUM, i(index) j(forecast)

*Normalize so that all forecasters have same real 2007Q4 forecasts
gen _RCbase = RCONSUM if forecast == 2
egen RCbase = max(_RCbase), by(index) 
sum RCONSUM if forecast == 2
replace RCONSUM = r(mean)*(RCONSUM/RCbase)


collapse (min) rpcmin = rpc_CP RCONS_min = RCONSUM (p25) rpc_p25 = rpc_CP (p50) rpc_p50 = rpc_CP RCONS_med = RCONSUM (p75) rpc_p75 = rpc_CP , by(forecast)

tsset forecast

gen qdate = `spfdate'
replace qdate = `spfdate' - 1 if forecast == 1
replace qdate = `spfdate' + 1 if forecast == 3
replace qdate = `spfdate' + 2 if forecast == 4
replace qdate = `spfdate' + 3 if forecast == 5
replace qdate = `spfdate' + 4 if forecast == 6

if "`year'" == "2001"{
	drop if forecast == 6
}

merge 1:1 qdate using `rpce', nogen keep(match)
format qdate %tq

	if "`year'" == "2008"{
		label define frombase 2 "2007 Q4" 3 "2008 Q1" 4 "2008 Q2" 5 "2008 Q3" 6 "2008 Q4"
		local end = "0704"
		local rebatestart = 4

	}
	else{
	label define frombase 1 "2001 Q1" 2 "2001 Q2" 3 "2001 Q3" 4 "2001 Q4" 5 "2002 Q1"
	local end = "0102"
	local rebatestart = 3

	}
label val forecast frombase



tw (scatter rpcmin forecast, c(l l) clp(l) clw(medthick) xline(`rebatestart', lp(-) lc(red))    ///
    ytitle("", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter rpc_p25 forecast, c(l l) clp(l) clw(medthick)  clc(green) mc(green)) ///
    (scatter rpc_p50 forecast, c(l l) clp(l) clw(medthick)  clc(blue) mc(blue)) ///
    (scatter rpc_p75 forecast, c(l l) clp(l) clw(medthick)  clc(brown) mc(brown)) ///
    , name(SPF_rcons, replace) scale(1.2) ///
    legend(order(1 "Minimum" 2 "P25" 3 "Median" 4 "p75" )) xlabel(, valuelabel) xtitle(Forecast) title(SPF: Real Consumption Growth) ytitle(Percent)


graph export ../output/SPF_dist_`end'.pdf, as(pdf) replace
cap graph export ../output/SPF_dist_`end'.png,  replace

tw (scatter rpcmin forecast, c(l l) clp(l) clw(medthick) xline(`rebatestart', lp(-) lc(red))    ///
    ytitle("", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter rpc_p25 forecast, c(l l) clp(l) clw(medthick)  clc(green) mc(green)) ///
    (scatter rpc_p50 forecast, c(l l) clp(l) clw(medthick)  clc(blue) mc(blue)) ///
    (scatter rpc_p75 forecast, c(l l) clp(l) clw(medthick)  clc(brown) mc(brown)) ///
    (scatter delta_rpce forecast, c(l l) clp(l) clw(medthick)  clc(black) mc(black)) ///
    , name(SPF_rcons, replace) scale(1.2) ///
    legend(order(1 "Minimum" 2 "P25" 3 "Median" 4 "p75" 5 "Actual") cols(1) ring(0) position(10)) xlabel(, valuelabel) xtitle(Forecast)  ytitle(Percent)


graph export ../output/SPF_dist_`end'_wactual.pdf, as(pdf) replace
cap graph export ../output/SPF_dist_`end'_wactual.png,  replace

keep qdate RCONS* pcec96

save ../output/levels_SPF`end'.dta, replace

	

}
/*
/*****************************************************************************
Minimum growth forecast from survey qdate to a year later
****************************************************************************/


collapse (min) rpc6min = rpc_CP6 (p25) rpc6_p25 = rpc_CP6 (p50) rpc6_p50 = rpc_CP6 (p75) rpc6_p75 = rpc_CP6 , by(qdate)

tsset qdate
format qdate %tq

merge 1:1 qdate using `rpce'





keep if qdate >=yq(2007,1) & qdate <= yq(2010,1)



tw (scatter rpc6min qdate, c(l l) clp(l) clw(medthick) xline(`=yq(2008,2)', lp(-) lc(red))   xline(`=yq(2008,3)', lp(-) lc(blue)) ///
    ytitle("", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter rpc6_p25 qdate, c(l l) clp(l) clw(medthick)  clc(green) mc(green)) ///
    (scatter rpc6_p50 qdate, c(l l) clp(l) clw(medthick)  clc(blue) mc(blue)) ///
    (scatter rpc6_p75 qdate, c(l l) clp(l) clw(medthick)  clc(brown) mc(brown)) ///
    , name(SPF_rcons, replace) scale(1.2) ///
    legend(order(1 "Minimum" 2 "P25" 3 "Median" 4 "p75" )) title(SPF: Real Consumption Growth Y/Y Forecast) xtitle(qdate) ytitle(Percent)
    
    
    graph export ../output/SPF_dist.pdf, as(pdf) replace
    cap graph export ../output/SPF_dist.png, replace

tw (scatter rpc6min qdate, c(l l) clp(l) clw(medthick) xline(`=yq(2008,2)', lp(-) lc(red))   xline(`=yq(2008,3)', lp(-) lc(blue)) ///
    ytitle("", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter rpc6_p25 qdate, c(l l) clp(l) clw(medthick)  clc(green) mc(green)) ///
    (scatter rpc6_p50 qdate, c(l l) clp(l) clw(medthick)  clc(blue) mc(blue)) ///
    (scatter rpc6_p75 qdate, c(l l) clp(l) clw(medthick)  clc(brown) mc(brown)) ///
    (scatter delta_yy qdate, c(l l) clp(l) clw(medthick)  clc(black) mc(black)) ///
    , name(SPF_rcons, replace) scale(1.2) ///
    legend(order(1 "Minimum" 2 "P25" 3 "Median" 4 "p75" 5 "Actual")) title(SPF: Real Consumption Growth Y/Y Forecast) xtitle(qdate) ytitle(Percent)
    
    
    graph export ../output/SPF_dist_wactual.eps, replace
    cap graph export ../output/SPF_dist_wactual.png, replace
*/
