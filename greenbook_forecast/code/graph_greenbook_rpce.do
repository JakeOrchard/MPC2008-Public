
set scheme s1color
set fredkey b2d9da12107485ec086c5691d41d626c



//Actual Q/Q AR RPCE spending

import fred PCE PCEPI,  clear 
gen idate = date(datestr,"YMD")
gen qdate = qofd(idate)
gen pcec96 = 100*PCE/PCEPI //(REAL CONSUMPTION)
collapse (mean) pcec96 , by(qdate)
format qdate %tq
tsset qdate
gen delta_qq_rpce = ((pcec96/l.pcec96)-1)
gen drpce_actual = 100*(((1+delta_qq_rpce)^4)-1)
gen drpce_notannualized =  100*(((1+delta_qq_rpce))-1)
keep qdate drpce*
tempfile rpce
save `rpce'





use ../output/greenbook_rpce.dta, clear
tostring Date, replace
gen qdate = quarterly(Date, "YQ")
format qdate %tq


merge 1:1 qdate using `rpce', nogen keep(match)

*Only keeps 2008/09 values
keep if qdate >= yq(2008,1) & qdate <= yq(2009,4)
tsset qdate

*Labels used surveys
lab var gRPCE_20071205 "December 05, 2007"
lab var gRPCE_20080123 "January 23, 2008"
lab var drpce_actual "Actual RPCE Growth"
	
	




/***************************************************************************
Graphs RPCE growth Forecast for2008 using 2007m12 and 2008m1 Greeenbooks
****************************************************************************/



tw (scatter gRPCE_20071205 qdate, c(l l) clp(l) clw(medthick) xline(4, lp(-) lc(red))   xline(5, lp(-) lc(blue)) ///
    ytitle("", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter gRPCE_20080123  qdate, c(l l) clp(l) clw(medthick)  clc(blue) mc(blue)) ///
    , name(greenbook_rpce08, replace) scale(1.2) ///
    xlabel(, valuelabel) xtitle(Forecast) title(Greenbook: Q/Q Real Consumption Growth) ytitle(Percent Growth)


graph export ../output/greenbook_fore_2008.pdf, as(pdf) replace
cap graph export ../output/greenbook_fore_2008.png,  replace




tw (scatter gRPCE_20071205 qdate, c(l l) clp(l) clw(medthick) xline(4, lp(-) lc(red))   xline(5, lp(-) lc(blue)) ///
    ytitle("", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter gRPCE_20080123  qdate, c(l l) clp(l) clw(medthick)  clc(blue) mc(blue)) ///
    (scatter drpce_actual qdate, c(l l) clp(l) clw(medthick)  clc(black) mc(black)) ///
    , name(greenbook_rpce08_wactual, replace) scale(1.2) ///
    xlabel(, valuelabel) xtitle(Forecast) title(Greenbook: Q/Q Real Consumption Growth) ytitle(Percent Growth)


graph export ../output/greenbook_fore_2008_wactual.pdf, as(pdf) replace
cap graph export ../output/greenbook_fore_2008_wactual.png,  replace


    
    
    graph export ../output/SPF_dist_wactual.pdf, as(pdf) replace
    cap graph export ../output/SPF_dist_wactual.png, replace
    
    
 keep gRPCE_20071205 gRPCE_20080123 qdate
 save ../output/greenbook08.dta, replace




