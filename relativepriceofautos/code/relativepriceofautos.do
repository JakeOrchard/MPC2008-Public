set fredkey b2d9da12107485ec086c5691d41d626c

* get experimental series
import excel using ../input/r-cpi-u-nv-data.xlsx, sheet("NewVehiclesUntaxed") case(lower) clear firstrow

keep cp nationnewvehiclesuntaxed
rename nationnewvehiclesuntaxed cpinewvehiclesresearchnsa

gen date = ym(floor(cp/100),cp - floor(cp/100)*100)

drop cp
tsset date, monthly

tempfile research
save `research'

* relative new car price
clear
import fred CUSR0000SETA01 CUUR0000SETA01 CPILFESL

rename CUSR0000SETA01 cpinewvehicles
rename CUUR0000SETA01 cpinewvehiclesnsa
rename CPILFESL cpicore

gen cpiseasonals = cpinewvehicles / cpinewvehiclesnsa


drop datestr
gen date = mofd(daten)
tsset date, monthly

merge 1:1 date using `research', nogenerate

* alternative SA
gen month = month(dofm(date))
gen lcpinewvehiclesresearchnsa = log(cpinewvehiclesresearchnsa)
reg lcpinewvehiclesresearchnsa i.month
predict lcpinewvehiclesresearchsadum, r
gen cpinewvehiclesresearchsad = exp(lcpinewvehiclesresearchsadum)

gen cpinewvehiclesresearch = cpinewvehiclesresearchnsa * cpiseasonals



foreach var in newvehicles newvehiclesresearch newvehiclesresearchsad {
	su cpi`var' if date == ym(2008,1)
	replace cpi`var' = cpi`var' / `=r(mean)' * 100
}



replace cpinewvehiclesresearch = cpinewvehicles if date<=ym(2007,12)
replace cpinewvehiclesresearchsad = cpinewvehicles if date<=ym(2007,12)

foreach var in newvehicles newvehiclesresearch newvehiclesresearchsad {	
	gen `var'relp = cpi`var' / cpicore
	su `var'relp if date == ym(2008,1)
	replace `var'relp = `var'relp / `=r(mean)' * 100
	
	* create simple forecast
	reg `var'relp date if date>=ym(2007,12) & date<=ym(2008,4)
	predict `var'relpfc
	
	* create simple forecast
	reg `var'relp date if date>=ym(2007,8) & date<=ym(2008,4)
	predict `var'relpfclong
}


tw scatter newvehiclesresearchrelp date if tin(2007m1,2009m12), ///
	xline(`=ym(2008,5)', lcolor(red) lpattern(dash)) ///
	title("Relative New Motor Vehicle Price") ///
	c(l ) clp(l ) ms(d ) clw(medthick ) ytitle("Index, Dec. 2007 = 100") ///
	clc(black) mc(black) xtitle("month")  scale(1.2) ysize(4) xsize(7)

graph export ../output/newvehiclerelpresearch_format.eps, replace	




tsline newvehiclesresearchrelp newvehiclesresearchrelpfc newvehiclesresearchrelpfclong  if tin(2007m8,2008m9), ///
			tscale(range(2007m8 2008m9)) ///
			tlabel(2007m8(3)2008m9) ///
			lwidth(medthick medthick medthick) ///
			lpattern(solid solid) lcolor(black red blue) ///
			ylabel(,nogrid tposition(outside) angle(horizontal))  ///
			legend(label(1 "Relative New Vehicle Price") label(2 "Trend: Dec 2007 through April 2008") label(3 "Trend: Aug 2007 through April 2008")  cols(1)) ///
			ytitle("Index, Dec. 2007 = 100")  ttitle("month") ///
			xline(`=ym(2008,5)', lcolor(red) lpattern(dash)) tlabel(, format(%tmMonCCYY) )


graph export ../output/newvehiclerelpresearch.eps, replace




tsline newvehiclesresearchrelp  newvehiclesresearchrelpfclong  if tin(2007m8,2008m9), ///
			tscale(range(2007m8 2008m9)) ///
			tlabel(2007m8(3)2008m9) ///
			lwidth(medthick medthick) ///
			lpattern(solid solid) lcolor(black blue) ///
			ylabel(,nogrid tposition(outside) angle(horizontal))  ///
			legend(label(1 "Relative New Vehicle Price") label(2 "Trend: Aug 2007 through April 2008")  cols(1)) ///
			ytitle("Index, Dec. 2007 = 100")  ttitle("month") ///
			xline(`=ym(2008,5)', lcolor(red) lpattern(dash)) tlabel(, format(%tmMonCCYY) )


graph export ../output/newvehiclerelpresearchlong.eps, replace

tsline newvehiclesresearchrelp   if tin(2007m8,2008m9), ///
			tscale(range(2007m8 2008m9)) ///
			tlabel(2007m8(3)2008m9) ///
			lwidth(medthick) ///
			lpattern(solid) lcolor(blue) ///
			ylabel(,nogrid tposition(outside) angle(horizontal))  ///
			legend(label(1 "Relative New Vehicle Price")  cols(1)) ///
			ytitle("Index, Dec. 2007 = 100")  ttitle("month") ///
			xline(`=ym(2008,5)', lcolor(red) lpattern(dash)) tlabel(, format(%tmMonCCYY) ) ///
			graphregion(color(white))

graph export ../output/newvehiclerelpresearchlongnotrend.eps, replace

* save output for later use
preserve
keep date newvehiclesresearchrelp
save ../output/motorvehicleprice.dta, replace
restore


// tsline newvehiclerelp newvehiclerelpfc if tin(`=ym(2008,1)',`=ym(2008,8)'), ///
// 			tscale(range(`=ym(2008,1)'(1)`=ym(2008,8)')) ///
// 			tlabel(`=ym(2008,1)'(3)`=ym(2008,8)') ///
// 			graphregion(color(white)) plotregion(style(none) margin(zero)) lwidth(medthick medthick) ///
// 			lpattern(solid solid) lcolor(black red) ///
// 			ylabel(,nogrid tposition(outside) angle(horizontal))  ///
// 			legend(label(1 "Relative New Vehicle Price") label(2 "Trend: January through April 2008")  cols(1) ring(0) position(7)) ///
// 			ytitle("Relative New Vehicle Price, Jan. 2008 = 100") yline(0, lcolor(black) lstyle(dot)) ttitle("Date") ///
// 			xline(`=ym(2008,5)', lcolor(red) lpattern(dash)) tlabel(, format(%tmMonCCYY) )

// graph export ../output/newvehiclerelp.eps, replace
