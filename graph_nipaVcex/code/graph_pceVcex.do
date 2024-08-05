clear all
set more off




/******************************************************************************
*1. Pull in NIPA data from FRED and convert to per-capita
******************************************************************************/
set fredkey b2d9da12107485ec086c5691d41d626c

import fred PCE PCEDG PCESV POPTHM PCEND  DFXARC1Q027SBEA DCLORC1Q027SBEA DGOERC1Q027SBEA DMOTRC1Q027SBEA DFDHRC1Q027SBEA DGDSRG3M086SBEA PCEPI, clear


rename DFXARC1Q027SBEA PCE_foodbev_home
rename DCLORC1Q027SBEA PCE_apparel
rename DGOERC1Q027SBEA PCE_gasoline
rename DMOTRC1Q027SBEA PCE_mv_parts
rename DFDHRC1Q027SBEA PCE_home_dur
rename DGDSRG3M086SBEA pcons_good
rename PCEDG PCE_durable
rename PCEND PCE_nondurable
rename PCESV PCE_service

*Convert to monthly expenditure

foreach var of varlist PCE PCE_durable PCE_service PCE_nondurable PCE_foodbev_home PCE_apparel PCE_gasoline PCE_mv_parts PCE_home_dur{

	gen `var'_month = `var'/12
	gen l`var'_month = log(`var'_month)

}

gen PCEG_month = PCE_durable_month + PCE_nondurable_month

gen date = mofd(daten)

/******************************************************************************
2. Merges with CEX and graphs
*****************************************************************************/

merge 1:1 date using ../output/sa_cex.dta

*Convert CEX data from per-capita to population expenditure (Billions of Dollars)

foreach var of varlist nipa*_pc {
	
	local newvar  "`=regexr("`var'","_pc$","")'"
	disp "`var':`newvar'"
	
	gen `newvar'_pop = POPTHM * `var'/10^6
	gen `newvar'_pop_sa = POPTHM * `var'_sa/10^6
	gen l`newvar'_pop_sa = ln(`newvar'_pop_sa)
	drop `var' `newvar'
}


tsset date
format date %tm

set scheme s1color

local timelist `"" "& date > ym(2005,12) & date < ym(2011,1)"'

foreach time in "`timelist'"{
	
	
	if "`time'" == "& date > ym(2005,12) & date < ym(2011,1)"{
		local timename "06to11"
		local timelabel "552(6)612"
	}
	else{
		local timename ""
		local timelabel "430(24)720"
	}
	*Durables
	tsline PCE_durable_month nipa_durable_nipa_pop nipa_durable_nipa_pop_sa if _merge == 3 `time' ,         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title(DURABLES)

	graph export ../output/durables_cexvnipa`timename'.pdf, as(pdf) replace

	*MOTOR VEHICLES AND PARTS
	tsline PCE_mv_parts_month nipa_mv_parts_pop nipa_mv_parts_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("MOTOR VEHICLES AND PARTS")

	graph export ../output/mvparts_cexvnipa`timename'.pdf, as(pdf) replace

	*FURNITURE AND OTHER HOUSEHOLD GOODS
	tsline PCE_home_dur_month nipa_home_dur_pop nipa_home_dur_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("HOUSEHOLD DURABLES INCLUDING FURNITURE")

	graph export ../output/furn_cexvnipa`timename'.pdf, as(pdf) replace



	*NONDURABLES
	tsline PCE_nondurable_month nipa_nondurable_nipa_pop nipa_nondurable_nipa_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("NON-DURABLES")

	graph export ../output/nondurables_cexvnipa`timename'.pdf, as(pdf) replace
	*FOOD and BEVERAGES

	tsline PCE_foodbev_home_month nipa_foodbev_home_pop nipa_foodbev_home_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("AT HOME FOOD AND BEVERAGES")
	graph export ../output/food_cexvnipa`timename'.pdf, as(pdf) replace

	*apparel

	tsline PCE_apparel_month nipa_clothing_pop nipa_clothing_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("CLOTHING AND FOOTWEAR")
	graph export ../output/clothing_cexvnipa`timename'.pdf, as(pdf) replace

	*GAS and Other fuel

	tsline PCE_gasoline_month nipa_gasoline_pop nipa_gasoline_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("GASOLINE AND OTHER FUEL GOODS")
	graph export ../output/gas_cexvnipa`timename'.pdf, as(pdf) replace

	*Services
	tsline PCE_service_month nipa_service_nipa_pop nipa_service_nipa_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("SERVICES")

	graph export ../output/services_cexvnipa`timename'.pdf, as(pdf) replace

	*PCE
	tsline PCE_month nipa_pce_pop nipa_pce_pop_sa if _merge == 3 `time',         ///
	legend(label(1 "PCE")  label(2 "CEX") label(3 "CEX Seasonally Adjusted")) ///                                        ///                                     ///
	xtitle("") ytitle("billions of $, monthly rate")           ///
	tlabel("`timelabel'", format(%tmCCYY)) title("PERSONAL CONSUMPTION EXPENDITURES")
	graph export ../output/pce_cexvnipa`timename'.pdf, as(pdf) replace

}






/********************************************************************************
*Graphs for November 2021 presentation
******************************************************************************/



format date %tm


*Goods
gen rPCEG = 100*PCEG_month/pcons_good
gen rpceg = 100*nipa_total_good_pop_sa/pcons_good

label var rpceg "Real CEX Goods"
label var rPCEG "Real PCE for goods"

tw (scatter rpceg date if date>m(2007m1) & date<=m(2009m12), c(l l) clp(l) clw(medthick) xline(`=ym(2008,5)', lp(-) lc(red)) ///
    ytitle("billions of $, monthly rate", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter rPCEG date if date>m(2007m1) & date<=m(2009m12), c(l l) clp(l) clw(medthick) yaxis(2) clc(blue) mc(blue))


graph export ../output/pres_cexvnipa_goods.pdf, as(pdf) replace


*Total PCE 

gen rPCE = 100*PCE_month/PCEPI
gen rpce = 100*nipa_pce_pop_sa/PCEPI

label var rpce "CEX Real``PCE'' Series"
label var rPCE "Real PCE"

tw (scatter rpce date if date>m(2007m1) & date<=m(2009m12), c(l l) clp(l) clw(medthick) xline(`=ym(2008,5)', lp(-) lc(red)) ///
    ytitle("billions of $, monthly rate", size(medsmall)) clc(purple) mc(purple))  ///
    (scatter rPCE date if date>m(2007m1) & date<=m(2009m12), c(l l) clp(l) clw(medthick) yaxis(2) clc(blue) mc(blue))


graph export ../output/pres_cexvnipa_total.pdf, as(pdf) replace


*********************************************************************************
* Save data
*********************************************************************************

gen mdate = date
save ../output/pce_cex_series.dta, replace

