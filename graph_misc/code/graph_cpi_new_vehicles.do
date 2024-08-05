
set scheme s1color


freduse CUSR0000SETA01, clear
gen mdate = mofd(daten)

format mdate %tm
rename CUSR0000SETA01 cpi_newcars

*Base Year (2007,1)
gen _cpi2007 = cpi_newcars if mdate == ym(2007,1)
egen cpi2007 = max(_cpi2007)
replace cpi_newcars =  100*(cpi_newcars/cpi2007)
label var cpi_newcars "CPI New Cars"

*Log CPI
gen lcpi_nc = log(cpi_newcars)
label var lcpi_nc "Log(CPI New Cars)"




tw (scatter lcpi_nc mdate if mdate < ym(2010,1) & mdate > ym(2006,12), c(l l) clp(l) clw(mdthick) xline(580, lp(-) lc(red)) clc(forest_green) mc(forest_green)) ///
    , name(lcpi_newcars, replace) scale(1.2) ///
    xlabel(, valuelabel)  yscale(titlegap(*10))


graph export ../output/mv_lCPI.pdf, as(pdf) replace
cap graph export ../output/mv_lCPI.png,  replace



tw (scatter cpi_newcars mdate if mdate < ym(2010,1) & mdate > ym(2006,12), c(l l) clp(l) clw(mdthick) xline(580, lp(-) lc(red)) clc(forest_green) mc(forest_green)) ///
    , name(cpi_newcars, replace) scale(1.2) ///
    xlabel(, valuelabel)  yscale(titlegap(*10))


graph export ../output/mv_CPI.pdf, as(pdf) replace
cap graph export ../output/mv_CPI.png,  replace

