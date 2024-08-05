clear all
set more off



use ../input/nipavariablesmonthly.dta, clear


*******************************************************************************
*1. Collapses data and creates per-capita measure
*******************************************************************************


*drops values with missing weight or family size information
drop if fam_size == . | finlwt21 == .



collapse (sum) fam_size nipa* rbtamt [w=finlwt21], by(date)

*Total good expenditure 
gen nipa_total_good = nipa_durable_nipa + nipa_nondurable_nipa

foreach var of varlist nipa* {
	
	gen `var'_pc = `var'/fam_size
	
}





sort date


*Drops last two months since not fully covered by interviews
drop if date >= ym(2019,11)

*******************************************************************************
*2. Seasonally adjusts data and creates per-capita measure
*******************************************************************************_

gen monthnum = month(dofm(date))


foreach var of varlist *_pc{
	
	reg `var' i.monthnum
	predict `var'_sa, resid
	replace `var'_sa =`var'_sa + _b[_cons] 
	forval i = 2/12{

	replace `var'_sa = _b[`i'.monthnum]/12 + `var'_sa
	
	
}

}

foreach var of varlist *_pc*{

	gen l`var' = log(`var')

}


save ../output/sa_cex.dta, replace
