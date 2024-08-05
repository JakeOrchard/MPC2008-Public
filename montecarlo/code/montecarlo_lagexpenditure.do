*Monte Carlo with Lagged Expenditure control

// keep if intdate>=ym(2008,6) & intdate<=ym(2008,11)
// gen counter = 1
// collapse (max) rbtind (sum) counter, by(cuid)
// keep if counter==2
// su rbtind

set seed 3213799

clear all

global  nreps = 10000
global samplesize_hh = 8000
global samplesize_T  = 4
global burnin = 6

*Sets values of the parameters
global rebate_coefficient = 300

global rebateonly = 0

//0.69 report, 85\% of households should have gotten it


scalar rho = 0.68 //From CE data in 2008
scalar mu = 0 //3629
scalar sigma_x = 7689.687
scalar sigma_fe = 0
scalar gammamax = 0.00006
scalar gammadelta = 42

// scalar sigma_dx = 8117

global versionlist "A" "B" "C"


*Degree of measurement error
scalar frac_rebate = 0.85

*Test for significance
scalar alpha = 0.05

cap program drop _all
program define simulate_expenditure, rclass
	clear
	*Create Samplesize_hh households for $samplesize_T time periods
	set obs $samplesize_hh
	gen cuid = _n 
	gen x_fe = rnormal(0,sigma_fe)
	expand  $samplesize_T + $burnin
	bys cuid: gen t = _n
	xtset cuid t 
	

	*Rebate dummies last 85 percent are treated in three waves
	gen rebate = 0
	replace rebate = 1 if t == $burnin + 2 & (cuid > 0 & cuid < 0.5*frac_rebate*$samplesize_hh)
	replace rebate = 1 if t == $burnin + 3 & (cuid > 0.5*frac_rebate*$samplesize_hh & cuid < 0.83*frac_rebate*$samplesize_hh)
	replace rebate = 1 if t == $burnin + 4 & (cuid > 0.83*frac_rebate*$samplesize_hh & cuid < frac_rebate*$samplesize_hh)
	
	if $rebateonly == 1 {
		drop if cuid >= frac_rebate*$samplesize_hh
	}
	egen everrebate = max(rebate), by(cuid)
	
	*Expenditure process, starts 6 periods prior to available dataset 
	gen epsilon_x = rnormal(0,sigma_x)
	gen lx = rnormal(mu/(1-rho),sigma_x/sqrt(1-rho^2)) if t == 1
	replace lx = mu + rho*L.lx  +  epsilon_x if t != 1 
	
	keep if t > $burnin
// 	gen x = exp(lx) + x_fe
	gen x = lx + x_fe
	
	*Rebate Effect
	gen x_without_rebate = x 
	gen dx_without_rebate = S.x 
	replace x = x + $rebate_coefficient * rebate
	
	keep if t > $burnin
	
	gen _rebatedate = t if rebate 
	egen rebatedate = max(_rebatedate), by(cuid)
	
	* difference
	gen dx = S.x 
	

	
	****** No Measurement Error ********
	reg dx rebate l.rebate, robust
	return scalar b_1 = _b[rebate]

	reg dx rebate l.rebate l.x, robust
	return scalar b_2 = _b[rebate]

	
	****** Correlated Measurement Error - Linear FF ********
	

	
	*Reported rebate with measurement error (correlated with absolute change in expenditure)
	local sigma_u = 1 
	gen u = rnormal(0 , `sigma_u')
	
	foreach version in "$versionlist"{
		
		if "`version'" == "A"{
			scalar type_i = 0.00000001
			scalar type_ii = 1-0.7/frac_rebate //From CE Survey portion of recipients in summer 2008
			
		}
		if "`version'" == "B"{
			scalar type_i = 0.05
			scalar type_ii = 0.05
		}
		if "`version'" == "C"{
			scalar type_i = 0.15 //From Rebate interview schedule table panel A) 5 pp. avg over reporting divded by 33 
			scalar type_ii = 0.0000001
		}
		
	
		
		local i = 1
		forvalues i = 1(1)`=gammadelta' {

	// 		local gamma = 0.0001 * `i'
			local gamma = (`i' - 1) / (gammadelta - 1) * gammamax
			
			gen rcl_latent = rebate * (invnormal(1-type_ii)) - `gamma'*(L.x_without_rebate) + u 

			replace rcl_latent = 0 if t == 7


			if "`version'" == "A"{
				*False Negatives 
				gen rr_cl = 0
				_pctile rcl_latent if rebate == 1, p(`=type_ii*100')
				replace rr_cl = rcl_latent  > r(r1)  if rebate == 1
				}
				
			if "`version'" == "B"{
				*False Negatives 
				gen rr_cl = 0
				_pctile rcl_latent if rebate == 1, p(`=type_ii*100')
				replace rr_cl = rcl_latent  > r(r1)  if rebate == 1
				*False Positives
				_pctile rcl_latent if rebate == 0 & t>7, p(`=100-type_i*100')
				replace rr_cl = rcl_latent  > r(r1)  if rebate == 0  & t>7
				}
		
			if "`version'" == "C"{
				* False Positives -> Misreport timing by recipients only
				gen rr_cl = rebate 
				_pctile rcl_latent if rebate == 0 & t> 7 & t>rebatedate & everrebate , p(`=100-type_i*100')
				replace rr_cl = rcl_latent  > r(r1)  if rebate == 0  & t>7 & t>rebatedate & everrebate
				egen sumrr_cl = sum(rr_cl),by(cuid)
				replace rr_cl = 0 if sumrr_cl > 1 & rebate == 1 //Shifts rebate report 
				
			}
				
			
			
			sum x if rebate & rr_cl == 0
			return scalar x_noreport`i'`version' = r(mean)
			sum x if rebate & rr_cl == 1
			return scalar x_report`i'`version' = r(mean)
			
				
			su rr_cl if L.x_without_rebate<.
			return scalar rr_cl`i'`version' = r(mean)
			su rebate if L.x_without_rebate<.
			return scalar rebate`i'`version' = r(mean)
			
			*Bias regression
			reg x f.rr_cl rr_cl, robust
			return scalar b_forward_lin`i'`version' = _b[f.rr_cl]
			
			reg dx rr_cl l.rr_cl, robust
			return scalar b_clin_nolag`i'`version' = _b[rr_cl]
		
			reg dx rr_cl l.rr_cl l.x, robust
			return scalar b_clin_lag`i'`version' = _b[rr_cl]
			
			return scalar gamma`i'`version' = `gamma'
			
			drop rcl_latent rr_cl 
			cap drop sumrr_cl
			
		}
	}
end


local simulstr = ""
foreach version in "$versionlist"{
	forvalues i = 1(1)`=gammadelta' {
		local simulstr = "`simulstr' b_forward_lin`i'`version' = r(b_forward_lin`i'`version') b_clin_nolag`i'`version' = r(b_clin_nolag`i'`version') b_clin_lag`i'`version' = r(b_clin_lag`i'`version') gammas`i'`version' = r(gamma`i'`version') rr_cl`i'`version' = r(rr_cl`i'`version') rebate`i'`version' = r(rebate`i'`version') x_noreport`i'`version' = r(x_noreport`i'`version') x_report`i'`version' = r(x_report`i'`version')"
	}

}

simulate b_1 = r(b_1)  b_2 = r(b_2) `simulstr' , reps( $nreps ) : simulate_expenditure
gen simulid = _n

save ../output/simulation_data.dta, replace
use ../output/simulation_data.dta, clear

*Distribution of b_1 
hist b_1, kdensity xtitle("Beta-Baseline") xlabel(-100(200)700) xline(300, lwidth(thick) lcolor(red))   color(blue%30)
	graph export ../output/hist_baseline.png, replace

qui foreach version in "$versionlist"{

hist b_clin_nolag1`version' , kdensity xtitle("Beta Model `version'", size(huge)) xlabel(-100(200)700, labsize(huge)) xline(300, lwidth(thick) lcolor(red))   color(blue%30) ylabel(0(0.001)0.003, labsize(huge)) ytitle("Density", size(huge))
	graph export ../output/hist_baseline`version'.png, replace
}

reshape long b_forward_lin b_clin_nolag b_clin_lag rr_cl rebate gammas x_report x_noreport, i(b_1 b_2 simulid) j(inum) string
gen version = substr(inum,-1,. )

collapse (mean) b* rr_cl rebate x_report x_noreport, by(gammas version)



foreach var of varlist b_1 b_2 b_clin* {
	
	gen bias_`var' = `var' - $rebate_coefficient
}


// label var bias_b_clog_lag "Lagged Expenditure"
// label var bias_b_clog_nolag "No Lagged Expenditure"
// label var b_bias "Consumption Response to Lead of Reported Rebate "
label var bias_b_clin_lag "Lagged Expenditure"
label var bias_b_clin_nolag "No Lagged Expenditure"
label var b_forward_lin "zeta"
label var x_report "Reports Rebate"
label var x_noreport "Does Not Report Rebate"

qui foreach version in "$versionlist"{

	if "`version'" == "A"{
			local longversion "Model A-False Negatives"
			local longversion2 "Model A-False Negatives"

			local extralegend ""
		}
		if "`version'" == "B"{
			local longversion "Model B-False Positives and Negatives"
			local longversion1 "Model B-False Positives"
			local longversion2 "Model B-False Negatives"
			local extralegend "legend(off)"
		}
		if "`version'" == "C"{
			local longversion "Model C-False Timing"
			local longversion1 "Model C-False Timing"

			local extralegend "legend(off)"
		}
		di "`longversion'"

	preserve 
	keep if version == "`version'"


	tw (connected  bias_b_clin_lag gamma) (connected bias_b_clin_nolag gamma) (connected b_forward_lin gamma)  ,   ytitle(Bias) xtitle(Gamma) title("`longversion'") legend(pos(5)) `extralegend' name(linear_gamma`version')
	graph export ../output/linear_gamma`version'.png, replace


	// label var bias_b_clog_lag "Logarithmic: Lagged Expenditure"
	// label var bias_b_clog_nolag "Logarithmic: No Lagged Expenditure"
	label var bias_b_clin_lag "Lagged Expenditure Control"
	label var bias_b_clin_nolag "No Control"

	tw (line bias_b_clin_nolag b_forward_lin if b_forward_lin > -1000 & b_forward_lin < 0, lwidth(1) xline(-575, lwidth(thick)))  (line  bias_b_clin_lag b_forward_lin  if b_forward_lin > -1000 & b_forward_lin < 0, lwidth(1) lpattern(longdash)), ytitle(Bias of Rebate Coefficient--Beta) xtitle(Coefficient on Lead of Rebate--Zeta)  legend(  cols(1) ring(0) position(1)) ylabel(, nogrid) yline(0, lwidth(0.25) lpattern(dot) lcolor(black)) graphregion(color(white))  name(all_ff`version') `extralegend' xlabel(0(-200)-800) 
	graph export ../output/all_ff`version'.eps, replace



	tw (line bias_b_clin_nolag gamma, lwidth(1)) , ytitle(Bias, size(huge)) xtitle(gamma, size(huge))  legend(  cols(1) ring(0) position(10)) ylabel(, nogrid) yline(0, lwidth(0.25) lpattern(dot) lcolor(black)) graphregion(color(white))  name(bias_gamma`version') `extralegend' ylabel(0(500)1500,labsize(huge)) xlabel(0(0.00002)0.00006, labsize(huge))

	graph export ../output/bias_gamma`version'.eps, replace

	tw (line b_forward_lin gamma, lwidth(1)), ytitle(zeta, size(huge)) xtitle(gamma, size(huge))  legend(  cols(1) ring(0) position(1)) ylabel(, nogrid) yline(0, lwidth(0.25) lpattern(dot) lcolor(black)) graphregion(color(white)) name(forward_gamma`version') `extralegend' ylabel(0(-1000)-3000, labsize(huge))

	graph export ../output/forward_gamma`version'.eps, replace
	
	if "`version'" == "C"{
		tw (line x_noreport gamma if gamma < 0.00004, lwidth(1)) (line x_report gamma if gamma < 0.00004, lwidth(1)), ytitle(Expenditure) xtitle(gamma)  legend(  cols(1) ring(0) position(3)) ylabel(, nogrid) yline(0, lwidth(0.25) lpattern(dot) lcolor(black)) graphregion(color(white)) title("`longversion'") name(x_report`version') 
	graph export ../output/x_report`version'.eps, replace

	}


	gen joinvar = 1

	gen disttotarget = abs(b_forward_lin + 574.9)
	su disttotarget
	replace disttotarget = disttotarget - r(min)

	tempfile biases
	save `biases'

	// 	su L.x_without_rebate, detail
	clear
	local plotpoints = 49
	set obs `plotpoints'
	gen lagx = .
	forvalues jj = 0(1)`plotpoints' {
		replace lagx = (4 * `jj' / (`=(`plotpoints'+1)') - 2)  * sigma_x/sqrt(1-rho^2) if _n == `jj'
	}
	replace lagx = lagx 
	label var lagx "Lagged Expenditure Rel. to Mean"
	gen exgl = (rho - 1) * lagx
	label var exgl "Expected Change in Expenditure"
	gen joinvar = 1


	joinby joinvar using `biases'

	drop joinvar

	sort gammas lagx
	
	if "`version'" == "A"{
		local type_ii = 1-0.7/frac_rebate
		gen probreport_cond_rebate = normal( invnormal(1-`type_ii') - gammas* lagx )
	}
	if "`version'" == "B"{
		local type_i = 0.05 
		local type_ii = 0.05
		gen probreport_cond_rebate = normal( invnormal(1-`type_ii')  - gammas* lagx )
		gen probreport_cond_norebate = normal(  - gammas* lagx - invnormal(1-`type_i') )
	}
	if "`version'" == "C"{
		local type_i =0.13
		local type_ii = 0.0
		gen probreport_cond_norebate = normal(  - gammas* lagx - invnormal(1-`type_i') )
	}
	
	if "`version'" != "C"{
	su b_forward_lin if gamma == 0
	local gamma0forward = r(mean)
	noi di "`version'"
	noi su gammas if disttotarget == 0
	local gammaminforward = round(r(mean) * 10^6 ) / 10^6
	disp "`gammaminforward'"

	tw (line probreport_cond_rebate exgl if gammas==0, lwidth(1) lpattern(longdash)) (line probreport_cond_rebate exgl if disttotarget == 0, lwidth(1) ),  ytitle(Prob(Report Rebate Conditional on Receipt))  legend(order(1 "gamma = 0" 2 "gamma = `gammaminforward'")  cols(1) ring(0) position(11)) ylabel(, nogrid)    graphregion(color(white)) xlabel(-6000(3000)6000) title("`longversion2'") name(probability_of_reporting`version') 

	graph export ../output/probability_of_reporting`version'.eps, replace
	}
	if "`version'" == "B" | "`version'" == "C"{
		su b_forward_lin if gamma == 0
		local gamma0forward = r(mean)
		noi di "`version'"
		noi su gammas if disttotarget == 0
		local gammaminforward = round(r(mean) * 10^6 ) / 10^6
		disp "`gammaminforward'"

		tw (line probreport_cond_norebate exgl if gammas==0, lwidth(1) lpattern(longdash)) (line probreport_cond_norebate exgl if disttotarget == 0, lwidth(1) ),  ytitle(Prob(Report Rebate Conditional on No Rebate Receipt))  legend(order(1 "gamma = 0" 2 "gamma = `gammaminforward'")  cols(1) ring(0) position(11)) ylabel(, nogrid)    graphregion(color(white)) xlabel(-6000(3000)6000) title("`longversion1'") name(probability_of_reporting_nr`version') 

		graph export ../output/probability_of_reporting_norebate`version'.eps, replace
		
	}
	
	restore

}


/**************************************************************************
Combined Graphs
**************************************************************************/

graph combine linear_gammaA linear_gammaB linear_gammaC, ycommon xcommon cols(3)
graph export ../output/linear_gamma_combined.eps, replace
graph combine all_ffA all_ffB all_ffC, ycommon xcommon cols(3)
graph export ../output/all_ff_combined.eps, replace

graph combine bias_gammaA bias_gammaB bias_gammaC, ycommon xcommon cols(3)
graph export ../output/bias_gamma_combined.eps, replace

graph combine forward_gammaA forward_gammaB forward_gammaC, ycommon xcommon cols(3)
graph export ../output/forward_gamma_combined.eps, replace


