clear all
set more off
version 16

cap cd ../psmjregressions/code


/*******************************************************************************

Step 0: Preliminaries: Formatting program and Locals for Regressions

*******************************************************************************/


*Formatting for tabulation
cap program drop format_lag
program format_lag ,eclass
	args b V N
	
    
    ereturn post b V 
    ereturn scalar N = N

   
end




local dependentvars_table "d_pce d_mv d_other_pce rbtamt" 

local dependentvars "rbtamt d_pce d_mv d_other_pce"

local absorbvars "intdate"

local controlvars "age d_num_adults d_perslt18"
local extravar "_Iincdinte* lag_pce lag_mv"

local versionlist `""  "lag" "hetconlag" "hetconlag_lc"' 
local exclude2rebate `"v1" ""'
local exclude2rebate `""'

local samplelist `"" "& everrbtindicator"'





/*******************************************************************************

Step 1: Pull in data and label variables

*******************************************************************************/
foreach vtype in "`exclude2rebate'"{ 
foreach rebatesample in "`samplelist'" {

foreach depvar in `dependentvars_table'{
	estimates clear

local regnumber = 1

foreach version in "`versionlist'"{
		

use ../output/psmjsampleinterview_wlabels.dta, clear


*Separate cohort variable for heterogenous TWFE specications
gen firstrbtintdate_esi = firstrbtintdate
gen rel_interview = intdate-firstrbtintdate
replace rel_interview = -9 if rel_interview < -7

egen numrbts = sum(rbtindicator), by(cuid)
tsset cuid intdate

	if "`vtype'" == "v1"{
		di "If households have two rebates, drop observations after first rebate"
		drop if numrbts >1 & rel_interview > 0
		gen _wrongint = rel_interview != 0 & rbtindicator == 1
		egen wrongint = max(_wrongint), by(cuid)
		replace firstrbtintdate_esi = intdate if wrongint & rbtindicator
		replace rel_interview = intdate-firstrbtintdate_esi
		local vtypename "singlerebatev1"		
		}
		
	else if "`vtype'" == ""{
		di "If households have two rebates, drop all household observations"
		drop if numrbts >1 
	}


keep if insample 


xtset cuid intdate

*Label Dependent Variables
label var d_pce "PCE"
label var lag_mv "Lag Motor Vehicle"
label var d_mv "Motor Vehicle and Parts Spending"
label var d_other_pce "Other Spending"
label var lag_pce "Lag Total Expenditure"
label var lag1rbtindicator "Lag Rebate Indicator"

local variablename "`: variable label `depvar''"

xi gen i.incdinterview 


/*******************************************************************************

Step 2: Regressions

*******************************************************************************/

		
		
		di "`version'"		
 				

	if "`version'" == ""{ // Column 1: TWFE with no lag
			
	local weighttype "OLS"
				
	*Estimates coefficients 
	noi di "reg  `depvar' rbtindicator  `controlvars'  i.intdate [w=finlwt21] if insample `rebatesample', cluster(cuid) "
	reg  `depvar' rbtindicator  `controlvars'  i.intdate [w=finlwt21] if insample `rebatesample', cluster(cuid)
			
	
	*Saves estimate	
	qui eststo col`regnumber'
	
	*Estimates MPC (3 month) 
	foreach sudep in `dependentvars' {
	
		reg  `sudep' rbtindicator  `controlvars'  i.intdate [w=finlwt21] if insample `rebatesample' 
	estimates store b_`version'`sudep'
	
	}
	
	suest b_`version'`depvar' b_`version'rbtamt b_`version'd_pce   b_`version'd_mv, vce(cluster cuid)
	if "`depvar'" != "rbtamt"{
		nlcom (_b[b_`version'`depvar'_mean: rbtindicator] )/(_b[b_`version'rbtamt_mean: rbtindicator] )
		local temp_mpc: di %3.2f r(b)[1,1]		
	}
	         }


        if "`version'" == "lag"{ // Column 2: TWFE with lag
	
	local weighttype "OLS"

				
        noi di "reg  `depvar' rbtindicator lag1rbtindicator `controlvars' i.intdate [w=finlwt21] if insample `rebatesample', cluster(cuid) "
	reg  `depvar' rbtindicator lag1rbtindicator  `controlvars'  i.intdate [w=finlwt21] if insample `rebatesample' , cluster(cuid)
			
	
	*Saves estimate
		qui eststo col`regnumber'

	
	
	*Estimates MPC (3 and 6 months) 
	foreach sudep in `dependentvars' {
	
		reg  `sudep' rbtindicator lag1rbtindicator  `controlvars'  i.intdate [w=finlwt21] if insample `rebatesample' 
	estimates store b_`version'`sudep'
	
	}
	
	

	suest b_`version'`depvar' b_`version'rbtamt b_`version'd_pce   b_`version'd_mv, vce(cluster cuid)
	if "`depvar'" != "rbtamt"{
		
		* 3-month MPC
		nlcom (_b[b_`version'`depvar'_mean: rbtindicator] )/(_b[b_`version'rbtamt_mean: rbtindicator] )
		local temp_mpc: di %3.2f r(b)[1,1]
	
	*6 month coefficient incorporates information from lagged spending

		nlcom (2*_b[b_`version'`depvar'_mean: rbtindicator] + _b[b_`version'`depvar'_mean: lag1rbtindicator])/(_b[b_`version'rbtamt_mean: rbtindicator] + _b[b_`version'rbtamt_mean: lag1rbtindicator])
		local temp_mpc6mo: di %3.2f r(b)[1,1]
		local temp_mpc6mo_se: di %3.2f sqrt(r(V)[1,1])

	}
	
        }
	

		
*Columns 3 and 4 TWFE with Heterogeneous Contemporaneous and Lagged effects
if "`version'" == "hetconlag" |  "`version'" == "hetconlag_lc"{

	if "`version'" == "hetconlag"{
		local add = ""
	}
	if "`version'" == "hetconlag_lc"{
		local add = "`extravar'"
	}

	 noi di "ivreghdfe  	`depvar' c.rbtindicator#firstrbtintdate_esi c.lag1rbtindicator#firstrbtintdate_esi `controlvars' `add' [weight=finlwt21] if insample `rebatesample' , cluster(cuid) absorb(`absorbvars' `hhfe')  "
        ivreghdfe  `depvar' c.rbtindicator#firstrbtintdate_esi c.lag1rbtindicator#firstrbtintdate_esi `controlvars' `add' [w=finlwt21] if insample `rebatesample' , cluster(cuid) absorb(`absorbvars' `hhfe')  
	
	*Save estimates prior to computing weights
	mat b = e(b)
	mat V = e(V)
	scalar N = e(N)

	*Compute weights
	if "`version'" == "hetconlag" | "`version'" == "hetconlag_lc"{
	local rebatevars "rbtindicator lag1rbtindicator"
		local weighttype "OLS"
		include compute_weights 
		local weight1 wlagOLS
		local weight2 wconOLS
	}
	
	format_lag b V N
		
	*Pull in OLS period-weights from graph weights
	frame create weights
	frame change weights
	use `weights', clear
	
	
	
	*Lag Coefficient
	sum `weight1' if intdate == 584 //Cohort lag = cohort date + 3
	scalar cw581 = _b[581b.firstrbtintdate_esi#c.lag1rbtindicator]*r(mean)
	local lagcoef = cw581
	scalar V581 = _se[581b.firstrbtintdate_esi#c.lag1rbtindicator]^2*r(mean)^2
	local lagVar =   V581 
	
	*Contemporaneous Coefficient
	sum `weight2' if intdate == 581
	scalar rbtw581 = _b[581b.firstrbtintdate_esi#c.rbtindicator]*r(mean)
	local coef = rbtw581
	scalar rbtV581 = _se[581b.firstrbtintdate_esi#c.rbtindicator]^2*r(mean)^2
	local rbtV =   rbtV581 
	
	
	qui forval intdatenum=582/590{

	*Lag Coef
	sum `weight1' if intdate == `intdatenum' + 3
	scalar cw`intdatenum' = _b[`intdatenum'.firstrbtintdate_esi#c.lag1rbtindicator]*r(mean)
	local addlag = cw`intdatenum'
	if missing(`addlag'){
		local addlag = 0 
	}
	local lagcoef = `addlag' + `lagcoef'
	
	*Lag S.E.
	scalar V`intdatenum' = _se[`intdatenum'.firstrbtintdate_esi#c.lag1rbtindicator]^2*r(mean)^2
	local addlagV = V`intdatenum'
	if missing(`addlagV'){
		local addlagV = 0 
	}
	local lagVar =   `addlagV' + `lagVar'
	*Contemporaneous Coef 
	sum `weight2' if intdate == `intdatenum'
	scalar rbtw`intdatenum' = _b[`intdatenum'.firstrbtintdate_esi#c.rbtindicator]*r(mean)
	local coef = rbtw`intdatenum' + `coef'
	scalar rbtV`intdatenum' = _se[`intdatenum'.firstrbtintdate_esi#c.rbtindicator]^2*r(mean)^2
	local rbtV =   rbtV`intdatenum' + `rbtV'
	}
	
	
	local covrbtind = 0
	local combinedV = `lagVar' + `rbtV'

	if "`rebatesample'" == ""{
		local endmonth = 590
		local endmonthshort = 587
		local begincon = 580
		local beginlag = 569
		
		
		
	}
	else{
		local endmonth = 590
		local endmonthshort = 587
		local begincon = 580
		local beginlag = 570	
		
	}
	
	qui forval ii=581/`endmonth'{ //Covariances
	
	qui forval jj=581/`endmonth'{
		
		local row1 = `ii' - `beginlag'
		local col1 = `jj' - `beginlag'
		local row2 = `ii' - `begincon'
		local col2 = `jj' - `begincon'
		
		scalar cv`ii'`jj' = e(V)[`row1',`col1']
		
		scalar crbtv`ii'`jj' = e(V)[`row2',`col2']
		
		
		*Lagged weights 
		
		
		local rintdate_corrected = `ii' + 3
		local cintdate_corrected = `jj' + 3
		sum `weight1' if intdate == `rintdate_corrected' 
		local wrow`ii'_1 = r(mean)
		sum `weight1' if intdate == `cintdate_corrected' 
		local wcol1 = r(mean)
		
		*Contemporaneous Weights
		sum `weight2' if intdate == `ii'
		local wrow`ii'_2 = r(mean)
		sum `weight2' if intdate == `jj'
		local wcol2 = r(mean)
		
		if `ii' ~= `jj'{
			local rbtV = `rbtV' +  `wcol2'*`wrow`ii'_2'*crbtv`ii'`jj'
			if `ii' <= `endmonthshort'  & `jj' <= `endmonthshort' {
				local lagVar = `lagVar' + `wcol1'*`wrow`ii'_1'*cv`ii'`jj'
				local combinedV = `combinedV' + `wcol1'*`wrow`ii'_1'*cv`ii'`jj' +  `wcol2'*`wrow`ii'_2'*crbtv`ii'`jj'
			}
		}
		
		}
	}
	
	
	frame change default
	frame drop weights
	
	
	*Post new coefficients
	mat b = e(b)
	scalar lagrbtind = `lagcoef'
	scalar rbtind = `coef'
	if "`rebatesample'" == ""{
		mat b[1,12] = lagrbtind
	}
	else{
		mat b[1,11] = lagrbtind
	}
	
	mat b[1,1] = rbtind
	local colnames: colnames b
	local colnames : subinstr local colnames "581b.firstrbtintdate_esi#c.lag1rbtindicator" "lag1rbtindicator"
	local colnames : subinstr local colnames "581b.firstrbtintdate_esi#c.rbtindicator" "rbtindicator"
	mat colnames b = `colnames'
	
	
	mat V = e(V)
	scalar lagrbtindV = `lagVar'
	scalar rbtindV = `rbtV'
	*scalar covrbtind = `covrbtind'
	if "`rebatesample'" == ""{
		mat V[12,12] = lagrbtindV
	}
	else{
		mat V[11,11] = lagrbtindV
	}
	mat V[1,1] =  rbtindV
	local colnames: colnames V
	local colnames : subinstr local colnames "581b.firstrbtintdate_esi#c.lag1rbtindicator" "lag1rbtindicator"
	local colnames : subinstr local colnames "581b.firstrbtintdate_esi#c.rbtindicator" "rbtindicator"	
	mat colnames V = `colnames'
	local rownames: rownames V
	local rownames : subinstr local rownames "581b.firstrbtintdate_esi#c.lag1rbtindicator" "lag1rbtindicator"
	local rownames : subinstr local rownames "581b.firstrbtintdate_esi#c.rbtindicator" "rbtindicator"	
	mat rownames V = `rownames'
	format_lag b V N

	*Saves estimate 
	
	qui eststo col`regnumber'
	
	*Estimate MPCs (3- and 6-month)
	*Estimate same reg model for other dependent variables and estimate S.E. via Delta method
	foreach sudep in `dependentvars' {
		reg  `sudep' c.rbtindicator#firstrbtintdate_esi c.lag1rbtindicator#firstrbtintdate_esi  `controlvars' `add' i.intdate [w=finlwt21] if insample `rebatesample' 
	estimates store b_`version'`sudep'
	}
		suest b_`version'`depvar' b_`version'rbtamt b_`version'd_pce   b_`version'd_mv, vce(cluster cuid)
	

	foreach sudep in `dependentvars' {
	local len = length("`sudep'")
	local level = substr("`sudep'",3,`len'-2)
	local lagv = "lag_" + "`level'"	
	

	forval ii = 581/`endmonth'{
		if `ii' == 581{
			local wcoef_lag`ii'_`sudep' = "`wrow581_1'*_b[b_`version'`sudep'_mean: 581b.firstrbtintdate_esi#c.lag1rbtindicator]"
			local wcoef_con`ii'_`sudep' = "`wrow581_2'*_b[b_`version'`sudep'_mean: 581b.firstrbtintdate_esi#c.rbtindicator]"
			local wcoef_con_lags`ii'_`sudep' = "`wrow581_2'*_b[b_`version'`sudep'_mean: 581b.firstrbtintdate_esi#c.rbtindicator]*_b[b_`version'`depvar'_mean: `lagv']"
			local wcoef_clr`ii'_`sudep' = "`wrow581_2'*_b[b_`version'`sudep'_mean: 581b.firstrbtintdate_esi#c.rbtindicator]*_b[b_`version'rbtamt_mean: `lagv']"

			}
			else{
				local wcoef_lag`ii'_`sudep' = "`wrow`ii'_1'*_b[b_`version'`sudep'_mean: `ii'.firstrbtintdate_esi#c.lag1rbtindicator]"
			local wcoef_con`ii'_`sudep' = "`wrow`ii'_2'*_b[b_`version'`sudep'_mean: `ii'.firstrbtintdate_esi#c.rbtindicator]"
			local wcoef_con_lags`ii'_`sudep' = "`wrow`ii'_2'*_b[b_`version'`sudep'_mean: `ii'.firstrbtintdate_esi#c.rbtindicator]*_b[b_`version'`depvar'_mean: `lagv']"
			local wcoef_clr`ii'_`sudep' = "`wrow`ii'_2'*_b[b_`version'`sudep'_mean: `ii'.firstrbtintdate_esi#c.rbtindicator]*_b[b_`version'rbtamt_mean: `lagv']"
			}
	}
		local hetcoefs_lag_`sudep' = "`wcoef_lag581_`sudep'' + `wcoef_lag582_`sudep'' + `wcoef_lag583_`sudep'' + `wcoef_lag584_`sudep'' + `wcoef_lag585_`sudep'' + `wcoef_lag586_`sudep'' + `wcoef_lag587_`sudep''"	
		local hetcoefs_con_`sudep' = "`wcoef_con581_`sudep'' + `wcoef_con582_`sudep'' + `wcoef_con583_`sudep'' + `wcoef_con584_`sudep'' + `wcoef_con585_`sudep'' + `wcoef_con586_`sudep'' + `wcoef_con587_`sudep'' + `wcoef_con588_`sudep'' + `wcoef_con589_`sudep''+  `wcoef_con590_`sudep''"	
		local hetcoefs_con_lags_`sudep' = "`wcoef_con_lags581_`sudep'' + `wcoef_con_lags582_`sudep'' + `wcoef_con_lags583_`sudep'' + `wcoef_con_lags584_`sudep'' + `wcoef_con_lags585_`sudep'' + `wcoef_con_lags586_`sudep'' + `wcoef_con_lags587_`sudep'' + `wcoef_con_lags588_`sudep'' + `wcoef_con_lags589_`sudep''+ `wcoef_con_lags590_`sudep''"	
		local hetcoefs_clr`sudep' = "`wcoef_clr581_`sudep'' + `wcoef_clr582_`sudep'' + `wcoef_clr583_`sudep'' + `wcoef_clr584_`sudep'' + `wcoef_clr585_`sudep'' + `wcoef_clr586_`sudep'' + `wcoef_clr587_`sudep'' + `wcoef_clr588_`sudep'' + `wcoef_clr589_`sudep''+ `wcoef_clr590_`sudep''"
		
	}
	
		if  "`depvar'" ~= "rbtamt"{
		* 3-month MPC
		nlcom ("`hetcoefs_con_`depvar''" )
		local coef = r(b)[1,1]
		nlcom ("`hetcoefs_con_rbtamt'" )
		local rebate = r(b)[1,1]
		local temp_mpc: di %3.2f `coef'/`rebate'
		
		*6 month coefficient incorporates information from lagged spending
		if "`version'" == "hetconlag_lc"{
		local hetcoefs = "`hetcoefs_con_`depvar'' + `hetcoefs_con_`depvar''+ `hetcoefs_lag_`depvar''  + `hetcoefs_con_lags_d_pce' + `hetcoefs_con_lags_d_mv'"
		nlcom  `hetcoefs'
	
		local coef6mo = r(b)[1,1]
		local coef6moV = r(V)[1,1]
		local rbtcoefs = "`hetcoefs_con_rbtamt'  + `hetcoefs_lag_rbtamt' + `hetcoefs_clrd_pce' + `hetcoefs_clrd_mv'"
		nlcom `rbtcoefs'
		
		local rebate6mo = r(b)[1,1]
		}
		else if "`version'" == "hetconlag"{
		local hetcoefs = "`hetcoefs_con_`depvar'' + `hetcoefs_con_`depvar''+ `hetcoefs_lag_`depvar''"
		nlcom  `hetcoefs'
		
		local coef6mo = r(b)[1,1]
		local coef6moV = r(V)[1,1]
		local rbtcoefs = "`hetcoefs_con_rbtamt'  + `hetcoefs_lag_rbtamt'" 
		nlcom `rbtcoefs'
		local rebate6mo = r(b)[1,1]
		}
		local temp_mpc6mo: di %3.2f `coef6mo'/`rebate6mo'
		local temp_mpc6mo_se: di %3.2f sqrt(((1/`rebate6mo')^2)*`coef6moV')
		}
	

		}
		
		
		
			

		
		
		
		
		
	
	if "`rebatesample'" == ""{
	local rebatename ""
	
	}
	else{
		
	local rebatename "_rbtonly"
			}
	
	
	
		 
		
		
		*Add implied MPC 
		
		

		qui estadd local implied_mpc "`temp_mpc'" : col`regnumber'
		if "`depvar'" ~= "rbtamt" & ("`version'" == "lag" | "`version'" == "hetconlag" | "`version'" =="hetconlag_lc"){
				qui estadd local implied_mpc6mo "`temp_mpc6mo'" : col`regnumber'
				qui estadd local implied_mpc6mo_se "(`temp_mpc6mo_se')" : col`regnumber'
				
		}
		if "`version'" == ""{
 			qui estadd local implied_mpc6mo "" : col`regnumber'
				qui estadd local implied_mpc6mo_se "" : col`regnumber'
				
		}
		
		*Add line about income decile FE
		if "`version'" =="hetconlag_lc"{
			qui estadd local decilefe `"Yes"' : col`regnumber'
		}
		else{
			qui estadd local decilefe `"No"' : col`regnumber'
		}
		
		local regnumber = `regnumber' + 1

	}
	
	
	estimates drop b_*
		
	*Labels
	if "`depvar'" == "d_pce"{
		local newlabel "Personal Consumption Expenditure (PCE)"
		local tablename `"Table_E1_`depvar'"'
		local filename `"Table_E1_`depvar'`vtype'`rebatename'"'
		local filename2 `"Table_E1pres_`depvar'`vtype'`rebatename'"'


	}
	
	if "`depvar'" == "rbtamt"{
		local tablename `"Table_AE1"'
		local filename `"Table_AE1`vtype'`rebatename'"'
		local filename2 `"Table_AE1pres`vtype'`rebatename'"'

	}
	if "`depvar'" == "d_mv"{
		local newlabel "motor vehicles and parts expenditure"
		local tablename `"Table_AE2"'
		local filename `"Table_AE2`vtype'`rebatename'"'
		local filename2 `"Table_AE2pres`vtype'`rebatename'"'

	}
	if "`depvar'" == "d_other_pce"{
		local newlabel "PCE excluding motor vehicles and parts"
		local tablename `"Table_AE3"'
		local filename `"Table_AE3`vtype'`rebatename'"'
		local filename2 `"Table_AE3pres`vtype'`rebatename'"'

	}

	* this code produces the regression table
	local nregs = 4
	
	local stats_label `" `"Implied 3-month MPC"' `"Implied 6-month MPC"' `"6-Month MPC S.E."' `"Income Decile FE"' `"Observations"' "'
	
	
	local varorder "rbtindicator lag1rbtindicator lag_pce lag_mv"
	local order `"`varorder'"' 
	local midrules1 `" \cmidrule(l{.75em}){2-3} \cmidrule(l{.75em}){4-5}"'
	local groups2 `" &   \multicolumn{2}{c}{ Homogeneous Treatment} &  \multicolumn{2}{c}{Heterogeneous Treatment} \\ "'  	 
	local groups `" "`groups2'" "`midrules1'"  "'

	local stats " implied_mpc implied_mpc6mo implied_mpc6mo_se decilefe N"
	
	local stats_fmt " %3.2f %3s  %12.0fc"
	
	
	if "`depvar'" == "rbtamt" | "`depvar'" == "lag_pce"{
	
		local stats "decilefe N"
		local stats_fmt " %3s  %12.0fc"
		local stats_label `"`"Income Decile FE"' `"Observations"' "'
	}
	local num_stats: word count `stats' 
	local layout
	forvalues l = 1/`num_stats' {
		local layout `"`layout' "\multicolumn{1}{c}{@}" "'
	}
	local keepvars `"`varorder'"' 
	local dropvars 
	
	local title `"Household `variablename'  Response to Rebate `extratitle'"'
	
	local notes `"Notes: The dependent variable is the change in `newlabel'.    Regressions include interview (time) fixed effects, and household level controls for age, change in number of adults, and change in number of children. Standard errors for the 6-month MPC are estimated via Delta-method.  The rebate coefficients in  columns (3) and (4) are the weighted average of the interaction between rebate cohort and the (lagged) rebate indicator with weights computed following Sun and Abraham (2021). Standard errors, in parentheses, are clustered at the household level: $ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 $."'
	

	if "`depvar'" == "rbtamt"{
		local title `" First Stage: Rebate Amount Conditional on Rebate Receipt `extratitle'"'
	
		local notes `"Notes: The dependent variable is the dollar value of Econoimic Stimulus Payments (ESP) received by the household.  Standard errors, in parentheses, are clustered at the household level. Regressions include interview (time) fixed effects, and household level controls for age, change in number of adults, and change in number of children. The rebate coefficients in  columns (3) and (4) are the weighted average of the interaction between rebate cohort and the (lagged) rebate indicator with weights computed following Sun and Abraham (2021). Significance is indicated by: $ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 $. "'
	
	}
	
	
	
	if "`rebatesample'" == ""{

	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=1.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"  "\small \caption{`title'}" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}"  "\toprule     \multicolumn{`nregs'}{l}{\textbf{Panel A: Full Sample}} \\ \midrule"  "'
	local postfoot `"postfoot(`"\hline \end{tabularx} "')"'
	}
	else if "`rebatesample'" ~= ""{
		local table_preamble `" " \vspace{2mm} \\ \begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}" "\toprule     \multicolumn{`nregs'}{l}{\textbf{Panel B: Rebate Recipients Only}} \\ \midrule"  "'
		
		
			local postfoot `"postfoot(`"\hline\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \footnotesize `notes'  \end{minipage} \label{tab:`tablename'} \end{table}"')"'
	}
	local prehead `"prehead(`table_preamble' `groups')"'			
	local posthead `"posthead(`"\hline"' `"\\"')"'
	
	local prefoot(" ")
	


	esttab * using "../output/`filename'.tex", replace cells(b(star fmt(%5.2f)) se(par fmt(%5.2f) abs)) starlevels(\$^{*}$ 0.1 \$^{**}$ 0.05 \$^{***}$ 0.01)  ///
		keep(`keepvars') `prehead' `posthead' `postfoot' order(`order') label varlabel(`varlabels') stats(`stats', layout(`layout') fmt(`stats_fmt') ////
		labels(`stats_label')) collabels(,none) numbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize) eqlabels("" "")
		
		*Presentation (Separate tables)
		
	if "`rebatesample'" == ""{

	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=1.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"  "\small \caption{`title'}" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}"  "\toprule     \multicolumn{`nregs'}{l}{\textbf{Full Sample}} \\ \midrule"  "'
	}
	else if "`rebatesample'" ~= ""{
		local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=1.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"  "\scriptsize" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}"  "\toprule     \multicolumn{`nregs'}{l}{\textbf{Rebate Only Sample}} \\ \midrule"  "'
	}
		
local postfoot `"postfoot(`"\hline\hline \end{tabularx} \end{table}"')"'
	
	local prehead `"prehead(`table_preamble' `groups')"'			
	local posthead `"posthead(`"\hline"' `"\\"')"'
	
	local prefoot(" ")
	


	esttab * using "../output/`filename2'.tex", replace cells(b(star fmt(%5.2f)) se(par fmt(%5.2f) abs)) starlevels(\$^{*}$ 0.1 \$^{**}$ 0.05 \$^{***}$ 0.01)  keep(`keepvars') `prehead' `posthead' `postfoot' order(`order') label varlabel(`varlabels') stats(`stats', layout(`layout') fmt(`stats_fmt') labels(`stats_label')) collabels(,none) numbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize) eqlabels("" "")
		
		
	
	eststo clear
	
	
	}
	
	
}

}	
	
	
	
	
	
	
	
	
	
	
