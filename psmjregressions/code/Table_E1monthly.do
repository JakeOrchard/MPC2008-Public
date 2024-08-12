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




local dependentvars_table "d_mv" 

local dependentvars "rbtamt d_pce d_mv d_other_pce"

local absorbvars "intdate cuid"

local controlvars "age d_num_adults d_perslt18"
local extravar "_Iincdinte* lag_pce lag_mv"

local versionlist `"hetconlag_lc"' 
local exclude2rebate `"v1" ""'
local exclude2rebate `""'

local samplelist `"" "& everrbtindicator"'





/*******************************************************************************

Step 1: Pull in data and label variables

*******************************************************************************/
foreach vtype in "`exclude2rebate'"{ 

foreach depvar in `dependentvars_table'{
	estimates clear

local regnumber = 1
foreach rebatesample in "`samplelist'" {

foreach version in "`versionlist'"{
	

use ../output/psmjsamplemonthly_wlabels.dta, clear

*Separate cohort variable for heterogenous TWFE specications
gen firstrbtdate_esi = firstrbtdate
gen rel_month = date-firstrbtdate
replace rel_month = -12 if rel_month < -11

egen numrbts = sum(rbtindicator), by(cuid)
xtset cuid date

	
		
	if "`vtype'" == ""{
		di "If households have two rebates, drop all household observations"
		drop if numrbts >1 
	}


keep if insample 



*Label Dependent Variables
label var pce "PCE"
label var lag_mv "Lag Motor Vehicle"
label var d_mv "Motor Vehicle and Parts Spending"
label var d_other_pce "Other Spending"
label var lag_pce "Lag Total Expenditure"
label var lag3rbtindicator "Lag 3 Rebate Indicator"

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
	         }


        if "`version'" == "lag"{ // Column 2: TWFE with lag
	
	local weighttype "OLS"

				
        noi di "reg  `depvar' rbtindicator lead1rbtindicator lag1rbtindicator lag2rbtindicator lag3rbtindicator `controlvars' i.intdate [w=finlwt21] if insample `rebatesample', cluster(cuid) "
	reg  `depvar' rbtindicator lead1rbtindicator lag1rbtindicator lag2rbtindicator lag3rbtindicator `controlvars'  i.intdate [w=finlwt21] if insample `rebatesample' , cluster(cuid)
			
	
	*Saves estimate
		qui eststo col`regnumber'
	
        }
	

		
*Columns 3 and 4 TWFE with Heterogeneous Contemporaneous and Lagged effects
if "`version'" == "hetconlag" |  "`version'" == "hetconlag_lc"{

	if "`version'" == "hetconlag"{
		local add = ""
	}
	if "`version'" == "hetconlag_lc"{
		local add = "`extravar'"
	}

	 noi di "ivreghdfe  	`depvar' c.rbtindicator#firstrbtdate_esi c.lag1rbtindicator#firstrbtdate_esi  c.lag2rbtindicator#firstrbtdate_esi c.lag3rbtindicator#firstrbtdate_esi c.lead1rbtindicator#firstrbtdate_esi `controlvars' `add' [weight=finlwt21] if insample `rebatesample' , cluster(cuid) absorb(`absorbvars' `hhfe')  "
        ivreghdfe  `depvar'  c.rbtindicator#firstrbtdate_esi c.lag1rbtindicator#firstrbtdate_esi  c.lag2rbtindicator#firstrbtdate_esi c.lag3rbtindicator#firstrbtdate_esi c.lead1rbtindicator#firstrbtdate_esi  `controlvars' `add' [w=finlwt21] if insample `rebatesample' , cluster(cuid) absorb(`absorbvars' `hhfe')  
	
	*Save estimates prior to computing weights
	mat b = e(b)
	mat V = e(V)
	scalar N = e(N)

	*Compute weights
	if "`version'" == "hetconlag" | "`version'" == "hetconlag_lc"{
	local rebatevars "lead1rbtindicator rbtindicator lag1rbtindicator lag2rbtindicator lag3rbtindicator"
		local weighttype "OLS"
		include compute_weights_monthly 
		local weightp1 wOLSlead1
		local weight0 wOLS0
		local weight1 wOLS1
		local weight2 wOLS2
	}
	
	format_lag b V N
		
	*Pull in OLS period-weights from graph weights
	frame create weights
	frame change weights
	use `weights', clear
	
	*Beginning of each rebate lags entries
	if "`rebatesample'" == ""{
		local endmonth = 589
		local beginp1 = 530 //49
		local begin0 = 578 //1
		local begin1 = 566 //13
		local begin2 = 554 //25
		local begin3 = 542 //37
	}
	else{
		local endmonth = 589
		local beginp1 = 534 //45
		local begin0 = 578 //1
		local begin1 = 567 //12
		local begin2 = 556 //23
		local begin3 = 545 //34

	}
	
	local h = -1
	foreach rbt in  "lead1rbtindicator" "rbtindicator" "lag1rbtindicator" "lag2rbtindicator" "lag3rbtindicator"{
		local hname "`h'"
		if "`rbt'" == "lead1rbtindicator"{
			local hname "p1"
		}
		if `h' > 0{
			local endmonthshort = 586
		}
		else{
			local endmonthshort = `endmonth'
			
		}
	
	sum `weight`hname'' if date == 579 + `h' //Cohort lag = cohort date + 3
	scalar cw579 = _b[579b.firstrbtdate_esi#c.`rbt']*r(mean)
	local coef`hname' = cw579
	
	scalar V579 = _se[579b.firstrbtdate_esi#c.`rbt']^2*r(mean)^2
	local V`hname' =   V579 
	
	
	qui forval datenum=580/`endmonth'{
	*Coef	
	if `datenum' <= `endmonthshort'{
		sum `weight`hname'' if date == `datenum' + `h'
		scalar cw`datenum' = _b[`datenum'.firstrbtdate_esi#c.`rbt']*r(mean)
		local addcoef = cw`datenum'
		local coef`hname' = `addcoef' + `coef`hname''

		*S.E.
		scalar V`datenum' = _se[`datenum'.firstrbtdate_esi#c.lag1rbtindicator]^2*r(mean)^2
		local addV = V`datenum'
		
		local V`hname' =   `addV' + `V`hname''
	}
	}
	

	
	
	qui forval ii=579/`endmonth'{ //Covariances
	
	qui forval jj=579/`endmonth'{
		
		local row1 = `ii' - `begin`hname''
		local col1 = `jj' - `begin`hname''
		scalar cv`ii'`jj' = e(V)[`row1',`col1']		
		
		*Weights 
		local rdate_corrected = `ii' + `h'
		local cdate_corrected = `jj' + `h'
		sum `weight`h'' if date == `rdate_corrected' 
		local wrow`ii'_1 = r(mean)
		sum `weight`h'' if date == `cdate_corrected' 
		local wcol1 = r(mean)
		
		if `ii' ~= `jj'{
			if `ii' <= `endmonthshort'  & `jj' <= `endmonthshort' {
				local V`hname' = `V`hname'' +  `wcol1'*`wrow`ii'_1'*cv`ii'`jj'
			}
		}
		
		}
	}
	local h = `h' + 1
	}
	
		

	frame change default
	frame drop weights
	
	
	*Post new coefficients
	mat b = e(b)
	scalar lead1rbtind = `coefp1'
	scalar rbtind = `coef0'
	scalar lag1rbtind = `coef1'
	scalar lag2rbtind = `coef2'
	scalar lag3rbtind = `coef3'

	if "`rebatesample'" == ""{
		mat b[1,49] = lead1rbtind
		mat b[1,13] = lag1rbtind
		mat b[1,25] = lag2rbtind
		mat b[1,37] = lag3rbtind
	}
	else{
		mat b[1,45] = lead1rbtind
		mat b[1,12] = lag1rbtind
		mat b[1,23] = lag2rbtind
		mat b[1,34] = lag3rbtind

	}
	
	mat b[1,1] = rbtind
	local colnames: colnames b
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lead1rbtindicator" "lead1rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.rbtindicator" "rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lag1rbtindicator" "lag1rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lag2rbtindicator" "lag2rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lag3rbtindicator" "lag3rbtindicator"

	mat colnames b = `colnames'
	
	
	mat V = e(V)
	scalar lead1rbtindV = `Vp1'
	scalar rbtindV = `V0'
	scalar lag1rbtindV = `V1'
	scalar lag2rbtindV = `V2'
	scalar lag3rbtindV = `V3'

	if "`rebatesample'" == ""{
		mat V[49,49] = lead1rbtindV
		mat V[13,13] = lag1rbtindV
		mat V[25,25] = lag2rbtindV
		mat V[37,37] = lag3rbtindV
	}
	else{
		mat V[45,45] = lead1rbtindV
		mat V[12,12] = lag1rbtindV
		mat V[23,23] = lag2rbtindV
		mat V[34,34] = lag3rbtindV
	}
	mat V[1,1] =  rbtindV
	local colnames: colnames V
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lead1rbtindicator" "lead1rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.rbtindicator" "rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lag1rbtindicator" "lag1rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lag2rbtindicator" "lag2rbtindicator"
	local colnames : subinstr local colnames "579b.firstrbtdate_esi#c.lag3rbtindicator" "lag3rbtindicator"

	mat colnames V = `colnames'
	local rownames: rownames V
	local rownames : subinstr local rownames "579b.firstrbtdate_esi#c.lead1rbtindicator" "lead1rbtindicator"
	local rownames : subinstr local rownames "579b.firstrbtdate_esi#c.rbtindicator" "rbtindicator"
	local rownames : subinstr local rownames "579b.firstrbtdate_esi#c.lag1rbtindicator" "lag1rbtindicator"
	local rownames : subinstr local rownames "579b.firstrbtdate_esi#c.lag2rbtindicator" "lag2rbtindicator"
	local rownames : subinstr local rownames "579b.firstrbtdate_esi#c.lag3rbtindicator" "lag3rbtindicator"

	mat rownames V = `rownames'
	format_lag b V N

	*Saves estimate 
	
	qui eststo col`regnumber'
	
	

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
}	
	
		
	*Labels
	
	
	if "`depvar'" == "d_mv"{
		local newlabel "motor vehicles and parts expenditure"
		local tablename `"Table_AEMonthly"'
		local filename `"Table_AEMonthly"'
	}
	
	* this code produces the regression table
	local nregs = 2
	
	local stats_label `"  `"Income Decile FE"' `"Observations"' "'
	
	
	local varorder "lead1rbtindicator rbtindicator lag1rbtindicator lag2rbtindicator lag3rbtindicator lag_pce lag_mv"
	local order `"`varorder'"' 
	local midrules1 `" \cmidrule(l{.75em}){2-2} \cmidrule(l{.75em}){3-3}"'
	local groups1 `"\toprule &  \multicolumn{1}{c}{Full Sample} &  \multicolumn{1}{c}{Rebate Only Sample} \\ "' 

	local groups `" "`groups1'" "`midrules1'"  "'

	local stats "decilefe N"
	
	local stats_fmt " %3s  %12.0fc"
	
	
	if "`depvar'" == "rbtamt" | "`depvar'" == "lag_pce"{
	
		local stats "decilefe N"
		local stats_fmt " %12.0fc"
		local stats_label `"`"Income Decile FE"' `"Observations"' "'
	}
	local num_stats: word count `stats' 
	local layout
	forvalues l = 1/`num_stats' {
		local layout `"`layout' "\multicolumn{1}{c}{@}" "'
	}
	local keepvars `"`varorder'"' 
	local dropvars 
	
	local title `"Monthly: Household `variablename' Response to Rebate `extratitle'"'
	
	local notes `"Notes: The dependent variable is the change in `newlabel'.    Regressions include month fixed effects, and household level controls for age, change in number of adults, and change in number of children. The rebate coefficients are the weighted average of the interaction between rebate month cohort and the rebate indicators with weights computed following Sun and Abraham (2021). Standard errors, in parentheses, are clustered at the household level: $ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 $."'
	

	if "`depvar'" == "rbtamt"{
		local title `" First Stage: Rebate Amount Conditional on Rebate Receipt `extratitle'"'
	
		local notes `"Notes: The dependent variable is the dollar value of Econoimic Stimulus Payments (ESP) received by the household.  Standard errors, in parentheses, are clustered at the household level. Regressions include interview (time) fixed effects, and household level controls for age, change in number of adults, and change in number of children. The rebate coefficients in  columns (3) and (4) are the weighted average of the interaction between rebate cohort and the (lagged) rebate indicator with weights computed following Sun and Abraham (2021). Significance is indicated by: $ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 $. "'
	
	}
	
	
	

	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=1.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"  "\small \caption{`title'}" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}}" "'
	local postfoot `"postfoot(`"\hline \end{tabularx} "')"'
	
			local postfoot `"postfoot(`"\hline\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \footnotesize `notes'  \end{minipage} \label{tab:`tablename'} \end{table}"')"'
	
	local prehead `"prehead(`table_preamble' `groups')"'			
	local posthead `"posthead(`"\hline"' `"\\"')"'
	
	local prefoot(" ")
	


	esttab * using "../output/`filename'.tex", replace cells(b(star fmt(%5.2f)) se(par fmt(%5.2f) abs)) starlevels(\$^{*}$ 0.1 \$^{**}$ 0.05 \$^{***}$ 0.01)  ///
		keep(`keepvars') `prehead' `posthead' `postfoot' order(`order') label varlabel(`varlabels') stats(`stats', layout(`layout') fmt(`stats_fmt') ////
		labels(`stats_label')) collabels(,none) numbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize) eqlabels("" "")
		
		
	
	eststo clear
	
	
	}
	
	
}

	
	
	
	
	
	
	
	
	
	
	
