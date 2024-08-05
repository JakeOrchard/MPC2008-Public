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




local dependentvars_table "lag_pce" 

local absorbvars "intdate"

local controlvars "age d_num_adults d_perslt18"
local extravar "_Iincdinte* lag_pce lag_mv"

local versionlist `"lag"' 
local exclude2rebate `"v1" ""'
local exclude2rebate `""'

local samplelist `"" "& everrbtindicator"'





/*******************************************************************************

Step 1: Pull in data and label variables

*******************************************************************************/

foreach depvar in `dependentvars_table'{
		estimates clear
local regnumber = 1
	foreach rebatesample in "`samplelist'" {



foreach version in "`versionlist'"{
		

use ../output/psmjsampleinterview_wlabels.dta, clear


*Separate cohort variable for heterogenous TWFE specications
gen firstrbtintdate_esi = firstrbtintdate
gen rel_interview = intdate-firstrbtintdate
replace rel_interview = -9 if rel_interview < -7

egen numrbts = sum(rbtindicator), by(cuid)
tsset cuid intdate

	
		di "If households have two rebates, drop all household observations"
		drop if numrbts >1 
	


keep if insample 


xtset cuid intdate

*Label Dependent Variables
label var d_pce "PCE"
label var lag_mv "Lag Motor Vehicle"
label var lag_pce "PCE"
local variablename "`: variable label `depvar''"
label var rbtindicator "Lead Rebate Indicator"
label var lag1rbtindicator "Rebate Indicator"

xi gen i.incdinterview 


/*******************************************************************************

Step 2: Regressions

*******************************************************************************/

		
		
		di "`version'"		
 				


        if "`version'" == "lag"{ // Column 2: TWFE with lag
	
	local weighttype "OLS"

				
        noi di "reg  `depvar' rbtindicator lag1rbtindicator `controlvars' i.intdate [w=finlwt21] if insample `rebatesample', cluster(cuid)"
	reg  `depvar' rbtindicator lag1rbtindicator  `controlvars'  i.intdate [w=finlwt21] if insample `rebatesample', cluster(cuid)
	
			
	*Saves estimate
		qui eststo col`regnumber'
	
	
	
        }
	
	 	
	*Add implied MPC 
		
		

		

	}
	local regnumber = `regnumber' + 1
}
		
	
	if "`depvar'" == "lag_pce"{
		local newlabel "Lag Personal Consumption Expenditure (PCE)"
		local tablename `"Table_bias"'
		local filename `"Table_bias"'


	}
	

	* this code produces the regression table
	local nregs = 2
	
	local stats_label `"`"Observations"' "'
	
	
	local varorder "rbtindicator lag1rbtindicator"
	local order `"`varorder'"' 
	local groups2 `" &   \multicolumn{1}{c}{ Full Sample} &  \multicolumn{1}{c}{Rebate Recipients Only} \\ "'  
	local groups `" "`groups2'" "`midrules1'"  "'	
	local stats "N"
	
	local stats_fmt " %12.0fc"
	
	
	local num_stats: word count `stats' 
	local layout
	forvalues l = 1/`num_stats' {
		local layout `"`layout' "\multicolumn{1}{c}{@}" "'
	}
	local keepvars `"`varorder'"' 
	local dropvars 
	
	local title `"Negative effect of future rebate receipt on current expenditure"'
	
	local notes `"Notes: The dependent variable is the Level of PCE.    Regressions include interview (time) fixed effects, and household level controls for age, change in number of adults, and change in number of children.Standard errors, in parentheses, are clustered at the household level: $ \:^{*}\:p<0.1,\:\:^{**}\:p<0.05,\:\:^{***}\:p<0.01 $."'
		
	

	local table_preamble `" "\begin{table}[!t] \centering \sisetup{table-format=1.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"  "\small \caption{`title'}" "\begin{tabularx}{\hsize}{@{\hskip\tabcolsep\extracolsep\fill}l*{`nregs'}{S}} \midrule""'
	
			local postfoot `"postfoot(`"\hline\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \footnotesize `notes'  \end{minipage} \label{tab:`tablename'} \end{table}"')"'
	
	local prehead `"prehead(`table_preamble' `groups')"'			
	local posthead `"posthead(`"\hline"' `"\\"')"'
	
	local prefoot(" ")
	


	esttab * using "../output/`filename'.tex", rename(tau0 rbtindicator tau3 lag1rbtindicator 0b.lag1rbtindicator#581b.firstrbtintdate_esi lag1rbtindicator) replace cells(b(star fmt(a2)) se(par fmt(a2) abs)) starlevels(\$^{*}$ 0.1 \$^{**}$ 0.05 \$^{***}$ 0.01) drop(`dropvars', relax) ///
		keep(`keepvars') `prehead' `posthead' `postfoot' order(`order') label varlabel(`varlabels') stats(`stats', layout(`layout') fmt(`stats_fmt') ////
		labels(`stats_label')) collabels(,none) numbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize) eqlabels("" "")
		
		
	
	eststo clear
	
	
	}
	
	


	
	
	
	
	
	
	
	
	
	
	
