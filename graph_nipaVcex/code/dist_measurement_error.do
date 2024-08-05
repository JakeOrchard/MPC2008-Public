

use ../output/pce_cex_series.dta, replace
 
cap ssc install stripplot

/***************************************************************************
1. Create Variables.
*************************************************************************/

*Log variables
gen lpce_pop_sa = log(nipa_pce_pop_sa)
gen lservice_nipa_pop_sa = log(nipa_service_nipa_pop_sa)
gen  ldurable_nipa_pop_sa = log(nipa_durable_nipa_pop_sa)
gen  lnondurable_nipa_pop_sa = log( nipa_nondurable_nipa_pop_sa)


*Log Difference
gen diff_lpce = lPCE_month-lnipa_pce_pop_sa
gen diff_lservice = lPCE_service_month - lservice_nipa_pop_sa
gen diff_ldurable = lPCE_durable_month -  ldurable_nipa_pop_sa
gen diff_lnondurable = lPCE_nondurable_month -  lnondurable_nipa_pop_sa

*Billions difference
gen diff_pce = PCE_month-nipa_pce_pop_sa
gen diff_service = PCE_service_month - nipa_service_nipa_pop_sa
gen diff_durable = PCE_durable_month -  nipa_durable_nipa_pop_sa
gen diff_nondurable = PCE_nondurable_month -  nipa_nondurable_nipa_pop_sa


*
gen timetrend = date-420
tsset date

*Month fixed effects
gen rebatemonth = (date == ym(2008,5)) | (date == ym(2008,6)) | (date == ym(2008,7)) | (date == ym(2008,8))
gen four_before = (date == ym(2008,1)) | (date == ym(2008,2)) | (date == ym(2008,3)) | (date == ym(2008,4))
gen four_after = (date == ym(2008,9)) | (date == ym(2008,10)) | (date == ym(2008,11)) | (date == ym(2008,12))


/***************************************************************************
2. Remove linear time trend, extract FE
*************************************************************************/



*Regressions

local dependentvars "diff_lpce diff_lnondurable diff_ldurable diff_lservice diff_pce diff_nondurable diff_durable diff_service"
keep if date > 420

foreach depvar in `dependentvars' { 

	reg `depvar'  timetrend
	predict `depvar'_notrend , resid

}

/***************************************************************************
3. Graph Residuals
*************************************************************************/


*Label variables
lab var diff_lpce_notrend  "PCE"
lab var diff_lnondurable_notrend  "Non-durable"
lab var diff_ldurable_notrend  "Durable"
lab var diff_lservice_notrend  "Service"

lab var diff_pce_notrend  "PCE"
lab var diff_nondurable_notrend  "Non-durable"
lab var diff_durable_notrend  "Durable"
lab var diff_service_notrend  "Service"


*Why Remove linear time trend?

set scheme s1color
tsline diff_nondurable, ytitle("PCE-CEX: billions of $, monthly rate") 
graph export ../output/nondurable_difference.pdf, as(pdf) replace

tsline diff_lnondurable, ytitle("log(PCE/CEX)") 
graph export ../output/lnondurable_difference.pdf, as(pdf) replace

*Truncate diff_ldurable_notrend
sum diff_ldurable_notrend, detail

replace diff_ldurable_notrend = . if diff_ldurable_notrend > r(p95) | diff_ldurable_notrend < r(p5)



stripplot  diff_lservice_notrend   diff_ldurable_notrend  diff_lnondurable_notrend diff_lpce_notrend /*
*/ , stack width(.01)  separate(rebatemonth) legend(order(2 "May-Aug. 2008"))/*
*/ variablelabels xtitle("Log Difference: PCE v CEX")  text(2.8 -.15 "{it:CEX Higher}") /*
*/text(2.8 -.2 "{&lArr}", size(large)) text( 2.8 .15 "{it:PCE Higher}") text(2.8 .2 "{&rArr}", size(large))


graph export ../output/log_difference.pdf, as(pdf) replace

stripplot  diff_service_notrend   diff_durable_notrend  diff_nondurable_notrend diff_pce_notrend /*
*/ , stack width(1)  separate(rebatemonth) legend(order(2 "May-Aug. 2008"))/*
*/ variablelabels xtitle("PCE-CEX: billions of $, monthly rate") text(2.5 -30 "{it:CEX Higher}") /*
*/text(2.5 -40 "{&lArr}", size(large)) text( 2.5 30 "{it:PCE Higher}") text(2.5 40 "{&rArr}", size(large))

graph export ../output/billion_difference.pdf, as(pdf) replace


*2006-2011


stripplot  diff_lservice_notrend   diff_ldurable_notrend  diff_lnondurable_notrend diff_lpce_notrend /*
*/ if date > ym(2005,12) & date < ym(2011,1), stack width(.01)  separate(rebatemonth) legend(order(2 "May-Aug. 2008"))/*
*/ variablelabels xtitle("Log Difference: PCE v CEX")  text(2.8 -.15 "{it:CEX Higher}") /*
*/text(2.8 -.2 "{&lArr}", size(large)) text( 2.8 .15 "{it:PCE Higher}") text(2.8 .2 "{&rArr}", size(large))


graph export ../output/log_difference0611.pdf, as(pdf) replace

stripplot  diff_service_notrend   diff_durable_notrend  diff_nondurable_notrend diff_pce_notrend /*
*/ if date > ym(2005,12) & date < ym(2011,1), stack width(1)  separate(rebatemonth) legend(order(2 "May-Aug. 2008"))/*
*/ variablelabels xtitle("PCE-CEX: billions of $, monthly rate") text(2.5 -24 "{it:CEX Higher}") /*
*/text(2.5 -30 "{&lArr}", size(large)) text( 2.5 14 "{it:PCE Higher}") text(2.5 20 "{&rArr}", size(large))

graph export ../output/billion_difference0611.pdf, as(pdf) replace



/*



**Export Table**

local nregs = 8
local nregs1 = `nregs'+1
local keepvars `"rebatemonth four_before four_after"'
local order `"rebatemonth four_before four_after"'
local indicateca

local stype `" & \multicolumn{2}{c}{$ \log\left( \frac{PCE^{NIPA}}{PCE^{CEX}} \right) $}  & \multicolumn{2}{c}{$ \log\left( \frac{ND.^{NIPA}}{ND.^{CEX}} \right)$}  & \multicolumn{2}{c}{$ \log\left( \frac{Dur.^{NIPA}}{Dur.^{CEX}} \right)$}  & \multicolumn{2}{c}{$ \log\left( \frac{Sv.^{NIPA}}{Sv.^{CEX}} \right) $}"'
local lower `" `stype'   \\  "'
local groups `" "`groups1'   `lower' "  "'

local stats "timetrend N"
local stats_fmt " %3s %12.0fc"
local stats_label `" `"Time Trend"'   `"Observations"' "'
local num_stats: word count `stats'

local layout
forvalues l = 1/`num_stats' {
	local layout `"`layout' "\multicolumn{1}{c}{@}" "'
}
local dropvars 
local table_preamble `" "\begin{table}[!t] \caption{Consumption Spending NIPA v. CEX} \centering \sisetup{tableformat=1.2} \def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi}"    "\label{table:pceVcex}""\begin{tabularx}{\textwidth}{lXXXXXXXX}" "\\" "\hline\hline" "'
local prehead `"prehead(`table_preamble' `groups')"'			
local posthead `"posthead(`"\hline"' `"\multicolumn{`=`nregs'+1'}{l}{Right hand side variables:}\\"' `"\\"')"'

local notes `"Notes: NeweyWest HAC Standard errors are in parentheses. Significance at the 1, 5, and 10 percent levels indicated by ***,**, and *."'
local prefoot(" ")
local postfoot `"postfoot(`"\hline\hline \end{tabularx} \begin{minipage}{\hsize} \rule{0pt}{9pt} \tiny \centering `notes'  \end{minipage} \label{tab:`filename'} \end{table}"')"'
*This first one produces tables for the presentation (panel A and b are self contained)

esttab * using "../output/PCEvCEX.tex",  replace cells(b(fmt(%9.3f)) se(par star fmt(%9.3f) abs)) starlevels( \$^{*}$ 0.1 \$^{**}$ 0.05 \$^{***}$ 0.01  ) drop(`dropvars', relax) keep(`keepvars') indicate(`indicate') `prehead' `posthead' `postfoot' order(`order') label coeflabels( rebatemonth "MayAugust 2008" four_before "JanApril 2008" four_after "SeptemberDecember 2008"  ) stats(`stats', layout(`layout') fmt(`stats_fmt') labels(`stats_label')) collabels(,none) numbers nomtitles substitute(# `" X "' tabular* tabularx `"{1}{c}{("' `"{1}{L}{("') width(\hsize)



 
