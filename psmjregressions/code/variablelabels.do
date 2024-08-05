version 16
cap ssc install did_imputation


foreach frequency in "interview" "monthly"     { //

    /*******************************************************************************

    Step 1: Set time and X-sectional variable

    *******************************************************************************/

    use ../input/psmjsample`frequency'.dta, clear
	
	
    *generates some aggregates of nipa variables 
    gen pce = nipa_pce 
    gen dpce = nipa_durable_nipa 
    gen mv = nipa_mv_parts_ser
    gen spce = nipa_service_nipa
    gen ndpce = nipa_nondurable_nipa
    gen mv_parts = nipa_mv_parts
    gen net_used  = nipa_net_used
    gen cars_nu = nipa_cars_nu
    
    
	*Financed car purchases
	gen car_finance = cars_nu - nipa_cars_downpayment

	*OTHER NIPA Categories
	gen other_nipa = pce -  mv
	
	*Parker Total 
	gen pmv = cartkn + cartku

    if "`frequency'"=="interview" {
	xtset cuid intdate
	foreach var of varlist pce dpce mv spce ndpce cars_nu mv_parts net_used car_finance other_nipa nipa* pmv{
		gen d_`var' = `var' -l3.`var' //First difference of NIPA variables
		gen lag_`var' = l3.`var' // Lag of NIPA Variables

	}
    }
    else if "`frequency'"=="monthly" {
        xtset cuid date
	foreach var of varlist pce dpce mv spce ndpce mv_parts cars_nu net_used car_finance other_nipa nipa* pmv{
		
		cap gen d_`var' = `var' -l.`var' //First difference of NIPA variables
		cap gen lag_`var' = l.`var' // Lag of NIPA Variables
		gen lag2_`var' = l2.`var' // Lag of NIPA Variables
		gen lag3_`var' = l3.`var' // Lag of NIPA Variables

	}
    }


    /*******************************************************************************

    Step 2: Labels for variables

    *******************************************************************************/

    // Label variables
    label variable rbtamt "Rebate Amount"
    label variable rbtindicator "Rebate Indicator"
    label variable rbt01amt "Rebate Amount"
    label variable rbt01indicator "Rebate Indicator"
	
	forvalues jj=1(1)2 {
		label variable lag`jj'rbtamt "Lag `jj' Rebate Amount"
		label variable lag`jj'rbtindicator "Lag `jj' Rebate Indicator"
		label variable lag`jj'rbt01amt "Lag `jj' Rebate Amount"
		label variable lag`jj'rbt01indicator "Lag `jj' Rebate Indicator"
	}
	forvalues jj=1(1)1 {
		label variable lead`jj'rbtamt "Lead `jj' Rebate Amount"
		label variable lead`jj'rbtindicator "Lead `jj' Rebate Indicator"
		label variable lead`jj'rbt01amt "Lead `jj' Rebate Amount"
		label variable lead`jj'rbt01indicator "Lead `jj' Rebate Indicator"
	}
	


    label variable foodbevs "FOOD"
    label variable sndexp "SND"
    label variable ndexp "ND"
    label variable totexp2 "EXP"
    label variable cartkn "CAR"
    label variable l_foodbevs "ln(FOOD)"
    label variable l_sndexp "ln(SND)"
    label variable l_ndexp "ln(ND)"
    label variable l_totexp2 "ln(EXP)"
    label variable pce  "CEX PCE Equivalent"
    label variable spce "CEX NIPA service equivalent"
    label variable dpce "CEX NIPA durable equivalent"
    label variable ndpce "CEX NIPA nondurable equivalent"
    
    

    label variable d_foodbevs "$\Delta$ FOOD"
    label variable d_sndexp "$\Delta$ SND"
    label variable d_ndexp "$\Delta$ ND"
     label variable d_dexp "$\Delta$ D"
      label variable d_sexp "$\Delta$ S"
    label variable d_totexp2 "$\Delta$ EXP"
    label variable d_cartkn "$\Delta$ New Vehicle (PSJM)"
    label variable d_newpcars "$\Delta$ New Car (PSJM)"
    label variable d_newtrucks "$\Delta$ New Truck (PSJM)"
     label variable d_cartku "$\Delta$ Used Vehicle (PSJM)"
    label variable d_l_foodbevs "$\Delta$ ln(FOOD)"
    label variable d_l_sndexp "$\Delta$ ln(SND)"
    label variable d_l_ndexp "$\Delta$ ln(ND)"
    label variable d_l_totexp2 "$\Delta$ ln(EXP)"
    label variable d_pce  " $\Delta$ PCE"
    label variable d_dpce "$\Delta$ Service"
    label variable d_spce "$\Delta$ Durable"
    label variable d_ndpce "$\Delta$ Nondurable"
    label variable d_net_used "$\Delta$ Net Used Vehicles"
    label variable d_cars_nu "$\Delta$ Vehicles (BEA)"
     label variable d_cars_n "$\Delta$ New Vehicle (BEA)"
     label variable d_nipa_pcars_n "$\Delta$ New Car (BEA)"
     label variable d_nipa_trucks_n "$\Delta$ New Truck (BEA)"
     label variable d_nipa_usedautos "$\Delta$ Net Used Car (BEA)"
     label variable d_nipa_usedlighttrucks "$\Delta$ Net Used Truck (BEA)"

    rename cat_eft_x cat_eft
    lab var cat_eft "CU Ever Received Rebate Via"
    lab def catvalues 0 "Non-Recipient" 1 "CHECK only recipient" 2 "EFT only recipient" 3 "Both" 9 "Recipient with missing EFT/Check information"
    lab val cat_eft catvalues
    
    lab var lincome "Log of Pre-Tax Income"
    
    *In correct income decile?
    bys incdpool: sum lincome
    

    
    /*******************************************************************************

    Step 3: Save

    *******************************************************************************/

    save ../output/psmjsample`frequency'_wlabels.dta, replace

}
