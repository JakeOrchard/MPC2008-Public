


import delimited "../input/sca-tableall-on-2021-Dec-02.csv", varnames(2) rowrange(1) case(upper) clear




/*****************************************************************************
Labels and renames Variablesa
*****************************************************************************/


label var ICS_ALL "Michigan Index of Consumer Sentiment"
label var ICE_ALL "Michigan Index of Consumer Expectations"
label var ICC_ALL "Michigan Index of Economic Conditions"
label var PAGO_R "Current Personal Financial Situation (100 = Neutral)"
label var PAGORN_NY "Percent Higher - Lower Income"
label var RINC_R "Expected 1-2 Year Change in Real Income (100 = Neutral)"
label var NEWS_R "Favorable v. Unfavorable News (100 = Neutral) "
label var NEWSRN_F_GOVT "Favorable Government News"
label var NEWSRN_F_DEM "Favorable Consumer Demand News"
label var NEWSRN_U_GOVT "Unfavorable Government News"
label var NEWSRN_U_DEM "Unfavorable Consumer Demand News"
label var NEWSRN_NP "Net Price News"
label var NEWSRN_NG "Net Government News"
label var NEWSRN_NE "Net Employment News"
label var RATEX_R "Expected Change in Interest Rates (Down-UP + 100)"
label var PX1_MED "Median Inflation Expectations 12 Months"
label var PX1_MEAN "Mean Inflation Expectations 12 Months"
label var PX5_MED "Median Inflation Expectations 5 Years"
label var PX5_MEAN "Mean Inflation Expectations 5 Years"
label var GOVT_R "Current Gov. Economic Policy (100 = Neutral)"
label var DUR_R "Good time to buy Durables? (100 = Neutral)"

*Specific Vehicle Questions

label var VEH_R "Net Good-Bad Time to Buy a Vehicle (100 = Neutral)"
label var VEHRN_LP "Vehicle: Good Time Low Prices"
label var VEHRN_BIAP "Vehicle: Good time, price will rise"
label var VEHRN_LR "Vehicle: Good time, low interest rates"
label var VEHRN_BIAR "Vehicle: Good time, interest rates will rise"
label var VEHRN_GT "Vehicle: Good time, prosperity"
label var VEHRN_MPG "Vehicle: Good time, new low mpg"
label var VEHRN_HP "Vehicle: Bad time, high prices"
label var VEHRN_HR "Vehicle: Bad time, high interst rates"
label var VEHRN_TB "Vehicle: Bad time, times are bad/cannot afford"
label var VEHRN_FB "Vehcile: Bad time, uncertainty/bad times ahead"
label var VEHRN_GAS "Vehicle: Bad time, High gas prices"
label var VEHRN_SEL "Vehicle: Bad time, poor selction"




*Follow naming convention in Barsky and Sims AER 2012

rename BUS5_R E5Y
label var E5Y "5-Year Expected Economic Conditions (100 = Neutral)"

rename BUS12_R E12M
label var E12M "12-Month Expected Economic Conditions (100 = Neutral)"

rename PEXP_R PFE
label var PFE "12-Month Expected Personal Financial Situation (100 = Neutral)"




/*****************************************************************************
Saves
*****************************************************************************/

gen date = mofd(mdy(MONTH,1,YYYY))
tsset date

save ../output/clean_michigan.dta, replace

*Smaller sample of only Barsky and Sims used variables

keep date ICS_ALL ICE_ALL ICC_ALL E5Y E12M PFE PAGO_R VEH* E5Y E12M PFE PX1_MEAN PX1_MED PX5_MEAN PX5_MED RINC_R NEWS_R  NEWS* RATEX_R DUR_R

format date %tm

save ../output/clean_michigan_small.dta, replace


