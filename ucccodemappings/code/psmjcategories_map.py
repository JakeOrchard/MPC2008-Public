import yaml

# dictionary with mappings
psmjcat = {}

psmjcat['FOODBEVS'] =   {                   # (1) FOOD+BEV
                            'FDHOME',       # Food consumed away from home ---- FDAWAY
                            'FDAWAY',       # Food consumed at home ----------- FDHOME
                            'ALCBEV',       # Purchases of alcholic beverages - ALCBEV
                        }

psmjcat['SNDEXP']   =   {                   # (2) Strictly nondurable edxpenditures
                            'FOODBEVS',     # food+bev ------------------- foodbev
                            'TOBACC',       # Tobaccoo --------------- tobacc
                            'UTIL',         # utilities -------------- util
                            'PERSCA',       # personal care ---------- persca
                            'HOUSOP',       # household operations --- housop
                            'PUBTRA',       # public transportation -- pubtra
                            'GASMO',        # gas and motor oil ------ gasmo
                            'MISC',         # miscellaneous expenses - misc
                        }   

psmjcat['NDEXP']    =   {                   # (3) Nondurable expenditures
                            'SNDEXP',       # snd --------------- snd
                            'APPAR',        # apparel ----------- appar
                            'HEALTH',       # health ------------ health
                            'READ',         # reading materials - read
                        }   

psmjcat['USED_CAR_PCE'] =  { #Net Used Cars according to BLS staff CEX->PCE crosswalk

                        'usedautos', 
                        'usedlighttrucks',
}                                                     

psmjcat['DEXP'] =  {                   # (4) Durable expenditures
                            'CARTKN',       # vehicle purchases (new) -- cartkn
                            'CARTKU',       # vehicle purchases (used) -- cartku
                            'OTHVEH',       # other vehicle (used) -- othveh
                            'HOUSEQ',       # furniture + equipment ---- houseq
                        } 

psmjcat['SEXP'] =  {                        # (4) Service Expenditure 
                            'EDUCA',        # education ----------- educa
                            'MAINRP',       # maintenance ----------- mainrp
                            'VEHINS',       # insurance ----------- vehins 
                            'ENTERT',       # entertainment-------- entert
                            'VRNTLO',      # vehicle lease (used) -- vrntlo
                            'VEHFIN',       # vehicle finance charges -- vehfin
                            'SHELT',        # shelter -------------- shelt
                        } 

                        
psmjcat['TOTEXP2']  =   {                   # (4) Total expenditures
                            'NDEXP',        # ndexp --------------- ndexp
                            'EDUCA',        # education ----------- educa
                            'HOUSEQ',       # furniture + equipment ---- houseq
                            'SHELT',        # shelter -------------- shelt
                            'CARTKN',       # vehicle purchases (new) -- cartkn
                            'CARTKU',       # vehicle purchases (used) -- cartku
                            'VRNTLO',       # vehicle lease (used) -- vrntlo
                            'VEHFIN',       # vehicle finance charges -- vehfin
                            'OTHVEH',       # other vehicle (used) -- othveh
                            'MAINRP',       # maintenance ----------- mainrp
                            'VEHINS',       # insurance ----------- vehins 
                            'ENTERT',       # entertainment-------- entert
                        } 

psmjcat['TOTEXP_NMV']  =   {                   # (4) Total expenditures minus new and used cars
                            'NDEXP',        # ndexp --------------- ndexp
                            'EDUCA',        # education ----------- educa
                            'HOUSEQ',       # furniture + equipment ---- houseq
                            'SHELT',        # shelter -------------- shelt
                            'VRNTLO',       # vehicle lease (used) -- vrntlo
                            'VEHFIN',       # vehicle finance charges -- vehfin
                            'OTHVEH',       # other vehicle (used) -- othveh
                            'MAINRP',       # maintenance ----------- mainrp
                            'VEHINS',       # insurance ----------- vehins 
                            'ENTERT',       # entertainment-------- entert
                        } 

#psmjcat['finance_psmj']     =   {                   # (4) Total expenditures minus new and used cars   
 #                           'VEHFIN',       # vehicle finance charges -- vehfin
  #                          'FINCHRG',
   #                         'MORTGAGEINT',
   #                     }      Can reconstruct in STATA              


# save output to yaml file
# need sort_keys=False to preserve ordering since construction depends on previous
# entries
with open('../output/psmjcategories_map.yml', 'w') as yamlfile:
    yaml.dump(psmjcat, yamlfile, default_flow_style=False, sort_keys=False) #, f, default_flow_style=False                    