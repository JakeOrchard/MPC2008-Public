import yaml

# dictionary with mappings
level3to2 = {}

level3to2['ALCBEVER']   =   {                   # Total alcohol expenditure  
                                'ALCHOME',      # Alcohol at hom
                                'ALCAWAY',      # Alcohol away
                            }


level3to2['UTILS']      =   {                   # Total utilities expenditure  
                                'ENERGY',       # Energy
                                'PHONE',        # Phone 
                                'WATER',        # Water
                            }  

# level3to2['UTILS_NIPA']      =   {                   # Total utilities expenditure  
#                                 'ENERGY_SERVICE',       # Energy
#                                 'PHONE',        # Phone 
#                                 'WATER',        # Water
#                             } 

level3to2['PERSCARE']   =   {                   # Personal Care
                                'PERSCARED',    # Personal Care Durables
                                'PERSCARES',    # Personal Care Services
                            }   

level3to2['HHOPER']     =   {                   # Household operations
                                'BABYSIT',      # Babysitter expenses
                                'ELDERLY',      # Elderly care expenses
                                'HOUSEHOLD',    # Other household operating expenses
                            }    

level3to2['MISCEL']     =   {                   # Miscellaneous expenses
                                'OCCUPEXP',     # Occupational expenses
                                'FINCHRG',      # Finance charges
                                'FEENCHRG',     # Other fees and charges
                            }  

level3to2['APPAREL']    =   {                   # Apparel
                                'CLOTHD',       # Clothing Durables
                                'CLOTHS',       # Clothing Services
                                'JEWELRY',      # Jewelery
                            }                          

level3to2['HEALTHCR']   =   {                   # Helath expenditures
                                'HEALTHINS',    # Health insurance
                                'HEALTHEXPS',   # Health services
                                'HEALTHEXPD',   # Health goods
                            } 

level3to2['SHELTER']    =   {                   # Shelter expenses
                                'MORTGAGEINT',  # Mortgage interest rate
                                'PROPTAX',      # Property taxes
                                'HOUSEEXPS',    # House expenditures services
                                'HOUSEEXPD',    # House expenditures durables
                                'HOUSEAWAY',    # Housing away 
                                'RENTPAID',     # Rent paid
                                'RENTASPAY',    # Rent as pay 
                                'RENTALEXPS',   # Rental expenditures services
                                'RENTALEXPD',   # Rental expenditures durables 
                            }   

level3to2['VEHRNT']     =   {                   # Vehicle rent, lease, other expenditure
                                'VEHLEASE',     # Vehicle lease payments
                                'VEHRENT',      # Vehiclee rental payments
                                'VEHOTHCH',     # Vehicle other payments
                            }

level3to2['ENTERTAIN']  =   {                   # Entertainment expenditures
                                'FEESADM',      # Fees and admissions
                                'TVAUDIO',      # TV and Audio
                                'PETSPLAY',     # Pets, toys
                                'ENTEROTH',     # Other entertainment
                            }  



# save output to yaml file
with open('../output/level3to2_map.yml', 'w') as yamlfile:
    yaml.dump(level3to2, yamlfile, default_flow_style=False) #, f, default_flow_style=False                    