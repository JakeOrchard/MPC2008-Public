import yaml

# dictionary with mappings
# Used CEX UCC hierachy and BEA Personal Consumption Expenditure underlying detail Table 2.4.5U
nipacat = {}

#DURABLES

nipacat['MV_PARTS'] = {

                        'newautos',      
                        'newlighttrucks',
                        'usedautos',
                        'usedlighttrucks',            
                        'tires',
                        'acces_parts',
                                    }

nipacat['NET_USED'] = {

                        'usedautos',
                        'usedlighttrucks',            
                                    }

nipacat['CARS_NU'] = {

                        'newautos',      
                        'newlighttrucks',
                        'usedautos',
                        'usedlighttrucks',            
                                    }

nipacat['CARS_N'] = {

                        'newautos',      
                        'newlighttrucks',          
                                    }

nipacat['CARS_DOWNPAYMENT'] = {

                        'usedcardown',      
                        'newcardown',          
                                    }

nipacat['PCARS_N'] = {

                        'newautos',      
                                    }

nipacat['TRUCKS_N'] = {

                        'newlighttrucks',          
                                    }

nipacat['HOME_DUR'] = {
                        'furniture',
                        'clock',
                        'carpets',
                        'window_cover',
                        'major_app',
                        'small_app',
                        'dishes',
                        'nelec_cook',
                        'tools',
}

nipacat['REC_DUR'] = {
                        'tv' ,
                        'other_video',
                        'audio_equip',
                        'audio_record',
                        'video_record',
                        'photo_equip',
                        'computer_equip',
                        'software',
                        'calculator',
                        'sport_equip',
                        'motorcycles',
                        'bicycles',
                        'boats',
                        'planes',
                        'other_recv',
                        'rec_book',
                        'music_instrument',
                    }

nipacat['OTHER_DUR'] = {
                        'jewelry',
                        'watches',
                        'med_equip_t',
                        'eyeglass',
                        'educ_book',
                        'luggage',
                        'telephone',
                        }

nipacat['DURABLE_NIPA'] = {

                        'MV_PARTS',
                        'HOME_DUR',
                        'REC_DUR',
                        'OTHER_DUR',
}

#NON-DURABLES
nipacat['FOODBEV_HOME'] = {
                        'food_out',       
                        'food_stamps',        
                        'foodhome',
                        'alchome',
}

nipacat['CLOTHING'] = {
                     'women_cloth',
                     'men_cloth',
                     'child_cloth',
                     'material_cloth',
                     'shoes',
}

nipacat['GASOLINE_cat'] = {

                    'gasoline', 
                    'motor_oil',
                    'fuel_oil',
                    'other_fuel',
}

nipacat['OTHER_NONDUR'] = {  
                    'prescription',
                    'non_prescription',
                    'other_medical',
                    'games',
                    'pet_food',
                    'flower',
                    'photo_supply',
                    'house_clean',
                    'house_paper',
                    'house_linen',
                    'sewing',
                    'house_misc',
                    'pers_care',
                    'cosmetic',
                    'pers_elec',
                    'tobacco',
                    'newspaper',
                    'stationary',
}

nipacat['NONDURABLE_NIPA'] = {

                        'FOODBEV_HOME',
                        'CLOTHING',
                        'GASOLINE_cat',
                        'OTHER_NONDUR',
}

#SERVICES

nipacat['HOUSE_NIPA'] ={
                        'rent',
                        'own_home',        
                        'group_home',     
                        'water',  
                        'trash',
                        'electric',
                        'nat_gas',

}

nipacat['HEALTH_NIPA'] = {  #HEALTH SERVICES: Very unlikely to match, as CEX does not include most direct payments
                        'dr_service',
                        'dental_service', 
                        'home_health',
                        'med_lab',
                        'med_outpatient',
                        'other_med',
                        'hospital',
                        'nurse_home',
                        }

nipacat['TRANSPORTATION'] ={ 
                
                        'auto_repair',
                        'car_lease',
                        'car_rental',
                        'parking_fee',
                        'train',
                        'bus',
                        'taxi',
                        'intracity',
                        'carpool',
                        'air',
                        'ship',
}

nipacat['REC_SERVICE'] = {
                        'member_club',
                        'amusement',
                        'theater',
                        'live_entertain',
                        'spec_sport',
                        'cable',
                        'film_process',
                        'photo_service',
                        'repair_tvaudio',
                        'video_rental',
                        'gambling',
                        'vet_service',
                        'repair_rec',
}

nipacat['ACC_SERVICE'] = {
                        'food_school',
                        'foodaway',
                        'alcaway',
                        'hotels',
                        'school_housing',
}

nipacat['FIN_SERVICE_cat'] = {   # Also unlikely to match since PCE is net finance charges
                        'fin_service',
                        'life_insur',
                        'home_insur',
                        'health_insur',
                        'auto_insur',
}

nipacat['OTHER_SERVICE'] = {
                        'tel_service',
                        'cell_service',
                        'post',
                        'internet',
                        'college_tuition',
                        'elem_tuition',
                        'preschool',
                        'other_school',
                        'legal',
                        'accounting',
                        'other_bus_service',
                        'funeral',
                        'hair',
                        'other_pers_service',
                        'laundry',
                        'cloth_repair',
                        'shoe_repair',
                        'childcare',
                        'soc_assist',
                        'charity',
                        'dom_service',
                        'moving',
                        'furn_repair',
                        'app_repair',
                        'other_house_service',
                        'foreign_travel',
}

nipacat['SERVICE_NIPA'] = {

                        'HOUSE_NIPA',
                        'HEALTH_NIPA',
                        'TRANSPORTATION',
                        'REC_SERVICE',
                        'ACC_SERVICE',
                        'FIN_SERVICE_cat',
                        'OTHER_SERVICE',
}

nipacat['PCE'] = {

                        'SERVICE_NIPA',
                        'NONDURABLE_NIPA',
                        'DURABLE_NIPA',
}


#NEW Categories (for paper)

nipacat['SNDEXP_NIPA']   =   {                  
                            'food_school',
                            'foodaway',
                            'alcaway', 
                            'food_out',       
                            'foodhome',
                            'alchome',
                            'tobacco',  
                            'gasoline', 
                            'motor_oil',
                            'fuel_oil',
                            'other_fuel',  
                            'water',  
                            'trash',
                            'electric',
                            'nat_gas',
                            'tel_service',
                            'cell_service',
                            'pers_elec' ,
                            'pers_care', #Many UCC not in PSMJ       
                            'hair',
                            'other_pers_service',       # personal care ---------- persca
                            'childcare',
                            'preschool',
                            'soc_assist',
                            'other_house_service',
                            'dom_service',
                            'internet',
                            'moving',
                            'train',
                            'bus',
                            'taxi',
                            'intracity',
                            'carpool',
                            'air',
                            'ship',
                            'fin_service',  
                            'gambling',
                            'legal',
                            'funeral',
                            'accounting',
}   


nipacat['NDEXP_NIPA']   =   {  
                            'SNDEXP_NIPA',
                            'women_cloth',
                            'men_cloth',
                            'child_cloth',
                            'material_cloth',
                            'shoes',
                            'jewelry',
                            'watches',
                            'laundry',
                            'cloth_repair',
                            'shoe_repair',
                            'health_insur',
                            'dr_service',
                            'dental_service', 
                            'home_health',
                            'med_lab',
                            'med_outpatient',
                            'other_med',
                            'hospital',
                            'nurse_home',
                            'med_equip_t',
                            'eyeglass',
                            'newspaper',
                            'rec_book',                       

}

                        
nipacat['MV_PARTS_SER'] = {

                        'newautos',      
                        'newlighttrucks',
                        'usedautos',
                        'usedlighttrucks',            
                        'tires',
                        'acces_parts',
                        'auto_repair',
}
        
nipacat['TUITION'] = {
                        'college_tuition',
                        'elem_tuition',
                        'preschool',
                        'other_school',
}


nipacat['ALCBEVER']   =   {                   # Total alcohol expenditure  
                                'ALCHOME',      # Alcohol at hom
                                'ALCAWAY',      # Alcohol away
}

nipacat['FOODBEVS'] =   {                   # (1) FOOD+BEV
                            'FOODHOME',       # Food consumed away from home ---- FDAWAY
                            'FOODAWAY',       # Food consumed at home ----------- FDHOME
                            'ALCBEVER',       # Purchases of alcholic beverages - ALCBEV
}

# save output to yaml file
# need sort_keys=False to preserve ordering since construction depends on previous
# entries
with open('../output/nipacategories_map.yml', 'w') as yamlfile:
    yaml.dump(nipacat, yamlfile, default_flow_style=False, sort_keys=False) #, f, default_flow_style=False                    