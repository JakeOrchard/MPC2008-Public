import yaml

# dictionary with mappings
nameconversion = {}

nameconversion['TOBACCO']     = 'TOBACC'
nameconversion['FOODHOME']    = 'FDHOME'
nameconversion['FOODAWAY']    = 'FDAWAY'
nameconversion['ALCBEVER']    = 'ALCBEV'
nameconversion['UTILS']       = 'UTIL'
nameconversion['PERSCARE']    = 'PERSCA'
nameconversion['HHOPER']      = 'HOUSOP'
nameconversion['PUBTRANS']    = 'PUBTRA'
nameconversion['GASOLINE']    = 'GASMO'
nameconversion['MISCEL']      = 'MISC'        
nameconversion['HEALTHCR']    = 'HEALTH'
nameconversion['APPAREL']     = 'APPAR'       
nameconversion['READING']     = 'READ'
nameconversion['EDUCATION']   = 'EDUCA'
nameconversion['FURNITURE']   = 'HOUSEQ'
nameconversion['SHELTER']     = 'SHELT'
nameconversion['NEWCARS']     = 'CARTKN'
nameconversion['USEDCARS']    = 'CARTKU'
nameconversion['CAREPAIR']    = 'MAINRP'
nameconversion['VEHRNT']      = 'VRNTLO'
nameconversion['VEHFINCH']    = 'VEHFIN'
nameconversion['OTHVEHCL']    = 'OTHVEH'
nameconversion['VEHINSUR']    = 'VEHINS'
nameconversion['ENTERTAIN']   = 'ENTERT'
nameconversion['FEESADM']     = 'FEEADM'
nameconversion['TVAUDIO']     = 'TVRDIO'
nameconversion['PETSPLAY']    = 'PETTOY'
nameconversion['ENTEROTH']    = 'OTHENT'
nameconversion['LIFEINSUR']   = 'LIFINS'
nameconversion['PENSIONS']    = 'RETPEN'      
# nameconversion['TOTEXP2']     = 'TOTEXP'
# nameconversion['TOTEXP3']     = 'TOTEXP'

# save output to yaml file
with open('../output/cexname_conversion.yml', 'w') as yamlfile:
    yaml.dump(nameconversion, yamlfile, default_flow_style=False)