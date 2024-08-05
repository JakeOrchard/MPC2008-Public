#%%
import pandas as pd



#%%
# ------------------------------------------------------------------------
# LOAD INPUT FILE
# ------------------------------------------------------------------------

# reads mtbi parquet file
df = pd.read_parquet('../input/mtbi.parquet')

#%%  ------------------------------------------------------------------------
# Some UCC codes should be split across NIPA categories
# this creates new ucc codes that are appropriate fractions of the originals.
# New UCC codes taken from BLS PCE-CEX concordance
# ------------------------------------------------------------------------

#UCC 450199	I	Used cars adjustment – new car trade-in

dfnew = df[df['UCC']==450116]
dfnew['COST'] = -1*dfnew['COST']
dfnew['UCC'] = 450199

df = df.append(dfnew)

#450299	I	Used trucks adjustment – new truck trade-in

dfnew = df[df['UCC']==450216]
dfnew['COST'] = -1*dfnew['COST']
dfnew['UCC'] = 450299

df = df.append(dfnew)

#322130	D	Infants' equipment, .1 of 320130

dfnew = df[df['UCC']==320130]
dfnew['COST'] = .1*dfnew['COST']
dfnew['UCC'] = 322130

df = df.append(dfnew)

#322904  Closet and storage 0.15 of 320904

dfnew = df[df['UCC']==320904]
dfnew['COST'] = 0.15*dfnew['COST']
dfnew['UCC'] = 322904

df = df.append(dfnew)

#321150  Outdoor Equipment 0.81 of 320150

dfnew = df[df['UCC']==320150]
dfnew['COST'] = 0.81*dfnew['COST']
dfnew['UCC'] = 321150

df = df.append(dfnew)

#Some appliances receive greater weight, but uccs do not change
#Dishwashers Built in
df.loc[df['UCC']== 230118,'COST']= 1.13*df.loc[df['UCC']== 230118,'COST']
#Dishwashers 
df.loc[df['UCC']== 220612,'COST']= 1.13*df.loc[df['UCC']== 220612,'COST']
#Fridge (owned home)
df.loc[df['UCC']== 300112,'COST']= 1.13*df.loc[df['UCC']== 300112,'COST']
#300212	I	Washing machines (owned home) (thru 2013 Q1)
df.loc[df['UCC']== 300212,'COST']= 1.13*df.loc[df['UCC']== 300212,'COST']
#300222	I	Clothes dryers (owned home) (thru 2013 Q1)
df.loc[df['UCC']== 300222,'COST']= 1.13*df.loc[df['UCC']== 300222,'COST']
#300312	I	Cooking stoves, ovens (owned home)
df.loc[df['UCC']== 300312,'COST']= 1.13*df.loc[df['UCC']== 300312,'COST']
#300322	I	Microwave ovens (owned home)
df.loc[df['UCC']== 300322,'COST']= 1.13*df.loc[df['UCC']== 300322,'COST']
#300332	I	Portable dishwasher (owned home)
df.loc[df['UCC']== 300332,'COST']= 1.13*df.loc[df['UCC']== 300332,'COST']


#322150 Outdoor Equipment 0.19 of 320150

dfnew = df[df['UCC']==320150]
dfnew['COST'] = 0.19*dfnew['COST']
dfnew['UCC'] = 322150 

df = df.append(dfnew)

#321905 Miscellaneous household equipment and parts 0.09 of 320905

dfnew = df[df['UCC']==320905]
dfnew['COST'] = 0.09*dfnew['COST']
dfnew['UCC'] = 321905

df = df.append(dfnew)

#333510	Miscellaneous household products .15 of 330510

dfnew = df[df['UCC']==330510]
dfnew['COST'] = 0.15*dfnew['COST']
dfnew['UCC'] = 333510

df = df.append(dfnew)

#323904		Closet and storage items .45 of 320904

dfnew = df[df['UCC']==320904]
dfnew['COST'] = 0.45*dfnew['COST']
dfnew['UCC'] = 323904

df = df.append(dfnew)

# 322905		Miscellaneous household equipment and parts .88 of  320905

dfnew = df[df['UCC']==320905]
dfnew['COST'] = 0.88*dfnew['COST']
dfnew['UCC'] = 322905

df = df.append(dfnew)

#311900	D	Accessories for electronic equipment 0.05 of UCC 310900


dfnew = df[df['UCC']==310900]
dfnew['COST'] = 0.05*dfnew['COST']
dfnew['UCC'] = 311900

df = df.append(dfnew)

# 692230	I	Business equipment for home use .6 of 690230

dfnew = df[df['UCC']==690230]
dfnew['COST'] = 0.6*dfnew['COST']
dfnew['UCC'] = 692230
df = df.append(dfnew)

#311230	D	Video game hardware/software (thru 2011 Q1) .5 of UCC 310230 

dfnew = df[df['UCC']==310230 ]
dfnew['COST'] = 0.5*dfnew['COST']
dfnew['UCC'] = 311230
df = df.append(dfnew)

#691230	I	Business equipment for home use, .4 of UCC 690230

dfnew = df[df['UCC']== 690230 ]
dfnew['COST'] = 0.4*dfnew['COST']
dfnew['UCC'] = 691230
df = df.append(dfnew)	

#662000	D	School supplies, etc. - unspecified .2 of UCC 660000


dfnew = df[df['UCC']== 660000 ]
dfnew['COST'] = 0.2*dfnew['COST']
dfnew['UCC'] = 662000
df = df.append(dfnew)

# 551320	I	Medical equip. for general use .5 of UCC 550320

dfnew = df[df['UCC']== 550320]
dfnew['COST'] = 0.5*dfnew['COST']
dfnew['UCC'] = 551320
df = df.append(dfnew)

#661900	I	"School books, supplies, & equip for day care, nursery, other .25 of UCC 660900 

dfnew = df[df['UCC']== 660900]
dfnew['COST'] = 0.25*dfnew['COST']
dfnew['UCC'] = 661900
df = df.append(dfnew)

#661901	I	"School books, supplies, & equip for day care, nursery .25 of UCC 660901

dfnew = df[df['UCC']== 660901]
dfnew['COST'] = 0.25*dfnew['COST']
dfnew['UCC'] = 661901	
df = df.append(dfnew)

#661210	I	School books/supplies & equip for elementary/high school .75 of UCC 660210
dfnew = df[df['UCC']== 660210]
dfnew['COST'] = 0.75*dfnew['COST']
dfnew['UCC'] = 661210	
df = df.append(dfnew)

#323130	D	Infants' equipment .1 of UCC 320130
dfnew = df[df['UCC']== 320130]
dfnew['COST'] = 0.1*dfnew['COST']
dfnew['UCC'] = 323130	
df = df.append(dfnew)

#313900	D	Accessories for electronic equipment .02 of UCC 310900
dfnew = df[df['UCC']== 310900]
dfnew['COST'] = 0.02*dfnew['COST']
dfnew['UCC'] = 313900	
df = df.append(dfnew)

#411130	D	Infants' underwear .62 of UCC 410130
dfnew = df[df['UCC']== 410130]
dfnew['COST'] = .62*dfnew['COST']
dfnew['UCC'] = 411130	
df = df.append(dfnew)

#251911	I	Coal, wood, other fuels (renter) (introduced 2005 Q2) .24 of UCC 250911
dfnew = df[df['UCC']== 250911]
dfnew['COST'] = .24*dfnew['COST']
dfnew['UCC'] = 251911	
df = df.append(dfnew)

#251912	I	Coal, wood, other fuels (introduced 2005 Q2) .24 of UCC 250912
dfnew = df[df['UCC']== 250912]
dfnew['COST'] = .24*dfnew['COST']
dfnew['UCC'] = 251912	
df = df.append(dfnew)

#251913	I	Coal, wood, other fuels (owned vacation) (introduced 2005 Q2) .24 of UCC 250913
dfnew = df[df['UCC']== 250913]
dfnew['COST'] = .24*dfnew['COST']
dfnew['UCC'] = 251913	
df = df.append(dfnew)

#251914	I	Coal, wood, other fuels (rented vacation) (introduced 2005 Q2) .24 of UCC 250914
dfnew = df[df['UCC']== 250914]
dfnew['COST'] = .24*dfnew['COST']
dfnew['UCC'] = 251914	
df = df.append(dfnew)

#252911	I	Coal, wood, other fuels (renter) (introduced 2005 Q2) 0.76 of UCC 250911
dfnew = df[df['UCC']== 250911]
dfnew['COST'] = .76*dfnew['COST']
dfnew['UCC'] = 252911	
df = df.append(dfnew)

#252912	I	Coal, wood, other fuels  0.76 of UCC 250912
dfnew = df[df['UCC']== 250912]
dfnew['COST'] = .76*dfnew['COST']
dfnew['UCC'] = 252912	
df = df.append(dfnew)

#252913	I	Coal, wood, other fuels  0.76 of UCC 250913
dfnew = df[df['UCC']== 250913]
dfnew['COST'] = .76*dfnew['COST']
dfnew['UCC'] = 252913	
df = df.append(dfnew)

#252914	I	Coal, wood, other fuels 0.76 of UCC 250914
dfnew = df[df['UCC']== 250914]
dfnew['COST'] = .76*dfnew['COST']
dfnew['UCC'] = 252914	
df = df.append(dfnew)

#552320	I	Medical equip. for general use 0.5 of UCC 550320
dfnew = df[df['UCC']== 550320]
dfnew['COST'] = .5*dfnew['COST']
dfnew['UCC'] = 552320	
df = df.append(dfnew)

#312230	D	Video game hardware/software (thru 2011 Q1) .5 of UCC 310230
dfnew = df[df['UCC']== 310230]
dfnew['COST'] = .5*dfnew['COST']
dfnew['UCC'] = 312230	
df = df.append(dfnew)

#331510	D	Miscellaneous household products .5 of UCC 330510
dfnew = df[df['UCC']== 330510]
dfnew['COST'] = .5*dfnew['COST']
dfnew['UCC'] = 331510	
df = df.append(dfnew)

#321904	D	Closet and storage items .4 of UCC 320904
dfnew = df[df['UCC']== 320904]
dfnew['COST'] = .4*dfnew['COST']
dfnew['UCC'] = 321904		
df = df.append(dfnew)

#332510	D	Miscellaneous household products .35 of UCC 330510
dfnew = df[df['UCC']== 330510]
dfnew['COST'] = .35*dfnew['COST']
dfnew['UCC'] = 332510		
df = df.append(dfnew)

#321130	D	Infants' equipment .8 of UCC 320130
dfnew = df[df['UCC']== 320130]
dfnew['COST'] = .8*dfnew['COST']
dfnew['UCC'] = 321130		
df = df.append(dfnew)

#412130	D	Infants' underwear .38 of UCC 410130
dfnew = df[df['UCC']== 410130]
dfnew['COST'] = .38*dfnew['COST']
dfnew['UCC'] = 412130		
df = df.append(dfnew)

#662900	I	"School books, supplies, & equip for day care, nursery, other .75 of UCC 660900
dfnew = df[df['UCC']== 660900]
dfnew['COST'] = .75*dfnew['COST']
dfnew['UCC'] = 662900		
df = df.append(dfnew)

#662901	I	"School books,supplies, & equip for day care, nursery(introduced 2006 Q2)" .75 of UCC 660901
dfnew = df[df['UCC']== 660901]
dfnew['COST'] = .75*dfnew['COST']
dfnew['UCC'] = 662901		
df = df.append(dfnew)

#662210	I	School books/supplies & equip for elementary/high school .25 of UCC 660210
dfnew = df[df['UCC']== 660210]
dfnew['COST'] = .25*dfnew['COST']
dfnew['UCC'] = 662210		
df = df.append(dfnew)

#661000	D	School supplies, etc. - unspecified .8 of UCC 660000
dfnew = df[df['UCC']== 660000]
dfnew['COST'] = .8*dfnew['COST']
dfnew['UCC'] = 661000		
df = df.append(dfnew)

#Owned Equivalent rent multiplied by 12 
df.loc[df['UCC']== 910050,'COST']= 12*df.loc[df['UCC']== 910050,'COST']
df.loc[df['UCC']== 910101,'COST']= 12*df.loc[df['UCC']== 910101,'COST']

#Vacation home equivalent rent multiplied by 6
df.loc[df['UCC']== 910100,'COST']= 6*df.loc[df['UCC']== 910100,'COST']
df.loc[df['UCC']== 910102,'COST']= 6*df.loc[df['UCC']== 910102,'COST']

#212310	I	Housing while attending school .05 of UCC 210310
dfnew = df[df['UCC']== 210310]
dfnew['COST'] = .05*dfnew['COST']
dfnew['UCC'] = 212310		
df = df.append(dfnew)

#571230	I	Other medical care services .5 of UCC 570230
dfnew = df[df['UCC']== 570230]
dfnew['COST'] = .5*dfnew['COST']
dfnew['UCC'] = 571230		
df = df.append(dfnew)

#572230	I	Other medical care services .5 of UCC 570230
dfnew = df[df['UCC']== 570230]
dfnew['COST'] = .5*dfnew['COST']
dfnew['UCC'] = 572230	
df = df.append(dfnew)

#531110	I	Airline fares .56 of  UCC 530110
dfnew = df[df['UCC']== 530110]
dfnew['COST'] = .56*dfnew['COST']
dfnew['UCC'] = 531110	
df = df.append(dfnew)

#531901	I	Ship fares .31 of UCC 530901
dfnew = df[df['UCC']== 530901]
dfnew['COST'] = .31*dfnew['COST']
dfnew['UCC'] = 531901	
df = df.append(dfnew)

#621310	I	Fees for recreational lessons .5 of UCC 620310
dfnew = df[df['UCC']== 620310]
dfnew['COST'] = .5*dfnew['COST']
dfnew['UCC'] = 621310		
df = df.append(dfnew)

#621904	I	Rent/repair music instruments .8 of UCC 620904
dfnew = df[df['UCC']== 620904]
dfnew['COST'] = .8*dfnew['COST']
dfnew['UCC'] = 621904	
df = df.append(dfnew)

#211310	I	Housing while attending school .95 of UCC 210310
dfnew = df[df['UCC']== 210310]
dfnew['COST'] = .95*dfnew['COST']
dfnew['UCC'] = 211310
df = df.append(dfnew)

#622310	I	Fees for recreational lessons .5 of UCC 620310
dfnew = df[df['UCC']== 620310]
dfnew['COST'] = .5*dfnew['COST']
dfnew['UCC'] = 622310
df = df.append(dfnew)

#532110	I	Airline fares .44 of UCC 530110
dfnew = df[df['UCC']== 530110]
dfnew['COST'] = .44*dfnew['COST']
dfnew['UCC'] = 532110
df = df.append(dfnew)



#532901	I	Ship fares .69 of UCC 530901
dfnew = df[df['UCC']== 530901]
dfnew['COST'] = .69*dfnew['COST']
dfnew['UCC'] = 532901
df = df.append(dfnew)


#Save as new file
df.to_parquet('../output/mtbi_wnewUCC.parquet')








# %%
