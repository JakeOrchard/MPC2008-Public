import yaml

# dictionary with mappings
uccmapping = {}

# JW: perhaps these mappings should be sets rather than lists?
#JO: creates some alternative mappings to match with NIPA current list of variables can overlap

uccmapping['FOODHOME']  = [         #	3	FOODHOME	Food at home	
                            190904 	,	#	4	190904	Food prepared by consumer unit on out-of-town trips	
                                        #		790210	Total purchases at grocery stores (Note: This is not a primary expenditure)
                            790220	,	#	4	790220	Grocery stores (thru Q20071)	
                            790230	,	#	4	790230	Convenience stores (thru Q20071)	
                            790240]		#	4	790240	Food and non alcoholic beverages (new UCC Q20072)	


uccmapping['FOODAWAY']	= [    	    #	3	FOODAWAY	Food away from home	
                            190901,	#	4	190901	Food or board at school	
                            190902,	#	4	190902	Catered affairs	
                            190903,	#	4	190903	Food on out-of-town trips	
                            790410,	#	4	790410	Meals at restaurants carry-outs and other	
                            790430,		#	4	790430	School lunches
                            800700]		#	4	800700	Meals as pay	
                            
                            
uccmapping['ALCHOME']   = [    	        #	3	ALCHOME	Alcoholic beverages consumed at home	
                            790310, 	#	4	790310	Beer and wine (thru Q20071)	
                            790320,	    #	4	790320	Other alcoholic beverages (thru Q20071)	
                            790330]		#	4	790330	Beer wine other alcohol (new UCC Q20072) (Note:790330=790310+790320)
                            
uccmapping['ALCAWAY']   = [    	        #	3	ALCAWAY	Away from home	
                            790420,	    #	4	790420	Alcoholic beverages at restaurants taverns	
                            200900]		#	4	200900	Alcoholic beverages purchased on trips	
                                    
                                        
                                        
                                        #	3	SHELTER	Shelter	
                                        #	4	OWNDWELL	Owned dwellings	 
uccmapping['MORTGAGEINT']   = [    	    #		MORTGAGEINT	Mortgage interest expenses	
                            220311,	    #	6	220311	Mortgage interest	
                            220312,	    #	7	220312	Mortgage interest (vacation homes) 	
                            220313,	    #	6	220313	Interest paid; home equity loan	
                            220314,	    #	7	220314	Interest paid; home equity loan (vacation homes)	
                            220321,	    #	6	220321	Prepayment penalty charges (thru Q20071)	
                            220322,	    #	7	220322	Prepayment penalty charge (vacation homes) (thru Q20071)	
                            880110,	    #	6	880110	Interest paid; home equity line of credit		
                            880310]		#	7	880310	Interest paid; home equity line of credit (vacation homes)	

uccmapping['PROPTAX']   = [  	#		PROPTAX	Property taxes	
                            220211,	    #	5	220211	Property taxes	
                            220212]		#	6	220212	Property taxes (vacation homes) 	
                        
uccmapping['HOUSEEXPS'] = [    	        #		HOUSEEXPS	House expenditures (services)	
                #	5	OWNEXPEN	Maintenance; repairs; insurance; other expenses	
                            210901,	    #	6	210901	Ground rent	
                            210902,	    #	7	210902	Ground rent (vacation homes)	
                            220111,	    #		220111	Fire and extended coverage for owned homes is now included in UCC 220121
                            220112,	    #		220112	Fire and extended coverage for owned vacation homes now included in 220122
                            220121,	    #	6	220121	Homeowners insurance	
                            220122,	    #	7	220122	Homeowners insurance (vacation homes)	
                #	6	OWNREPSV	Maintenance and repair services	
                            230112,	    #	7	230112	Painting and papering	
                            230113,	    #	7	230113	Plumbing and water heating	
                            230114,	    #	7	230114	Heat; a/c; electrical work	
                            230115,	    #	7	230115	Roofing and gutters	
                            230116,	    #		230116	Contractors' labor and material costs, owned home (now mapped to 230151)
                            230119,	    #		230119	Same as 230116 - owned vacation home, vacation condos and co-ops (to 230152)
                            230122,	    #	7	230122	Repair and replacement of hard surface flooring	
                            230123,	    #	8	230123	Repair and replacement of hard surface flooring (vacation homes)	
                            230142,	    #	7	230142	Repair of built-in appliances	
                            230151,	    #	7	230151	Other repair and maintenance services	
                            230152,	    #	8	230152	Repair and remodeling services (vacation homes)	
                #	6	OWNMNAGE	Property management and security	
                            220901,	    #	6	220901	Parking	
                            220902,	    #	7	220902	Parking (vacation homes)	
                            230901,	    #	7	230901	Property management	
                            230902,	    #	8	230902	Property management (vacation homes)	
                            340911,	    #	7	340911	Management and upkeep services for security	
                            340912,	    #	8	340912	Management and upkeep services for security (vacation homes)
                            680905]	    #	3	680905	Vacation clubs (new UCC Q20072)	
                                        #               [Hierarchical listing puts 680905 into MISC dictionary places it here 
                                        #                and the latter is how FMLI files are constructed]	
                            	
                        
uccmapping['HOUSEEXPD']   = [    	        #		HOUSEEXPD	House expenditures (durables)	
                #	6	OWNREPSP	Maintenance and repair commodities	
                            240112,	    #	7	240112	Paints; wallpaper and supplies	
                            240113,	    #	8	240113	Paints; wallpaper; supplies (vacation homes)	
                            240122,	    #	7	240122	Tools and equipment for painting and wallpapering	
                            240123,	    #	8	240123	Tools and equipment for painting and wallpapering (vacation homes)	
                            240212,	    #	7	240212	Materials for plaster.; panel.; siding; windows; doors; screens; awnings
                            240213,	    #	7	240213	Materials and equipment for roof and gutters	
                            240214,	    #	8	240214	Materials for plastering; paneling; roofing; gutters; ...  (vacation homes)
                            240222,	    #	7	240222	Materials for patio; walk; fence; driveway;masonry; brick and stucco work
                            240223,	    #	8	240223	Material for patio; walk; fence; drive;masonry; brick; stucco (vac. homes)
                            240312,	    #	7	240312	Plumbing supplies and equipment	
                            240313,	    #	8	240313	Plumbing supplies and equipment (vacation homes)	
                            240322,	    #	7	240322	Electrical supplies; heating and cooling equipment	
                            240323,	    #	8	240323	Electrical supplies; heating and cooling equipment (vacation homes)	
                            320622,	    #	7	320622	Materials for hard surface flooring; repair and replacement	
                            320623,	    #	8	320623	Materials for hard surface flooring (vacation homes)	
                            320632,	    #	7	320632	Materials for landscaping maintenance	
                            320633,	    #	8	320633	Materials for landscaping maintenance (vacation homes)	
                #	7	OWNMISC	Miscellaneous supplies and equipment	
                            320612,	    #	8	320612	Material for insulation; other maintenance and repair	
                            320613,	    #	9	320613	Material for insulation; other maintenance and repair (vacation homes)	
                            990930,	    #	8	990930	Materials to finish basement; remodel rooms or build patios; walks; etc. 
                            990940]		#	9	990940	Material for finishing basements and remodeling rooms (vacation homes)	
                        
uccmapping['HOUSEAWAY'] = [    	        #		HOUSEAWAY	Housing away	
                            210210,	    #	5	210210	Lodging on out-of-town trips	
                            210310]		#	5	210310	Housing while attending school	
                        
uccmapping['RENTPAID']  = [    	        #		RENTPAID	Rent paid	
                            210110]		#	5	210110	Rent	
                        
uccmapping['RENTASPAY'] = [    	        #		RENTASPAY	Rent as pay  	
                            800710]		#	5	800710	Rent as pay (Note: I assume that the employer pays the CUs rent;         
                #                           but this is probably NOT income from renting an appartment!)
                        
uccmapping['RENTALEXPS']    = [    	        #		RENTALEXPS	rental expenditures (services)	
                #	5	RNTEXPEN	Maintenance; insurance and other expenses	
                            350110,	    #	6	350110	Tenant's insurance	
                #	6	RNTREPSV	Maintenance and repair services	
                            230111,	    #	7	230111	Contractors' labor and material costs, renter (mapped to 230150) (deleted 94)
                            230121,	    #	7	230121	Repair and replacement of hard surface flooring	
                            230141,	    #	7	230141	Repair of built-in appliances	
                            230150]		#	7	230150	Repair or maintenance services	
                        
uccmapping['RENTALEXPD']    = [    	        #		RENTALEXPD	rental expenditures (durables)	
                #	6	RNTREPSP	Maintenance and repair commodities	
                            240111,	    #	7	240111	Paint; wallpaper; and supplies	
                            240121,	    #	7	240121	Tools and equipment for painting and wallpapering	
                            240211,	    #	7	240211	Materials for plastering; panels; roofing;gutters; etc.	
                            240221,	    #	7	240221	Materials for patio; walk; fence; driveway;masonry; brick and stucco work
                            240311,	    #	F7	240311	Plumbing supplies and equipment	
                            240321,	    #	7	240321	Electrical supplies; heating and cooling equipment	
                #	7	RNTMISC	Miscellaneous supplies and equipment	
                            320611,	    #	8	320611	Material for insulation; other maintenance and repair	
                            320621,	    #	7	320621	Material for hard surface flooring	
                            320631,	    #	7	320631	Material for landscape maintenance	
                            790690,	    #	8	790690	Construction materials for jobs not started	
                            990920]		#	8	990920	Materials for additions; finishing basements;remodeling rooms	
                        
# uccmapping['FUELOIL'] = [ #Fuel Oil and other fuels (goods)

#  #	5	FUELOIL	Fuel oil	
#                             250111,	    #	6	250111	Fuel oil (renter)	
#                             250112,	    #	6	250112	Fuel oil (owned home)	
#                             250113,	    #	6	250113	Fuel oil (owned vacation)	
#                             250114,    #	6	250114	Fuel oil (rented vacation)	]
#                             #	5	BOTTLGAS	Bottled gas	
#                             250211,	    #	6	250211	Gas; btld/tank (renter)	
#                             250212,	    #	6	250212	Gas; btld/tank (owned home)	
#                             250213,	    #	6	250213	Gas; btld/tank (owned vacation)	
#                             250214,	    #	6	250214	Gas; btld/tank (rented vacation)	
#                 #	5	CLWDOTFL	Coal; wood; and other fuels	
#                             250221,	    #	6	250221	Coal - renter (deleted 05)	
#                             250222,	    #	6	250222	Coal - owned dwelling (deleted 05)	
#                             250223,	    #	6	250223	Coal - owned vacation property (deleted 05)	
#                             250224,	    #	6	250224	Coal - rented vacation property (deleted 05)	
#                             250901,	    #	6	250901	Wood/ Kerosene/ Other fuels - renter (deleted 05)	
#                             250902,	    #	6	250902	Wood/ Kerosene/ Other fuels - owned dwelling (deleted 05)	
#                             250903,	    #	6	250903	Wood/ Kerosene/ Other fuels - owned vacation property (deleted 05)	
#                             250904,	    #	6	250904	Wood/ Kerosene/ Other fuels - rented vacation property (deleted 05)	
#                             250911,	    #	6	250911	Coal; wood; other fuels (renter) [2005Q2] 	
#                             250912,	    #	6	250912	Coal; wood; other fuels (owned home) [2005Q2]	
#                             250913,	    #	6	250913	Coal; wood; other fuels (owned vacation) [2005Q2]	
#                             250914]		#	6	250914	Coal; wood; other fuels (rented vacation) [2005Q2]


# uccmapping['ENERGY_SERVICE']    = [    	        #		ENERGY	Energy	
#                 #	4	NATRLGAS	Natural gas	
#                             260211,	    #	5	260211	Utility--natural gas (renter)	
#                             260212,	    #	5	260212	Utility--natural gas (owned home)	
#                             260213,	    #	5	260213	Utility--natural gas (owned vacation)	
#                             260214,	    #	5	260214	Utility--natural gas (rented vacation)	
#                 #	4	ELECTRIC	Electricity	
#                             260111,	    #	5	260111	Electricity (renter)	
#                             260112,	    #	5	260112	Electricity (owned home)	
#                             260113,	    #	5	260113	Electricity (owned vacation)	
#                             260114]	    #	5	260114	Electricity (rented vacation)	
               
                	
                        


uccmapping['ENERGY']    = [    	        #		ENERGY	Energy	
                #	4	NATRLGAS	Natural gas	
                            260211,	    #	5	260211	Utility--natural gas (renter)	
                            260212,	    #	5	260212	Utility--natural gas (owned home)	
                            260213,	    #	5	260213	Utility--natural gas (owned vacation)	
                            260214,	    #	5	260214	Utility--natural gas (rented vacation)	
                #	4	ELECTRIC	Electricity	
                            260111,	    #	5	260111	Electricity (renter)	
                            260112,	    #	5	260112	Electricity (owned home)	
                            260113,	    #	5	260113	Electricity (owned vacation)	
                            260114,	    #	5	260114	Electricity (rented vacation)	
                #	4	OTHRFUEL	Fuel oil and other fuels	
                #	5	FUELOIL	Fuel oil	
                            250111,	    #	6	250111	Fuel oil (renter)	
                            250112,	    #	6	250112	Fuel oil (owned home)	
                            250113,	    #	6	250113	Fuel oil (owned vacation)	
                            250114,	    #	6	250114	Fuel oil (rented vacation)	
                #	5	BOTTLGAS	Bottled gas	
                            250211,	    #	6	250211	Gas; btld/tank (renter)	
                            250212,	    #	6	250212	Gas; btld/tank (owned home)	
                            250213,	    #	6	250213	Gas; btld/tank (owned vacation)	
                            250214,	    #	6	250214	Gas; btld/tank (rented vacation)	
                #	5	CLWDOTFL	Coal; wood; and other fuels	
                            250221,	    #	6	250221	Coal - renter (deleted 05)	
                            250222,	    #	6	250222	Coal - owned dwelling (deleted 05)	
                            250223,	    #	6	250223	Coal - owned vacation property (deleted 05)	
                            250224,	    #	6	250224	Coal - rented vacation property (deleted 05)	
                            250901,	    #	6	250901	Wood/ Kerosene/ Other fuels - renter (deleted 05)	
                            250902,	    #	6	250902	Wood/ Kerosene/ Other fuels - owned dwelling (deleted 05)	
                            250903,	    #	6	250903	Wood/ Kerosene/ Other fuels - owned vacation property (deleted 05)	
                            250904,	    #	6	250904	Wood/ Kerosene/ Other fuels - rented vacation property (deleted 05)	
                            250911,	    #	6	250911	Coal; wood; other fuels (renter) [2005Q2] 	
                            250912,	    #	6	250912	Coal; wood; other fuels (owned home) [2005Q2]	
                            250913,	    #	6	250913	Coal; wood; other fuels (owned vacation) [2005Q2]	
                            250914]		#	6	250914	Coal; wood; other fuels (rented vacation) [2005Q2]	
                        
uccmapping['PHONE']     = [    	        #	4	PHONE	Telephone services	
                            270000,	    #	5	270000	Telephone service, incl. public pay phones (deleted 91)	
                            270101,	    #	5	270101	Residential telephone/pay phones [1991Q2] 	
                            270102,	    #	5	270102	Cellular phone service  [1993Q2]	
                            270103,	    #	5	270103	Pager service (2001-2006)	
                            270104,	    #	5	270104	Phone cards (2001Q2)	
                            270105]		#	5	270105	Voice over IP service (new UCC Q20072)	


uccmapping['PHONED']     = [    	        #	4	PHONE	Telephone Durables
                            320232,	    #	5	320232	Telephones and accessories	
                            690210]    #	5	690210	Telephone answering devices	
                        
                        
uccmapping['WATER']     = [    	        #	4	WATER	Water and other public services	
                #	5	SEWER	Water and sewerage maintenance	
                            270211,	    #	6	270211	Water/sewer maint. (renter)	
                            270212,	    #	6	270212	Water/sewer maint. (owned home)	
                            270213,	    #	6	270213	Water/sewer maint. (owned vacation)	
                            270214,	    #	6	270214	Water/sewer maint. (rented vacation)	
                #	5	TRASH	Trash and garbage collection	
                            270411,	    #	6	270411	Trash/garb. coll. (renter)	
                            270412,	    #	6	270412	Trash/garb. coll. (owned home)	
                            270413,	    #	6	270413	Trash/garb. coll. (owned vacation)	
                            270414,	    #	6	270414	Trash/garb. coll. (rented vacation)	
                #	5	SEPTANK	Septic tank cleaning	
                            270901,	    #	6	270901	Septic tank clean. (renter)	
                            270902,	    #	6	270902	Septic tank clean. (owned home)	
                            270903,	    #	6	270903	Septic tank clean. (owned vacation)	
                            270904]		#	6	270904	Septic tank clean. (rented vacation)	
                        
uccmapping['BABYSIT']   = [    	        #		BABYSIT	Babysitting  	
                            340210,	    #		340210	Babysitting (deleted 93)	
                            340211,	    #	5	340211	Babysitting and child care in your own home (1993Q2)	
                            340212,	    #	5	340212	Babysitting and child care in someone else's home (1993Q2)	
                            670310]		#	5	670310	Day care centers; nursery; and preschools	
                        
uccmapping['ELDERLY']   = [    	        #		ELDERLY	Elderly care	
                            340906,	    #	5	340906	Care for elderly; invalids; handicapped; etc.	
                            340910]		#	5	340910	Adult day care centers	


# uccmapping['HOUSEHOLDS'] = [    	        #		HOUSEHOLD	Household related services	
#                 #	4	HHOTHXPN	Other household expenses	
#                             340310,	    #	5	340310	Housekeeping services	
#                             340410,	    #	5	340410	Gardening; lawn care service	
#                             340420,	    #	5	340420	Water softening service	
#                             340510,	    #	5	340510	Moving; storage; freight
#                             340620,	    #	5	340620	Appliance repair; including service center	
#                             340630,	    #	5	340630	Reupholstering; furniture repair	
#                             340901,	    #	5	340901	Repairs/rentals of lawn and garden equipment; other household equipment
#                             340903,	    #	5	340903	Other home services	
#                             340907,	    #	5	340907	Appliance rental	
#                             340908,	    #	5	340908	Rental of office equipment for non-business use	
#                             340914,	    #	5	340914	Services for termite/pest control	
#                             340915,	    #	5	340915	Home security system service fee	
#                             690113,	    #	5	690113	Repair of computer systems for non-business use	
#                             990900]	    #	5	990900	Rental and installation of dishwashers; range hoods;and garbage disposals


uccmapping['HOUSEHOLD'] = [    	        #		HOUSEHOLD	Household related expenses	
                #	4	HHOTHXPN	Other household expenses	
                            330511,	    #	5	330511	Termite/pest control products	
                            340310,	    #	5	340310	Housekeeping services	
                            340410,	    #	5	340410	Gardening; lawn care service	
                            340420,	    #	5	340420	Water softening service	
                            340510,	    #	5	340510	Moving; storage; freight	
                            340520,	    #	5	340520	Household laundry and dry cleaning; sent out (non-clothing) not coin-operated
                            340530,	    #	5	340530	Coin-operated household laundry and dry cleaning (non-clothing)	
                            340620,	    #	5	340620	Appliance repair; including service center	
                            340630,	    #	5	340630	Reupholstering; furniture repair	
                            340901,	    #	5	340901	Repairs/rentals of lawn and garden equipment; other household equipment
                            340903,	    #	5	340903	Other home services	
                            340907,	    #	5	340907	Appliance rental	
                            340908,	    #	5	340908	Rental of office equipment for non-business use	
                            340914,	    #	5	340914	Services for termite/pest control	
                            340915,	    #	5	340915	Home security system service fee	
                            690113,	    #	5	690113	Repair of computer systems for non-business use	
                            690114,	    #	5	690114	Computer information services	
                            690116,	    #	5	690116	Internet services away from home [2005]	
                                        #       [690116: Per dictionary included here but hierarchical listing says it belongs in household equipment]
                            690310,	    #	5	690310	Installation of computer (new UCC Q20072)	
                            990900,	    #	5	990900	Rental and installation of dishwashers; range hoods;and garbage disposals
                            990910]		#		990910	Materials for pest control are no longer collected as cap. improvement

uccmapping['FURNITURE'] = [    	        #		FURNITURE	Furniture 	
                #	4	HHTXTILE	Household textiles	
                            280110,	    #	5	280110	Bathroom linens	
                            280120,	    #	5	280120	Bedroom linens	
                            280130,	    #	5	280130	Kitchen and dining room linens	
                            280210,	    #	5	280210	Curtains and draperies	
                            280220,	    #	5	280220	Slipcovers; decorative pillows	
                            280230,	    #	5	280230	Sewing materials for slipcovers; curtains; other sewing materials for home
                            280900,	    #	5	280900	Other linens	
                #	4	FURNITUR	Furniture	
                            290110,	    #	5	290110	Mattress and springs	
                            290120,	    #	5	290120	Other bedroom furniture	
                            290210,	    #	5	290210	Sofas	
                            290310,	    #	5	290310	Living room chairs	
                            290320,	    #	5	290320	Living room tables	
                            290410,	    #	5	290410	Kitchen; dining room furniture	
                            290420,	    #	5	290420	Infants' furniture	
                            290430,	    #	5	290430	Outdoor furniture	
                            290440,	    #	5	290440	Wall units; cabinets and other occasional furniture	
                #	4	FLOORCOV	Floor coverings	
                            230131,	    #	5	230131	Installed Wall-to-Wall Carpeting - Renter  (deleted 99)	
                            230132,	    #	5	230132	Wall-to-wall carpet installed replacement owned (newly added UCC 230133)
                            230133,	    #	5	230133	Wall-to-wall carpet (replacement) (owned home)	
                            230134,	    #	6	230134	Wall-to-wall carpet (renter)	
                            320110,	    #	5	320110	Non-permanent floor coverings is now included in UCC 320111 (deleted 99)
                            320111,	    #	5	320111	Floor coverings; non-permanent	
                            320161,	    #	5	320161	Non-Installed Wall to Wall Carpeting and Carpet Squares - Renter 
                            320162,	    #	5	320162	Non-Installed Wall-to-Wall Carpeting and Carpet Squares - Homeowner 
                            320163,	    #	6	320163	Wall-to-wall carpet (replacement)(renter) [1999]	
                #	4	MAJAPPL	Major appliances	
                            230117,	    #	5	230117	Dishwashers (built-in); garbage disposals; range hoods; (renter)
                            230118,	    #	5	230118	Dishwashers (built-in); garbage disposals; range hoods; (owned home)
                            300111,	    #	5	300111	Refrigerators; freezers (renter)	
                            300112,	    #	5	300112	Refrigerators; freezers (owned home)	
                            300211,	    #	5	300211	Washing machines (renter)	
                            300212,	    #	5	300212	Washing machines (owned home)	
                            300221,	    #	5	300221	Clothes dryers (renter)	
                            300222,	    #	5	300222	Clothes dryers (owned home)	
                            300311,	    #	5	300311	Cooking stoves; ovens (renter)	
                            300312,	    #	5	300312	Cooking stoves; ovens (owned home)	
                            300321,	    #	5	300321	Microwave ovens (renter)	
                            300322,	    #	5	300322	Microwave ovens (owned home)	
                            300331,	    #	5	300331	Portable dishwasher (renter)	
                            300332,	    #	5	300332	Portable dishwasher (owned home)	
                            300411,	    #	5	300411	Window air conditioners (renter)	
                            300412,	    #	5	300412	Window air conditioners (owned home)	
                            320511,	    #	5	320511	Electric floor cleaning equipment	
                            320512,	    #	5	320512	Sewing machines	
                #	4	SMAPPHWR	Small appliances; miscellaneous housewares	
                            320310,	    #	6	320310	Plastic dinnerware	
                            320320,	    #	6	320320	China and other dinnerware	
                            320330,	    #	6	320330	Flatware	
                            320340,	    #	6	320340	Glassware	
                            320350,	    #	6	320350	Silver serving pieces	
                            320360,	    #	6	320360	Other serving pieces	
                            320370,	    #	6	320370	Nonelectric cookware	
                            320521,	    #	6	320521	Small electric kitchen appliances	
                            320522,	    #	6	320522	Portable heating and cooling equipment	
                #	4	MISCHHEQ	Miscellaneous household equipment	
                            320120,	    #	5	320120	Window coverings	
                            320130,	    #	5	320130	Infants' equipment	
                            320150,	    #	5	320150	Outdoor equipment	
                            320210,	    #	5	320210	Clocks (thru Q20071)	
                            320220,	    #	5	320220	Lamps and lighting fixtures	
                            320230,	    #	5	320230	Other Household Decorative Items (deleted 81)	
                            320231,	    #	5	320231	Other household decorative items (thru Q20071)	
                            320232,	    #	5	320232	Telephones and accessories	
                            320233,	    #	5	320233	Clocks and other household decorative items (new UCC Q20072)
                            320410,	    #	5	320410	Lawn and garden equipment	
                            320420,	    #	5	320420	Power tools	
                            320901,	    #	5	320901	Office furniture for home use	
                            320902,	    #	5	320902	Hand tools	
                            320903,	    #	5	320903	Indoor plants; fresh flowers	
                            320904,	    #	5	320904	Closet and storage items	
                            340904,	    #	5	340904	Rental of furniture	
                            430130,	    #	5	430130	Luggage	
                            690110,	    #	5	690110	Computers for Non-Business Use (Not in 90,91,92)  (1982-1987)
                            690111,	    #	5	690111	Computers and computer hardware for non-business use [1988]	
                            690112,	    #	5	690112	Computer software and accessories for non-business use [1988]
                            690115,	    #	5	690115	Personal digital assistants [2005]	
                            690210,	    #	5	690210	Telephone answering devices	
                            690220,	    #	5	690220	Calculators (deleted 05)	
                            690230,	    #	5	690230	Business equipment for home use	
                            690241,	    #	5	690241	Smoke alarms (renter)	
                            690242,	    #	5	690242	Smoke alarms (owned home)	
                            690243,	    #	5	690243	Smoke alarms (owned vacation)	
                            690244,	    #	5	690244	Other household appliances (renter)	
                            690245]		#	5	690245	Other household appliances (owned home)	


# uccmapping['HOUSE_TEXTILE'] = [ #Household textiles
#                 #	4	HHTXTILE	Household textiles	
#                             280110,	    #	5	280110	Bathroom linens	
#                             280120,	    #	5	280120	Bedroom linens	
#                             280130,	    #	5	280130	Kitchen and dining room linens	
#                             280210,	    #	5	280210	Curtains and draperies	
#                             280220,	    #	5	280220	Slipcovers; decorative pillows	
#                             280230,	    #	5	280230	Sewing materials for slipcovers; curtains; other sewing materials for home
#                             280900]	    #	5	280900	Other linens


# uccmapping['FURNITURE_NIPA'] = [    	        #		FURNITURE	Furniture 	
	
#                 #	4	FURNITUR	Furniture	
#                             290110,	    #	5	290110	Mattress and springs	
#                             290120,	    #	5	290120	Other bedroom furniture	
#                             290210,	    #	5	290210	Sofas	
#                             290310,	    #	5	290310	Living room chairs	
#                             290320,	    #	5	290320	Living room tables	
#                             290410,	    #	5	290410	Kitchen; dining room furniture	
#                             290420,	    #	5	290420	Infants' furniture	
#                             290430,	    #	5	290430	Outdoor furniture	
#                             290440,	    #	5	290440	Wall units; cabinets and other occasional furniture	
#                 #	4	FLOORCOV	Floor coverings	
#                             230131,	    #	5	230131	Installed Wall-to-Wall Carpeting - Renter  (deleted 99)	
#                             230132,	    #	5	230132	Wall-to-wall carpet installed replacement owned (newly added UCC 230133)
#                             230133,	    #	5	230133	Wall-to-wall carpet (replacement) (owned home)	
#                             230134,	    #	6	230134	Wall-to-wall carpet (renter)	
#                             320110,	    #	5	320110	Non-permanent floor coverings is now included in UCC 320111 (deleted 99)
#                             320111,	    #	5	320111	Floor coverings; non-permanent	
#                             320161,	    #	5	320161	Non-Installed Wall to Wall Carpeting and Carpet Squares - Renter 
#                             320162,	    #	5	320162	Non-Installed Wall-to-Wall Carpeting and Carpet Squares - Homeowner 
#                             320163,	    #	6	320163	Wall-to-wall carpet (replacement)(renter) [1999]	
#                 #	4	MAJAPPL	Major appliances	
#                             230117,	    #	5	230117	Dishwashers (built-in); garbage disposals; range hoods; (renter)
#                             230118,	    #	5	230118	Dishwashers (built-in); garbage disposals; range hoods; (owned home)
#                             300111,	    #	5	300111	Refrigerators; freezers (renter)	
#                             300112,	    #	5	300112	Refrigerators; freezers (owned home)	
#                             300211,	    #	5	300211	Washing machines (renter)	
#                             300212,	    #	5	300212	Washing machines (owned home)	
#                             300221,	    #	5	300221	Clothes dryers (renter)	
#                             300222,	    #	5	300222	Clothes dryers (owned home)	
#                             300311,	    #	5	300311	Cooking stoves; ovens (renter)	
#                             300312,	    #	5	300312	Cooking stoves; ovens (owned home)	
#                             300321,	    #	5	300321	Microwave ovens (renter)	
#                             300322,	    #	5	300322	Microwave ovens (owned home)	
#                             300331,	    #	5	300331	Portable dishwasher (renter)	
#                             300332,	    #	5	300332	Portable dishwasher (owned home)	
#                             300411,	    #	5	300411	Window air conditioners (renter)	
#                             300412,	    #	5	300412	Window air conditioners (owned home)	
#                             320511,	    #	5	320511	Electric floor cleaning equipment	
#                             320512,	    #	5	320512	Sewing machines	
#                 #	4	SMAPPHWR	Small appliances; miscellaneous housewares	
#                             320310,	    #	6	320310	Plastic dinnerware	
#                             320320,	    #	6	320320	China and other dinnerware	
#                             320330,	    #	6	320330	Flatware	
#                             320340,	    #	6	320340	Glassware	
#                             320350,	    #	6	320350	Silver serving pieces	
#                             320360,	    #	6	320360	Other serving pieces	
#                             320370,	    #	6	320370	Nonelectric cookware	
#                             320521,	    #	6	320521	Small electric kitchen appliances	
#                             320522,	    #	6	320522	Portable heating and cooling equipment	
#                 #	4	MISCHHEQ	Miscellaneous household equipment	
#                             320120,	    #	5	320120	Window coverings	
#                             320130,	    #	5	320130	Infants' equipment	
#                             320150,	    #	5	320150	Outdoor equipment	
#                             320210,	    #	5	320210	Clocks (thru Q20071)	
#                             320220,	    #	5	320220	Lamps and lighting fixtures	
#                             320230,	    #	5	320230	Other Household Decorative Items (deleted 81)	
#                             320231,	    #	5	320231	Other household decorative items (thru Q20071)	
#                             320233,	    #	5	320233	Clocks and other household decorative items (new UCC Q20072)
#                             320410,	    #	5	320410	Lawn and garden equipment	
#                             320420,	    #	5	320420	Power tools	
#                             320901,	    #	5	320901	Office furniture for home use	
#                             320902,	    #	5	320902	Hand tools	
#                             320903,	    #	5	320903	Indoor plants; fresh flowers	
#                             320904,	    #	5	320904	Closet and storage items	
#                             340904,	    #	5	340904	Rental of furniture	
#                             690241,	    #	5	690241	Smoke alarms (renter)	
#                             690242,	    #	5	690242	Smoke alarms (owned home)	
#                             690243,	    #	5	690243	Smoke alarms (owned vacation)	
#                             690244,	    #	5	690244	Other household appliances (renter)	
#                             690245]		#	5	690245	Other household appliances (owned home)	

# uccmapping['INFOPROCESS'] = [        #Information processing equipment as defined in NIPA

#                             690110,	    #	5	690110	Computers for Non-Business Use (Not in 90,91,92)  (1982-1987)
#                             690111,	    #	5	690111	Computers and computer hardware for non-business use [1988]	
#                             690112,	    #	5	690112	Computer software and accessories for non-business use [1988]
#                             690115,	    #	5	690115	Personal digital assistants [2005]	
#                             690220,	    #	5	690220	Calculators (deleted 05)	
#                             690230,	    #	5	690230	Business equipment for home use	

#                             ]


# uccmapping['INTERNET'] = [        #Computer Information Services
#                             690114,	    #	5	690114	Computer information services	
#                             690116,	    #	5	690116	Internet services away from home [2005]	
#                          ]


uccmapping['CLOTHD']    = [    	        #		CLOTH	Clothing (durables)	
                #	4	MENS	Men; 16 and over	
                            360110,	    #	5	360110	Men's suits	
                            360120,	    #	5	360120	Men's sportcoats; tailored jackets	
                            360210,	    #	5	360210	Men's coats and jackets	
                            360311,	    #	5	360311	Men's underwear	
                            360312,	    #	5	360312	Men's hosiery	
                            360320,	    #	5	360320	Men's nightwear	
                            360330,	    #	5	360330	Men's accessories	
                            360340,	    #	5	360340	Men's sweaters and vests	
                            360350,	    #	5	360350	Men's active sportswear	
                            360410,	    #	5	360410	Men's shirts	
                            360511,	    #	5	360511	Men's pants (thru Q20071)	
                            360512,	    #	5	360512	Men's shorts; shorts sets (thru Q20071)	
                            360513,	    #	5	360513	Men's pants and shorts (new UCC Q20072)	
                            360901,	    #	5	360901	Men's uniforms	
                            360902,	    #	5	360902	Men's costumes	
                #	4	BOYS	Boys; 2 to 15	
                            370110,	    #	5	370110	Boys' coats and jackets	
                            370120,	    #	5	370120	Boys' sweaters	
                            370130,	    #	5	370130	Boys' shirts	
                            370211,	    #	5	370211	Boys' underwear	
                            370212,	    #	5	370212	Boys' nightwear	
                            370213,	    #	5	370213	Boys' hosiery	
                            370220,	    #	5	370220	Boys' accessories	
                            370311,	    #	5	370311	Boys' suits; sportcoats; vests	
                            370312,	    #	5	370312	Boys' pants (thru Q20071)	
                            370313,	    #	5	370313	Boys' shorts; shorts sets (thru Q20071)	
                            370314,	    #	5	370314	Boys' pants and shorts (new UCC Q20072)	
                            370901,	    #		370901	Boy's uniforms, active sports wear (split into UCC 370903 and 370904) 
                            370902,	    #	5	370902	Boys' costumes	
                            370903,	    #	5	370903	Boys' uniforms	
                            370904,	    #	5	370904	Boys' active sportswear	
                #	4	WOMENS	Women; 16 and over	
                            380110,	    #	5	380110	Women's coats and jackets	
                            380210,	    #	5	380210	Women's dresses	
                            380311,	    #	5	380311	Women's sportcoats; tailored jackets	
                            380312,	    #	5	380312	Women's vests and sweaters	
                            380313,	    #	5	380313	Women's shirts; tops; blouses	
                            380320,	    #	5	380320	Women's skirts	
                            380331,	    #	5	380331	Women's pants (thru Q20071)	
                            380332,	    #	5	380332	Women's shorts; shorts sets (thru Q20071)	
                            380333,	    #	5	380333	Women's pants and shorts (new UCC Q20072)	
                            380340,	    #	5	380340	Women's active sportswear	
                            380410,	    #	5	380410	Women's sleepwear	
                            380420,	    #	5	380420	Women's undergarments	
                            380430,	    #	5	380430	Women's hosiery	
                            380510,	    #	5	380510	Women's suits	
                            380901,	    #	5	380901	Women's accessories	
                            380902,	    #	5	380902	Women's uniforms	
                            380903,	    #	5	380903	Women's costumes	
                #	4	GIRLS	Girls; 2 to 15	
                            390110,	    #	5	390110	Girls' coats and jackets	
                            390120,	    #	5	390120	Girls' dresses and suits	
                            390210,	    #	5	390210	Girls' shirts; blouses; sweaters	
                            390221,	    #	5	390221	Girls' skirts and pants (thru Q20071)	
                            390222,	    #	5	390222	Girls' shorts; shorts sets (thru Q20071)	
                            390223,	    #	5	390223	Girls' skirts; pants; and shorts (new UCC Q20072)	
                            390230,	    #	5	390230	Girls' active sportswear	
                            390310,	    #	5	390310	Girls' underwear and sleepwear	
                            390321,	    #	5	390321	Girls' hosiery	
                            390322,	    #	5	390322	Girls' accessories	
                            390901,	    #	5	390901	Girls' uniforms	
                            390902,	    #	5	390902	Girls' costumes	
                #	3	INFANT	Children under 2	
                            410110,	    #	4	410110	Infant coat; jacket; snowsuit	
                            410111,	    #	4	410111	Infant Coats, Jackets, and Snowsuits 9B (deleted 91Q2)	
                            410112,	    #	4	410112	Infant Coats, Jackets, and Snowsuits 9A (deleted 91Q2)	
                            410120,	    #	4	410120	Infant dresses; outerwear	
                            410121,	    #	4	410121	Infant Dresses and Outerwear 9B (deleted 91Q2)	
                            410122,	    #	4	410122	Infant Dresses and Outerwear 9A (deleted 91Q2)	
                            410130,	    #	4	410130	Infant underwear	
                            410131,	    #	4	410131	Infant Undergarments 9B, Including Diapers (deleted 91Q2)	
                            410132,	    #	4	410132	Infant Undergarments 9A, Including Diapers (deleted 91Q2)	
                            410140,	    #	4	410140	Infant nightwear; loungewear	
                            410141,	    #	4	410141	Infant Sleeping Garments 9B (deleted 91Q2)	
                            410142,	    #	4	410142	Infant Sleeping Garments 9A (deleted 91Q2)	
                            410901,	    #	4	410901	Infant accessories	
                            410902,	    #	4	410902	Infants' other clothing (deleted 91)	
                            410903,	    #	4	410903	Infant Accessories 9A (deleted 91Q2)	
                            410904,	    #	4	410904	Infant Hosiery, Footwear, and Other Clothing (deleted 91Q2)	
                #	3	FOOTWEAR	Footwear	
                            400110,	    #	4	400110	Men's footwear	
                            400210,	    #	4	400210	Boys' footwear	
                            400220,	    #	4	400220	Girls' footwear	
                            400310,	    #	4	400310	Women's footwear	
                #	3	OTHAPPRL	Other apparel products	
                            420110,	    #	4	420110	Material for making clothes	
                            420120]		#	4	420120	Sewing patterns and notions	
                        
uccmapping['CLOTHS']    = [    	        #		CLOTH	Clothing (services)	
                #	3	OTHAPPRL	Other apparel products and services	
                            440110,	    #	4	440110	Shoe repair and other shoe service	
                            440120,	    #	4	440120	Coin-operated apparel laundry and dry cleaning	
                            440130,	    #	4	440130	Alteration; repair and tailoring of apparel and accessories	
                            440140,	    #	4	440140	Clothing rental	
                            440150,	    #	4	440150	Watch and jewelry repair	
                            440210,	    #	4	440210	Apparel laundry and dry cleaning not coin-operated	
                            440900]		#	4	440900	Clothing storage	
                        
uccmapping['JEWELRY']   = [    	        #		JEWELRY	Jewelry	
                            430110,	    #	4	430110	Watches	
                            430120]		#	4	430120	Jewelry	

uccmapping['NEWCARS']  = [    	        #	3	VEHPURCH	Vehicle purchases (net outlay)	
                                        #	4	NEWCARS	Cars and trucks; new	
                            450110,	    #	5	450110	New cars	
                                        #	5	450116	Trade-In Allowance For New Cars (deleted) 

                #  Note: 
                #   (1) "trade-in allowance"=A reduction in the price of a new item when an old item is given as part of the deal. 
                #   (2) Since the vehicle purchases are already 'net outlays' we do not have to substract the trade-in allowances.

                            450210]     #	5	450210	New trucks	
                                        #	5	450216	Trade-In Allowance For New Trucks or Vans (deleted)	  

uccmapping['NEWPCARS']  = [    	       	
                            450110]	    #	5	450110	New cars	
                                       

uccmapping['NEWTRUCKS']  = [    	       
                          450210]     #	5	450210	New trucks		

uccmapping['USEDCARS']  = [             #	4	USEDCARS	Cars and trucks; used	
                            460110,	    #	5	460110	Used cars	
                                        #	5	460116	Trade-In Allowance For Used Cars (deleted)	
                            460901]	    #	5	460901	Used trucks	
                                        #	5	460907	Trade-In Allowance For Used Trucks or Vans (deleted)    

uccmapping['OTHVEHCL']  = [             #	4	OTHVEHCL	Other vehicles	
                            450220,	    #	5	450220	New motorcycles	
                                        #	5	450226	Trade-In Allowance For New Motorcycles, Motor Scooters, or Mopeds (deleted)
                            450900,	    #	5	450900	New aircraft (deleted 82/3)	
                                        #	5	450906	Trade-In allowance new airplanes (deleted 82/3)	
                            460902,	    #	5	460902	Used motorcycles	
                                        #	5	460908	Trade-In Allowance For Used Motorcycles, Motor Scooters, or Mopeds (deleted)
                            460903]		#	5	460903	Used aircraft (deleted 82/3)	
                                        #	5	460909	Trade-In allowance used airplanes (deleted 82/3)                                                                                                         
                        
uccmapping['VEHPURCH']  = [    	        #	3	VEHPURCH	Vehicle purchases (net outlay)	
                                        #	4	NEWCARS	Cars and trucks; new	
                            450110,	    #	5	450110	New cars	
                                        #	5	450116	Trade-In Allowance For New Cars (deleted) 

                #  Note: 
                #   (1) "trade-in allowance"=A reduction in the price of a new item when an old item is given as part of the deal. 
                #   (2) Since the vehicle purchases are already 'net outlays' we do not have to substract the trade-in allowances.

                            450210,	    #	5	450210	New trucks	
                                        #	5	450216	Trade-In Allowance For New Trucks or Vans (deleted)	
                                        #	4	USEDCARS	Cars and trucks; used	
                            460110,	    #	5	460110	Used cars	
                                        #	5	460116	Trade-In Allowance For Used Cars (deleted)	
                            460901,	    #	5	460901	Used trucks	
                                        #	5	460907	Trade-In Allowance For Used Trucks or Vans (deleted)	
                                        #	4	OTHVEHCL	Other vehicles	
                            450220,	    #	5	450220	New motorcycles	
                                        #	5	450226	Trade-In Allowance For New Motorcycles, Motor Scooters, or Mopeds (deleted)
                            450900,	    #	5	450900	New aircraft (deleted 82/3)	
                                        #	5	450906	Trade-In allowance new airplanes (deleted 82/3)	
                            460902,	    #	5	460902	Used motorcycles	
                                        #	5	460908	Trade-In Allowance For Used Motorcycles, Motor Scooters, or Mopeds (deleted)
                            460903]		#	5	460903	Used aircraft (deleted 82/3)	
                                        #	5	460909	Trade-In allowance used airplanes (deleted 82/3)	
                
                
# uccmapping['VEHPURCHBEA']   = [    	        #	3	VEHPURCHBEA	Vehicle purchases (net outlay) to align with BEA definitions	
#                 #	4	NEWCARS	Cars and trucks; new	
#                             450110,	    #	5	450110	New cars	
#                 #	5	450116	Trade-In Allowance For New Cars (deleted) 

#                 #  Note: 
#                 #   (1) "trade-in allowance"=A reduction in the price of a new item when an old item is given as part of the deal. 
#                 #   (2) Since the vehicle purchases are already 'net outlays' we do not have to substract the trade-in allowances.


#                             450210,	    #	5	450210	New trucks	
#                 #	5	450216	Trade-In Allowance For New Trucks or Vans (deleted)	
#                 #	4	USEDCARS	Cars and trucks; used	
#                             460110,	    #	5	460110	Used cars	
#                 #	5	460116	Trade-In Allowance For Used Cars (deleted)	
#                             460901]	 	#	5	460901	Used trucks	
#                 #	5	460907	Trade-In Allowance For Used Trucks or Vans (deleted)	
                
uccmapping['VEHLEASE']  = [    	        #	6	VEHLEASE	Leased vehicles	
                            450310,	    #	7	450310	Car lease payments	
                            450311,	    #	7	450311	Charges Other Than Basic Lease, Such as Insurance or Maintenance (Car Lease)
                #	7	450312	Trade-In Allowance (Car Lease) (deleted)	
                            450313,	    #	7	450313	Cash down-payment (car lease)	
                            450314,	    #	7	450314	Termination fee (car lease)	
                            450410,	    #	7	450410	Truck lease payments	
                            450411,	    #	7	450411	Charges Other Than Basic Lease (Truck/Van Lease) 
                #	7	450412	Trade-In Allowance (Truck/Van Lease) (deleted)	
                            450413,	    #	7	450413	Cash down-payment (truck lease)	
                            450414]	 	#	7	450414	Termination fee (truck lease)	

uccmapping['VEHRENT']   = [    	        #	6	RENTVEH	Rented vehicles	
                            520511,	    #	7	520511	Auto rental	
                            520512,	    #	7	520512	Auto rental; out-of-town trips	
                            520521,	    #	7	520521	Truck rental	
                            520522,	    #	7	520522	Truck rental; out-of-town trips	
                            520902,	    #	7	520902	Motorcycle rental	
                            520903,	    #	7	520903	Aircraft rental	
                            520905,	    #	7	520905	Motorcycle rental; out-of-town trips	
                            520906]	    #	7	520906	Aircraft rental, out of town trips (deleted 06)	

uccmapping['VEHINSUR']  = [             #		VEHINS		
                            500110]	    #	4	500110	Vehicle insurance	


uccmapping['VEHFINCH']  = [             #	4	VEHFINCH	Vehicle finance charges	
                            510110,	    #	5	510110	Automobile finance charges	
                            510901,	    #	5	510901	Truck finance charges	
                            510902,	    #	5	510902	Motorcycle and plane finance charges	
                            850300]	    #	5	850300	Other vehicle finance charges	                                     

uccmapping['VEHOTHCH']  = [   
                            520110,	    #	5	520110	State or local vehicle registration	
                            520111,	    #	5	520111	State vehicle registration (deleted 06)	
                            520112,	    #	5	520112	Local vehicle registration (deleted 06)	
                            520310,	    #	5	520310	Drivers' license	
                            520410,	    #	5	520410	Vehicle inspection	
     	                                #	5	PARKING	Parking fees	
                            520530,	    #	6	520530	Parking fees, incl. garages, meters, excl. costs of property ownership
                            520531,	    #	6	520531	Parking fees in home city; excluding residence	
                            520532,	    #	6	520532	Parking fees; out-of-town trips	
                            520541,	    #	5	520541	Tolls or electronic toll passes	
                            520542,	    #	5	520542	Tolls on out-of-town trips	
                            520550,	    #	5	520550	Towing charges	
                            520560,	    #	5	520560	Global positioning services	
                            620113]		#	5	620113	Automobile service clubs  

uccmapping['CAREPAIR']  = [    	        #	4	CAREPAIR	Maintenance and repairs	
                            470220,	    #	5	470220	Coolant; brake fluid; transmission fluid; and other additives	
                            480110,	    #	5	480110	Tires - purchased; replaced; installed	
                            480211,	    #	5	480211	Vehicle equipment installed by CU (split into UCC's 480213 and 480214)	
                            480212,	    #	5	480212	Vehicle products and cleaning services	
                            480213,	    #	5	480213	Parts; equipment; and accessories	
                            480214,	    #	5	480214	Vehicle audio equipment	
                            480215,		#	5	480215	Vehicle video equipment	
                            490110,	    #	5	490110	Body work and painting	
                            490211,	    #	5	490211	Clutch; transmission repair	
                            490212,	    #	5	490212	Drive shaft and rear-end repair	
                            490220,	    #	5	490220	Brake work (490220 and 490315 are now collapsed and collected under 490221)
                            490221,	    #	5	490221	Brake work; including adjustments	
                            490231,	    #	5	490231	Repair to steering or front-end	
                            490232,	    #	5	490232	Repair to engine cooling system	
                            490311,	    #	5	490311	Motor tune-up	
                            490312,	    #	5	490312	Lube; oil change; and oil filters	
                            490313,	    #	5	490313	Front-end alignment; wheel balance and rotation	
                            490314,	    #	5	490314	Shock absorber replacement	
                            490315,	    #	5	490315	Brake adjustment (490220 and 490315 now collapsed and collected in 490221)
                            490317,	    #	5	490317	Minor vehicle repairs and service out-of-town trips (deleted 91)
                            490318,	    #	5	490318	Repair tires and other repair work	
                            490319,	    #	5	490319	Vehicle air conditioning repair	
                            490411,	    #	5	490411	Exhaust system repair	
                            490412,	    #	5	490412	Electrical system repair	
                            490413,	    #	5	490413	Motor repair; replacement	
                            490500,	    #	5	490500	Purchase and installation of vehicle accessories (split 490501 and 490502)
                            490501,	    #	5	490501	Vehicle accessories including labor	
                            490502,	    #	5	490502	Vehicle audio equipment including labor (thru Q20051)	
                            490900]	    #	5	490900	Auto repair service policy                                                                                             

# uccmapping['VEHEXPD']   = [             #	4	CAREPAIR	Maintenance and repairs	
#                             470220,	    #	5	470220	Coolant; brake fluid; transmission fluid; and other additives	
#                             480110,	    #	5	480110	Tires - purchased; replaced; installed	
#                             480211,	    #	5	480211	Vehicle equipment installed by CU (split into UCC's 480213 and 480214)	
#                             480212,	    #	5	480212	Vehicle products and cleaning services	
#                             480213,	    #	5	480213	Parts; equipment; and accessories	
#                             480214,	    #	5	480214	Vehicle audio equipment	
#                             480215]		#	5	480215	Vehicle video equipment	
        
# uccmapping['VEHEXPS']   = [    	        #		VEHEXPS	Vehicle expenditures (services)	
#                 #	4	CAREPAIR	Maintenance and repairs	
#                             490110,	    #	5	490110	Body work and painting	
#                             490211,	    #	5	490211	Clutch; transmission repair	
#                             490212,	    #	5	490212	Drive shaft and rear-end repair	
#                             490220,	    #	5	490220	Brake work (490220 and 490315 are now collapsed and collected under 490221)
#                             490221,	    #	5	490221	Brake work; including adjustments	
#                             490231,	    #	5	490231	Repair to steering or front-end	
#                             490232,	    #	5	490232	Repair to engine cooling system	
#                             490311,	    #	5	490311	Motor tune-up	
#                             490312,	    #	5	490312	Lube; oil change; and oil filters	
#                             490313,	    #	5	490313	Front-end alignment; wheel balance and rotation	
#                             490314,	    #	5	490314	Shock absorber replacement	
#                             490315,	    #	5	490315	Brake adjustment (490220 and 490315 now collapsed and collected in 490221)
#                             490317,	    #	5	490317	Minor vehicle repairs and service out-of-town trips (deleted 91)
#                             490318,	    #	5	490318	Repair tires and other repair work	
#                             490319,	    #	5	490319	Vehicle air conditioning repair	
#                             490411,	    #	5	490411	Exhaust system repair	
#                             490412,	    #	5	490412	Electrical system repair	
#                             490413,	    #	5	490413	Motor repair; replacement	
#                             490500,	    #	5	490500	Purchase and installation of vehicle accessories (split 490501 and 490502)
#                             490501,	    #	5	490501	Vehicle accessories including labor	
#                             490502,	    #	5	490502	Vehicle audio equipment including labor (thru Q20051)	
#                             490900,	    #	5	490900	Auto repair service policy	
#                 #		VEHINS		
#                             500110,	    #	4	500110	Vehicle insurance	
#                 #	4	VEHFINCH	Vehicle finance charges	
#                             510110,	    #	5	510110	Automobile finance charges	
#                             510901,	    #	5	510901	Truck finance charges	
#                             510902,	    #	5	510902	Motorcycle and plane finance charges	
#                             850300,	    #	5	850300	Other vehicle finance charges	
#                 #	6	RENTVEH	Rented vehicles	
#                             520511,	    #	7	520511	Auto rental	
#                             520512,	    #	7	520512	Auto rental; out-of-town trips	
#                             520521,	    #	7	520521	Truck rental	
#                             520522,	    #	7	520522	Truck rental; out-of-town trips	
#                             520902,	    #	7	520902	Motorcycle rental	
#                             520903,	    #	7	520903	Aircraft rental	
#                             520905,	    #	7	520905	Motorcycle rental; out-of-town trips	
#                             520906,	    #	7	520906	Aircraft rental, out of town trips (deleted 06)	
#                 #	6	LEASVEH	Leased vehicles	
#                             450310,	    #	7	450310	Car lease payments	
#                             450311,	    #	7	450311	Charges Other Than Basic Lease, Such as Insurance or Maintenance (Car Lease)
#                 #	7	450312	Trade-In Allowance (Car Lease) (deleted)	
#                             450313,	    #	7	450313	Cash down-payment (car lease)	
#                             450314,	    #	7	450314	Termination fee (car lease)	
#                             450410,	    #	7	450410	Truck lease payments	
#                             450411,	    #	7	450411	Charges Other Than Basic Lease (Truck/Van Lease) 
#                 #	7	450412	Trade-In Allowance (Truck/Van Lease) (deleted)	
#                             450413,	    #	7	450413	Cash down-payment (truck lease)	
#                             450414,	    #	7	450414	Termination fee (truck lease)	
#                             520110,	    #	5	520110	State or local vehicle registration	
#                             520111,	    #	5	520111	State vehicle registration (deleted 06)	
#                             520112,	    #	5	520112	Local vehicle registration (deleted 06)	
#                             520310,	    #	5	520310	Drivers' license	
#                             520410,	    #	5	520410	Vehicle inspection	
#                 #	5	PARKING	Parking fees	
#                             520530,	    #	6	520530	Parking fees, incl. garages, meters, excl. costs of property ownership
#                             520531,	    #	6	520531	Parking fees in home city; excluding residence	
#                             520532,	    #	6	520532	Parking fees; out-of-town trips	
#                             520541,	    #	5	520541	Tolls or electronic toll passes	
#                             520542,	    #	5	520542	Tolls on out-of-town trips	
#                             520550,	    #	5	520550	Towing charges	
#                             520560,	    #	5	520560	Global positioning services	
#                             620113]		#	5	620113	Automobile service clubs                            
                        
uccmapping['GASOLINE']  = [    	        #		GASOLINE	Gasoline	
                #	3	GASOIL	Gasoline and motor oil	
                            470111,	    #	4	470111	Gasoline	
                            470112,	    #	4	470112	Diesel fuel	
                            470113,	    #	4	470113	Gasoline on out-of-town trips	
                            470211,	    #	4	470211	Motor oil	
                            470212]		#	4	470212	Motor oil on out-of-town trips	
                        
	
                        
uccmapping['PUBTRANS']  = [    	        #	3	PUBTRANS	Public transportation	
                            530110,	    #	4	530110	Airline fares	
                            530210,	    #	4	530210	Intercity bus fares	
                            530311,	    #	4	530311	Intercity mass transit fares	
                            530312,	    #	4	530312	Local trans. on out-of-town trips	
                            530411,	    #	4	530411	Taxi fares and limousine services on trips	
                            530412,	    #	4	530412	Taxi fares and limousine services	
                            530510,	    #	4	530510	Intercity train fares	
                            530901,	    #	4	530901	Ship fares	
                            530902]	    #	4	530902	School bus	
                                        #	2	HEALTH	Healthcare	
                        
uccmapping['HEALTHINS'] = [    	        #		HEALTHINS	Health insurance	
                                        #	4	COMHLTIN	Commercial health insurance	
                            580110,	    #	5	580110	Commercial health insurance (split into two new UCCs - 580111 and 580113)
                            580111,	    #	5	580111	Traditional fee for service health plan (not BCBS)	
                            580113,	    #	5	580113	Preferred provider health plan (not BCBS)	
                                        #	4	BCBS	Blue Cross; Blue Shield	
                            580112,	    #	5	580112	Traditional fee for service health plan (BCBS)	
                            580114,	    #	5	580114	Preferred provider health plan (BCBS)	
                            580210,	    #	5	580210	Blue Cross/Blue Shield (divided into 580112, 580114, 580312, 580904,580906)
                            580310,	    #	5	580310	Health maintenance plans (now collected under 580311) (deleted 95)
                            580311,	    #	4	580311	Health maintenance organization (not BCBS)	
                            580312,	    #	5	580312	Health maintenance organization (BCBS)	
                            580901,	    #	4	580901	Medicare payments	
                            580904,	    #	5	580904	Commercial medicare supplement (BCBS)	
                            580906,	    #	5	580906	Other health insurance (BCBS)	
                            580907,	    #	4	580907	Medicare prescription drug premium (new UCC Q20062)	
                                        #	4	COMEDOTH	Commercial medicare supplements and other health insurance	
                            580400,	    #	4	580400	Long term care insurance	
                            580902,	    #	5	580902	Commercial Medicare supplement and other health insurance (580903 and 580905)
                            580903,	    #	5	580903	Commercial medicare supplement (not BCBS)	
                            580905]		#	5	580905	Other health insurance (not BCBS)	
                        
                                        #		HEALTHEXP		
uccmapping['HEALTHEXPD']    = [    	    #		HEALTHEXPD	Health expenditures (durables)	
                                        #	3	MEDSUPPL	Medical supplies	
                            550110,	    #	4	550110	Eyeglasses and contact lenses	
                            550320,	    #	4	550320	Medical equipment for general use	
                            550330,	    #	4	550330	Supportive and convalescent medical equipment	
                            550340]		#	4	550340	Hearing aids	

uccmapping['DRUGS'] = [ #Presecription and non-prescription drugs
                            540000,	    #	3	540000	Prescription drugs	
                            580907]    #	4	580907	Medicare prescription drug premium (new UCC Q20062)	


uccmapping['HEALTHEXPS']    = [    	        #		HEALTHEXPS	Health expenditures (services)	
                #	3	MEDSUPPL	Medical supplies	
                            570901,	    #	4	570901	Rental of medical equipment	
                            570903,	    #	4	570903	Rental of supportive; convalescent medical equipment	
                #	3	MEDSERVS	Medical services	
                            540000,	    #	3	540000	Prescription drugs	
                            560110,	    #	4	560110	Physician's services	
                            560210,	    #	4	560210	Dental services	
                            560310,	    #	4	560310	Eye-care services	
                            560320,	    #	4	560320	Medical services by practitioner other than physician (deleted 91)
                            560330,	    #	4	560330	Lab tests; x-rays	
                            560400,	    #	4	560400	Service by professionals other than physician	
                            560900,	    #	4	560900	Nursing services and therapeutic treatments (deleted 91)	
                            570110,	    #	4	570110	Hospital room (deleted 05)	
                            570111,	    #	4	570111	Hospital room and services	
                            570210,	    #	4	570210	Hospital service other than room (deleted 05)	
                            570220,	    #	4	570220	Care in convalescent or nursing home	
                            570230,	    #	4	570230	Other medical care services	
                            570240]		#	4	570240	Medical care in retirement community	



uccmapping['FEESADM']   = [    	        #		ENTRTAINS	Entertainment (services)	
                                        #	3	FEESADM	Fees and admissions	
                            610900,	    #	4	610900	Recreation expenses; out-of-town trips	
                            620110,	    #	4	620110	Country clubs, swimming pools, fraternal org. (deleted 92)	
                            620111,	    #	4	620111	Social; recreation; health club membership	
                            620121,	    #	4	620121	Fees for participant sports	
                            620122,	    #	4	620122	Participant sports; out-of-town trips	
                            620211,	    #	4	620211	Movie; theater; amusement parks; and other	
                            620212,	    #	4	620212	Movie; other admissions; out-of-town trips	
                            620221,	    #	4	620221	Admission to sporting events	
                            620222,	    #	4	620222	Admission to sports events; out-of-town trips	
                            620310,	    #	4	620310	Fees for recreational lessons	
                            620903]	    #	4	620903	Other entertainment services; out-of-town trips	





# uccmapping['TVAUDIOD']   = [             #	3	TVAUDIO	Audio and visual equipment 
#                             310311,	    #	4	310311	Radios	
#                             310312,	    #	4	310312	PHONOGRAPHS This UCC was deleted as an interview UCC (now UCC 310333) 
#                             310313,	    #	4	310313	Tape recorders and players	
#                             310316,     #   4   310316  Stereos, radios, speakers, and sound components including those in vehicles 
#                             310110,	    #	4	310110	Black and white TV (deleted 05)	
#                             310120,	    #	4	310120	Color RV console (deleted 05)	
#                             310130,	    #	4	310130	Color TV portable/ table model (deleted 05)	
#                             310140,	    #	4	310140	Televisions	
#                             310331,     #   4   310331  Miscellaneous sound equipment                            
#                             310332,     #   4   310332  Sound equipment accessories	
#                             310210,	    #	4	310210	VCR's and video disc players	
#                             310315,     #   4   310315  Digital media players and recorders (new UCC Q20071)
#                             310230,	    #	4	310230	Video game hardware and software	
#                             310231,	    #	4	310231  Video game software 
#                             310232,	    #	4   310232  Video game hardware and accessories  	
#                             310220,	    #	4	310220	Video cassettes; tapes; and discs	
#                             310240,	    #	4	310240	Streaming; downloading video	
#                             310314,	    #	4	310314	Personal digital audio players	
#                             310320,	    #	4	310320	Sound components and component systems	
#                             310330,	    #	4	310330	Accessories and other sound equipment (deleted 96)	
#                             310333,	    #	4	310333	Accessories and other sound equipment	
#                             310334,	    #	4	310334	Satellite dishes	
#                             310340,	    #	4	310340	CDs; records; audio tapes	
#                             310341,	    #	4	310341	Record/tape/CD/video mail order club (deleted 05)
#                             310342,	    #	4	310342	Records, CDs, audio tapes, needles (deleted 05)	
#                             310350,	    #	4	310350	Streaming; downloading audio (2005)	
#                             310400,     #   4   310400  Applications, games, ringtones for handheld devices
#                             610130,	    #	4	610130	Musical instruments and accessories	  	
#                             ]   

# uccmapping['TVAUDIOS']   = [             #	3	TVAUDIOS	Audio and visual services		
#                             270310,	    #	4	270310	Cable and satellite television services	                   
#                             270311,	    #	4	270311	Satellite radio service (new UCC Q20072)	
#                             310243,     #   4   310243  Rental, streaming, downloading video 
#                             340610,	    #	4	340610	Repair of TV; radio; and sound equipment	
#                             340902,	    #	4	340902	Rental of televisions	
#                             340905,	    #	4	340905	Rental of VCR; radio; and sound equipment
#                             620904,	    #	4	620904	Rental and repair of musical instruments	
#                             620912,	    #	4	620912	Rental of video cassettes; tapes; films; and discs (82)	
#                             620916,	    #	4	620916	Rental of computer and video game hardware and software (new UCC Q20072)     
#                             620917,     #   4   620917  Rental of video hardware/accessories 
#                             620918,     #   4   620918  Rental of video software    
#                             690320,	    #	4	690320	Installation of televisions (new UCC Q20072)	
#                             690330,	    #	4	690330	Installation of satellite television equipment (new UCC Q20072)	
#                             690340,	    #	4	690340	Installation of sound systems (new UCC Q20072)	
#                             690350,	    #	4	690350	Installation of other video equipment of sound systems (new UCC Q20072)]  
#                             ]   

uccmapping['TVAUDIO']   = [             #	3	TVAUDIO	Audio and visual equipment and services	(= two previous categories)
                            310311,	    #	4	310311	Radios	
                            310312,	    #	4	310312	PHONOGRAPHS This UCC was deleted as an interview UCC (now UCC 310333) 
                            310313,	    #	4	310313	Tape recorders and players	
                            310316,     #   4   310316  Stereos, radios, speakers, and sound components including those in vehicles 
                            310110,	    #	4	310110	Black and white TV (deleted 05)	
                            310120,	    #	4	310120	Color RV console (deleted 05)	
                            310130,	    #	4	310130	Color TV portable/ table model (deleted 05)	
                            310140,	    #	4	310140	Televisions	
                            270310,	    #	4	270310	Cable and satellite television services	
                            310331,     #   4   310331  Miscellaneous sound equipment                            
                            270311,	    #	4	270311	Satellite radio service (new UCC Q20072)	
                            310332,     #   4   310332  Sound equipment accessories	
                            310210,	    #	4	310210	VCR's and video disc players	
                            310315,     #   4   310315  Digital media players and recorders (new UCC Q20071)
                            310230,	    #	4	310230	Video game hardware and software	
                            310231,	    #	4	310231  Video game software 
                            310232,	    #	4   310232  Video game hardware and accessories  	
                            310220,	    #	4	310220	Video cassettes; tapes; and discs	
                            310240,	    #	4	310240	Streaming; downloading video	
                            310243,     #   4   310243  Rental, streaming, downloading video 
                            340610,	    #	4	340610	Repair of TV; radio; and sound equipment	
                            340902,	    #	4	340902	Rental of televisions	
                            310314,	    #	4	310314	Personal digital audio players	
                            310320,	    #	4	310320	Sound components and component systems	
                            310330,	    #	4	310330	Accessories and other sound equipment (deleted 96)	
                            310333,	    #	4	310333	Accessories and other sound equipment	
                            310334,	    #	4	310334	Satellite dishes	
                            310340,	    #	4	310340	CDs; records; audio tapes	
                            310341,	    #	4	310341	Record/tape/CD/video mail order club (deleted 05)
                            310342,	    #	4	310342	Records, CDs, audio tapes, needles (deleted 05)	
                            310350,	    #	4	310350	Streaming; downloading audio (2005)	
                            310400,     #   4   310400  Applications, games, ringtones for handheld devices
                            340905,	    #	4	340905	Rental of VCR; radio; and sound equipment
                            610130,	    #	4	610130	Musical instruments and accessories	  	
                            620904,	    #	4	620904	Rental and repair of musical instruments	
                            620912,	    #	4	620912	Rental of video cassettes; tapes; films; and discs (82)	
                            620916,	    #	4	620916	Rental of computer and video game hardware and software (new UCC Q20072)     
                            620917,     #   4   620917  Rental of video hardware/accessories 
                            620918,     #   4   620918  Rental of video software    
                            690320,	    #	4	690320	Installation of televisions (new UCC Q20072)	
                            690330,	    #	4	690330	Installation of satellite television equipment (new UCC Q20072)	
                            690340,	    #	4	690340	Installation of sound systems (new UCC Q20072)	
                            690350,	    #	4	690350	Installation of other video equipment of sound systems (new UCC Q20072)]  
                            ]   

uccmapping['PETSPLAY']  = [             #	3	PETSPLAY	Pets; toys; hobbies; and playground equipment	
                            610110,	    #	4	610110	Toys; games; arts and crafts; and tricycles	
                            610120,	    #	4	610120	Playground equipment	
                            610140,	    #	4	610140	Stamp and coin collecting [2004]	
                            610320,	    #	5	610320	Pet purchase; supplies; medicine
                            620410,	    #	5	620410	Pet services	
                            620420,	    #	5	620420	Vet services		
                            ]  


# uccmapping['PET_SERVICE']  = [             #	3	PETS	Pet services	
#                             620410,	    #	5	620410	Pet services	
#                             620420,	    #	5	620420	Vet services		
#                             ]  

# uccmapping['ENTEROTHD']  = [            #   Other Entertainment Equipment
#                                         #	4	UNMTRBOT	Unmotored recreational vehicles	
#                             600121,	    #	5	600121	Boat without motor and boat trailers	
#                                         #	5	600127	Trade-In Allowance For Boat Without Motor or Non Camper-Type Trailer 
#                             600122,	    #	5	600122	Trailer and other attachable campers	
#                                         #	5	600128	Trade-In Allowance For Trailer-Type or Other Attachable-Type Camper
#                                         #	4	PWRSPVEH	Motorized recreational vehicles	
#                             600131,	    #	5	600131	Motorized camper-coach or other vehicles (net outlay) (600141 and 600142)
#                                         #	5	600137	Trade-in allowance for motorized camper-coach (now 600143 and 600144)	
#                             600132,	    #	5	600132	Purchase of boat with motor	
#                                         #	5	600138	Trade-In Allowance For Boat With Motor (deleted)	
#                             600141,	    #	5	600141	Purchase of motorized camper	
#                                         #	5	600143	Trade-In allowance, motorized camper (deleted)	
#                             600142,	    #	5	600142	Purchase of other vehicle	
#                                         #	5	600144	Trade-In allowance, other vehicle (deleted)	
#                                         #	4	RECEQUIP	Sports; recreation and exercise equipment	
#                             600210,	    #	5	600210	Athletic gear; game tables; and exercise equipment	
#                             600310,	    #	5	600310	Bicycles	
#                             600410,	    #	5	600410	Camping equipment	
#                             600420,	    #	5	600420	Hunting and fishing equipment	
#                             600430,	    #	5	600430	Winter sports equipment	
#                             600900,	    #	5	600900	Water and miscellaneous sports equipment (split into 600901 and 600902)
#                             600901,	    #	5	600901	Water sports equipment	
#                             600902,	    #	5	600902	Other sports equipment                                                   
#                                         #	4	PHOTOEQ	Photographic equipment; supplies and services	
#                             610210,	    #	5	610210	Film	
#                             610230]		#	5	610230	Photographic equipment

# uccmapping['ENTEROTHS']  = [            #   Other Entertainment Services
#                                         #	4	RNTSPVEH	Rental of recreational vehicles	
#                             520901,	    #	4	520901	Docking and landing fees	
#                             520904,	    #	5	520904	Rental non-camper trailer	
#                             520907,	    #	5	520907	Boat and trailer rental out-of-town trips	
#                             600110,	    #	4	600110	Outboard motors	
#                             620902,	    #	4	620902	Rental of campers and other vehicles on trips (deleted 91)	
#                             620906,	    #	5	620906	Rental of boat	
#                             620907,	    #	5	620907	Rental of all campers and other recreational vehicles (620921 and 620922)
#                             620909,	    #	5	620909	Rental of campers on out-of-town trips	
#                             620919,	    #	5	620919	Rental of other vehicles on out-of-town trips	
#                             620921,	    #	5	620921	Rental of motorized camper	
#                             620922,	    #	5	620922	Rental of other RV"s	
#                             600311,     #   5   600311  Bike sharing, E-scooters (new UCC Q20192)        
#                             620908,	    #	5	620908	Rental and repair of miscellaneous sports equipment	
#                             620930,	    #	4	620930	Online gaming services (2005)
#                                         #       [ONLINE GAMING: per dictionary included in other entertainment, but hierarchical listing says its in TVAUDIO]
#                                         #	4	PHOTOEQ	Photographic equipment; supplies and services	
#                             620320,	    #	5	620320	Photographer fees	
#                             620330,	    #	5	620330	Photo processing	
#                             620905,	    #	5	620905	Repair and rental of photographic equipment	
#                             680310,	    #	4	680310	Live entertainment for catered affairs (new UCC Q20072)	
#                             680320]		#	4	680320	Rental of party supplies for catered affairs (new UCC Q20072)	   	 	       

uccmapping['ENTEROTH']  = [            #   Other Entertainment (= to previous two categories)
                                        #	4	UNMTRBOT	Unmotored recreational vehicles	
                            600121,	    #	5	600121	Boat without motor and boat trailers	
                                        #	5	600127	Trade-In Allowance For Boat Without Motor or Non Camper-Type Trailer 
                            600122,	    #	5	600122	Trailer and other attachable campers	
                                        #	5	600128	Trade-In Allowance For Trailer-Type or Other Attachable-Type Camper
                                        #	4	PWRSPVEH	Motorized recreational vehicles	
                            600131,	    #	5	600131	Motorized camper-coach or other vehicles (net outlay) (600141 and 600142)
                                        #	5	600137	Trade-in allowance for motorized camper-coach (now 600143 and 600144)	
                            600132,	    #	5	600132	Purchase of boat with motor	
                                        #	5	600138	Trade-In Allowance For Boat With Motor (deleted)	
                            600141,	    #	5	600141	Purchase of motorized camper	
                                        #	5	600143	Trade-In allowance, motorized camper (deleted)	
                            600142,	    #	5	600142	Purchase of other vehicle	
                                        #	5	600144	Trade-In allowance, other vehicle (deleted)	
                                        #	4	RNTSPVEH	Rental of recreational vehicles	
                            520901,	    #	4	520901	Docking and landing fees	
                            520904,	    #	5	520904	Rental non-camper trailer	
                            520907,	    #	5	520907	Boat and trailer rental out-of-town trips	
                            600110,	    #	4	600110	Outboard motors	
                            620902,	    #	4	620902	Rental of campers and other vehicles on trips (deleted 91)	
                            620906,	    #	5	620906	Rental of boat	
                            620907,	    #	5	620907	Rental of all campers and other recreational vehicles (620921 and 620922)
                            620909,	    #	5	620909	Rental of campers on out-of-town trips	
                            620919,	    #	5	620919	Rental of other vehicles on out-of-town trips	
                            620921,	    #	5	620921	Rental of motorized camper	
                            620922,	    #	5	620922	Rental of other RV"s	
                                        #	4	RECEQUIP	Sports; recreation and exercise equipment	
                            600210,	    #	5	600210	Athletic gear; game tables; and exercise equipment	
                            600310,	    #	5	600310	Bicycles	
                            600311,     #   5   600311  Bike sharing, E-scooters (new UCC Q20192)
                            600410,	    #	5	600410	Camping equipment	
                            600420,	    #	5	600420	Hunting and fishing equipment	
                            600430,	    #	5	600430	Winter sports equipment	
                            600900,	    #	5	600900	Water and miscellaneous sports equipment (split into 600901 and 600902)
                            600901,	    #	5	600901	Water sports equipment	
                            600902,	    #	5	600902	Other sports equipment            
                            620908,	    #	5	620908	Rental and repair of miscellaneous sports equipment	
                            620930,	    #	4	620930	Online gaming services (2005)
                                        #       [ONLINE GAMING: per dictionary included in other entertainment, but hierarchical listing says its in TVAUDIO]
                                        #	4	PHOTOEQ	Photographic equipment; supplies and services	
                            610210,	    #	5	610210	Film	
                            610230,		#	5	610230	Photographic equipment	 
                            620320,	    #	5	620320	Photographer fees	
                            620330,	    #	5	620330	Photo processing	
                            620905,	    #	5	620905	Repair and rental of photographic equipment	
                            680310,	    #	4	680310	Live entertainment for catered affairs (new UCC Q20072)	
                            680320]		#	4	680320	Rental of party supplies for catered affairs (new UCC Q20072)	                                                              

# uccmapping['ENTRTAINS'] = [    	        #		ENTRTAINS	Entertainment (services)	
#                                         #	3	FEESADM	Fees and admissions	
#                             610900,	    #	4	610900	Recreation expenses; out-of-town trips	
#                             620110,	    #	4	620110	Country clubs, swimming pools, fraternal org. (deleted 92)	
#                             620111,	    #	4	620111	Social; recreation; health club membership	
#                             620121,	    #	4	620121	Fees for participant sports	
#                             620122,	    #	4	620122	Participant sports; out-of-town trips	
#                             620211,	    #	4	620211	Movie; theater; amusement parks; and other	
#                             620212,	    #	4	620212	Movie; other admissions; out-of-town trips	
#                             620221,	    #	4	620221	Admission to sporting events	
#                             620222,	    #	4	620222	Admission to sports events; out-of-town trips	
#                             620310,	    #	4	620310	Fees for recreational lessons	
#                             620903,	    #	4	620903	Other entertainment services; out-of-town trips	
#                                         #	3	TVAUDIO	Audio and visual equipment and services	
#                             270310,	    #	4	270310	Cable and satellite television services	
#                             270311,	    #	4	270311	Satellite radio service (new UCC Q20072)	
#                             310240,	    #	4	310240	Streaming; downloading video	
#                             310341,	    #	4	310341	Record/tape/CD/video mail order club (deleted 05)	
#                             310350,	    #	4	310350	Streaming; downloading audio (2005)	
#                             340610,	    #	4	340610	Repair of TV; radio; and sound equipment	
#                             340902,	    #	4	340902	Rental of televisions	
#                             340905,	    #	4	340905	Rental of VCR; radio; and sound equipment	
#                             620904,	    #	4	620904	Rental and repair of musical instruments	
#                             620912,	    #	4	620912	Rental of video cassettes; tapes; films; and discs (82)	
#                             620916,	    #	4	620916	Rental of computer and video game hardware and software (new UCC Q20072)
#                             620930,	    #	4	620930	Online gaming services (2005)	
#                             690320,	    #	4	690320	Installation of televisions (new UCC Q20072)	
#                             690330,	    #	4	690330	Installation of satellite television equipment (new UCC Q20072)	
#                             690340,	    #	4	690340	Installation of sound systems (new UCC Q20072)	
#                             690350,	    #	4	690350	Installation of other video equipment of sound systems (new UCC Q20072)
#                                         #	3	PETSPLAY	Pets; toys; hobbies; and playground equipment	
#                             620410,	    #	5	620410	Pet services	
#                             620420,	    #	5	620420	Vet services	
#                                         #	4	RNTSPVEH	Rental of recreational vehicles	
#                             520901,	    #	4	520901	Docking and landing fees	
#                             520904,	    #	5	520904	Rental non-camper trailer	
#                             520907,	    #	5	520907	Boat and trailer rental out-of-town trips	
#                             600110,	    #	4	600110	Outboard motors	
#                             620902,	    #	4	620902	Rental of campers and other vehicles on trips (deleted 91)	
#                             620906,	    #	5	620906	Rental of boat	
#                             620907,	    #	5	620907	Rental of all campers and other recreational vehicles (620921 and 620922)
#                             620909,	    #	5	620909	Rental of campers on out-of-town trips	
#                             620919,	    #	5	620919	Rental of other vehicles on out-of-town trips	
#                             620921,	    #	5	620921	Rental of motorized camper	
#                             620922,	    #	5	620922	Rental of other RV"s	
#                 #	4	RECEQUIP	Sports; recreation and exercise equipment	
#                             620908,	    #	5	620908	Rental and repair of miscellaneous sports equipment	
#                 #	4	PHOTOEQ	Photographic equipment; supplies and services	
#                             620320,	    #	5	620320	Photographer fees	
#                             620330,	    #	5	620330	Photo processing	
#                             620905,	    #	5	620905	Repair and rental of photographic equipment	
#                             680310,	    #	4	680310	Live entertainment for catered affairs (new UCC Q20072)	
#                             680320]		#	4	680320	Rental of party supplies for catered affairs (new UCC Q20072)	
                          
                      
# uccmapping['ENTRTAIND']    = [    	        #		ENTRTAIND	Entertainment (durables)	
#                 #	3	TVAUDIO	Audio and visual equipment and services	
#                             310110,	    #	4	310110	Black and white TV (deleted 05)	
#                             310120,	    #	4	310120	Color RV console (deleted 05)	
#                             310130,	    #	4	310130	Color TV portable/ table model (deleted 05)	
#                             310140,	    #	4	310140	Televisions	
#                             310210,	    #	4	310210	VCR's and video disc players	
#                             310220,	    #	4	310220	Video cassettes; tapes; and discs	
#                             310230,	    #	4	310230	Video game hardware and software	
#                             310311,	    #	4	310311	Radios	
#                             310312,	    #	4	310312	PHONOGRAPHS This UCC was deleted as an interview UCC (now UCC 310333) 
#                             310313,	    #	4	310313	Tape recorders and players	
#                             310314,	    #	4	310314	Personal digital audio players	
#                             310320,	    #	4	310320	Sound components and component systems	
#                             310330,	    #	4	310330	Accessories and other sound equipment (deleted 96)	
#                             310333,	    #	4	310333	Accessories and other sound equipment	
#                             310334,	    #	4	310334	Satellite dishes	
#                             310340,	    #	4	310340	CDs; records; audio tapes	
#                             310342,	    #	4	310342	Records, CDs, audio tapes, needles (deleted 05)	
#                             610130,	    #	4	610130	Musical instruments and accessories	
#                 #	3	PETSPLAY	Pets; toys; hobbies; and playground equipment	
#                             610110,	    #	4	610110	Toys; games; arts and crafts; and tricycles	
#                             610120,	    #	4	610120	Playground equipment	
#                             610140,	    #	4	610140	Stamp and coin collecting [2004]	
#                             610320,	    #	5	610320	Pet purchase; supplies; medicine	
#                 #	4	UNMTRBOT	Unmotored recreational vehicles	
#                             600121,	    #	5	600121	Boat without motor and boat trailers	
#                 #	5	600127	Trade-In Allowance For Boat Without Motor or Non Camper-Type Trailer 
#                             600122,	    #	5	600122	Trailer and other attachable campers	
#                 #	5	600128	Trade-In Allowance For Trailer-Type or Other Attachable-Type Camper
#                 #	4	PWRSPVEH	Motorized recreational vehicles	
#                             600131,	    #	5	600131	Motorized camper-coach or other vehicles (net outlay) (600141 and 600142)
#                 #	5	600137	Trade-in allowance for motorized camper-coach (now 600143 and 600144)	
#                             600132,	    #	5	600132	Purchase of boat with motor	
#                 #	5	600138	Trade-In Allowance For Boat With Motor (deleted)	
#                             600141,	    #	5	600141	Purchase of motorized camper	
#                 #	5	600143	Trade-In allowance, motorized camper (deleted)	
#                             600142,	    #	5	600142	Purchase of other vehicle	
#                 #	5	600144	Trade-In allowance, other vehicle (deleted)	
#                 #	4	RECEQUIP	Sports; recreation and exercise equipment	
#                             600210,	    #	5	600210	Athletic gear; game tables; and exercise equipment	
#                             600310,	    #	5	600310	Bicycles	
#                             600410,	    #	5	600410	Camping equipment	
#                             600420,	    #	5	600420	Hunting and fishing equipment	
#                             600430,	    #	5	600430	Winter sports equipment	
#                             600900,	    #	5	600900	Water and miscellaneous sports equipment (split into 600901 and 600902)
#                             600901,	    #	5	600901	Water sports equipment	
#                             600902,	    #	5	600902	Other sports equipment	
#                 #	4	PHOTOEQ	Photographic equipment; supplies and services	
#                             610210,	    #	5	610210	Film	
#                             610230]		#	5	610230	Photographic equipment	

                #	2	PERSCARE	Personal care products and services								
uccmapping['PERSCARED'] = [    	        #		PERSCARED	Personal care products
                            640130,	    #	3	640130	Wigs and hairpieces	
                            640420]		#	3	640420	Electric personal care appliances	
                        
uccmapping['PERSCARES'] = [    	        #		PERSCARES	personal care services	
                            650110,	    #	3	650110	Personal care services for females is now included in newly added UCC 650310
                            650210,	    #	3	650210	Personal care services for males is now included in newly added UCC 650310
                            650310,	    #	3	650310	Personal care services	
                            650900]		#	3	650900	The repair and servicing of personal care appliances no longer collected 
                        
uccmapping['READING']   = [    	        #	2	READING	Reading	
                            590110,	    #	3	590110	Newspapers-single copy and subscription (split into 590111 and 590112)
                            590111,	    #	3	590111	Newspaper subscriptions (deleted 05)	
                            590112,	    #	3	590112	Newspapers, non-subscription (deleted 05)	
                            590210,	    #	3	590210	Magazines and periodicals-single copy and subscription (now 590211 590212)
                            590211,	    #	3	590211	Magazine subscriptions (deleted 05)	
                            590212,	    #	3	590212	Magazines, non-subscription (deleted 05)	
                            590220,	    #	3	590220	Books thru book clubs	
                            590230,	    #	3	590230	Books not thru book clubs	
                            590310,	    #	3	590310	Newspaper; magazine by subscription (2005)	
                            590410,	    #	3	590410	Newspaper; magazine non-subscription (2005)	
                            660310]		#	3	660310	Encyclopedia and other sets of reference books	


uccmapping['BOOKS']   = [    	        #		Durable Books
                            
                            590220,	    #	3	590220	Books thru book clubs	
                            590230,	    #	3	590230	Books not thru book clubs	
                            660310]		#	3	660310	Encyclopedia and other sets of reference books	

uccmapping['MAGAZINES']   = [    	        #		MAGAZINES
                            590110,	    #	3	590110	Newspapers-single copy and subscription (split into 590111 and 590112)
                            590111,	    #	3	590111	Newspaper subscriptions (deleted 05)	
                            590112,	    #	3	590112	Newspapers, non-subscription (deleted 05)	
                            590210,	    #	3	590210	Magazines and periodicals-single copy and subscription (now 590211 590212)
                            590211,	    #	3	590211	Magazine subscriptions (deleted 05)	
                            590212,	    #	3	590212	Magazines, non-subscription (deleted 05)		
                            590310,	    #	3	590310	Newspaper; magazine by subscription (2005)	
                            590410]    #	3	590410	Newspaper; magazine non-subscription (2005)	
                        
uccmapping['EDUCATION'] = [    	        #		EDUCATION	
                            660110,	    #	3	660110	School books; supplies; equipment for college	
                            660210,	    #	3	660210	School books; supplies; equipment for elementary; high school
                            660410,	    #	3	660410	School books; supplies; equipment for vocational and technical schools 
                            660900,	    #	3	660900	School books; supplies; equipment for day care; nursery;preschool; other 
                            660901,	    #	3	660901	School books; supplies; equipment for day care; nursery (new UCC Q20062)
                            660902,	    #	3	660902	School books; supplies; equipment for other schools (new UCC Q20062)
                            670110,	    #	3	670110	College tuition	
                            670210,	    #	3	670210	Elementary and high school tuition	
                            670410,	    #	3	670410	Vocational and technical school tuition (new UCC Q20072)	
                            670901,	    #	3	670901	Other schools tuition	
                            670902]		#	3	670902	Other school expenses including rentals	

# uccmapping['EDUCATION_BOOKS'] = [    	        #		EDUCATIONAL BOOKS	
#                             660110,	    #	3	660110	School books; supplies; equipment for college	
#                             660210,	    #	3	660210	School books; supplies; equipment for elementary; high school
#                             660410,	    #	3	660410	School books; supplies; equipment for vocational and technical schools 
#                             660900,	    #	3	660900	School books; supplies; equipment for day care; nursery;preschool; other 
#                             660901,	    #	3	660901	School books; supplies; equipment for day care; nursery (new UCC Q20062)
#                             660902,	    #	3	660902	School books; supplies; equipment for other schools (new UCC Q20062)
#                             ]	


# uccmapping['EDUCATIONS'] = [    	        #		EDUCATION	Services
#                             670110,	    #	3	670110	College tuition	
#                             670210,	    #	3	670210	Elementary and high school tuition	
#                             670410,	    #	3	670410	Vocational and technical school tuition (new UCC Q20072)	
#                             670901,	    #	3	670901	Other schools tuition	
#                             670902]		#	3	670902	Other school expenses including rentals	

# uccmapping['LUGGAGE'] = [    	        #		LUGGAGE and SIMLAR EQUIPMENT	
#                             430130    #	5	430130	Luggage	

#                             ]

uccmapping['TOBACCO']   = [    	        #	2	TOBACCO	Tobacco products and smoking supplies	
                            630110,	    #	3	630110	Cigarettes	
                            630210]		#	3	630210	Other tobacco products

                                      
                                        #	2	MISC	Miscellaneous	


uccmapping['OCCUPEXP']  = [    	        #		OCCUPEXP	Occupational expenses	 (part of MISC1 in BLS accounting)
                            900001,	    #	3	900001	Occupational expenses(ITAB)	
                            900002]		#	3	900002	Occupational expenses	

uccmapping['FINCHRG']   = [    	        #	    FINCHRG 	Finance (MISC2 in BLS accounting)	
                              5420,     #   3   005420  Finance, late, interest charges for credit cards
                              5520,     #   3   005420  Finance, late, interest charges for student loans
                              5620,     #   3   005620  Finance, late, interest charges for other loans 
                            710110]	    #	3	710110	Finance charges excluding mortgage and vehicle	                           

uccmapping['FEENCHRG']  = [    	        #		FEENCHRG	Finance and other charges	 (part of MISC1 in BLS accounting)
                            620112,	    #	3	620112	Credit card memberships	
                            620115,	    #	3	620115	Shopping club membership fees	
                            620926,     #   3   620926  Lotteries and pari-mutuel losses
                            680110,	    #	3	680110	Legal fees	
                            680210,	    #	3	680210	Safe deposit box rental	
                            680220,	    #	3	680220	Checking accounts; other bank service charges	
                            680901,	    #	3	680901	Cemetery lots; vaults; maintenance fees	
                            680140,	    #	3	680140	Funeral expenses	
                            680902,	    #	3	680902	Accounting fees	
                            680904,	    #	3	680904	Dating services	
                            790840,		#	3	790840	Other charges in sale of other properties (deleted 91)
                            790600,		#	3	790600	Expenses for other properties (vacation homes)	
                            880210]	    #	3	880210	Interest paid; home equity line of credit (other property)
                        
uccmapping['LIFEINSUR'] = [    	        #		LIFEINSUR	Life and other personal insurance	
                            2120,	    #	4	2120	Other non-health insurance	
                            700110]		#	4	700110	Life; endowment; annuity; other personal insurance	
                        
uccmapping['PENSIONS']  = [    	        #		PENSIONS	Pensions and Social Security	
                            800910,	    #	4	800910	Deductions for government retirement	
                            800920,	    #	4	800920	Deductions for railroad retirement	
                            800931,	    #	4	800931	Deductions for private pensions	
                            800932,	    #	4	800932	Non-payroll deposit to retirement plans	
                            800940]		#	4	800940	Deductions for Social Security	
                        
                        
        #		Non-Consumption Expenditures (i.e. not part of total expenditures)	
                        
uccmapping['CASHCONT']  = [    	        #		CASHCONT	Cash contributions	
                #	2	CASHCONT	(replaces contribution variables on FMLY files starting in 2001Q2)
                            800111,	    #	3	800111	Alimony expenditures (Sec 19)	
                            800121,	    #	3	800121	Child support expenditures (Sec. 19)	
                            800800,	    #	3	800800	Cash gifts and contributions to non-CU members (deleted 92)	
                            800803,	    #	3	800803	Cash gifts to non-CU members contributions to organizations (deleted 02)
                            800804,	    #	3	800804	Support for college students (Sec. 19)	
                            800811,	    #	3	800811	Gift to non-CU members of stocks; bonds; and mutual funds	
                            800821,	    #	3	800821	Cash contributions to charities and other organizations	
                            800831,	    #	3	800831	Cash contributions to church; religious organizations	
                            800841,	    #	3	800841	Cash contribution to educational institutions	
                            800851,	    #	3	800851	Cash contribution to political organizations	
                            800861]		#	3	800861	Other cash gifts

uccmapping['usedautos'] = [  
                            860100,     #   860100 Used Car Sales 
                            460110,    #	460110	Used Cars	
                            450199]     #	450199 = 450116*(-1)	Trade-In Allowance For New Cars  (netted out by 450116)

uccmapping['usedlighttrucks'] = [  
                            860200,      #  860200  Used Truck sales 	       	
                            460901,      #	460901	Used Trucks	
                            450299]     #	450299 = 450216*(-1)	Trade-In Allowance For New Trucks  (netted out by 450216)
	

# total expenditure is the union of all these codes but excludes cash contributions
# totexpkeys = []
# for key in uccmapping:
#     totexpkeys = totexpkeys + uccmapping[key]

# uccmapping['TOTEXP3'] = totexpkeys

# save output to yaml file
with open('../output/ucc_category_map.yml', 'w') as yamlfile:
    yaml.dump(uccmapping, yamlfile, default_flow_style=False) #, f, default_flow_style=False