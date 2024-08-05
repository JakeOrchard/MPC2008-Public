import yaml

# dictionary with mappings
uccmapping = {}

#JO follows BLS 2017 PCE-CEX concordance
#all variables names are lower-cased to distinguish from UCC-FMLI hierarchy 

#DURABLE


uccmapping['newautos']  = [    	       	
                            450110,	    #	5	450110	New cars	
                            450116]     #	5	450116	Trade-In Allowance For New Cars  (netted out by 450199)
 
uccmapping['newlighttrucks']  = [
                            450210,    #	5	450210	New trucks	
                            450216]     #	5	450216	Trade-In Allowance For New Trucks or Vans (netted out by )
uccmapping['newcardown'] = [
                            870102,      # Cash downpayment for new cars, trucks, or vans, purchase financed
                            870101]       # New cars, trucks, or vans (net outlay), purchase not financed

uccmapping['usedautos'] = [  
                            860100,     #   860100 Used Car Sales 
                            460110,    #	460110	Used Cars	
                            450199]     #	450199 = 450116*(-1)	Trade-In Allowance For New Cars  (netted out by 450116)

uccmapping['usedlighttrucks'] = [  
                            860200,      #  860200  Used Truck sales 	       	
                            460901,      #	460901	Used Trucks	
                            450299]     #	450299 = 450116*(-1)	Trade-In Allowance For New Trucks  (netted out by 450216)

uccmapping['usedcardown'] = [
                            870202,      # Cash downpayment for used cars, trucks, or vans, purchase financed
                            870201]       # Used cars, trucks, or vans (net outlay), purchase not financed

uccmapping['tires'] =   [
                            480110]      # 480110 Tires

uccmapping['acces_parts'] = [

                            480212,      #	D	Vehicle products
                            480213]	     #  I	Parts/equip/accessories


uccmapping['furniture'] = [

                        290110,#	I	Mattress and springs
                        290120,#	I	Other bedroom furniture	
                        290210,#	I	Sofas	
                        290310,#	I	Living room chairs	
                        290320,#	I	Living room tables	
                        290410,#	I	Kitchen, dining room furniture	
                        290420,#	I	Infants' furniture	
                        290440,#	I	Wall units, cabinets and other occasional furniture	
                        320901,#	I	Office furn for home use	
                        290430,#	D	Outdoor furniture	
                        340904,#	I	Rental of furniture	
                        680320,#	I	Rental of party supplies for catered affairs	
                        322130] #	D	Infants' equipment	

uccmapping['clock'] = [
                        320210,#	D	Clocks and other household decorative items (thru 2006 Q4)	DCLFRC
                        320220,#	D	Lamps and lighting fixtures (thru 2012 Q4) 	DCLFRC
                        320221,#	D	Lamps, lighting fixtures, and ceiling fans (introduced 2013 Q1)	DCLFRC
                        320231,#	D	Other household decorative items (thru 2006 Q4)	DCLFRC
                        320233,#	D	Clocks and other household decorative items (introduced 2007 Q1)	DCLFRC
                        322904,#	D	Closet and storage items	DCLFRC
                        321150] #	D	Outdoor equipment	DCLFRC



uccmapping['carpets'] = [
                        320626,#	I	Flooring installation, replacement, and repair (owned vacation) (introduced 2013 Q2)
                        320625,#	I	Flooring installation, replacement, and repair (owned home) (introduced 2013 Q2)
                        230133,#	I	Wall-to-wall carpet (replacement) (owned home)
                        320111,#	I	Floor coverings, nonpermanent
                        220616] #	I	Wall-to-wall carpeting

uccmapping['window_cover'] = [

                        320120,#	I	Window coverings
                        280210]#	I	Curtains and drapes


uccmapping['major_app'] = [
                        230118,#	I	Dishwasher (built-in), garbage disposals, range hoods (owned home)
                        220612,#	I	Dishwasher, disposal, range hood, capital improvement, owned
                        300112,#	I	Refrigerators, freezers (owned home)
                        300217,#	I	Washing Machine/Clothes Dryers (owned home) (introduced 2013 Q2)
                        300212,#	I	Washing machines (owned home) (thru 2013 Q1)
                        300222,#	I	Clothes dryers (owned home) (thru 2013 Q1)
                        300312,#	I	Cooking stoves, ovens (owned home)
                        300322,#	I	Microwave ovens (owned home)
                        300332,#	I	Portable dishwasher (owned home)
                        300412,#	I	Window air conditioners (owned home)
                        320511,#	I	Electric floor cleaning equipment
                        320512,#	I	Sewing machines
                        300900,#	D	Miscellaneous household appliances
                        320522,#	I	Portable heating/cooling equip
                        690245,#	I	Other household appliances (owned home)
                        322150,#	D	Outdoor equipment
                        321905]#	D	Miscellaneous household equipment and parts

uccmapping['small_app'] = [
                        320521] #	I	Small electric kitchen appliances	DSEARC	     

uccmapping['dishes'] = [
                        320310,#	I	Plastic dinnerware (thru 2012 Q4)
                        320320,#	D	China and other dinnerware (thru 2012 Q4)
                        320330,#	D	Flatware (thru 2012 Q4)
                        320340,#	D	Glassware (thru 2012 Q4)
                        320350,#	D	Silver serving pieces (thru 2012 Q4)
                        320345,#	D	Dinnerware, glassware, serving pieces (introduced 2013 Q1)
                        333510,#	D	Miscellaneous household products
                        320360]#	I	Other serving pieces

    
uccmapping['nelec_cook'] = [ 
                        320370,#	D	Nonelectric cookware
                        320380,#	D	Tableware, nonelectric kitchenware
                        323904,#	D	Closet and storage items
                        322905]#	D	Miscellaneous household equipment and parts


uccmapping['tools'] = [ 
                        320420,#	D	Power tools
                        320902,#	D	Hand tools
                        320610,#	D	Misc supplies & equipment
                        320430]#	D	Other hardware

uccmapping['tv'] = [ 

                        310140,#	I	Televisions (introduced 2005 Q2)
                        310110,#	I	Black and white tv (thru 2005 Q1)
                        310120,#	I	Color tv - console (thru 2005 Q1)
                        310130]#	I	Color tv - portable/table model (thru 2005 Q1)

uccmapping['other_video'] = [
                        310334,#	I	Satellite dishes 
                        310210,#	I	VCR's/video disc players
                        480215,#	I	Vehicle video equipment (thru 2013 Q1)
                        310335,#	D	Misc video equipment
                        311900,#	D	Accessories for electronic equipment
                        310315] #	D	Digital media players and recorders


uccmapping['audio_equip'] = [ 
                        480214,#	I	Vehicle audio equipment, excluding labor (thru 2013 Q1)
                        310316,#	I	Radios, speakers, and sound computer systems (introduced 2013 Q2)
                        310311,#	I	Radios (thru 2013 Q1)
                        310312,#	D	Phonographs
                        310313,#	D	Tape recorders and players
                        310314,#	I	Personal digital audio players
                        310320,#	I	Sound components/component systems
                        310331,#	D	Misc sound equipment
                        310332,#	D	Sound equip accessories
                        490502]#	I	Vehicle audio equipment including labor

uccmapping['audio_record'] = [
                        310341,#	I	Record/tape/cd/video mail ord club (thru 2005 Q1)
                        310342,#	I	Records, cd's, audio tapes, needles (thru 2005 Q1)
                        310340,#	I	CDs, records, audio tapes (introduced 2005 Q2)
                        310350]#	I	Streaming/downloading audio

uccmapping['video_record'] = [
                        310220,#	I	Video cassettes/tapes/discs
                        310240,#	I	Streaming/downloading video (thru 2017 Q1)
                        310243]#	I	Rental, streaming, downloading video (introduced 2017 Q2)

uccmapping['photo_equip'] = [
                        610230,#	I	Photographic equipment
                        610903]#	D	Visual goods
                                

uccmapping['computer_equip'] = [
                        312900,#	D	Accessories for electronic equipment
                        690111,#	I	Computer/computer hardware nonbus use
                        690115,#	I	Personal digital assistants
                        690117,#	I	Portable memory
                        690118,#	I	Digital book readers (introduced 2011 Q2)
                        692230]#	I	Business equipment for home use

uccmapping['software'] = [
                        690112,#	I	Computer software/accessories nonbus use (thru 2011 Q1)
                        690119,#	I	Computer software (introduced 2011 Q2)
                        690120,#	D	Computer accessories (introduced 2011 Q2)
                        311230,#	D	Video game hardware/software (thru 2011 Q1)
                        310231]#	I	Video game software (introduced 2011 Q2)

uccmapping['calculator'] = [
                        690220,#	I	Calculators (thru 2005 Q1)
                        691230]#	I	Business equipment for home use

uccmapping['sport_equip'] = [
                        600142,#	I	Purchase of other vehicle
                        860302,#	I	Amount other vehicle sold or reimbursed
                        600210,#	D	Athletic gear/game tables/ex. equip
                        600410,#	D	Camping equipment
                        600420,#	D	Hunting, fishing equipment
                        600430,#	I	Winter sport equipment
                        600901,#	I	Water sports equipment
                        600902,#	I	Other sports equipment
                        610120,#	I	Playground equipment
                        610901]#	D	Fireworks

uccmapping['motorcycles'] = [
                        450220,#	I	New motorcycles
                        460902,#	I	Used motorcycles
                        860500]#	I	Amount motorcycle, motor scooter, or moped sold or reimbursed

uccmapping['bicycles'] = [
                        600310]#	I	Bicycles

uccmapping['boats'] = [
                        600121,#	I	Boat w/o motor/boat trailers
                        600132,#	I	Purchase of boat with motor
                        600110,#	I	Outboard motors
                        860600,#	I	Amount boat with motor sold or reimbursed
                        860700]#	I	Amount boat without motor or non camper-type trailer such as for motorcycle sold or reimbursed

uccmapping['planes'] = [
                        460903,#	I	Used aircraft
                        450900]#	I	New aircraft

uccmapping['other_recv'] = [
                        600122,#	I	Trailer/other attachable campers
                        600141,#	I	Purchase of motorized camper
                        860301,#	I	Amount motor home sold or reimbursed
                        860400]#	I	Amount trailer-type or other attchable-type camper sold or reimbursed

uccmapping['rec_book'] = [
                        590220,#	I	Books thru book clubs
                        590230,#	I	Books not thru book clubs
                        660310,#	I	Encyclopedia and other reference book sets
                        662000]#	D	School supplies, etc. - unspecified


uccmapping['music_instrument'] = [
                        610130]#	I	Music instruments/accessories

uccmapping['jewelry'] = [
                        430120]#	I	Jewelry

uccmapping['watches'] = [
                        430110]#	D	Watches

uccmapping['med_equip_t'] = [
                        550330,#	I	Supportive/conval med. equip.
                        550340,#	I	Hearing aids
                        551320]#	I	Medical equip. for general use (0.5 of 550320)

uccmapping['eyeglass'] = [
                        550110,#	I	Eyeglasses
                        560310]#	I	Eyecare services

uccmapping['educ_book'] = [
                        660110,#	I	School books/supplies& equipment for college
                        661900,#	I	"School books, supplies, & equip for day care, nursery, other(thru 2006 Q1)"
                        661901,#	I	"School books, supplies, & equip for day care, nursery(introduced 2006 Q2)"
                        661210,#	I	School books/supplies & equip for elementary/high school
                        660902,#	I	School books, supplies, & equip for other (introduced 2006 Q2) 
                        660410]#	I	"School books/supplies & equip for vocational and technical schools(introduced 2007 Q2)"

uccmapping['luggage'] = [
                        430130,#	D	Luggage
                        323130,#	D	Infants' equipment
                        313900]#	D	Accessories for electronic equipment

uccmapping['telephone'] = [
                        320232,#	D	Telephones and accessories
                        690210,#	I	Telephone answering devices
                        600903]#	D	Global positioning equipment devices

############### NONDURABLE
uccmapping['food_out'] = [
                        190904]#	I	Food prepared by CU on out-of-town trips

uccmapping['food_stamps'] = [
                        900150]#	I	Food Stamps]

uccmapping['foodhome']  = [             #BLS concordance uses diary here	
                                        #		790210	Total purchases at grocery stores (Note: This is not a primary expenditure)
                            790220	,	#	4	790220	Grocery stores (thru Q20071)	
                            790230	,	#	4	790230	Convenience stores (thru Q20071)	
                            790240]		#	4	790240	Food and non alcoholic beverages (new UCC Q20072)	

uccmapping['alchome']   = [    	       #BLS concordance uses diary here	
                            790310, 	#	4	790310	Beer and wine (thru Q20071)	
                            790320,	    #	4	790320	Other alcoholic beverages (thru Q20071)	
                            790330]		#	4	790330	Beer wine other alcohol (new UCC Q20072) (Note:790330=790310+790320)
                          
uccmapping['women_cloth'] = [
                            380110,#	D	Women's coats and jackets
                            380210,#	D	Women's dresses
                            380311,#	I	Women's sportcoats, tailored jackets
                            380312,#	D	Women's vests and sweaters (thru 2012 Q4)
                            380313,#	D	Women's shirts, tops,blouses (thru 2012 Q4)
                            380320,#	D	Women's skirts
                            380331,#	D	Women's pants (thru 2006)
                            380332,#	D	Women's shorts,shorts sets (thru 2006)
                            380333,#	D	Women's pants and shorts (introduced 2007)
                            380340,#	D	Women's active sportswear
                            380410,#	D	Women's sleepwear
                            380420,#	D	Women's undergarments
                            380430,#	D	Women's hosiery
                            380510,#	I	Women's suits
                            380901,#	D	Women's accessories
                            380902,#	I	Women's uniforms
                            380903,#	I	Women's costumes
                            380315,#	D	Women's shirts, tops, and vests (introduced 2013 Q1)
                            390110,#	D	Girls' coats and jackets
                            390120,#	D	Girls' dresses, suits
                            390210,#	D	Girls' shirts/blouses/sweaters
                            390221,#	I	Girls' skirts and pants (thru 2007 Q1)
                            390222,#	I	Girls' shorts, shorts sets (thru 2007 Q1)
                            390223,#	I	Girls' skirts, pants and shorts (introduced 2007 Q2)
                            390230,#	D	Girls' active sportswear
                            390310,#	D	Girls' underwear and sleepwear
                            390321,#	D	Girls' hosiery
                            390322,#	D	Girls' accessories
                            390901,#	I	Girls' uniforms
                            390902]#	I	Girls' costumes

uccmapping['men_cloth'] = [
                            360110,#	I	Men's suits
                            360120,#	I	Men's sportcoats/ tailored jackets
                            360210,#	D	Men's coats and jackets
                            360311,#	D	Men's underwear
                            360312,#	D	Men's hosiery
                            360320,#	I	Men's nightwear/loungewear
                            360330,#	D	Men's accessories
                            360340,#	D	Men's sweaters and vests (thru 2012 Q4)
                            360350,#	D	Men's active sportswear
                            360410,#	D	Men's shirts
                            360511,#	D	Men's pants (thru 2006)
                            360512,#	D	Men's shorts,shorts sets (thru 2006)
                            360513,#	D	Men's pants and shorts (introduced 2007)
                            360901,#	I	Men's uniforms
                            360902,#	I	Men's costumes
                            360420,#	D	Men's shirts, sweaters, and vests (introduced 2013 Q1)
                            370110,#	I	Boys' coats and jackets
                            370120,#	I	Boys' sweaters (thru 2012 Q4)
                            370130,#	D	Boys' shirts (thru 2012 Q4)
                            370211,#	D	Boys' underwear
                            370212,#	I	Boys' nightwear
                            370213,#	D	Boys' hosiery
                            370220,#	D	Boys' accessories
                            370311,#	I	Boys' suits, sportcoats,vests
                            370312,#	I	Boys' pants (thru 2007 Q1)
                            370313,#	I	Boys' shorts, shorts sets (thru 2007 Q1)
                            370314,#	I	Boys' pants and shorts (introduced 2007 Q2)
                            370902,#	I	Boys' costumes
                            370903,#	I	Boys' uniforms
                            370904,#	I	Boys' active sportswear
                            370125]#	D	Boys' shirts and sweaters (introduced 2013 Q1)

uccmapping['child_cloth'] = [
                            410110,#	I	Infant coat/jacket/snowsuit
                            410120,#	D	Infant dresses/outerwear
                            411130,#	D	Infants' underwear
                            410140,#	I	Infant nightwear/loungewear
                            410901]#	D	Infant accessories

uccmapping['material_cloth'] = [
                            420110,#	D	Materials for making clothes (thru 2012 Q4)
                            420120,#	D	Sewing patterns and notions (thru 2012 Q4)
                            420115]#	D	Material and supplies for sewing, needlework, and quilting (introduced 2013 Q1)

uccmapping['shoes'] = [
                            400110,#	D	Men's  footwear
                            400210,#	D	Boys' footwear
                            400310,#	D	Women's footwear
                            400220]#	D	Girls' footwear

uccmapping['gasoline'] =[
                            470111,#	I	Gasoline
                            470112,#	I	Diesel fuel
                            470113,#	I	Gasoline on out of town trips
                            470114,#	D	Gasohol
                            470311]#	I	Electric vehicle charging (introduced 2017 Q2)

uccmapping['motor_oil'] = [
                            470211,#	I	Motor oil
                            470212,#	I	Motor oil on out-of-town trip
                            470220]#	I	Coolant/additives/brake/transmission fluids

uccmapping['fuel_oil'] = [
                            250111,#	I	Fuel oil (renter)
                            250112,#	I	Fuel oil (owned home)
                            250113,#	I	Fuel oil (owned vacation)
                            250114,#	I	Fuel oil (rented vacation)
                            250211,#	I	Gas, bottled/tank (renter)
                            250212,#	I	Gas, bottled/tank (owned home)
                            250213,#	I	Gas, bottled/tank (owned vacation & rv''s)
                            250214,#	I	Gas, bottled/tank (rented vacation)
                            251911,#	I	Coal, wood, other fuels (renter) (introduced 2005 Q2)
                            251912,#	I	Coal, wood, other fuels (introduced 2005 Q2)
                            251913,#	I	Coal, wood, other fuels (owned vacation) (introduced 2005 Q2)
                            251914]#	I	Coal, wood, other fuels (rented vacation) (introduced 2005 Q2)

uccmapping['other_fuel'] = [
                            250901,#	I	Wood/other fuels (renter) (thru 2005 Q1)
                            250902,#	I	Wood/other fuels (owned home) (thru 2005 Q1)
                            250903,#	I	Wood/other fuels (owned vacation) (thru 2005 Q1)
                            250904,#	I	Wood/other fuels (rented vacation) (thru 2005 Q1)
                            252911,#	I	Coal, wood, other fuels (renter) (introduced 2005 Q2)
                            252912,#	I	Coal, wood, other fuels (introduced 2005 Q2)
                            252913,#	I	Coal, wood, other fuels (owned vacation) (introduced 2005 Q2)
                            252914,#	I	Coal, wood, other fuels (rented vacation) (introduced 2005 Q2)
                            250221,#	I	Coal (renter) (thru 2005 Q1)
                            250222,#	I	Coal (owned home) (thru 2005 Q1)
                            250223,#	I	Coal (owned vacation & rv''s) (thru 2005 Q1)
                            250224]#	I	Coal (rented vacation) (thru 2005 Q1)

uccmapping['prescription'] = [
                            540000]#	I	Prescription drugs

uccmapping['non_prescription'] = [
                            550210,#	D	Non-prescription drugs
                            550410,#	D	Non-prescription vitamins
                            180720]#	D	Vitamin supplements

uccmapping['other_medical'] = [
                            550310,#	D	Topicals and dressings
                            552320]#	I	Medical equip. for general use


uccmapping['games'] = [
                            312230,#	D	Video game hardware/software (thru 2011 Q1)
                            310232,#	D	Video game hardware (introduced 2011 Q2)
                            610110,#	D	Toys, games, hobbies, and tricycles
                            610140,#	I	Stamp and coin collecting
                            620930]#	D	Online gaming services

uccmapping['pet_food'] = [
                            610310,#	D	Pet food
                            610320]#	I	Pet-purchase/supplies/medicine

uccmapping['flower'] = [
                            320903]#	D	Indoor plants, fresh flowers

uccmapping['photo_supply'] = [
                            610210,#	I	Film
                            610220]#	D	Other photographic supplies

uccmapping['house_clean'] = [
                            320140,#	D	Laundry and cleaning equip.
                            330110,#	D	Soaps and detergents
                            330210,#	D	Other laundry and cleaning products
                            330511]#	I	Termite/pest control products

uccmapping['house_paper'] = [
                            330310,#	D	Cleansing and toilet tissue, paper towels and napkins
                            331510]#	D	Miscellaneous household products

uccmapping['house_linen'] = [
                            280110,#	D	Bathroom linens
                            280120,#	D	Bedroom linens
                            280130,#	D	Kitchen and dining room linens (thru 2012 Q4)
                            280140,#	D	Kitchen and dining room other linens (introduced 2013 Q1)
                            280220,#	I	Slipcovers, decorative pillows
                            280900,#	I	Other linens (thru 2012 Q4)
                            321904]#	D	Closet and storage items

uccmapping['sewing'] = [
                            420115,#	I	Sewing, Knitting, and Quilting materials (introduced 2013 Q2)
                            280230]#	I	Sewing materials for slipcovers, curtains, other sewing materials for the home (thru 2013 Q1)


uccmapping['house_misc'] =[
                            332510]#	D	Miscellaneous household products

uccmapping['pers_care'] = [
                            321130,#	D	Infants' equipment
                            640110,#	D	Hair care products
                            640120,#	D	Non-electric articles for the hair
                            640130,#	I	Wigs and hairpieces
                            640210,#	D	Oral hygiene products,articles
                            640220,#	D	Shaving needs
                            640410,#	D	Deodorants, feminine hygiene, misc. personal care
                            412130,#	D	Infants' underwear
                            640430]#	I	Adult diapers (introduced 2013 Q2)

uccmapping['cosmetic'] = [
                            640310]#	D	Cosmetics, perfume, bath prep

uccmapping['pers_elec'] = [
                            640420]#	I	Electric personal care appliances

uccmapping['tobacco'] = [
                            630110,#	I	Cigarettes
                            630210,#	I	Other tobacco products
                            630220]#	D	Smoking accessories


uccmapping['newspaper'] = [
                            590310,#	I	Newspaper, magazine by subscription (introduced 2005 Q2)
                            590410,#	I	Newspaper, magazine non-subscription (introduced 2005 Q2)
                            590211,#	I	Magazine subscriptions (thru 2005 Q1)
                            590212,#	I	Magazine, non-subscriptions (thru 2005 Q1)
                            590900,#	D	Newsletters
                            590111,#	I	Newspaper subscriptions (thru 2005 Q1)
                            590112]#	I	Newspaper, non-subscriptions (thru 2005 Q1)

uccmapping['stationary'] =[
                            330410,#	D	Stationery, stationery supplies, giftwraps
                            662900,#	I	"School books, supplies, & equip for day care, nursery, other(thru 2006 Q1)"
                            662901,#	I	"School books,supplies, & equip for day care, nursery(introduced 2006 Q2)"
                            662210,#	I	School books/supplies & equip for elementary/high school
                            661000,#	D	School supplies, etc. - unspecified
                            610902,#	D	Souvenirs
                            610900]#	I	Recreation expenses, out-of-town trips


##Services

uccmapping['rent'] =[
                            210110,#	I	Rent
                            230117,#	I	Dishwasher (built-in), garbage disposals, range hoods (renter)
                            300111,#	I	Refrigerators, freezers (renter)
                            300216,#	I	Washing Machine/Clothes Dryer (renter) (introduced 2013 Q2)
                            300211,#	I	Washing machines (renter) (thru 2013 Q1)
                            300221,#	I	Clothes dryers (renter) (thru 2013 Q1)
                            300311,#	I	Cooking stoves,ovens (renter)
                            300321,#	I	Microwave ovens (renter)
                            300331,#	I	Portable dishwasher (renter)
                            300411,#	I	Window air conditioners (renter)
                            690244,#	I	Other household appliances (renter)
                            320624,#	I	Flooring Installation, Replacement, and Repair (renter) (introduced 2013 Q2)
                            230134,#	I	Wall-to-wall carpet (renter)
                            320163,#	I	Wall-to-wall carpet (replacement) (renter)
                            230121,#	I	Repair/repl. hard surface flooring (renter) (thru 2013 Q1)
                            230141,#	I	Repair of disposal, b.i.d., range hood, renter
                            230150,#	I	Other repair and maintenance services, renter
                            240111,#	I	Paint, wallpaper, and supplies, renter
                            240121,#	I	Tools/equip. for painting and wallpapering, renter
                            240211,#	I	Materials for plaster., panel., roofing, gutters, renter
                            240221,#	I	Materials for patio, walk, fence, drive., masonry, brick, stucco, renter
                            240311,#	I	Plumbing supplies and equipment, renter
                            240321,#	I	Electrical supplies, heating/cooling equip., renter
                            320611,#	I	Materials for insulation, other maint./repair, renter
                            320621,#	I	Materials for hard surface flooring, renter
                            270901]#	I	Septic tank cleaning (renter)


uccmapping['own_home'] = [
                            910050,#	I	Rental equivalence of owned home
                            910060,#	I	Estimated monthly rental value of time share (thru 1999 Q1)
                            910070,#	I	Estimated monthly rental value of own. vac. Home (thru 1999 Q1)
                            910100,#	I	Rental equivalence of own. vac. home (thru 2007 Q1)
                            910101,#	I	"Rental equivalence of vacation home not available for rent(introduced 2007 Q2)"
                            910102,#	I	"Rental equivalence of vacation home available for rent(introduced 2007 Q2)"
                            910103]#	I	Rental equivalence of timeshares

uccmapping['group_home'] = [
                            800710,#	I	Rent as pay
                            212310]#	I	Housing while attending school

uccmapping['water'] = [
                            270211,#	I	Water/sewer maintenance (renter)
                            270212,#	I	Water/sewer maintenance (owned)
                            270213,#	I	Water/sewer maintenance (owned vac.)
                            270214]#	I	Water/sewer maintenance (rented vac.)

uccmapping['trash'] = [
                            270411,#	I	Trash/garbage coll. (renter)
                            270412,#	I	Trash/garbage coll. (owned)
                            270413,#	I	Trash/garbage coll. (owned vac.)
                            270414]#	I	Trash/garbage coll. (rented vac.)

uccmapping['electric'] = [
                            260111,#	I	Electricity (renter)
                            260112,#	I	Electricity (owned home)
                            260113,#	I	Electricity (owned vacation)
                            260114]#	I	Electricity (rented vacation)


uccmapping['nat_gas'] = [
                            260211,#	I	Util.--Natural gas (renter)
                            260212,#	I	Util.--Natural gas (own. samp.)
                            260213,#	I	Util.--Natural gas (own. vac.)
                            260214,#	I	Util.--Natural gas (rent. vac.)
                            270905]#	D	Steam heat

uccmapping['dr_service'] =[
                            560110]#	I	Physicians' services

uccmapping['dental_service'] =[
                            560210]#	I	Dental services

uccmapping['home_health'] =[
                           570240,#	I	Medical care incl. in homeowners expenses
                           340906]#	I	Care for elderly, invalids, handicapped, etc. (in the home)

uccmapping['med_lab'] =[
                           560330]#	I	Lab tests, x-rays

uccmapping['med_outpatient'] =[
                            560400,#	I	Service by professionals other than physicians (thru 2017 Q1)
                            560410,#	I	Non physician services inside home (introduced 2017 Q2)
                            560420,#	I	Non physician services outside home (introduced 2017 Q2)
                            571230]#	I	Other medical care services

uccmapping['other_med'] = [
                            572230,#	I	Other medical care services (0.5 of 570230)
                            570903,#	I	Rental of supportive/convalescent equipment
                            570901]#	I	Rental of medical equipment
                     
uccmapping['hospital'] = [
                            570110,#	I	Hospital room (thru 2005 Q1)
                            570210,#	I	Hospital service other than room (thru 2005 Q1)
                            570111]#	I	Hospital room and service (introduced 2005 Q2)

uccmapping['nurse_home'] = [
                            570220]#	I	Care in convalescent or nursing home

uccmapping['auto_repair'] = [
                            490000,#	D	Misc. auto repair/servicing
                            490110,#	I	Body work and painting
                            490211,#	I	Clutch, transmission repair (thru 2013 Q1)
                            490212,#	I	Drive shaft and rear end repair (thru 2013 Q1)
                            490221,#	I	Brake work, including adjustments (thru 2013 Q1)
                            490231,#	I	Steering or front end repair (thru 2013 Q1)
                            490232,#	I	Engine cooling system repair (thru 2013 Q1)
                            490311,#	I	Motor tune-up
                            490312,#	I	Lubrication, oil change and oil filters
                            490313,#	I	Front-end alignment,wheel balance/rotation
                            490314,#	I	Shock absorber replacement
                            490316,#	D	Gas tank repair, replacement
                            490318,#	I	Repair tire and other repair work
                            490319,#	I	Vehicle air conditioning repair (thru 2013 Q1)
                            490411,#	I	Exhaust system repair (thru 2013 Q1)
                            490412,#	I	Electrical system repair (thru 2013 Q1)
                            490413,#	I	Motor repair/replacement (thru 2013 Q1)
                            520550,#	I	Towing charges
                            490501,#	I	Vehicle accessories including labor
                            480216,#	I	Vehicle cleaning services
                            490300]#	I	Vehicle or engine repairs (introduced 2013 Q2)

uccmapping['car_lease'] =[
                            450310,#	I	Car lease payments (thru 2013 Q1)
                            450351,#	I	Extra fees for car/truck lease (introduced 2013 Q2)
                            452351,#	I	Extra fees for car lease (introduced 2014 Q1)
                            451451,#	I	Extra fees for truck lease (introduced 2014 Q1)
                            450313,#	I	Cash downpayment (car lease) (thru 2013 Q1)
                            450314,#	I	Termination fee (car lease) (thru 2013 Q1)
                            450312,#	I	Car lease, trade in (thru 2013 Q1)
                            450352,#	I	Car/Truck lease, trade in (introduced 2013 Q2)
                            452352,#	I	Car lease, trade in (introduced 2014 Q1)
                            451452,#	I	Truck lease, trade in (introduced 2014 Q1)
                            450350,#	I	Car/Truck lease payment (introduced 2013 Q2)
                            452350,#	I	Car lease payment (introduced 2014 Q1)
                            451450,#	I	Truck lease payment (introduced 2014 Q1)
                            450353,#	I	Car/Truck lease, down payment (introduced 2013 Q2)
                            452353,#	I	Car lease, down payment (introduced 2014 Q1)
                            451453,#	I	Truck lease, down payment (introduced 2014 Q1)
                            450354,#	I	Car/Truck lease, termination fee (introduced 2013 Q2)
                            452354,#	I	Car lease, termination fee (introduced 2014 Q1)
                            451454,#	I	Truck lease, termination fee (introduced 2014 Q1)
                            450410,#	I	Truck lease payments (thru 2013 Q1)
                            450413,#	I	Cash downpayment (truck lease) (thru 2013 Q1)
                            450414,#	I	Termination fee (truck lease) (thru 2013 Q1)
                            450412]#	I	Truck lease, trade in (thru 2013 Q1)


uccmapping['car_rental'] = [
                            520511,#	I	Auto rental (thru 2013 Q1)
                            520512,#	I	Auto rental, out-of-town trips (thru 2013 Q1)
                            520521,#	I	Truck rental (thru 2013 Q1)
                            520522,#	I	Truck rental, out-of-town trip (thru 2013 Q1)
                            520516,#	I	Auto/Truck rental (introduced 2013 Q2)
                            520517,#	I	Auto/Truck rental, out-of-town trip (introduced 2013 Q2)
                            520902,#	I	Motorcycle rental
                            520905,#	I	Motorcycle rental out-of-town trips
                            620909,#	I	Rental of camper on trips
                            620919,#	I	Rental of other vehicles on trips
                            520904,#	I	Rental non-camper trailer
                            620907,#	I	Rental of campers or other r.v.'s
                            620921,#	I	Rental of motorized camper
                            620922,#	I	Rental of other r.v.'s
                            620113,#	I	Auto service club membership (thru 2013 Q1)
                            620114]#	I	Auto service clubs and GPS services (introduced 2013 Q2)

uccmapping['parking_fee'] = [
                            520531,#	I	Parking fees excluding residence
                            520532,#	I	Parking fees while on trips
                            520541,#	I	Tolls
                            520542]# I	Tolls on out-of-town trips

uccmapping['train'] = [
                            530510]#	I	Intercity train fares

uccmapping['bus'] = [
                            530210]#	I	Intercity bus fares

uccmapping['taxi'] = [
                            530411,#	I	Taxi fares on trips
                            530412]#	D	Taxi fares and limo service

uccmapping['intracity'] = [
                            530311,#	I	Intracity mass transit fares
                            530312,#	I	Local transportation on out-of-town trips
                            530902]#	I	School bus

uccmapping['carpool'] = [
                            530903]#	D	Car/van pool & non-motor. trans.

uccmapping['air'] = [
                            531110]#	I	Airline fares

uccmapping['ship'] = [
                            531901]#	I	Ship fares

uccmapping['member_club'] = [
                            520901,#	I	Docking/landing fees
                            620111,#	I	Social/recreational/civic club membership
                            520903,#	I	Aircraft rental
                            520906,#	I	Aircraft rental/out-of-town trips
                            620121,#	D	Fees for participant sports
                            620122]#	I	Participant sports, out-of-town trips

uccmapping['amusement'] =[
                            621310,#	I	Fees for recreational lessons
                            620925,#	D	Miscellaneous fees
                            620913,#	D	Pinball, electronic video games
                            620906,#	I	Rental of boat
                            520907,#	I	Boat/trailer rent out of town
                            620710,#	D	Recreational camp fees
                            620903,#	I	Other entertnm. serv. on out-of-town trip
                            621904]#	I	Rent/repair music instruments


uccmapping['theater'] = [
                            620211,#	I	Movie, theater, opera, ballet (thru 2013 Q1)
                            620212]#	I	Movie, other admissions, out-of-town trips

uccmapping['live_entertain'] = [
                            620213,#	I	Play, theater, opera, concert (introduced 2013 Q2)
                            620214,#	I	Movies, parks, museums (introduced 2013 Q2) (thru 2017 Q1)
                            620215,#	I	Tickets to movies (introduced 2017 Q2)
                            620216,#	I	Tickets to parks or museums (introduced 2017 Q2)
                            680310]#	I	Live entertainment for catered affairs 

uccmapping['spec_sport'] =[
                            620221,#	I	Admission to sporting events
                            620222]#	I	Admission to sports events, out-of-town trips
                                
uccmapping['cable'] = [
                            270310,#	I	Community antenna or cable tv
                            270311,#	I	Satellite radio service 
                            520560,#	I	Global positioning service (thru 2013 Q1)
                            690330]#	I	Installation of satellite television equipment

uccmapping['film_process'] = [
                            620330]#	I	Film processing

uccmapping['photo_service'] = [ 
                            620320]#	D	Photographer fees

uccmapping['repair_tvaudio'] = [
                            340610,#	I	Repair of tv/radio/sound equipment
                            340902,#	I	Rental of televisions
                            340905,#	I	Rental of VCR/radio/sound equip
                            620905,#	I	Rent/rep photo equip
                            690320,#	I	Installation of televisions
                            690340,#	I	Installation of sound systems 
                            690350,#	I	Installation of other video equipment/sound systems
                            340908,#	I	Rental of office equipment for non-business use
                            620917,#	I	Rental video hardware/accessories (introduced 2011 Q2)
                            620916,#	I	Rental of computer and video game hardware and software (introduced 2007 Q1) (thru 2011 Q1)
                            620918]#	I	Rental video software (introduced 2011 Q2)


uccmapping['video_rental'] =[
                            620912]#	D	Rental video cassettes/tapes/discs/films

uccmapping['gambling'] =[
                            620926]#	D	Lotteries/pari-mutuel losses

uccmapping['vet_service'] = [
                            620410,#	I	Pet services
                            620420]#	D	Vet services

uccmapping['repair_rec'] =[
                            620908,#	I	Rent/repair of misc sports equipment
                            622904]#	I	Rent/repair music instruments

uccmapping['food_school'] = [
                            190901,#	I	Board (including at school)
                            790430]#	I	School lunches


uccmapping['foodaway']	= [    	    #	   BLS concordance uses diary inestead	
                            190902,	#	4	190902	Catered affairs	
                            190903,	#	4	190903	Food on out-of-town trips	
                            790410,	#	4	790410	Meals at restaurants carry-outs and other	
                            800700]		#	4	800700	Meals as pay

uccmapping['alcaway']   = [    	        #	BLS concordance uses diary instead	
                            790420,	    #	4	790420	Alcoholic beverages at restaurants taverns	
                            200900]		#	4	200900	Alcoholic beverages purchased on trips

uccmapping['hotels'] = [ 
                            210210,#	I	Lodging while out of town
                            680905]#	I	Vacation clubs 

uccmapping['school_housing'] = [
                            211310]#	I	Housing while attending school


uccmapping['fin_service'] =[
                            680210,#	I	Safe deposit box rental
                            680220,#	I	Checking accounts, other bank services
                            620112]#	I	Credit card membership

uccmapping['life_insur'] = [
                            700110]#	I	Life, endow., annuity,insurance

uccmapping['home_insur'] = [
                            220121,#	I	Homeowners insurance (ownd)
                            220122,#	I	Homeowners insurance (ownv)
                            350110]#	I	Tenant's insurance

uccmapping['health_insur'] = [
                            580111,#	I	Traditional fee for service health plan (not BCBS) (thru 2013 Q1)
                            580112,#	I	Traditional fee for service health plan (BCBS) (thru 2013 Q1)
                            580115,#	I	Fee for service without BCBS (introduced 2013 Q2)
                            580116,#	I	Fee for service with BCBS (introduced 2013 Q2)
                            580113,#	I	Preferred provider health plan (not BCBS) (thru 2013 Q1)
                            580114,#	I	Preferred provider health plan (BCBS) (thru 2013 Q1)
                            580311,#	I	Health maintenance organizaiton (not BCBS)
                            580312,#	I	Health maintenance organization (BCBS)
                            580400,#	I	Long term care insurance (thru 2017 Q1)
                            580401,#	I	Long term care insurance (not BCBS) (introduced 2017 Q2)
                            580402,#	I	Long term care insurance (BCBS) (introduced 2017 Q2)
                            580903,#	I	Commercial Medicare supplement (not BCBS)
                            580904,#	I	Commercial Medicare supplement (BCBS)
                            580905,#	I	Other health insurance (not BCBS)
                            580411,#	I	Dental care insurance (not BCBS) (introduced 2017 Q2)
                            580431,#	I	Vision care insurance (not BCBS) (introduced 2017 Q2)
                            580441,#	I	Other single service insurance (not BCBS) (introduced 2017 Q2)
                            580421,#	I	Prescription drug insurance (not BCBS) (introduced 2017 Q2)
                            580906,#	I	Other health insurance (BCBS) (thru 2017 Q1)
                            580412,#	I	Dental care insurance (BCBS) (introduced 2017 Q2)
                            580432,#	I	Vision care insurance (BCBS) (introduced 2017 Q2)
                            580442,#	I	Other single service insurance (BCBS) (introduced 2017 Q2)
                            580422,#	I	Prescription drug insurance (BCBS) (introduced 2017 Q2)
                            580908,#	I	Medicaid premiums (introduced 2017 Q2)
                            580909,#	I	Tricare/military premiums (introduced 2017 Q2)
                            580910]#	I	Children's Health Insurance Program (CHIP) premiums (introduced 2017 Q2)
                            
uccmapping['auto_insur'] = [
                            500110,#	D	Auto insurance 
                            490900]#	I	Auto repair service policy

uccmapping['tel_service'] = [
                            270105,#	I	Voice over IP service
                            270101,#	I	Residential telephone/pay phones (thru 2013 Q1)
                            270106]#	I	Residential telephone inc. VOIP (introduced 2013 Q2)

uccmapping['cell_service'] = [
                            270104,#	I	Prepaid phone cards
                            270102,#	I	Cellular phone service
                            270103,#	I	Pager service (thru 2006 Q1)
                            310400]#	I	Applications, games, ringtones for handheld devices (introduced 2011 Q2)

uccmapping['post'] = [  # Only in Diary
                            340110,#	D	Postage
                            340120]#	D	Delivery services

uccmapping['internet'] = [
                            690116,#	I	Internet services away from home
                            690114]#	I	Computer information services

uccmapping['college_tuition'] =[
                            670110,#	I	College tuition
                            670902]#	I	Other school expenses including rentals

uccmapping['elem_tuition'] = [
                            670210]#	I	Elementary/high school tuition

uccmapping['preschool'] = [
                            670310]#	I	Day-care centers, nursery, and preschools

uccmapping['other_school'] = [
                            670410,#	I	Vocational and technical school tuition 
                            670901,#	I	Other school tuition
                            670903,#	I	Test preparation, tutoring services 
                            622310]#	I	Fees for recreational lessons

uccmapping['legal'] = [
                            680110]#	I	Legal fees

uccmapping['accounting'] = [
                            680902]#	I	Accounting fees

uccmapping['other_bus_service'] = [ #Only in diary
                            680903]#	D	Miscellaneous personal services

uccmapping['funeral'] = [
                            680140,#	I	Funeral expenses
                            680901]#	I	Cemetery lots, vaults, & maintenance fees

uccmapping['hair'] = [
                            650310,#	I	Personal care services
                            650110,#	I	Personal care services for females* (introduced 2014 Q1)
                            650210]#	I	Personal care services for males* (introduced 2014 Q1)

uccmapping['other_pers_service'] =[
                            440150,#	I	Watch and jewelry repair
                            650900,#	D	Repair of personal care appliances
                            620115,#	I	Shopping club membership fees
                            680904]#	I	Dating Services

uccmapping['laundry'] =[
                            440120,#	I	Coin-operated apparel laundry/dry cleaning
                            440210,#	I	Apparel laundry/dry cleaning not coin operated
                            340520,#	I	"Household laundry/dry cleaning, sent out (non-clothing) not coin-operated"
                            340530]#	D	Coin-operated household laundry/dry cleaning (non-cloth)

uccmapping['cloth_repair'] = [
                            440130,#	I	Alteration, repair and tailoring of apparel
                            440140]#	I	Clothing rental

uccmapping['shoe_repair'] = [
                            440110]#	I	Shoe repair, other shoe service

uccmapping['childcare'] = [
                            340211,#	I	Babysitting and child care in own home (thru 2013 Q1)
                            340212,#	I	Babysitting and child care in someone else's home (thru 2013 Q1)
                            340210]#	I	Babysitting and child care (introduced 2013 Q2)

uccmapping['soc_assist'] = [
                            340910,#	I	Adult day care centers
                            800811]#	I	Gifts to non-CU members of stocks and bonds

uccmapping['charity'] =[
                            800821,#	I	Cash contributions to charity
                            800851,#	I	Cash contributions to political organizations
                            800831,#	I	Cash contributions to church
                            800841]#	I	Cash contributions to education organizations

uccmapping['dom_service'] =[
                            340310]#	I	Housekeeping services

uccmapping['moving'] = [
                            340510,#	I	Moving, storage, freight express
                            440900]#	I	Clothing storage

uccmapping['furn_repair'] = [
                            340630]#	I	Reupholstering/furniture repair

uccmapping['app_repair'] = [
                            340620,#	I	Appliance repair, including service center
                            230142,#	I	Repair of disposal, b.i.d., range hood, owned
                            340907,#	I	Appliance rental
                            340913,#	D	Repair of miscellaneous household equipment and furnishings
                            990900,#	I	"Rental and installation of dishwashers, range hoods,and garbage disposals"
                            570902]#	D	Repair of medical equipment

uccmapping['other_house_service'] = [
                            340420,#	I	Water softening service
                            340903,#	I	Other home services
                            340914,#	I	Services for termite/pest control
                            340911,#	I	Management and upkeep services for security, owner
                            340912,#	I	Management and upkeep services for security, own. vac.
                            790640,#	I	Management, security, parking, other props.
                            340915]#	I	Home security system service fee

uccmapping['foreign_travel'] = [
                            532110,#	I	Airline fares
                            532901]#	I	Ship fares


#These last categories are created in the NIPA file so that we can get the same observations as the PSMJ file
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
                                            


# save output to yaml file
with open('../output/ucc_nipa_map.yml', 'w') as yamlfile:
    yaml.dump(uccmapping, yamlfile, default_flow_style=False) #, f, default_flow_style=False