function maximize_product_data_dictionary(time_start,time_stop,time_step)

	# Get the default data_dictionary -
	data_dictionary = DataDictionary(time_start,time_stop,time_step)

	# setup the obj -
	objective_coefficient_array = data_dictionary["objective_coefficient_array"]
	objective_coefficient_array[2259] = -1.0 #Product
	objective_coefficient_array[2257] = -1.0 #Biomass

#=
	# Default flux bounds array -
	default_flux_bounds_array = data_dictionary["default_flux_bounds_array"]
	default_flux_bounds_array[21,2] = 0.0

	# ATP maintenance -
	# default_flux_bounds_array[20,1:2] = 7.6
=#
	#
	#= AEROBIC
	# setup exchange array -
	species_bounds_array = data_dictionary["species_bounds_array"]
	exchange_array = [
		0.0	1.0	;	# 73 M_ac_b
		0.0	0.0	;	# 74 M_acald_b
		0.0	0.0	;	# 75 M_akg_b
		0.0	100.0	;	# 76 M_co2_b
		0.0	1.0	;	# 77 M_etoh_b
		0.0	1.0	;	# 78 M_for_b
		0.0	0.0	;	# 79 M_fru_b
		0.0	0.0	;	# 80 M_fum_b
		-1.0	0.0	;	# 81 M_glc_D_b
		0.0	0.0	;	# 82 M_gln_L_b
		0.0	0.0	;	# 83 M_glu_L_b
		-10.0	10.0	;	# 84 M_h2o_b
		-100.0	100.0	;	# 85 M_h_b
		0.0	0.0	;	# 86 M_lac_D_b
		0.0	0.0	;	# 87 M_mal_L_b
		-10.0	10.0	;	# 88 M_nh4_b
		-1.7	0.0	;	# 89 M_o2_b
		-10.0	10.0	;	# 90 M_pi_b
		0.0	0.0	;	# 91 M_pyr_b
		0.0	0.0	;	# 92 M_succ_b
	]



	# how many unbalanced species do we have?
	offset = 72
	(number_of_exchange_species,number_of_bounds) = size(exchange_array)
	for exchange_index = 1:number_of_exchange_species

		bounds_row_index = offset+exchange_index

		# update the lower bound -
		species_bounds_array[bounds_row_index,1] = exchange_array[exchange_index,1]

		# update the upper bound -
		species_bounds_array[bounds_row_index,2] = exchange_array[exchange_index,2]
	end
=#
	return data_dictionary
end
