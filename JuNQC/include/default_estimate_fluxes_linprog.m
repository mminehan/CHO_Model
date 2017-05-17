function convex_flux_array = estimate_convex_fluxes(t,x,volume,kinetic_flux_array,data_dictionary)

	% What is my measurement tolerance?
	epsilon = 0.10;

	% Get the time step from the data_dictionary -
	Ts = data_dictionary.time_step_size;

	% Get some stuff from the data_dictionary -
	STM = data_dictionary.stoichiometric_matrix;
	[NM,NRATES] = size(STM);

	% Grab the experimental data -
	measurement_array_block = data_dictionary.experimental_data_array;

	% Call the dilution function -
	species_dilution_array = Dilution(t,x,volume,data_dictionary);

	% Formulate objective vector (default is to minimize fluxes)
	objective_vector = data_dictionary.objective_coefficient_array;

	% Get the index vectors from the data dictionary -
	index_vector_convex_species = data_dictionary.index_vector_convex_species;
	index_vector_convex_rates = data_dictionary.index_vector_convex_rates;
	index_vector_measured_species = data_dictionary.index_vector_measured_species;
	index_vector_kinetic_rates = data_dictionary.index_vector_kinetic_rates;
	index_vector_free_species = data_dictionary.index_vector_free_species;

	% = EQUALITY CONSTRAINTS ======================================================= %
	% Formulate the convex_array_block (eqality constrainta) -
	AEq = STM(index_vector_convex_species,index_vector_convex_rates);

	% Formulate the right-hand-side vector for the Equality constrainta -
	free_array_block = STM(index_vector_convex_species,index_vector_kinetic_rates);
	bEq = -1*free_array_block*kinetic_flux_array;
	% ============================================================================== %

	% = INEQUALITY CONSTRAINTS ===================================================== %
	% The inequality constraints are used to capture the measurements.
	% How many measurements do we have?
	number_of_measured_species = length(index_vector_measured_species);

	% Get the STM associated w/the measured species -
 	measured_stm_block = STM(index_vector_measured_species,index_vector_convex_rates);
	A = [
		measured_stm_block 			;	% upper bound
		-1*measured_stm_block		;	% lower bound
	];

	% Re-order the measursement array -
  measurement_selection_index_vector = data_dictionary.measurement_selection_index_vector;

	% We need to interpolate the measurement_array_block - (interpolate 1 step ahead)
	interpolated_measurements = interp1(60*measurement_array_block(:,1),measurement_array_block(:,measurement_selection_index_vector+1),t+Ts);

	% upper bound inequality constraints -
	number_of_free_species = length(index_vector_free_species);
	for measured_species_index = 1:number_of_measured_species

		% Get the measured and simulated species -
		measured_value = interpolated_measurements(measured_species_index);
		simulated_value = x(number_of_free_species+measured_species_index);
		feed_source_term = species_dilution_array(number_of_free_species+measured_species_index);

		% Calc the upper and lower bounds -
		upper_ineq(measured_species_index,1) = (1/Ts)*(measured_value*(1+epsilon)-simulated_value) - feed_source_term;
		lower_ineq(measured_species_index,1) = (1/Ts)*(measured_value*(1-epsilon)-simulated_value) - feed_source_term;
	end

	% create b -
	b = [
		upper_ineq		;	% upper bounds
		-1*lower_ineq	;	% lower bounds
	];
	% ============================================================================== %

	% = FLUX BOUNDS ================================================================ %
	default_bounds_array = data_dictionary.default_flux_bounds_array;
	LB = default_bounds_array(:,1);
	UB = default_bounds_array(:,2);
	% ============================================================================== %

	% What are my options?
	options = optimset('TolFun',1e-6,'Display','none');

  % Call the LP solver -
  [convex_flux_array,fVal,status,OUT,LAM] = linprog(objective_vector,A,b,AEq,bEq,LB,UB,[],options);

	% Shutdown -
	if (status ~= 1)
		number_of_convex_rates = length(index_vector_convex_rates);
		convex_flux_array = zeros(number_of_convex_rates,1);
	end

return;
