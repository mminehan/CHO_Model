function runSim(time_array,state_array,data_dictionary,feed_flag)

number_of_timesteps = length(time_array)

# initialize the problem -
number_of_external_states = length(state_array)

# Problem specific kinetic parameters -
vmax_glucose_uptake = 10.5
K_glucose_uptake = 36.0

#Distinguish between feed and no feed conditions
if feed_flag == true
  F_in = 0.005 #Feed rate (L/min)
  F_out = 0.005 #Removal rate (L/min)
  feed_array = zeros(number_of_external_states - 1) #Exclude Volume
  feed_array[3] = 115.0 #concentration of glucose in feed
else
  F_in = 0.0 #Feed rate (L/min)
  F_out = 0.0 #Removal rate (L/min)
  feed_array = zeros(number_of_external_states - 1) #Exclude Volume
  feed_array[3] = 0.0 #concentration of glucose in feed
end

# capture the exit flags -
exit_flag_array = Int[]

# main loop -
for time_step_index = 1:number_of_timesteps-1

  # make a deepcopy of the data_dictionary -
  copy_data_dictionary = deepcopy(data_dictionary)

  # grab the state -
  product = state_array[time_step_index,1]
  glucose = state_array[time_step_index,3]
  product_mRNA = state_array[time_step_index,4]
  volume = state_array[time_step_index,5]
  #acetate = state_array[time_step_index,2]
  cellmass = state_array[time_step_index,2]
  #formate = state_array[time_step_index,4]
  #ethanol = state_array[time_step_index,5]

  # Default flux bounds array -
	v_TX = 954.5*60/1629 #BP/min CHO cell, 1629 = mRNA length
	v_TL = 300*60/526 #AA/min Homo Sapiens, 526 = protein length
	km = 10.0^-3 #Polysome density?
	kd = 0.1 #mRNA degradation constant, this one's from Homo Sapiens
	default_bounds_array = data_dictionary["default_flux_bounds_array"]
	default_bounds_array[2258,2] = v_TX
	default_bounds_array[2259,2] = v_TL*(product_mRNA/kd)/(km+(product_mRNA/kd))
  data_dictionary["default_flux_bounds_array"] = default_bounds_array

  # calculate glucose uptake -
  qGlc = vmax_glucose_uptake*(glucose)/(K_glucose_uptake+glucose)

  @show qGlc

  # calculate volume increase (dilution)
  qV = F_in - F_out

  # setup the species bounds -
  # Beta-Glucose-6-Phosphate_c1 is s_0055
  # Alpha-D-Glucose_c1 is s_0548 THIS IS UPTAKEN GLUCOSE
  # Extracelluar Glucose is s_1248
  species_bounds_array = copy_data_dictionary["species_bounds_array"]
  #species_bounds_array[1222:1353] = replaceBounds(species_bounds_array)
  species_bounds_array[549,1] = -qGlc
  species_bounds_array[549,2] = -0.99*qGlc

  # calculate the fluxes using the LP -
  (objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(copy_data_dictionary)

  k_g = 5
  # grab the fluxes from the flux_array -
  qProd = flux_array[2259]
  qMRNA = flux_array[2258]
  mu = 0.03466*glucose/(k_g+glucose)
  k_death = 0.03466*0.005

  @show mu

  #flux_array[2257]
  #qAcetate_production = flux_array[35]
  #qFormate_production = flux_array[41]
  #qEthanol_production = flux_array[40]

  # update the external state -
  state_array[time_step_index+1,1] = product + qProd*cellmass + (F_in/volume)*feed_array[1]-(F_out/volume)*product - (qV/volume)*product
  state_array[time_step_index+1,2] = cellmass + mu*cellmass - k_death*cellmass + (F_in/volume)*feed_array[2]-(F_out/volume)*cellmass - (qV/volume)*cellmass
  state_array[time_step_index+1,3] = glucose -1.0*qGlc*cellmass + (F_in/volume)*feed_array[3]-(F_out/volume)*glucose - (qV/volume)*glucose
  state_array[time_step_index+1,4] = product_mRNA + qMRNA*cellmass + (F_in/volume)*feed_array[4]-(F_out/volume)*product_mRNA - (qV/volume)*product_mRNA
  state_array[time_step_index+1,5] = volume + qV

  #state_array[time_step_index+1,2] = acetate + qAcetate_production*cellmass
  #state_array[time_step_index+1,3] = cellmass + mu*cellmass
  #state_array[time_step_index+1,4] = formate + qFormate_production*cellmass
  #state_array[time_step_index+1,5] = ethanol + qEthanol_production*cellmass

  # correct negatives -
  idx_nz = find(state_array[time_step_index+1,:].<0)
  state_array[time_step_index+1,idx_nz] = 0.0

  # capture the exit flag -
  push!(exit_flag_array,exit_flag)
end

out = (time_array,state_array)
end

function make_atom_array(species_url)

#FUNCTION TO GENERATE ATOM MATRIX FROM URL

element_array = ["C","H","N","O","P","S"]
atom_array = zeros(length(species_url),length(element_array))

for num_species = 1:length(species_url)
  id = String[]
  url = species_url[num_species]
  for nn = 5:-1:1
    push!(id,string(url[end-nn+1]))
  end
  id_str = join(id)

  query_url = "http://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity?chebiId="*id_str

  a = get(query_url)
  main = parse_string(readstring(a))
  trunk = root(main)

  formula = String[]

  for c1 in child_nodes(trunk)
    header1 = XMLElement(c1) #Body
    println(name(header1))
    if name(header1) == "Body"
      for c2 in child_nodes(header1)
        header2 = XMLElement(c2) #getCompleteEntityResponse
        println(name(header2))
        if name(header2) == "getCompleteEntityResponse"
          for c3 in child_nodes(header2)
            header3 = XMLElement(c3) #return
            println(name(header3))
            if name(header3) == "return"
              for c4 in child_nodes(header3)
                header4 = XMLElement(c4) #Formulae
                if name(header4) == "Formulae"
                  println(name(header4))
                  for c5 in child_nodes(header4)
                    header5 = XMLElement(c5) #data
                    if name(header5) == "data"
                      println(name(header5))
                      push!(formula,string(header5))
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
  # sample formula: "<data>C9H16NO3</data>"
  #<data> + </data> = 6 + 7 = 13 chars

  #if length(formula) == 0
  #  break
  #end

  try #Found errors: length(formula) = 0, formula has characters like (C3H16)n
    split_form = split(formula[1],"")
    form_string = join(split_form[7:length(split_form)-7])

    num_ind = []

    for j = 1:length(element_array) #Iterate thru elements
      for n = 1:length(form_string) #Search formula for each element
        num_flag = isnumeric(string(form_string[n])) # True/false for number/nonnumber
        try
          num_flag_plus = isnumeric(string(form_string[n+1]))
        end
        push!(num_ind,num_flag) #Create array of indicies for number locations

        if !num_flag & !num_flag_plus #Capture the 1 case (N in C6H12NO2)
          if string(form_string[n]) == element_array[j]
            atom_array[num_species,j] = 1.0
          end
        end

        if num_flag & num_flag_plus #Capture the double-digit case
          if string(form_string[n-1]) == element_array[j]
            atom_array[num_species,j] = float(string(form_string[n])*string(form_string[n+1]))
          end

        elseif num_flag #Capture the single-digit case
          if string(form_string[n-1]) == element_array[j]
            atom_array[num_species,j] = float(string(form_string[n]))
          end
        end
        num_flag_plus = false
      end
    end
  end
end
return atom_array
#=
#Try-catch method--BROKEN
for n = 7:(length(split_form)-7)
  if split_form[n] == "C"
    temp_c1 = float(split_form[n+1])
    try
      temp_c1 = float(split_form[n+1])
      try
        temp_c2 = float(split_form[n+2])
        temp_c1 = float(split_form[n+1]*split_form[n+2])
      end
    catch
      temp_c1 = 1.0
    end
  end
  println(temp_c1)
end
=#
#
end

function isnumeric{T<:String}(s::T)
   isa(parse(s), Number)
end

function replaceBounds(species_bounds_array)
  for nn in 1222:1353
    species_bounds_array[nn,1]=-1.0
    species_bounds_array[nn,2]=1.0
  end
  return species_bounds_array
end

function findID(species_name,species_id,name)
  ##BROKEN
#assuming c1 compartment entries are near top, works backwords to ID
#cytosolic metabolite entries
for ind in 1:length(species_name)
  if species_name[(length(species_name) - ind +1)] == name
    a_id = species_id[(length(species_name) - ind +1)]
  end
end
return a_id

#=
c = 0
for ind in 1:length(species_name)
  if species_name[(length(species_name) - ind +1)] == "Adenine"
    a_id = species_id[(length(species_name) - ind +1)]
  #  c = 1
  end
  if species_name[(length(species_name) - ind +1)] == "Cytosine"
    c_id = species_id[(length(species_name) - ind +1)]
  end
  if species_name[(length(species_name) - ind +1)] == "Guanine"
    g_id = species_id[(length(species_name) - ind +1)]
  end
  if species_name[(length(species_name) - ind +1)] == "Uracil"
    u_id = species_id[(length(species_name) - ind +1)]
  end
end
=#

end

function searchall(s, t, overlap::Bool=false)
    idxfcn = overlap ? first : last
    r = search(s, t)
    idxs = Array(typeof(r), 0) # Or to only count: n = 0
    while last(r) > 0
        push!(idxs, r) # n += 1
        r = search(s, t, idxfcn(r) + 1)
    end
    idxs # return n
end

function getWebData(url)

request = get(url*".json")
request_struct = JSON.parse(request.data)

formula_array = String[]

for children in request_struct["data"]["children"]
  for child in children["data"]
    if first(child) == "Formula"
      push!(formula_array,last(child))
    end
  end
end
return formula_array
end
function find_missing_species()
# This function finds the species_id's of the metabolites without an included CHEBI page
id_list = String[]

for n in 1:986
  num = string(dec(n,4))
  push!(id_list,"s_"*num)
end

known_species_array = ParseXML()
known_species_id = known_species_array[2]

out = deleteat!(id_list, findin(id_list, known_species_id))

return out
end

function getFormula(url)



end

function ParseXML()
  #Gets all species with CHEBI pages from the XML file
  #Outputs [name,id,url]
main = parse_file("cho_cobra.xml")

trunk = root(main)
println(name(trunk))
url_id = String[]
known_species_id = String[]
known_species_name = String[]

for c1 in child_nodes(trunk)  # c is an instance of XMLNode
    println(nodetype(c1))
    if is_elementnode(c1)
        branch = XMLElement(c1)  # this makes an XMLElement instance
        #println(name(branch))
        for c2 in child_nodes(branch)
          if is_elementnode(c2)
              branch2 = XMLElement(c2)  # this makes an XMLElement instance
              #println(name(branch2))
              if name(branch2) == "listOfSpecies"
                #println(nodetype(branch2))
                for c3 in child_nodes(branch2) #ITERATE THRU SPECIES HERE
                  if is_elementnode(c3)
                  species = XMLElement(c3)
                  #println(name(species))
                  if name(species) == "species"
                    for c4 in child_nodes(species)
                      if is_elementnode(c4)
                        annotation = XMLElement(c4)
                        #println(name(annotation))
                        if name(annotation) == "annotation"
                          for c5 in child_nodes(annotation)
                            if is_elementnode(c5)
                              RDF = XMLElement(c5)
                              #println(name(RDF))
                              if name(RDF) == "RDF"
                                for c6 in child_nodes(RDF)
                                  if is_elementnode(c6)
                                    Description = XMLElement(c6)
                                    #println(name(RDF))
                                    if name(Description) == "Description"
                                      #println(Description["is"])
                                      for c7 in child_nodes(Description)
                                        if is_elementnode(c7)
                                          is = XMLElement(c7)
                                          #println(name(is))
                                          if name(is) == "is"
                                            for c8 in child_nodes(is)
                                              if is_elementnode(c8)
                                                Bag = XMLElement(c8)
                                                #println(name(Bag))
                                                #println(Bag["li"])
                                                for c9 in child_nodes(Bag)
                                                  if is_elementnode(c9)
                                                    li = XMLElement(c9)
                                                    println(attribute(li,"resource"))
                                                    push!(url_id,attribute(li,"resource"))
                                                    push!(known_species_id,attribute(species,"id"))
                                                    push!(known_species_name,attribute(species,"name"))
                                                  end
                                                end
                                              end
                                            end
                                          end
                                        end
                                      end
                                    end
                                  end
                                end
                              end
                            end
                          end
                        end
                      end
                      end
                    end
                  end
                end
              end
          end
        end
    end
end

#species_name = union(species_name)
#species_id = union(species_id)
#url_id = union(url_id)

out = [known_species_name,known_species_id,url_id]
return out

end

function show_eigenreaction_profile(eigenreaction_array_column::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # what species coefficients > epsilon?
  idx_cutoff = find(abs(eigenreaction_array_column).>epsilon)

  # get the list of species -
  list_of_species_symbols = data_dictionary["list_of_metabolite_symbols"]

  # create my species list -
  eigenspecies_list = []
  for species_index in idx_cutoff

      value = list_of_species_symbols[species_index]
      record = "$(value)"
      push!(eigenspecies_list,record)
  end

  return eigenspecies_list
end

function show_eigenconnectivity_profile(eigenconnection_array_column::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # what species coefficients > epsilon?
  idx_cutoff = find(abs(eigenconnection_array_column).>epsilon)

  # get the list of species -
  list_of_reaction_strings = data_dictionary["list_of_reaction_strings"]

  # create my species list -
  eigenmode_array = []
  for flux_index in idx_cutoff

    # key,value -
    key = list_of_reaction_strings[flux_index]
    record = "$(flux_index),$(key)"
    push!(eigenmode_array,record)
  end

  return eigenmode_array
end

function show_flux_profile(flux_array::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # what fluxes are > epsilon?
  idx_cutoff = find(flux_array.>epsilon)

  # what is the list of reaction strings?
  list_of_reaction_strings = data_dictionary["list_of_reaction_strings"]

  # create a list of reactions?
  list_of_flux_records = String[]
  for flux_index in idx_cutoff

    # key,value -
    key = list_of_reaction_strings[flux_index]
    value = flux_array[flux_index]
    record = "$(flux_index),$(key),$(value)"
    push!(list_of_flux_records,record)

  end

  return list_of_flux_records
end

function find_missing_metabolites_in_atom_dataset(path_to_atom_file::AbstractString,data_dictionary::Dict{AbstractString,Any})

   # how many metabolite symbols do we have in *the model*?
	list_of_metabolite_symbols_model = data_dictionary["list_of_metabolite_symbols"]
	number_of_metabolites = length(list_of_metabolite_symbols_model)

  # initialize -
	atom_names_array = AbstractString[]
	missing_species_array = AbstractString[]
	tmp_array::Array{AbstractString} = AbstractString[]

  # load the atom file -
  try

    open(path_to_atom_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    # build atom_names_array
    for record in tmp_array
	    @show record

      # split -
      split_array = split(record,",")

      # get my key -
      key = split_array[1]  # Metabolite symbol -
	    @show key
	     push!(atom_names_array, key)
    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  #check if specise are in atom_names_array
	for species in list_of_metabolite_symbols_model
		if(!in(species,atom_names_array))
			push!(missing_species_array, species)
		end
	end

  return missing_species_array
end

function generate_atom_matrix(path_to_atom_file::AbstractString,data_dictionary::Dict{AbstractString,Any})


  # how many metabolite symbols do we have in *the model*?
  list_of_metabolite_symbols_model = data_dictionary["list_of_metabolite_symbols"]
  number_of_metabolites = length(list_of_metabolite_symbols_model)

  # initialize -
  tmp_array::Array{AbstractString} = AbstractString[]
  atom_array = zeros(number_of_metabolites,6)


  # load the atom file -
  try

    open(path_to_atom_file,"r") do model_file
      for line in eachline(model_file)

          if (contains(line,"//") == false && search(line,"\n")[1] != 1)
            push!(tmp_array,chomp(line))
          end
      end
    end

    # ok, create a local dictionary w/the atom records -
    local_dictionary::Dict{AbstractString,Any} = Dict{AbstractString,Any}()
    for record in tmp_array

      # split -
      split_array = split(record,",")

      # get my key -
      key = split_array[1]  # Metabolite symbol -

      # local array -
      local_atom_array = zeros(6)
      local_atom_array[1] = parse(Float64,split_array[2]) # C
      local_atom_array[2] = parse(Float64,split_array[3]) # H
      local_atom_array[3] = parse(Float64,split_array[4]) # N
      local_atom_array[4] = parse(Float64,split_array[5]) # O
      local_atom_array[5] = parse(Float64,split_array[6]) # P
      local_atom_array[6] = parse(Float64,split_array[7]) # S

      # store -
      local_dictionary[key] = local_atom_array
    end

    # ok, so now we have the local dictionary, we can lookup (in order) the metabolites in the model -
    for (index,model_metabolite_symbol) in enumerate(list_of_metabolite_symbols_model)

      # what is the atom array for *this metabolite*?
      local_atom_array = local_dictionary[model_metabolite_symbol]

      for (atom_index,coefficient) in enumerate(local_atom_array)
        atom_array[index,atom_index] = local_atom_array[atom_index]
      end

    end

  catch err
    showerror(STDOUT, err, backtrace());println()
  end

  return atom_array
end

function generate_mode_file_buffer(uptake_archive::Array{Float64,2},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # Grab the list of metabolites -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]
  stoichiometric_matrix = data_dictionary["stoichiometric_matrix"]

  # Write the modes mapping file -
  buffer = ""
  (number_of_species,number_of_modes) = size(uptake_archive)
  for mode_index in 1:number_of_modes

    # record -
    buffer *= "M$(mode_index)"

    # ok, for this mode, find the pivot index (if we have multiple, choose the first)
    idx_non_zero = find(uptake_archive[:,mode_index].<0)
    if (length(idx_non_zero)>1)
      idx_non_zero = idx_non_zero[1]
    end

    @show idx_non_zero

    idx_pivot = find(stoichiometric_matrix[idx_non_zero,:].<0)[1]
    buffer *= ",$(idx_pivot)"

    # grab the flux -
    reaction_string = generate_net_reaction_string(uptake_archive[:,mode_index],epsilon,data_dictionary)

    buffer *=",$(reaction_string)"

    # what species are consumed by this mode?
    idx_reactants = find(uptake_archive[:,mode_index].<0.0)
    for reactant_index in idx_reactants
      metabolite_symbol = list_of_metabolite_symbols[reactant_index]
      buffer *=",$(metabolite_symbol)"
    end

    buffer *= "\n"
  end

  return buffer
end

function generate_net_reaction_string(uptake_array::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # get list of metabolite symbols -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]

  # check for smalls -
  idx_small = find(abs(uptake_array).<epsilon)
  uptake_array[idx_small] = 0.0

  # which elememts are positive (products)?
  idx_product_array = find(uptake_array.>0)

  # which elements are negative (reactants?)
  idx_reactant_array = find(uptake_array.<0)

  # build the string ...
  net_reaction_buffer = ""
  for idx_reactant in idx_reactant_array

    metabolite_symbol = list_of_metabolite_symbols[idx_reactant]
    st_coeff = round(abs(uptake_array[idx_reactant]),2)

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # add the arrow -
  net_reaction_buffer *= " --> "

  # write the trailing stuff -
  for idx_product in idx_product_array

    metabolite_symbol = list_of_metabolite_symbols[idx_product]
    st_coeff = round(abs(uptake_array[idx_product]),2)

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # return -
  return net_reaction_buffer
end
