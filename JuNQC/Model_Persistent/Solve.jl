include("Include.jl")

# load the data dictionary -
#data_dictionary = DataDictionary(0,0,0)

# solve the lp problem -
#(objective_value, flux_array, dual_array, uptake_array, exit_flag) = FluxDriver(data_dictionary)


known_species_array = ParseXML()
#unpack known_species_array
species_name = known_species_array[1]
species_id = known_species_array[2]
species_url = known_species_array[3]

atom_array = make_atom_array(species_url)

print(atom_array)

#=
wanted_species = ["Adenine","Cytosine","Guanine","Uracil","Glycine","L-Alanine","L-Arginine","L-Asparagine","L-Aspartate","L-Cysteine","L-Glutamate","L-Histidine","L-Isoleucine","L-Leucine","L-Lysine","L-Methionine","L-Phenylalanine","L-Proline","L-Serine","L-Tryptophan","L-Tyrosine","L-Valine","L-Glutamine"]
wanted_ids = zeros(1,length(wanted_species))

#for el in wanted_species
#  wanted_ids[el] = findID(species_name,species_id,"Adenine")
#end
=#
#missing_species = find_missing_species()


#url = "http://identifiers.org/chebi/CHEBI:12266"#Debugging
