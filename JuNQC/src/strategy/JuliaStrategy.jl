function build_function_header_buffer(comment_dictionary)

  # initialize -
  buffer = ""

  # get some data from the comment_dictionary -
  function_name = comment_dictionary["function_name"]
  function_description = comment_dictionary["function_description"]
  input_arg_array = comment_dictionary["input_args"]
  output_arg_array = comment_dictionary["output_args"]

  buffer *= "# ----------------------------------------------------------------------------------- #\n"
  buffer *= "# Function: $(function_name)\n"
  buffer *= "# Description: $(function_description)\n"
  buffer *= "# Generated on: $(now())\n"
  buffer *= "#\n"
  buffer *= "# Input arguments:\n"

  for argument_dictionary in input_arg_array

    arg_symbol = argument_dictionary["symbol"]
    arg_description = argument_dictionary["description"]

    # write the buffer -
    buffer *= "# $(arg_symbol) => $(arg_description) \n"
  end

  buffer *= "#\n"
  buffer *= "# Output arguments:\n"
  for argument_dictionary in output_arg_array

    arg_symbol = argument_dictionary["symbol"]
    arg_description = argument_dictionary["description"]

    # write the buffer -
    buffer *= "# $(arg_symbol) => $(arg_description) \n"
  end
  buffer *= "# ----------------------------------------------------------------------------------- #\n"

  # return the buffer -
  return buffer
end


function build_copyright_header_buffer(problem_object::ProblemObject)

  # What is the current year?
  current_year = string(Dates.year(now()))

  buffer = ""
  buffer*= "# ----------------------------------------------------------------------------------- #\n"
  buffer*= "# Copyright (c) $(current_year) Varnerlab\n"
  buffer*= "# Robert Frederick Smith School of Chemical and Biomolecular Engineering\n"
  buffer*= "# Cornell University, Ithaca NY 14850\n"
  buffer*= "#\n"
  buffer*= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
  buffer*= "# of this software and associated documentation files (the \"Software\"), to deal\n"
  buffer*= "# in the Software without restriction, including without limitation the rights\n"
  buffer*= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
  buffer*= "# copies of the Software, and to permit persons to whom the Software is\n"
  buffer*= "# furnished to do so, subject to the following conditions:\n"
  buffer*= "#\n"
  buffer*= "# The above copyright notice and this permission notice shall be included in\n"
  buffer*= "# all copies or substantial portions of the Software.\n"
  buffer*= "#\n"
  buffer*= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
  buffer*= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
  buffer*= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
  buffer*= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
  buffer*= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
  buffer*= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
  buffer*= "# THE SOFTWARE.\n"
  buffer*= "# ----------------------------------------------------------------------------------- #\n"

  # return -
  return buffer
end

function build_debug_buffer(problem_object::ProblemObject)
  filename = "Debug.txt"

  buffer = ""

  # write the list of reactions -
  list_of_reactions::Array{ReactionObject} = problem_object.list_of_reactions
  counter = 1
  for (index,reaction_object) in enumerate(list_of_reactions)

    reaction_string = reaction_object.reaction_name
    reaction_type = reaction_object.reaction_type

    # Build comment string -
    comment_string = build_reaction_comment_string(reaction_object)
    buffer *= "$(counter) $(reaction_string)::$(comment_string)\n"

    # update the counter -
    counter = counter + 1;
  end

  # type SpeciesObject
  #
  #   species_type::Symbol
  #   species_bound_type::Symbol
  #   species_symbol::AbstractString
  #   species_index::Int
  #   stoichiometric_coefficient::Float64
  #   species_compartment::Symbol
  #
  #   function SpeciesObject()
  #     this = new()
  #   end
  # end

  buffer *= "\n"

  # write the list of species -
  list_of_species::Array{SpeciesObject} = problem_object.list_of_species
  counter = 1
  for (index,species_object) in enumerate(list_of_species)

    species_symbol = species_object.species_symbol
    buffer *= "$(index) $(species_symbol)\n"
  end

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer
  program_component.component_type = :buffer

  # return -
  return (program_component)
end

function build_data_dictionary_buffer(problem_object::ProblemObject)

  filename = "DataDictionary.jl"

  # build the header -
  header_buffer = build_copyright_header_buffer(problem_object)

  # get the comment buffer -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["data_dictionary_function"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  debug_message = "DataDictionary - finished building the copyright and function headers..."
  println(debug_message)

  # Get the default -
  default_parameter_dictionary = problem_object.configuration_dictionary["default_parameter_dictionary"]
  enzyme_initial_condition = parse(Float64,default_parameter_dictionary["default_protein_initial_condition"])
  default_rate_constant = parse(Float64,default_parameter_dictionary["default_enzyme_kcat"])
  default_upper_bound = default_rate_constant*enzyme_initial_condition

  # initialize the buffer -
  buffer = ""
  buffer *= header_buffer
  buffer *= "#\n"
  buffer *= function_comment_buffer
  buffer *= "function DataDictionary(time_start,time_stop,time_step)\n"
  buffer *= "\n"
  buffer *= "\t# Load the stoichiometric network from disk - \n"
  buffer *= "\tstoichiometric_matrix = readdlm(\"Network.dat\");\n"

  # generate flux bounds array -
  buffer *= "\n"
  buffer *= "\t# Setup default flux bounds array - \n"
  default_bounds_array = build_default_flux_bounds(problem_object)
  buffer *= default_bounds_array
  buffer *= "\n"

  # generate species bounds constraints -
  buffer *= "\t# Setup default species bounds array - \n"
  buffer *= "\tspecies_bounds_array = [\n";

  list_of_species::Array{SpeciesObject} = problem_object.list_of_species
  number_of_species = length(list_of_species)
  counter = 1
  for species_object in list_of_species

    # Get the bound type, and species -
    species_bound_type = species_object.species_bound_type
    species_symbol = species_object.species_symbol

    debug_message = "Processing $(species_symbol) index $(counter) of $(number_of_species)"
    println(debug_message)

    if (species_bound_type == :unbalanced)

      # this species is unbounded -
      buffer *= "\t\t-1.0\t1.0\t;\t#"

    else

      # this species is bounded -
      buffer *= "\t\t0.0\t0.0\t;\t#"
    end

    buffer *= " $(counter) $(species_symbol)\n"
    counter = counter+1
  end

  buffer *= "\t];\n"

  debug_message = "DataDictionary - finished building the species bounds array ..."
  println(debug_message)

  buffer *= "\n"
  buffer *= "\t# Setup the objective coefficient array - \n"
  buffer *= "\tobjective_coefficient_array = [\n";

  # What species are calculated? ()
  counter = 1
  list_of_reactions::Array{ReactionObject} = problem_object.list_of_reactions
  number_of_reactions = length(list_of_reactions)
  for reaction_object in list_of_reactions

    reaction_string = reaction_object.reaction_name
    reaction_type = reaction_object.reaction_type

    # Build comment string -
    comment_string = build_reaction_comment_string(reaction_object)
    buffer *= "\t\t0.0\t;\t# $(counter) $(reaction_string)::$(comment_string)\n"

    debug_message = "Processing $(reaction_string) index $(counter) of $(number_of_reactions)"
    println(debug_message)

    # update the counter -
    counter = counter + 1;
  end

  buffer *= "\t];\n"

  debug_message = "DataDictionary - finished building the reactions bounds array ..."
  println(debug_message)

  # set min/max flag -
  buffer *= "\n"
  buffer *= "\t# Min/Max flag - default is minimum - \n"
  buffer *= "\tis_minimum_flag = true\n"

  # write the reaction string list -
  counter = 1
  buffer *= "\n"
  buffer *= "\t# List of reation strings - used to write flux report \n"
  buffer *= "\tlist_of_reaction_strings = [\n"
  for reaction_object in list_of_reactions

    reaction_string = reaction_object.reaction_name
    reaction_type = reaction_object.reaction_type

    debug_message = "Processing $(reaction_string) index $(counter) of $(number_of_reactions)"
    println(debug_message)

    # Build comment string -
    comment_string = build_reaction_comment_string(reaction_object)
    buffer *= "\t\t\"$(reaction_string)::$(comment_string)\"\n"

    # update the counter -
    counter = counter + 1;
  end

  buffer *= "\t];\n"

  debug_message = "DataDictionary - finished building the reaction string array ..."
  println(debug_message)

  # wtite list of metabolite symbols -
  counter = 1
  buffer *= "\n"
  buffer *= "\t# List of metabolite strings - used to write flux report \n"
  buffer *= "\tlist_of_metabolite_symbols = [\n"
  for species_object in list_of_species

    # Get the bound type, and species -
    species_bound_type = species_object.species_bound_type
    species_symbol = species_object.species_symbol
    buffer *= "\t\t\"$(species_symbol)\"\n"

    # debug -
    debug_message = "Processing species symbol $(species_symbol) (index $(counter) of $(number_of_species))"
    println(debug_message)
    counter = counter + 1
  end
  buffer *= "\t];\n"

  # return block -
  buffer *= "\n"
  buffer *= "\t# =============================== DO NOT EDIT BELOW THIS LINE ============================== #\n"
  buffer *= "\tdata_dictionary = Dict{AbstractString,Any}()\n"
  buffer *= "\tdata_dictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix\n"
  buffer *= "\tdata_dictionary[\"objective_coefficient_array\"] = objective_coefficient_array\n"
  buffer *= "\tdata_dictionary[\"default_flux_bounds_array\"] = default_bounds_array;\n"
  buffer *= "\tdata_dictionary[\"species_bounds_array\"] = species_bounds_array\n"
  buffer *= "\tdata_dictionary[\"list_of_reaction_strings\"] = list_of_reaction_strings\n"
  buffer *= "\tdata_dictionary[\"list_of_metabolite_symbols\"] = list_of_metabolite_symbols\n"
  buffer *= "\tdata_dictionary[\"is_minimum_flag\"] = is_minimum_flag\n"
  buffer *= "\t# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #\n"
  buffer *= "\treturn data_dictionary\n"
  buffer *= "end\n"

  # build the component -
  program_component::ProgramComponent = ProgramComponent()
  program_component.filename = filename
  program_component.buffer = buffer
  program_component.component_type = :buffer

  # return -
  return (program_component)
end

function build_default_flux_bounds(problem_object::ProblemObject)

  # get the comment buffer -
  comment_header_dictionary = problem_object.configuration_dictionary["function_comment_dictionary"]["data_dictionary_function"]
  function_comment_buffer = build_function_header_buffer(comment_header_dictionary)

  # Get the default -
  default_parameter_dictionary = problem_object.configuration_dictionary["default_parameter_dictionary"]
  enzyme_initial_condition = parse(Float64,default_parameter_dictionary["default_protein_initial_condition"])
  default_rate_constant = parse(Float64,default_parameter_dictionary["default_enzyme_kcat"])
  default_upper_bound = default_rate_constant*enzyme_initial_condition

  # initialize the buffer -
  buffer = ""

  # open -
  buffer *= "\tdefault_bounds_array = [\n"

  # Grab the list of reactions -
  list_of_reactions::Array{ReactionObject} = problem_object.list_of_reactions
  counter = 1
  for (index,reaction_object) in enumerate(list_of_reactions)

    reaction_string = reaction_object.reaction_name
    reaction_type = reaction_object.reaction_type

    # Generate the comment string -
    comment_string = build_reaction_comment_string(reaction_object)

    # build bounds record -
    buffer *= "\t\t0\t$(default_upper_bound)\t;\t# $(counter) $(comment_string)\n"

    # update counter -
    counter = counter + 1

  end

  # close -
  buffer *= "\t];\n"
end
