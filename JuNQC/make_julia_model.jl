using ArgParse
using JSON

include("./src/Types.jl")
include("./src/Macros.jl")
include("./src/Parser.jl")
include("./src/Problem.jl")
include("./src/strategy/JuliaStrategy.jl")
include("./src/Common.jl")

# Grab the required functions for code generation -
# const parser_function::Function = parse_vff_file

function parse_commandline()
    settings_object = ArgParseSettings()
    @add_arg_table settings_object begin
      "-o"
        help = "Directory where the Matlab model files will be written."
        arg_type = AbstractString
        default = "."

      "-m"
        help = "Path to the biochemical reaction file written in the vff format."
        arg_type = AbstractString
        required = true
    end

    # return a dictionary w/args -
    return parse_args(settings_object)
end


function main()

  # Build the arguement dictionary -
  parsed_args = parse_commandline()

  # Load the statement_vector -
  path_to_model_file = parsed_args["m"]
  metabolic_statement_vector::Array{VFFSentence} = parse_vff_metabolic_statements(path_to_model_file)

  # Load the JSON configuration file -
  config_dict = JSON.parsefile("./config/Configuration.json")

  # Generate the problem object -
  problem_object = generate_problem_object(metabolic_statement_vector,config_dict)
  problem_object.configuration_dictionary = config_dict

  # Write the DataDictionary -
  component_set = Set{ProgramComponent}()

  # write debug buffer -
  program_component_debug = build_debug_buffer(problem_object)
  push!(component_set,program_component_debug)
  progress_message = "constructed the debug buffer ..."
  println(progress_message)

  # data dictionary -
  program_component_data_dictionary = build_data_dictionary_buffer(problem_object)
  push!(component_set,program_component_data_dictionary)
  progress_message = "constructed the data_dictionary buffer ..."
  println(progress_message)

  # Write the stoichiometric_matrix --
  program_component_stoichiometric_matrix = generate_stoichiomteric_matrix_buffer(problem_object)
  push!(component_set,program_component_stoichiometric_matrix)

  # Dump the component_set to disk -
  path_to_output_file = parsed_args["o"]
  write_program_components_to_disk(path_to_output_file,component_set)

  # Transfer distrubtion files to the output -
  transfer_distribution_file("./distribution","FluxDriver.jl",path_to_output_file,"FluxDriver.jl")
  transfer_distribution_file("./distribution","Solve.jl",path_to_output_file,"Solve.jl")
  transfer_distribution_file("./distribution","Include.jl",path_to_output_file,"Include.jl")
  transfer_distribution_file("./distribution","Utility.jl",path_to_output_file,"Utility.jl")
  transfer_distribution_file("./distribution","README_INCLUDE.md",path_to_output_file,"README.md")
  transfer_distribution_file("./distribution","Atom.csv",path_to_output_file,"Atom.txt")

end

main()
