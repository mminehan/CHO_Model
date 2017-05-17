To run fed-batch simulation, navigate to Model_Persistent and run include("SimDriver.jl") from Julia

tPA.jl includes data and functions for turning sequence information into reactions for Network.net

Model_Working is meant to be overwritten with the generation of new models from Network.net, and the relevant files (DataDictionary.jl, Network.dat) be copied into Model_Persistent---The DataDictionary file will need to be updated with the maximize_product_data_dictionary function

The simulation and atom array generation files are located in Utility.jl, and are titled "runSim" and "make_atom_array" respectively.
