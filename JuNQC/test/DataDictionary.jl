# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #
#
# ----------------------------------------------------------------------------------- #
# Function: DataDictionary
# Description: Holds simulation and model parameters as key => value pairs in a Julia Dict()
# Generated on: 2017-02-02T10:48:46.296
#
# Input arguments:
# time_start::Float64 => Simulation start time value (scalar)
# time_stop::Float64 => Simulation stop time value (scalar)
# time_step::Float64 => Simulation time step (scalar)
#
# Output arguments:
# data_dictionary::Dict{AbstractString,Any} => Dictionary holding model and simulation parameters as key => value pairs
# ----------------------------------------------------------------------------------- #
function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("Network.dat");

	# Setup default flux bounds array -
	default_bounds_array = [
		0	100.0	;	# 1 M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c+M_h_c
		0	100.0	;	# 2 M_g6p_c --> M_f6p_c
		0	100.0	;	# 3 M_f6p_c --> M_g6p_c
		0	100.0	;	# 4 M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c+M_h_c
		0	100.0	;	# 5 M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0	100.0	;	# 6 M_fdp_c --> M_dhap_c+M_g3p_c
		0	100.0	;	# 7 M_dhap_c+M_g3p_c --> M_fdp_c
		0	100.0	;	# 8 M_dhap_c --> M_g3p_c
		0	100.0	;	# 9 M_g3p_c --> M_dhap_c
		0	100.0	;	# 10 M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0	100.0	;	# 11 M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0	100.0	;	# 12 M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0	100.0	;	# 13 M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0	100.0	;	# 14 M_3pg_c --> M_2pg_c
		0	100.0	;	# 15 M_2pg_c --> M_3pg_c
		0	100.0	;	# 16 M_2pg_c --> M_h2o_c+M_pep_c
		0	100.0	;	# 17 M_h2o_c+M_pep_c --> M_2pg_c
		0	100.0	;	# 18 M_adp_c+M_h_c+M_pep_c --> M_atp_c+M_pyr_c
		0	100.0	;	# 19 M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0	100.0	;	# 20 M_co2_c+M_h2o_c+M_pep_c --> M_h_c+M_oaa_c+M_pi_c
		0	100.0	;	# 21 M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c
		0	100.0	;	# 22 M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+2.0*M_h_c+M_pep_c+M_pi_c
		0	100.0	;	# 23 M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0	100.0	;	# 24 M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0	100.0	;	# 25 M_6pgl_c+M_h2o_c --> M_6pgc_c+M_h_c
		0	100.0	;	# 26 M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c
		0	100.0	;	# 27 M_ru5p_D_c --> M_xu5p_D_c
		0	100.0	;	# 28 M_xu5p_D_c --> M_ru5p_D_c
		0	100.0	;	# 29 M_r5p_c --> M_ru5p_D_c
		0	100.0	;	# 30 M_ru5p_D_c --> M_r5p_c
		0	100.0	;	# 31 M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0	100.0	;	# 32 M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 33 M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0	100.0	;	# 34 M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0	100.0	;	# 35 M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0	100.0	;	# 36 M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0	100.0	;	# 37 M_6pgc_c --> M_2ddg6p_c+M_h2o_c
		0	100.0	;	# 38 M_2ddg6p_c --> M_g3p_c+M_pyr_c
		0	100.0	;	# 39 M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c+M_h_c
		0	100.0	;	# 40 M_cit_c --> M_icit_c
		0	100.0	;	# 41 M_icit_c --> M_cit_c
		0	100.0	;	# 42 M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c
		0	100.0	;	# 43 M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c
		0	100.0	;	# 44 M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0	100.0	;	# 45 M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0	100.0	;	# 46 M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c
		0	100.0	;	# 47 M_fum_c+M_h2o_c --> M_mal_L_c
		0	100.0	;	# 48 M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c
		0	100.0	;	# 49 2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_h_e
		0	100.0	;	# 50 4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_h_e
		0	100.0	;	# 51 2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_h_e
		0	100.0	;	# 52 M_adp_c+M_pi_c+4.0*M_h_e --> M_atp_c+3.0*M_h_c+M_h2o_c
		0	100.0	;	# 53 3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_h_e
		0	100.0	;	# 54 M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0	100.0	;	# 55 M_nadh_c+M_nadp_c+2.0*M_h_e --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0	100.0	;	# 56 M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c
		0	100.0	;	# 57 M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c
		0	100.0	;	# 58 M_atp_c+M_pi_c --> M_adp_c+M_ppi_c
		0	100.0	;	# 59 M_ppi_c+M_h2o_c --> 2.0*M_pi_c+M_h_c
		0	100.0	;	# 60 M_icit_c --> M_glx_c+M_succ_c
		0	100.0	;	# 61 M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_h_c+M_mal_L_c
		0	100.0	;	# 62 M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c
		0	100.0	;	# 63 M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c
		0	100.0	;	# 64 M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0	100.0	;	# 65 M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0	100.0	;	# 66 M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0	100.0	;	# 67 M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0	100.0	;	# 68 M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c
		0	100.0	;	# 69 M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c
		0	100.0	;	# 70 M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c
		0	100.0	;	# 71 M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c
		0	100.0	;	# 72 M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c
		0	100.0	;	# 73 M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0	100.0	;	# 74 M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c
		0	100.0	;	# 75 M_glu_L_c+M_accoa_c+4.0*M_atp_c+M_nadph_c+3.0*M_h2o_c+M_gln_L_c+M_asp_L_c+M_co2_c --> M_arg_L_c+M_coa_c+5.0*M_h_c+3.0*M_adp_c+3.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c
		0	100.0	;	# 76 M_fum_c+M_nh4_c --> M_asp_L_c
		0	100.0	;	# 77 M_asp_L_c --> M_fum_c+M_nh4_c
		0	100.0	;	# 78 M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c
		0	100.0	;	# 79 M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_h_c+M_ppi_c+M_amp_c
		0	100.0	;	# 80 M_asp_L_c+M_atp_c+M_nh4_c --> M_asn_L_c+M_h_c+M_ppi_c+M_amp_c
		0	100.0	;	# 81 M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_h_c+M_ac_c
		0	100.0	;	# 82 M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c
		0	100.0	;	# 83 M_akg_c+M_nadph_c+M_nh4_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0	100.0	;	# 84 M_glu_L_c+M_atp_c+M_nh4_c --> M_gln_L_c+M_h_c+M_adp_c+M_pi_c
		0	100.0	;	# 85 M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c
		0	100.0	;	# 86 M_gln_L_c+M_r5p_c+3.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_adp_c+2.0*M_nadh_c+M_pi_c+2.0*M_ppi_c+6.0*M_h_c
		0	100.0	;	# 87 M_thr_L_c+2.0*M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh4_c+M_co2_c+M_nadp_c+M_akg_c
		0	100.0	;	# 88 2.0*M_pyr_c+M_glu_L_c+M_h_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
		0	100.0	;	# 89 M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
		0	100.0	;	# 90 M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh4_c+M_pyr_c
		0	100.0	;	# 91 M_chor_c+M_h_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
		0	100.0	;	# 92 M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c
		0	100.0	;	# 93 M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
		0	100.0	;	# 94 M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c
		0	100.0	;	# 95 M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+2.0*M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+2.0*M_adp_c+M_h_c
		0	100.0	;	# 96 M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c
		0	100.0	;	# 97 2.0*M_pyr_c+2.0*M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
		0	100.0	;	# 98 M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+2.0*M_nad_c+M_akg_c --> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+2.0*M_nadh_c+M_succ_c
		0	100.0	;	# 99 M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+M_nad_c+M_nadp_c+M_akg_c --> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+M_nadh_c+M_nadph_c+M_succ_c
		0	100.0	;	# 100 M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+2.0*M_nadp_c+M_akg_c --> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+2.0*M_nadph_c+M_succ_c
		0	100.0	;	# 101 M_arg_L_c+3.0*M_h2o_c+2.0*M_akg_c+2.0*M_nad_c --> 3.0*M_h_c+M_co2_c+M_urea_c+2.0*M_glu_L_c+2.0*M_nadh_c+M_succ_c
		0	100.0	;	# 102 M_arg_L_c+3.0*M_h2o_c+2.0*M_akg_c+M_nad_c+M_nadp_c --> 3.0*M_h_c+M_co2_c+M_urea_c+2.0*M_glu_L_c+M_nadh_c+M_nadph_c+M_succ_c
		0	100.0	;	# 103 M_arg_L_c+M_accoa_c+4.0*M_h2o_c+M_akg_c+M_nad_c --> M_coa_c+M_h_c+M_co2_c+2.0*M_nh4_c+2.0*M_glu_L_c+M_nadh_c+M_succ_c
		0	100.0	;	# 104 M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
		0	100.0	;	# 105 M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh4_c+M_mglx_c
		0	100.0	;	# 106 M_thr_L_c+M_coa_c+M_nad_c --> M_gly_L_c+M_accoa_c+M_nadh_c+M_h_c
		0	100.0	;	# 107 M_thr_L_c+M_pi_c+M_adp_c --> M_h_c+M_h2o_c+M_for_c+M_atp_c+M_prop_c
		0	100.0	;	# 108 M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> 2.0*M_h2o_c+M_co2_c+M_nadp_c+M_akg_c+M_ile_L_c
		0	100.0	;	# 109 M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh4_c+M_mglx_c
		0	100.0	;	# 110 M_mglx_c+M_nadp_c+M_h2o_c --> M_pyr_c+M_nadph_c+M_h_c
		0	100.0	;	# 111 M_pyr_c+M_nadph_c+M_h_c --> M_mglx_c+M_nadp_c+M_h2o_c
		0	100.0	;	# 112 M_ser_L_c --> M_nh4_c+M_pyr_c
		0	100.0	;	# 113 M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> 2.0*M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
		0	100.0	;	# 114 M_trp_L_c+M_h2o_c --> M_indole_c+M_nh4_c+M_pyr_c
		0	100.0	;	# 115 M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh4_c+M_pyr_c
		0	100.0	;	# 116 M_ala_L_c+M_h2o_c+M_q8_c --> M_q8h2_c+M_nh4_c+M_pyr_c
		0	100.0	;	# 117 M_lys_L_c --> M_co2_c+M_cadav_c
		0	100.0	;	# 118 M_gln_L_c+M_h2o_c --> M_nh4_c+M_glu_L_c
		0	100.0	;	# 119 M_glu_L_c+M_h_c --> M_co2_c+M_gaba_c
		0	100.0	;	# 120 M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+2.0*M_h_c+M_nadh_c
		0	100.0	;	# 121 M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+2.0*M_h_c+M_nadph_c
		0	100.0	;	# 122 M_asn_L_c+M_h2o_c+M_adp_c+M_pi_c --> M_nh4_c+M_asp_L_c+M_atp_c
		0	100.0	;	# 123 M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c
		0	100.0	;	# 124 M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c
		0	100.0	;	# 125 M_4adochor_c --> M_4abz_c+M_pyr_c
		0	100.0	;	# 126 M_4abz_c --> M_ppi_c+M_78dhf_c
		0	100.0	;	# 127 M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c
		0	100.0	;	# 128 M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c
		0	100.0	;	# 129 M_thf_c --> M_10fthf_c
		0	100.0	;	# 130 M_h2o_c+M_methf_c --> M_10fthf_c
		0	100.0	;	# 131 M_10fthf_c --> M_h2o_c+M_methf_c
		0	100.0	;	# 132 M_mlthf_c+M_nadp_c --> M_h_c+M_methf_c+M_nadph_c
		0	100.0	;	# 133 M_h_c+M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c
		0	100.0	;	# 134 M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c
		0	100.0	;	# 135 M_r5p_c+M_atp_c --> M_prpp_c+M_adp_c
		0	100.0	;	# 136 2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_atp_c+M_glu_L_c+M_pi_c+M_clasp_c
		0	100.0	;	# 137 M_clasp_c+M_q8_c+M_h2o_c --> M_or_c+M_q8h2_c
		0	100.0	;	# 138 M_prpp_c+M_or_c --> M_omp_c+M_ppi_c
		0	100.0	;	# 139 M_omp_c --> M_ump_c+M_co2_c
		0	100.0	;	# 140 M_utp_c+M_atp_c+M_nh4_c --> M_ctp_c+M_adp_c+M_pi_c
		0	100.0	;	# 141 M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c+3.0*M_h_c
		0	100.0	;	# 142 M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c
		0	100.0	;	# 143 M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c
		0	100.0	;	# 144 M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c
		0	100.0	;	# 145 M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c
		0	100.0	;	# 146 M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c
		0	100.0	;	# 147 M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c
		0	100.0	;	# 148 M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c
		0	100.0	;	# 149 M_saicar_c --> M_fum_c+M_aicar_c
		0	100.0	;	# 150 M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c
		0	100.0	;	# 151 M_faicar_c --> M_imp_c+M_h2o_c
		0	100.0	;	# 152 M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c
		0	100.0	;	# 153 M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c
		0	100.0	;	# 154 M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c
		0	100.0	;	# 155 M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c+M_h_c
		0	100.0	;	# 156 M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c+M_h_c
		0	100.0	;	# 157 M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c+M_h_c
		0	100.0	;	# 158 M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c+M_h_c
		0	100.0	;	# 159 M_atp_c+M_h2o_c --> M_adp_c+M_pi_c+M_h_c
		0	100.0	;	# 160 M_utp_c+M_h2o_c --> M_udp_c+M_pi_c+M_h_c
		0	100.0	;	# 161 M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c+M_h_c
		0	100.0	;	# 162 M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c+M_h_c
		0	100.0	;	# 163 M_udp_c+M_atp_c --> M_utp_c+M_adp_c
		0	100.0	;	# 164 M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c
		0	100.0	;	# 165 M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c
		0	100.0	;	# 166 M_atp_c+M_ump_c --> M_adp_c+M_udp_c
		0	100.0	;	# 167 M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c
		0	100.0	;	# 168 M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c
		0	100.0	;	# 169 M_amp_c+M_atp_c --> 2.0*M_adp_c
		0	100.0	;	# 170 M_amp_c+M_utp_c --> 2.0*M_udp_c
		0	100.0	;	# 171 M_amp_c+M_ctp_c --> 2.0*M_cdp_c
		0	100.0	;	# 172 M_amp_c+M_gtp_c --> 2.0*M_gdp_c
		0	100.0	;	# 173 GENE_CAT+RNAP --> OPEN_GENE_CAT
		0	100.0	;	# 174 OPEN_GENE_CAT+151.0*M_gtp_c+144.0*M_ctp_c+189.0*M_utp_c+176.0*M_atp_c --> mRNA_CAT+GENE_CAT+RNAP+1320.0*M_pi_c
		0	100.0	;	# 175 mRNA_CAT --> 151.0*M_gmp_c+144.0*M_cmp_c+189.0*M_ump_c+176.0*M_amp_c
		0	100.0	;	# 176 mRNA_CAT+RIBOSOME --> RIBOSOME_START_CAT
		0	100.0	;	# 177 RIBOSOME_START_CAT+438.0*M_gtp_c+15.0*M_ala_L_c_tRNA+5.0*M_arg_L_c_tRNA+10.0*M_asn_L_c_tRNA+12.0*M_asp_L_c_tRNA+5.0*M_cys_L_c_tRNA+12.0*M_glu_L_c_tRNA+13.0*M_gln_L_c_tRNA+10.0*M_gly_L_c_tRNA+12.0*M_his_L_c_tRNA+9.0*M_ile_L_c_tRNA+13.0*M_leu_L_c_tRNA+12.0*M_lys_L_c_tRNA+9.0*M_met_L_c_tRNA+20.0*M_phe_L_c_tRNA+7.0*M_pro_L_c_tRNA+10.0*M_ser_L_c_tRNA+13.0*M_thr_L_c_tRNA+5.0*M_trp_L_c_tRNA+11.0*M_tyr_L_c_tRNA+16.0*M_val_L_c_tRNA --> RIBOSOME+mRNA_CAT+PROTEIN_CAT+438.0*M_gdp_c+438.0*M_pi_c+219.0*tRNA
		0	100.0	;	# 178 15.0*M_ala_L_c+15.0*M_atp_c+15.0*tRNA --> 15.0*M_ala_L_c_tRNA+15.0*M_amp_c+30.0*M_pi_c
		0	100.0	;	# 179 5.0*M_arg_L_c+5.0*M_atp_c+5.0*tRNA --> 5.0*M_arg_L_c_tRNA+5.0*M_amp_c+10.0*M_pi_c
		0	100.0	;	# 180 10.0*M_asn_L_c+10.0*M_atp_c+10.0*tRNA --> 10.0*M_asn_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
		0	100.0	;	# 181 12.0*M_asp_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_asp_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0	100.0	;	# 182 5.0*M_cys_L_c+5.0*M_atp_c+5.0*tRNA --> 5.0*M_cys_L_c_tRNA+5.0*M_amp_c+10.0*M_pi_c
		0	100.0	;	# 183 12.0*M_glu_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_glu_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0	100.0	;	# 184 13.0*M_gln_L_c+13.0*M_atp_c+13.0*tRNA --> 13.0*M_gln_L_c_tRNA+13.0*M_amp_c+26.0*M_pi_c
		0	100.0	;	# 185 10.0*M_gly_L_c+10.0*M_atp_c+10.0*tRNA --> 10.0*M_gly_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
		0	100.0	;	# 186 12.0*M_his_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_his_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0	100.0	;	# 187 9.0*M_ile_L_c+9.0*M_atp_c+9.0*tRNA --> 9.0*M_ile_L_c_tRNA+9.0*M_amp_c+18.0*M_pi_c
		0	100.0	;	# 188 13.0*M_leu_L_c+13.0*M_atp_c+13.0*tRNA --> 13.0*M_leu_L_c_tRNA+13.0*M_amp_c+26.0*M_pi_c
		0	100.0	;	# 189 12.0*M_lys_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_lys_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0	100.0	;	# 190 9.0*M_met_L_c+9.0*M_atp_c+9.0*tRNA --> 9.0*M_met_L_c_tRNA+9.0*M_amp_c+18.0*M_pi_c
		0	100.0	;	# 191 20.0*M_phe_L_c+20.0*M_atp_c+20.0*tRNA --> 20.0*M_phe_L_c_tRNA+20.0*M_amp_c+40.0*M_pi_c
		0	100.0	;	# 192 7.0*M_pro_L_c+7.0*M_atp_c+7.0*tRNA --> 7.0*M_pro_L_c_tRNA+7.0*M_amp_c+14.0*M_pi_c
		0	100.0	;	# 193 10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA --> 10.0*M_ser_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
		0	100.0	;	# 194 13.0*M_thr_L_c+13.0*M_atp_c+13.0*tRNA --> 13.0*M_thr_L_c_tRNA+13.0*M_amp_c+26.0*M_pi_c
		0	100.0	;	# 195 5.0*M_trp_L_c+5.0*M_atp_c+5.0*tRNA --> 5.0*M_trp_L_c_tRNA+5.0*M_amp_c+10.0*M_pi_c
		0	100.0	;	# 196 11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA --> 11.0*M_tyr_L_c_tRNA+11.0*M_amp_c+22.0*M_pi_c
		0	100.0	;	# 197 16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA --> 16.0*M_val_L_c_tRNA+16.0*M_amp_c+32.0*M_pi_c
		0	100.0	;	# 198 [] --> tRNA
		0	100.0	;	# 199 tRNA --> []
		0	100.0	;	# 200 PROTEIN_CAT --> []
		0	100.0	;	# 201 [] --> M_o2_c
		0	100.0	;	# 202 M_co2_c --> []
		0	100.0	;	# 203 M_h_c --> []
		0	100.0	;	# 204 [] --> M_h_c
		0	100.0	;	# 205 [] --> M_h2s_c
		0	100.0	;	# 206 M_h2s_c --> []
		0	100.0	;	# 207 [] --> M_h2o_c
		0	100.0	;	# 208 M_h2o_c --> []
		0	100.0	;	# 209 [] --> M_pi_c
		0	100.0	;	# 210 M_pi_c --> []
		0	100.0	;	# 211 [] --> M_nh4_c
		0	100.0	;	# 212 M_nh4_c --> []
		0	100.0	;	# 213 [] --> M_glc_D_c
		0	100.0	;	# 214 M_glc_D_c --> []
		0	100.0	;	# 215 [] --> M_pyr_c
		0	0.0;	# 216 M_pyr_c --> []
		0	0.0;	# 217 [] --> M_ac_c
		0	0.0;	# 218 M_ac_c --> []
		0	0.0;	# 219 [] --> M_lac_D_c
		0	0.0;	# 220 M_lac_D_c --> []
		0	0.0;	# 221 [] --> M_hco3_c
		0	0.0;	# 222 M_hco3_c --> []
		0	0.0;	# 223 [] --> M_ala_L_c
		0	0.0;	# 224 M_ala_L_c --> []
		0	0.0;	# 225 [] --> M_arg_L_c
		0	0.0;	# 226 M_arg_L_c --> []
		0	0.0;	# 227 [] --> M_asn_L_c
		0	0.0;	# 228 M_asn_L_c --> []
		0	0.0;	# 229 [] --> M_asp_L_c
		0	0.0;	# 230 M_asp_L_c --> []
		0	0.0;	# 231 [] --> M_cys_L_c
		0	0.0;	# 232 M_cys_L_c --> []
		0	0.0;	# 233 [] --> M_glu_L_c
		0	0.0;	# 234 M_glu_L_c --> []
		0	0.0;	# 235 [] --> M_gln_L_c
		0	0.0;	# 236 M_gln_L_c --> []
		0	0.0;	# 237 [] --> M_gly_L_c
		0	0.0;	# 238 M_gly_L_c --> []
		0	0.0;	# 239 [] --> M_his_L_c
		0	0.0;	# 240 M_his_L_c --> []
		0	0.0;	# 241 [] --> M_ile_L_c
		0	0.0;	# 242 M_ile_L_c --> []
		0	0.0;	# 243 [] --> M_leu_L_c
		0	0.0;	# 244 M_leu_L_c --> []
		0	0.0;	# 245 [] --> M_lys_L_c
		0	0.0;	# 246 M_lys_L_c --> []
		0	0.0;	# 247 [] --> M_met_L_c
		0	0.0;	# 248 M_met_L_c --> []
		0	0.0;	# 249 [] --> M_phe_L_c
		0	0.0;	# 250 M_phe_L_c --> []
		0	0.0;	# 251 [] --> M_pro_L_c
		0	0.0;	# 252 M_pro_L_c --> []
		0	0.0;	# 253 [] --> M_ser_L_c
		0	0.0;	# 254 M_ser_L_c --> []
		0	0.0;	# 255 [] --> M_thr_L_c
		0	0.0;	# 256 M_thr_L_c --> []
		0	0.0;	# 257 [] --> M_trp_L_c
		0	0.0;	# 258 M_trp_L_c --> []
		0	0.0;	# 259 [] --> M_tyr_L_c
		0	0.0;	# 260 M_tyr_L_c --> []
		0	0.0;	# 261 [] --> M_val_L_c
		0	0.0;	# 262 M_val_L_c --> []
		0	100.0	;	# 263 M_atp_c --> []
	];

	# Setup default species bounds array -
	# Setup the objective coefficient array -
	objective_coefficient_array = [
		0.0	;	# 1 R_glk_atp::M_atp_c+M_glc_D_c --> M_adp_c+M_g6p_c+M_h_c
		0.0	;	# 2 R_pgi::M_g6p_c --> M_f6p_c
		0.0	;	# 3 R_pgi_reverse::M_f6p_c --> M_g6p_c
		0.0	;	# 4 R_pfk::M_atp_c+M_f6p_c --> M_adp_c+M_fdp_c+M_h_c
		0.0	;	# 5 R_fdp::M_fdp_c+M_h2o_c --> M_f6p_c+M_pi_c
		0.0	;	# 6 R_fbaA::M_fdp_c --> M_dhap_c+M_g3p_c
		0.0	;	# 7 R_fbaA_reverse::M_dhap_c+M_g3p_c --> M_fdp_c
		0.0	;	# 8 R_tpiA::M_dhap_c --> M_g3p_c
		0.0	;	# 9 R_tpiA_reverse::M_g3p_c --> M_dhap_c
		0.0	;	# 10 R_gapA::M_g3p_c+M_nad_c+M_pi_c --> M_13dpg_c+M_h_c+M_nadh_c
		0.0	;	# 11 R_gapA_reverse::M_13dpg_c+M_h_c+M_nadh_c --> M_g3p_c+M_nad_c+M_pi_c
		0.0	;	# 12 R_pgk::M_13dpg_c+M_adp_c --> M_3pg_c+M_atp_c
		0.0	;	# 13 R_pgk_reverse::M_3pg_c+M_atp_c --> M_13dpg_c+M_adp_c
		0.0	;	# 14 R_gpm::M_3pg_c --> M_2pg_c
		0.0	;	# 15 R_gpm_reverse::M_2pg_c --> M_3pg_c
		0.0	;	# 16 R_eno::M_2pg_c --> M_h2o_c+M_pep_c
		0.0	;	# 17 R_eno_reverse::M_h2o_c+M_pep_c --> M_2pg_c
		0.0	;	# 18 R_pyk::M_adp_c+M_h_c+M_pep_c --> M_atp_c+M_pyr_c
		0.0	;	# 19 R_pck::M_atp_c+M_oaa_c --> M_adp_c+M_co2_c+M_pep_c
		0.0	;	# 20 R_ppc::M_co2_c+M_h2o_c+M_pep_c --> M_h_c+M_oaa_c+M_pi_c
		0.0	;	# 21 R_pdh::M_coa_c+M_nad_c+M_pyr_c --> M_accoa_c+M_co2_c+M_nadh_c
		0.0	;	# 22 R_pps::M_atp_c+M_h2o_c+M_pyr_c --> M_amp_c+2.0*M_h_c+M_pep_c+M_pi_c
		0.0	;	# 23 R_zwf::M_g6p_c+M_nadp_c --> M_6pgl_c+M_h_c+M_nadph_c
		0.0	;	# 24 R_zwf_reverse::M_6pgl_c+M_h_c+M_nadph_c --> M_g6p_c+M_nadp_c
		0.0	;	# 25 R_pgl::M_6pgl_c+M_h2o_c --> M_6pgc_c+M_h_c
		0.0	;	# 26 R_gnd::M_6pgc_c+M_nadp_c --> M_co2_c+M_nadph_c+M_ru5p_D_c
		0.0	;	# 27 R_rpe::M_ru5p_D_c --> M_xu5p_D_c
		0.0	;	# 28 R_rpe_reverse::M_xu5p_D_c --> M_ru5p_D_c
		0.0	;	# 29 R_rpi::M_r5p_c --> M_ru5p_D_c
		0.0	;	# 30 R_rpi_reverse::M_ru5p_D_c --> M_r5p_c
		0.0	;	# 31 R_talAB::M_g3p_c+M_s7p_c --> M_e4p_c+M_f6p_c
		0.0	;	# 32 R_talAB_reverse::M_e4p_c+M_f6p_c --> M_g3p_c+M_s7p_c
		0.0	;	# 33 R_tkt1::M_r5p_c+M_xu5p_D_c --> M_g3p_c+M_s7p_c
		0.0	;	# 34 R_tkt1_reverse::M_g3p_c+M_s7p_c --> M_r5p_c+M_xu5p_D_c
		0.0	;	# 35 R_tkt2::M_e4p_c+M_xu5p_D_c --> M_f6p_c+M_g3p_c
		0.0	;	# 36 R_tkt2_reverse::M_f6p_c+M_g3p_c --> M_e4p_c+M_xu5p_D_c
		0.0	;	# 37 R_edd::M_6pgc_c --> M_2ddg6p_c+M_h2o_c
		0.0	;	# 38 R_eda::M_2ddg6p_c --> M_g3p_c+M_pyr_c
		0.0	;	# 39 R_gltA::M_accoa_c+M_h2o_c+M_oaa_c --> M_cit_c+M_coa_c+M_h_c
		0.0	;	# 40 R_acn::M_cit_c --> M_icit_c
		0.0	;	# 41 R_acn_reverse::M_icit_c --> M_cit_c
		0.0	;	# 42 R_icd::M_icit_c+M_nadp_c --> M_akg_c+M_co2_c+M_nadph_c
		0.0	;	# 43 R_sucAB::M_akg_c+M_coa_c+M_nad_c --> M_co2_c+M_nadh_c+M_succoa_c
		0.0	;	# 44 R_sucCD::M_adp_c+M_pi_c+M_succoa_c --> M_atp_c+M_coa_c+M_succ_c
		0.0	;	# 45 R_sdh::M_q8_c+M_succ_c --> M_fum_c+M_q8h2_c
		0.0	;	# 46 R_frd::M_fum_c+M_mql8_c --> M_mqn8_c+M_succ_c
		0.0	;	# 47 R_fum::M_fum_c+M_h2o_c --> M_mal_L_c
		0.0	;	# 48 R_mdh::M_mal_L_c+M_nad_c --> M_oaa_c+M_h_c+M_nadh_c
		0.0	;	# 49 R_cyd::2.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+2.0*M_h_e
		0.0	;	# 50 R_cyo::4.0*M_h_c+0.5*M_o2_c+M_q8h2_c --> M_h2o_c+M_q8_c+4.0*M_h_e
		0.0	;	# 51 R_app::2.0*M_h_c+M_mql8_c+0.5*M_o2_c --> M_h2o_c+M_mqn8_c+2.0*M_h_e
		0.0	;	# 52 R_atp::M_adp_c+M_pi_c+4.0*M_h_e --> M_atp_c+3.0*M_h_c+M_h2o_c
		0.0	;	# 53 R_nuo::3.0*M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c+2.0*M_h_e
		0.0	;	# 54 R_pnt1::M_nad_c+M_nadph_c --> M_nadh_c+M_nadp_c
		0.0	;	# 55 R_pnt2::M_nadh_c+M_nadp_c+2.0*M_h_e --> 2.0*M_h_c+M_nad_c+M_nadph_c
		0.0	;	# 56 R_ndh1::M_h_c+M_nadh_c+M_q8_c --> M_nad_c+M_q8h2_c
		0.0	;	# 57 R_ndh2::M_h_c+M_mqn8_c+M_nadh_c --> M_mql8_c+M_nad_c
		0.0	;	# 58 R_ppk::M_atp_c+M_pi_c --> M_adp_c+M_ppi_c
		0.0	;	# 59 R_ppa::M_ppi_c+M_h2o_c --> 2.0*M_pi_c+M_h_c
		0.0	;	# 60 R_aceA::M_icit_c --> M_glx_c+M_succ_c
		0.0	;	# 61 R_aceB::M_accoa_c+M_glx_c+M_h2o_c --> M_coa_c+M_h_c+M_mal_L_c
		0.0	;	# 62 R_maeA::M_mal_L_c+M_nad_c --> M_co2_c+M_nadh_c+M_pyr_c
		0.0	;	# 63 R_maeB::M_mal_L_c+M_nadp_c --> M_co2_c+M_nadph_c+M_pyr_c
		0.0	;	# 64 R_pta::M_accoa_c+M_pi_c --> M_actp_c+M_coa_c
		0.0	;	# 65 R_pta_reverse::M_actp_c+M_coa_c --> M_accoa_c+M_pi_c
		0.0	;	# 66 R_ackA::M_actp_c+M_adp_c --> M_ac_c+M_atp_c
		0.0	;	# 67 R_ackA_reverse::M_ac_c+M_atp_c --> M_actp_c+M_adp_c
		0.0	;	# 68 R_acs::M_ac_c+M_atp_c+M_coa_c --> M_accoa_c+M_amp_c+M_ppi_c
		0.0	;	# 69 R_adhE::M_accoa_c+2.0*M_h_c+2.0*M_nadh_c --> M_coa_c+M_etoh_c+2.0*M_nad_c
		0.0	;	# 70 R_adhE_reverse::M_coa_c+M_etoh_c+2.0*M_nad_c --> M_accoa_c+2.0*M_h_c+2.0*M_nadh_c
		0.0	;	# 71 R_ldh_f::M_pyr_c+M_nadh_c+M_h_c --> M_lac_D_c+M_nad_c
		0.0	;	# 72 R_ldh_r::M_lac_D_c+M_nad_c --> M_pyr_c+M_nadh_c+M_h_c
		0.0	;	# 73 R_pflAB::M_coa_c+M_pyr_c --> M_accoa_c+M_for_c
		0.0	;	# 74 R_alaAC::M_pyr_c+M_glu_L_c --> M_ala_L_c+M_akg_c
		0.0	;	# 75 R_arg::M_glu_L_c+M_accoa_c+4.0*M_atp_c+M_nadph_c+3.0*M_h2o_c+M_gln_L_c+M_asp_L_c+M_co2_c --> M_arg_L_c+M_coa_c+5.0*M_h_c+3.0*M_adp_c+3.0*M_pi_c+M_nadp_c+M_akg_c+M_ac_c+M_amp_c+M_ppi_c+M_fum_c
		0.0	;	# 76 R_aspA::M_fum_c+M_nh4_c --> M_asp_L_c
		0.0	;	# 77 R_aspA_reverse::M_asp_L_c --> M_fum_c+M_nh4_c
		0.0	;	# 78 R_aspC::M_glu_L_c+M_oaa_c --> M_asp_L_c+M_akg_c
		0.0	;	# 79 R_asnB::M_asp_L_c+M_gln_L_c+M_h2o_c+M_atp_c --> M_asn_L_c+M_glu_L_c+M_h_c+M_ppi_c+M_amp_c
		0.0	;	# 80 R_asnA::M_asp_L_c+M_atp_c+M_nh4_c --> M_asn_L_c+M_h_c+M_ppi_c+M_amp_c
		0.0	;	# 81 R_cysEMK::M_ser_L_c+M_accoa_c+M_h2s_c --> M_cys_L_c+M_coa_c+M_h_c+M_ac_c
		0.0	;	# 82 R_gltBD::M_gln_L_c+M_akg_c+M_nadph_c+M_h_c --> 2.0*M_glu_L_c+M_nadp_c
		0.0	;	# 83 R_gdhA::M_akg_c+M_nadph_c+M_nh4_c+M_h_c --> M_glu_L_c+M_h2o_c+M_nadp_c
		0.0	;	# 84 R_glnA::M_glu_L_c+M_atp_c+M_nh4_c --> M_gln_L_c+M_h_c+M_adp_c+M_pi_c
		0.0	;	# 85 R_glyA::M_ser_L_c+M_thf_c --> M_gly_L_c+M_h2o_c+M_mlthf_c
		0.0	;	# 86 R_his::M_gln_L_c+M_r5p_c+3.0*M_atp_c+2.0*M_nad_c+3.0*M_h2o_c --> M_his_L_c+M_akg_c+M_aicar_c+2.0*M_adp_c+2.0*M_nadh_c+M_pi_c+2.0*M_ppi_c+6.0*M_h_c
		0.0	;	# 87 R_ile::M_thr_L_c+2.0*M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> M_ile_L_c+M_h2o_c+M_nh4_c+M_co2_c+M_nadp_c+M_akg_c
		0.0	;	# 88 R_leu::2.0*M_pyr_c+M_glu_L_c+M_h_c+M_nad_c+M_nadph_c+M_accoa_c --> M_leu_L_c+2.0*M_co2_c+M_nadp_c+M_coa_c+M_nadh_c+M_akg_c
		0.0	;	# 89 R_lys::M_asp_L_c+M_atp_c+2.0*M_nadph_c+2.0*M_h_c+M_pyr_c+M_succoa_c+M_glu_L_c --> M_lys_L_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_coa_c+M_akg_c+M_succ_c+M_co2_c
		0.0	;	# 90 R_met::M_asp_L_c+M_cys_L_c+M_succoa_c+M_atp_c+2.0*M_nadph_c+M_5mthf_c --> M_met_L_c+M_coa_c+M_succ_c+M_adp_c+M_pi_c+2.0*M_nadp_c+M_thf_c+M_nh4_c+M_pyr_c
		0.0	;	# 91 R_phe::M_chor_c+M_h_c+M_glu_L_c --> M_phe_L_c+M_co2_c+M_h2o_c+M_akg_c
		0.0	;	# 92 R_pro::M_glu_L_c+M_atp_c+2.0*M_h_c+2.0*M_nadph_c --> M_pro_L_c+M_adp_c+2.0*M_nadp_c+M_pi_c+M_h2o_c
		0.0	;	# 93 R_serABC::M_3pg_c+M_nad_c+M_glu_L_c+M_h2o_c --> M_ser_L_c+M_nadh_c+M_h_c+M_akg_c+M_pi_c
		0.0	;	# 94 R_thr::M_asp_L_c+2.0*M_atp_c+2.0*M_nadph_c+M_h_c+M_h2o_c --> M_thr_L_c+2.0*M_adp_c+2.0*M_pi_c+2.0*M_nadp_c
		0.0	;	# 95 R_trp::M_chor_c+M_gln_L_c+M_ser_L_c+M_r5p_c+2.0*M_atp_c --> M_trp_L_c+M_glu_L_c+M_pyr_c+M_ppi_c+2.0*M_h2o_c+M_co2_c+M_g3p_c+2.0*M_adp_c+M_h_c
		0.0	;	# 96 R_tyr::M_chor_c+M_glu_L_c+M_nad_c --> M_tyr_L_c+M_akg_c+M_nadh_c+M_co2_c
		0.0	;	# 97 R_val::2.0*M_pyr_c+2.0*M_h_c+M_nadph_c+M_glu_L_c --> M_val_L_c+M_co2_c+M_nadp_c+M_h2o_c+M_akg_c
		0.0	;	# 98 R_arg_deg1::M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+2.0*M_nad_c+M_akg_c --> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+2.0*M_nadh_c+M_succ_c
		0.0	;	# 99 R_arg_deg2::M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+M_nad_c+M_nadp_c+M_akg_c --> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+M_nadh_c+M_nadph_c+M_succ_c
		0.0	;	# 100 R_arg_deg3::M_arg_L_c+5.0*M_h2o_c+M_atp_c+M_o2_c+2.0*M_nadp_c+M_akg_c --> 4.0*M_h_c+M_co2_c+M_urea_c+M_glu_L_c+M_pi_c+M_adp_c+M_nh4_c+M_h2o2_c+2.0*M_nadph_c+M_succ_c
		0.0	;	# 101 R_arg_deg4::M_arg_L_c+3.0*M_h2o_c+2.0*M_akg_c+2.0*M_nad_c --> 3.0*M_h_c+M_co2_c+M_urea_c+2.0*M_glu_L_c+2.0*M_nadh_c+M_succ_c
		0.0	;	# 102 R_arg_deg5::M_arg_L_c+3.0*M_h2o_c+2.0*M_akg_c+M_nad_c+M_nadp_c --> 3.0*M_h_c+M_co2_c+M_urea_c+2.0*M_glu_L_c+M_nadh_c+M_nadph_c+M_succ_c
		0.0	;	# 103 R_arg_deg6::M_arg_L_c+M_accoa_c+4.0*M_h2o_c+M_akg_c+M_nad_c --> M_coa_c+M_h_c+M_co2_c+2.0*M_nh4_c+2.0*M_glu_L_c+M_nadh_c+M_succ_c
		0.0	;	# 104 R_thr_deg1::M_thr_L_c+M_nad_c+M_coa_c --> M_nadh_c+M_h_c+M_accoa_c+M_gly_L_c
		0.0	;	# 105 R_thr_deg2::M_thr_L_c+M_nad_c+M_o2_c+M_h2o_c --> M_nadh_c+M_co2_c+M_h2o2_c+M_nh4_c+M_mglx_c
		0.0	;	# 106 R_thr_deg3::M_thr_L_c+M_coa_c+M_nad_c --> M_gly_L_c+M_accoa_c+M_nadh_c+M_h_c
		0.0	;	# 107 R_thr_deg4::M_thr_L_c+M_pi_c+M_adp_c --> M_h_c+M_h2o_c+M_for_c+M_atp_c+M_prop_c
		0.0	;	# 108 R_thr_deg5::M_thr_L_c+M_h_c+M_pyr_c+M_nadph_c+M_glu_L_c --> 2.0*M_h2o_c+M_co2_c+M_nadp_c+M_akg_c+M_ile_L_c
		0.0	;	# 109 R_gly_deg::M_gly_L_c+M_accoa_c+M_h_c+M_o2_c+M_h2o_c --> M_coa_c+M_co2_c+M_h2o2_c+M_nh4_c+M_mglx_c
		0.0	;	# 110 R_mglx_deg::M_mglx_c+M_nadp_c+M_h2o_c --> M_pyr_c+M_nadph_c+M_h_c
		0.0	;	# 111 R_mglx_deg_reverse::M_pyr_c+M_nadph_c+M_h_c --> M_mglx_c+M_nadp_c+M_h2o_c
		0.0	;	# 112 R_ser_deg::M_ser_L_c --> M_nh4_c+M_pyr_c
		0.0	;	# 113 R_pro_deg::M_pro_L_c+M_q8_c+2.0*M_h2o_c+M_nad_c --> 2.0*M_h_c+M_q8h2_c+M_nadh_c+M_glu_L_c
		0.0	;	# 114 R_trp_deg::M_trp_L_c+M_h2o_c --> M_indole_c+M_nh4_c+M_pyr_c
		0.0	;	# 115 R_cys_deg::M_cys_L_c+M_h2o_c --> M_h2s_c+M_nh4_c+M_pyr_c
		0.0	;	# 116 R_ala_deg::M_ala_L_c+M_h2o_c+M_q8_c --> M_q8h2_c+M_nh4_c+M_pyr_c
		0.0	;	# 117 R_lys_deg::M_lys_L_c --> M_co2_c+M_cadav_c
		0.0	;	# 118 R_gln_deg::M_gln_L_c+M_h2o_c --> M_nh4_c+M_glu_L_c
		0.0	;	# 119 R_glu_deg::M_glu_L_c+M_h_c --> M_co2_c+M_gaba_c
		0.0	;	# 120 R_gaba_deg1::M_gaba_c+M_akg_c+M_h2o_c+M_nad_c --> M_succ_c+M_glu_L_c+2.0*M_h_c+M_nadh_c
		0.0	;	# 121 R_gaba_deg2::M_gaba_c+M_akg_c+M_h2o_c+M_nadp_c --> M_succ_c+M_glu_L_c+2.0*M_h_c+M_nadph_c
		0.0	;	# 122 R_asn_deg::M_asn_L_c+M_h2o_c+M_adp_c+M_pi_c --> M_nh4_c+M_asp_L_c+M_atp_c
		0.0	;	# 123 R_chor::M_e4p_c+2.0*M_pep_c+M_nadph_c+M_atp_c --> M_chor_c+M_nadp_c+M_adp_c+4.0*M_pi_c
		0.0	;	# 124 R_fol_1::M_chor_c+M_gln_L_c --> M_4adochor_c+M_glu_L_c
		0.0	;	# 125 R_fol_2::M_4adochor_c --> M_4abz_c+M_pyr_c
		0.0	;	# 126 R_fol_2::M_4abz_c --> M_ppi_c+M_78dhf_c
		0.0	;	# 127 R_fol_3::M_78dhf_c+M_atp_c+M_glu_L_c --> M_adp_c+M_pi_c+M_dhf_c
		0.0	;	# 128 R_fol_4::M_dhf_c+M_nadph_c+M_h_c --> M_thf_c+M_nadp_c
		0.0	;	# 129 R_fol_5::M_thf_c --> M_10fthf_c
		0.0	;	# 130 R_mthfc::M_h2o_c+M_methf_c --> M_10fthf_c
		0.0	;	# 131 R_mthfc_reverse::M_10fthf_c --> M_h2o_c+M_methf_c
		0.0	;	# 132 R_mthfd::M_mlthf_c+M_nadp_c --> M_h_c+M_methf_c+M_nadph_c
		0.0	;	# 133 R_mthfd_reverse::M_h_c+M_methf_c+M_nadph_c --> M_mlthf_c+M_nadp_c
		0.0	;	# 134 R_mthfr2::M_mlthf_c+M_h_c+M_nadh_c --> M_5mthf_c+M_nad_c
		0.0	;	# 135 R_prpp_syn::M_r5p_c+M_atp_c --> M_prpp_c+M_adp_c
		0.0	;	# 136 R_or_syn_1::2.0*M_atp_c+M_gln_L_c+M_hco3_c+M_h2o_c --> 2.0*M_atp_c+M_glu_L_c+M_pi_c+M_clasp_c
		0.0	;	# 137 R_or_syn_2::M_clasp_c+M_q8_c+M_h2o_c --> M_or_c+M_q8h2_c
		0.0	;	# 138 R_omp_syn::M_prpp_c+M_or_c --> M_omp_c+M_ppi_c
		0.0	;	# 139 R_ump_syn::M_omp_c --> M_ump_c+M_co2_c
		0.0	;	# 140 R_ctp_1::M_utp_c+M_atp_c+M_nh4_c --> M_ctp_c+M_adp_c+M_pi_c
		0.0	;	# 141 R_ctp_2::M_utp_c+M_gln_L_c+M_atp_c+M_h2o_c --> M_ctp_c+M_glu_L_c+M_adp_c+M_pi_c+3.0*M_h_c
		0.0	;	# 142 R_A_syn_1::M_gln_L_c+M_prpp_c+M_h2o_c --> M_5pbdra+M_ppi_c+M_glu_L_c
		0.0	;	# 143 R_A_syn_2::M_5pbdra+M_gly_L_c --> M_adp_c+M_pi_c+M_gar_c
		0.0	;	# 144 R_A_syn_3::M_10fthf_c+M_gar_c --> M_thf_c+M_fgar_c
		0.0	;	# 145 R_A_syn_4::M_atp_c+M_fgar_c+M_gln_L_c+M_h2o_c --> M_adp_c+M_pi_c+M_fgam_c+M_glu_L_c
		0.0	;	# 146 R_A_syn_5::M_atp_c+M_fgam_c --> M_adp_c+M_pi_c+M_air_c
		0.0	;	# 147 R_A_syn_6::M_atp_c+M_air_c+M_hco3_c --> M_adp_c+M_pi_c+M_cair_c
		0.0	;	# 148 R_A_syn_7::M_atp_c+M_cair_c+M_asp_L_c --> M_adp_c+M_pi_c+M_saicar_c
		0.0	;	# 149 R_A_syn_8::M_saicar_c --> M_fum_c+M_aicar_c
		0.0	;	# 150 R_A_syn_9::M_aicar_c+M_10fthf_c --> M_thf_c+M_faicar_c
		0.0	;	# 151 R_A_syn_10::M_faicar_c --> M_imp_c+M_h2o_c
		0.0	;	# 152 R_A_syn_12::M_imp_c+M_gtp_c+M_asp_L_c --> M_gdp_c+M_pi_c+M_fum_c+M_amp_c
		0.0	;	# 153 R_xmp_syn::M_imp_c+M_nad_c+M_h2o_c --> M_xmp_c+M_nadh_c+M_h_c
		0.0	;	# 154 R_gmp_syn::M_atp_c+M_xmp_c+M_gln_L_c+M_h2o_c --> M_amp_c+M_ppi_c+M_gmp_c+M_glu_L_c
		0.0	;	# 155 R_atp_amp::M_atp_c+M_h2o_c --> M_amp_c+M_ppi_c+M_h_c
		0.0	;	# 156 R_utp_ump::M_utp_c+M_h2o_c --> M_ump_c+M_ppi_c+M_h_c
		0.0	;	# 157 R_ctp_cmp::M_ctp_c+M_h2o_c --> M_cmp_c+M_ppi_c+M_h_c
		0.0	;	# 158 R_gtp_gmp::M_gtp_c+M_h2o_c --> M_gmp_c+M_ppi_c+M_h_c
		0.0	;	# 159 R_atp_adp::M_atp_c+M_h2o_c --> M_adp_c+M_pi_c+M_h_c
		0.0	;	# 160 R_utp_adp::M_utp_c+M_h2o_c --> M_udp_c+M_pi_c+M_h_c
		0.0	;	# 161 R_ctp_adp::M_ctp_c+M_h2o_c --> M_cdp_c+M_pi_c+M_h_c
		0.0	;	# 162 R_gtp_adp::M_gtp_c+M_h2o_c --> M_gdp_c+M_pi_c+M_h_c
		0.0	;	# 163 R_udp_utp::M_udp_c+M_atp_c --> M_utp_c+M_adp_c
		0.0	;	# 164 R_cdp_ctp::M_cdp_c+M_atp_c --> M_ctp_c+M_adp_c
		0.0	;	# 165 R_gdp_gtp::M_gdp_c+M_atp_c --> M_gtp_c+M_adp_c
		0.0	;	# 166 R_atp_ump::M_atp_c+M_ump_c --> M_adp_c+M_udp_c
		0.0	;	# 167 R_atp_cmp::M_atp_c+M_cmp_c --> M_adp_c+M_cdp_c
		0.0	;	# 168 R_atp_gmp::M_atp_c+M_gmp_c --> M_adp_c+M_gdp_c
		0.0	;	# 169 R_adk_atp::M_amp_c+M_atp_c --> 2.0*M_adp_c
		0.0	;	# 170 R_adk_utp::M_amp_c+M_utp_c --> 2.0*M_udp_c
		0.0	;	# 171 R_adk_ctp::M_amp_c+M_ctp_c --> 2.0*M_cdp_c
		0.0	;	# 172 R_adk_gtp::M_amp_c+M_gtp_c --> 2.0*M_gdp_c
		0.0	;	# 173 transcriptional_initiation_CAT::GENE_CAT+RNAP --> OPEN_GENE_CAT
		0.0	;	# 174 transcription_CAT::OPEN_GENE_CAT+151.0*M_gtp_c+144.0*M_ctp_c+189.0*M_utp_c+176.0*M_atp_c --> mRNA_CAT+GENE_CAT+RNAP+1320.0*M_pi_c
		0.0	;	# 175 mRNA_degradation_CAT::mRNA_CAT --> 151.0*M_gmp_c+144.0*M_cmp_c+189.0*M_ump_c+176.0*M_amp_c
		0.0	;	# 176 translation_initiation_CAT::mRNA_CAT+RIBOSOME --> RIBOSOME_START_CAT
		0.0	;	# 177 translation_CAT::RIBOSOME_START_CAT+438.0*M_gtp_c+15.0*M_ala_L_c_tRNA+5.0*M_arg_L_c_tRNA+10.0*M_asn_L_c_tRNA+12.0*M_asp_L_c_tRNA+5.0*M_cys_L_c_tRNA+12.0*M_glu_L_c_tRNA+13.0*M_gln_L_c_tRNA+10.0*M_gly_L_c_tRNA+12.0*M_his_L_c_tRNA+9.0*M_ile_L_c_tRNA+13.0*M_leu_L_c_tRNA+12.0*M_lys_L_c_tRNA+9.0*M_met_L_c_tRNA+20.0*M_phe_L_c_tRNA+7.0*M_pro_L_c_tRNA+10.0*M_ser_L_c_tRNA+13.0*M_thr_L_c_tRNA+5.0*M_trp_L_c_tRNA+11.0*M_tyr_L_c_tRNA+16.0*M_val_L_c_tRNA --> RIBOSOME+mRNA_CAT+PROTEIN_CAT+438.0*M_gdp_c+438.0*M_pi_c+219.0*tRNA
		0.0	;	# 178 tRNA_charging_M_ala_L_c_CAT::15.0*M_ala_L_c+15.0*M_atp_c+15.0*tRNA --> 15.0*M_ala_L_c_tRNA+15.0*M_amp_c+30.0*M_pi_c
		0.0	;	# 179 tRNA_charging_M_arg_L_c_CAT::5.0*M_arg_L_c+5.0*M_atp_c+5.0*tRNA --> 5.0*M_arg_L_c_tRNA+5.0*M_amp_c+10.0*M_pi_c
		0.0	;	# 180 tRNA_charging_M_asn_L_c_CAT::10.0*M_asn_L_c+10.0*M_atp_c+10.0*tRNA --> 10.0*M_asn_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
		0.0	;	# 181 tRNA_charging_M_asp_L_c_CAT::12.0*M_asp_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_asp_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0.0	;	# 182 tRNA_charging_M_cys_L_c_CAT::5.0*M_cys_L_c+5.0*M_atp_c+5.0*tRNA --> 5.0*M_cys_L_c_tRNA+5.0*M_amp_c+10.0*M_pi_c
		0.0	;	# 183 tRNA_charging_M_glu_L_c_CAT::12.0*M_glu_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_glu_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0.0	;	# 184 tRNA_charging_M_gln_L_c_CAT::13.0*M_gln_L_c+13.0*M_atp_c+13.0*tRNA --> 13.0*M_gln_L_c_tRNA+13.0*M_amp_c+26.0*M_pi_c
		0.0	;	# 185 tRNA_charging_M_gly_L_c_CAT::10.0*M_gly_L_c+10.0*M_atp_c+10.0*tRNA --> 10.0*M_gly_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
		0.0	;	# 186 tRNA_charging_M_his_L_c_CAT::12.0*M_his_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_his_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0.0	;	# 187 tRNA_charging_M_ile_L_c_CAT::9.0*M_ile_L_c+9.0*M_atp_c+9.0*tRNA --> 9.0*M_ile_L_c_tRNA+9.0*M_amp_c+18.0*M_pi_c
		0.0	;	# 188 tRNA_charging_M_leu_L_c_CAT::13.0*M_leu_L_c+13.0*M_atp_c+13.0*tRNA --> 13.0*M_leu_L_c_tRNA+13.0*M_amp_c+26.0*M_pi_c
		0.0	;	# 189 tRNA_charging_M_lys_L_c_CAT::12.0*M_lys_L_c+12.0*M_atp_c+12.0*tRNA --> 12.0*M_lys_L_c_tRNA+12.0*M_amp_c+24.0*M_pi_c
		0.0	;	# 190 tRNA_charging_M_met_L_c_CAT::9.0*M_met_L_c+9.0*M_atp_c+9.0*tRNA --> 9.0*M_met_L_c_tRNA+9.0*M_amp_c+18.0*M_pi_c
		0.0	;	# 191 tRNA_charging_M_phe_L_c_CAT::20.0*M_phe_L_c+20.0*M_atp_c+20.0*tRNA --> 20.0*M_phe_L_c_tRNA+20.0*M_amp_c+40.0*M_pi_c
		0.0	;	# 192 tRNA_charging_M_pro_L_c_CAT::7.0*M_pro_L_c+7.0*M_atp_c+7.0*tRNA --> 7.0*M_pro_L_c_tRNA+7.0*M_amp_c+14.0*M_pi_c
		0.0	;	# 193 tRNA_charging_M_ser_L_c_CAT::10.0*M_ser_L_c+10.0*M_atp_c+10.0*tRNA --> 10.0*M_ser_L_c_tRNA+10.0*M_amp_c+20.0*M_pi_c
		0.0	;	# 194 tRNA_charging_M_thr_L_c_CAT::13.0*M_thr_L_c+13.0*M_atp_c+13.0*tRNA --> 13.0*M_thr_L_c_tRNA+13.0*M_amp_c+26.0*M_pi_c
		0.0	;	# 195 tRNA_charging_M_trp_L_c_CAT::5.0*M_trp_L_c+5.0*M_atp_c+5.0*tRNA --> 5.0*M_trp_L_c_tRNA+5.0*M_amp_c+10.0*M_pi_c
		0.0	;	# 196 tRNA_charging_M_tyr_L_c_CAT::11.0*M_tyr_L_c+11.0*M_atp_c+11.0*tRNA --> 11.0*M_tyr_L_c_tRNA+11.0*M_amp_c+22.0*M_pi_c
		0.0	;	# 197 tRNA_charging_M_val_L_c_CAT::16.0*M_val_L_c+16.0*M_atp_c+16.0*tRNA --> 16.0*M_val_L_c_tRNA+16.0*M_amp_c+32.0*M_pi_c
		0.0	;	# 198 tNRA_exchange::[] --> tRNA
		0.0	;	# 199 tNRA_exchange_reverse::tRNA --> []
		0.0	;	# 200 PROTEIN_export_CAT::PROTEIN_CAT --> []
		0.0	;	# 201 M_o2_c_exchange::[] --> M_o2_c
		0.0	;	# 202 M_co2_c_exchange::M_co2_c --> []
		0.0	;	# 203 M_h_c_exchange::M_h_c --> []
		0.0	;	# 204 M_h_c_exchange_reverse::[] --> M_h_c
		0.0	;	# 205 M_h2s_c_exchange::[] --> M_h2s_c
		0.0	;	# 206 M_h2s_c_exchange_reverse::M_h2s_c --> []
		0.0	;	# 207 M_h2o_c_exchange::[] --> M_h2o_c
		0.0	;	# 208 M_h2o_c_exchange_reverse::M_h2o_c --> []
		0.0	;	# 209 M_pi_c_exchange::[] --> M_pi_c
		0.0	;	# 210 M_pi_c_exchange_reverse::M_pi_c --> []
		0.0	;	# 211 M_nh4_c_exchange::[] --> M_nh4_c
		0.0	;	# 212 M_nh4_c_exchange_reverse::M_nh4_c --> []
		0.0	;	# 213 M_glc_D_c_exchange::[] --> M_glc_D_c
		0.0	;	# 214 M_glc_D_c_exchange_reverse::M_glc_D_c --> []
		0.0	;	# 215 M_pyr_c_exchange::[] --> M_pyr_c
		0.0	;	# 216 M_pyr_c_exchange_reverse::M_pyr_c --> []
		0.0	;	# 217 M_ac_c_exchange::[] --> M_ac_c
		0.0	;	# 218 M_ac_c_exchange_reverse::M_ac_c --> []
		0.0	;	# 219 M_lac_D_c_exchange::[] --> M_lac_D_c
		0.0	;	# 220 M_lac_D_c_exchange_reverse::M_lac_D_c --> []
		0.0	;	# 221 M_hco3_c_exchange::[] --> M_hco3_c
		0.0	;	# 222 M_hco3_c_exchange_reverse::M_hco3_c --> []
		0.0	;	# 223 M_ala_L_c_exchange::[] --> M_ala_L_c
		0.0	;	# 224 M_ala_L_c_exchange_reverse::M_ala_L_c --> []
		0.0	;	# 225 M_arg_L_c_exchange::[] --> M_arg_L_c
		0.0	;	# 226 M_arg_L_c_exchange_reverse::M_arg_L_c --> []
		0.0	;	# 227 M_asn_L_c_exchange::[] --> M_asn_L_c
		0.0	;	# 228 M_asn_L_c_exchange_reverse::M_asn_L_c --> []
		0.0	;	# 229 M_asp_L_c_exchange::[] --> M_asp_L_c
		0.0	;	# 230 M_asp_L_c_exchange_reverse::M_asp_L_c --> []
		0.0	;	# 231 M_cys_L_c_exchange::[] --> M_cys_L_c
		0.0	;	# 232 M_cys_L_c_exchange_reverse::M_cys_L_c --> []
		0.0	;	# 233 M_glu_L_c_exchange::[] --> M_glu_L_c
		0.0	;	# 234 M_glu_L_c_exchange_reverse::M_glu_L_c --> []
		0.0	;	# 235 M_gln_L_c_exchange::[] --> M_gln_L_c
		0.0	;	# 236 M_gln_L_c_exchange_reverse::M_gln_L_c --> []
		0.0	;	# 237 M_gly_L_c_exchange::[] --> M_gly_L_c
		0.0	;	# 238 M_gly_L_c_exchange_reverse::M_gly_L_c --> []
		0.0	;	# 239 M_his_L_c_exchange::[] --> M_his_L_c
		0.0	;	# 240 M_his_L_c_exchange_reverse::M_his_L_c --> []
		0.0	;	# 241 M_ile_L_c_exchange::[] --> M_ile_L_c
		0.0	;	# 242 M_ile_L_c_exchange_reverse::M_ile_L_c --> []
		0.0	;	# 243 M_leu_L_c_exchange::[] --> M_leu_L_c
		0.0	;	# 244 M_leu_L_c_exchange_reverse::M_leu_L_c --> []
		0.0	;	# 245 M_lys_L_c_exchange::[] --> M_lys_L_c
		0.0	;	# 246 M_lys_L_c_exchange_reverse::M_lys_L_c --> []
		0.0	;	# 247 M_met_L_c_exchange::[] --> M_met_L_c
		0.0	;	# 248 M_met_L_c_exchange_reverse::M_met_L_c --> []
		0.0	;	# 249 M_phe_L_c_exchange::[] --> M_phe_L_c
		0.0	;	# 250 M_phe_L_c_exchange_reverse::M_phe_L_c --> []
		0.0	;	# 251 M_pro_L_c_exchange::[] --> M_pro_L_c
		0.0	;	# 252 M_pro_L_c_exchange_reverse::M_pro_L_c --> []
		0.0	;	# 253 M_ser_L_c_exchange::[] --> M_ser_L_c
		0.0	;	# 254 M_ser_L_c_exchange_reverse::M_ser_L_c --> []
		0.0	;	# 255 M_thr_L_c_exchange::[] --> M_thr_L_c
		0.0	;	# 256 M_thr_L_c_exchange_reverse::M_thr_L_c --> []
		0.0	;	# 257 M_trp_L_c_exchange::[] --> M_trp_L_c
		0.0	;	# 258 M_trp_L_c_exchange_reverse::M_trp_L_c --> []
		0.0	;	# 259 M_tyr_L_c_exchange::[] --> M_tyr_L_c
		0.0	;	# 260 M_tyr_L_c_exchange_reverse::M_tyr_L_c --> []
		0.0	;	# 261 M_val_L_c_exchange::[] --> M_val_L_c
		0.0	;	# 262 M_val_L_c_exchange_reverse::M_val_L_c --> []
		-1.0	;	# 263 M_atp_exchange::M_atp_c --> []
	];

	# Min/Max flag - default is minimum -
	is_minimum_flag = true

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["is_minimum_flag"] = is_minimum_flag

	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
