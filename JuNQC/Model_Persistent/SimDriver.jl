# Script to estimate the acetate and celmass in W3110
# Varma A, Palsson BO (1994) Stoichiometric flux balance models quantitatively predict growth and metabolic by-product
# secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60: 3724-31.

# include -
include("Include.jl")

# setup the time-scale - Total simulation
time_start = 0.0
time_feed_start = 10.0
time_feed_stop = 225.0
time_stop = 250.0
time_step = 0.5

# setup the time-scale - Phase 1
time_start_phase_1 = time_start
time_stop_phase_1 = time_feed_start
time_array_1 = collect(time_start_phase_1:time_step:time_stop_phase_1)

# initialize the problem -
number_of_timesteps = length(time_array_1)
number_of_external_states = 5
state_array_1 = zeros(number_of_timesteps,number_of_external_states)

# set the ic -
state_array_1[1,1] = 0.0 # product
state_array_1[1,2] = 0.3 #cellmass (10^6 cell/mL)
state_array_1[1,3] = 0.5   # 3 glucose (mmol)
state_array_1[1,4] = 0.0 #product_mRNA
state_array_1[1,5] = 4.0 #Bioreactor Volume (L)
#state_array[1,2] = 0.5     # 2 acetate
#state_array[1,3] = 0.001   # 3 cellmass
#state_array[1,4] = 0.5     # 4 formate
#state_array[1,5] = 0.5     # 5 ethanol

# Fire up the max cellmass -
data_dictionary = maximize_product_data_dictionary(time_start,time_stop,time_step)

feed_flag = false

(T1,X1) = runSim(time_array_1,state_array_1,data_dictionary,feed_flag)

# Begin feed time
time_start_phase_2 = T1[end]+time_step
time_stop_phase_2 = time_feed_stop
time_array_2 = collect(time_start_phase_2:time_step:time_stop_phase_2)

number_of_timesteps = length(time_array_2)
state_array_2 = zeros(number_of_timesteps,number_of_external_states)
state_array_2[1,:] = X1[end,:]

feed_flag = true

(T2,X2) = runSim(time_array_2,state_array_2,data_dictionary,feed_flag)

# Stop feed
time_start_phase_3 = T2[end]+time_step
time_stop_phase_3 = time_stop
time_array_3 = collect(time_start_phase_3:time_step:time_stop_phase_3)

number_of_timesteps = length(time_array_3)
state_array_3 = zeros(number_of_timesteps,number_of_external_states)
state_array_3[1,:] = X2[end,:]

feed_flag = false

(T3,X3) = runSim(time_array_3,state_array_3,data_dictionary,feed_flag)

T = [T1 ; T2 ; T3];
X = [X1 ; X2 ; X3];

Titles = ["Product (mmol)","Cell Mass (10^6 Cell/mL)","Glucose (mmol)","Product mRNA (mmol)","Volume (L)"]
for ii in 1:length(Titles)
  figure(ii)
  plot(T,X[:,ii])
  xlabel("Time (hr)")
  ylabel(Titles[ii])
end
