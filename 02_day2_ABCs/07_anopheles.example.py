# looking to run anopheles ms prime scenario

import msprime
import demesdraw
from tqdm import tqdm
import numpy as np

from matplotlib import pyplot as plt

# add demography
demography = msprime.Demography()
demography.add_population(name="N1", initial_size=90_000) # pop1 30x larger than pop2
demography.add_population(name="N2", initial_size=3_000)

# add ancestral group 
demography.add_population(name="ANC", initial_size=7_000_000)
demography.add_population_split(time=6000, derived=["N1", "N2"], ancestral="ANC") # split 9k generations ago - as default for now

# Plot a schematic of the model
# plt.clf()
# demesdraw.tubes(demography.to_demes(), ax=plt.gca(), seed=1, log_time=True)
# plt.show()

ts = msprime.sim_ancestry(
        {"N1": 50, "N2": 50}, # sample sizes 
        demography=demography, # defined above
        recombination_rate=3.5e-9, # as defined in pdf 
        sequence_length=1_000, # set sequence length
        random_seed=1234)
# print(ts)

# then add mutations
mts = msprime.sim_mutations(ts, rate=3.5e-9, random_seed=1234)
# print(mts.tables.sites)

# Define the samples between which Fst will be calculated
pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
sample_sets=[mts.samples(pop_id["N1"]), mts.samples(pop_id["N2"])]

print(mts.Fst(sample_sets))

## Hurrah- this works!  Now adding into a loop 

def repeat_simulations(n1_start, n2_start, anc_start, split_t, seq_length, num_simulations, seed=None):
    results = []
    for i in tqdm(range(num_simulations), desc="Running simulations"): 
        if seed is not None:
            np.random.seed(seed + i) 

        # add demography
        demography = msprime.Demography()
        demography.add_population(name="N1", initial_size= n1_start) # pop1 30x larger than pop2
        demography.add_population(name="N2", initial_size= n2_start)
        
        # add ancestral group 
        demography.add_population(name="ANC", initial_size=anc_start)
        demography.add_population_split(time=split_t, derived=["N1", "N2"], ancestral="ANC") # split 9k generations ago - as default for now
        
        # Plot a schematic of the model
        # plt.clf()
        # demesdraw.tubes(demography.to_demes(), ax=plt.gca(), seed=1, log_time=True)
        # plt.show()
        
        ts = msprime.sim_ancestry(
                {"N1": 50, "N2": 50}, # sample sizes 
                demography=demography, # defined above
                recombination_rate=3.5e-9, # as defined in pdf 
                sequence_length=seq_length, # set sequence length
                random_seed=1234)
        # print(ts)
        
        # then add mutations
        mts = msprime.sim_mutations(ts, rate=3.5e-9, random_seed=1234)
        # print(mts.tables.sites)

        diversity = mts.diversity()
        tajimas_d = mts.Tajimas_D()
        allele_frequency_spectrum = mts.allele_frequency_spectrum(polarised=True)
        
        pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
        sample_sets=[mts.samples(pop_id["N1"]), mts.samples(pop_id["N2"])]
        Fst = mts.Fst(sample_sets)
        
        results.append((mts, None, diversity, tajimas_d, allele_frequency_spectrum, Fst))
    return results
  
# define variables:
n1_start = 90_000
n2_start = 3_000
anc_start = 7_000_000 
split_t = 6000
seq_length = 1_000
num_simulations = 1
seed = 1234

results = repeat_simulations(n1_start , n2_start, anc_start, split_t, seq_length, num_simulations, seed=seed)

diversities = [result[2] for result in results]
tajimas_ds = [result[3] for result in results]
allele_frequency_spectra = [result[4] for result in results]

diversity_values = []
tajima_values = []
allele_freq = []

times_to_split = [5000,6000]

# trying to make a loop 
for time in times_to_split: 
   
   # define variables:
    n1_start = 90_000
    n2_start = 3_000
    anc_start = 7_000_000 
    split_t = time
    seq_length = 1_000
    num_simulations = 1
    seed = 1234
   
    results = repeat_simulations(n1_start , n2_start, anc_start, split_t, seq_length, num_simulations, seed=seed)
    
    return(results)
  

# Output directory
output_directory = "."
# Output file name
output_file = os.path.join(output_directory, "mosquito-task2.csv")
    
with open(output_file, "w", newline="") as csvfile:
writer = csv.writer(csvfile, delimiter=",")

# Write header
writer.writerow(["N1", "N2", "T_split", "MigRate", "Diversity", "Tajima", "AlleleFreq", "Fst"])
  
# Write data to file or print data
writer.writerow([params["N1"], params["N2"], params["T_split"], params["mig"],Diversity, Tajima, AlleleFreq, Fst)
    
diversity_values.append([result[2] for result in results])
tajima_values.append([result[3] for result in results])
allele_freq.append([result[4] for result in results]) 
Fst.append([result[5] for result in results])
