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
plt.clf()
demesdraw.tubes(demography.to_demes(), ax=plt.gca(), seed=1, log_time=True)
plt.show()

ts = msprime.sim_ancestry(
        {"N1": 50, "N2": 50}, # sample sizes 
        demography=demography, # defined above
        recombination_rate=3.5e-9, # as defined in pdf 
        sequence_length=1_000, # set sequence length
        random_seed=1234)
print(ts)

# then add mutations
mts = msprime.sim_mutations(ts, rate=3.5e-9, random_seed=1234)
print(mts.tables.sites)

# Define the samples between which Fst will be calculated
pop_id = {p.metadata["name"]: p.id for p in mts.populations()}
sample_sets=[mts.samples(pop_id["N1"]), mts.samples(pop_id["N2"])]

print(mts.Fst(sample_sets))

## Hurrah- this works!  Now adding into a loop 

def repeat_simulations(mut, sample_sizes, length, reco, pop_size, num_simulations, seed=None):
    results = []
    for i in tqdm(range(num_simulations), desc="Running simulations"): 
        if seed is not None:
            np.random.seed(seed + i) 
        # Simulate 10 diploid samples under the coalescent with recombination on a 10kb region.
        ts = msprime.sim_ancestry(
            samples=sum(sample_sizes),
            recombination_rate=reco,
            sequence_length=length,
            population_size=pop_size,
            random_seed=np.random.randint(99999999))
        
        # we can add mutations
        mutated_ts = msprime.sim_mutations(ts, rate=mut, random_seed=np.random.randint(99999999))

        diversity = mutated_ts.diversity()
        tajimas_d = mutated_ts.Tajimas_D()
        allele_frequency_spectrum = mutated_ts.allele_frequency_spectrum(polarised=True)
        results.append((mutated_ts, None, diversity, tajimas_d, allele_frequency_spectrum))
    return results


