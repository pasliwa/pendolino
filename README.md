# pendolino

generate_chromosome.py
Arguments:
* -n number of parts of chromosome (default = 1000)
* -r required average distance of a part from the center (default = 5.1)
* -m prefix of the names of generated files
* -d delta (the larger, the less likely to make distance-incresing moves) (default = 3)

The script creates three files:
* chainA.pdb - chromosome-chain before condensation
* chainB.pdb - chromosome-chain after condensation
* chain - python-pickled chromosome-chain data after condensation

The script prints some self-explanatory information during the process.


position_the_chromosomes.py
Arguments:
* -n number of steps for decondensation (default = 50000000)

The script creates multiple files:
* final_chain_i.pdb - decondensed i-th chromosome-chain
* positioned_chain_i.pdb - positioned condensed i-th chromosome-chain
* lamin_laminfin.pdb - lamin
* chain_i_step_k.pdb - i-th chromosome-chain after k steps of decondesation (after every 10 million steps)
* final_chains - python-pickled chromosomes-chain data after decondensation

The script prints the progress (current step [out of n steps]).
