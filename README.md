# pendolino

generate_chromosome.py
Arguments:
-n number of parts of chromosome (default = 1000)
-r required average distance of a part from the center (default = 5.1)
-m prefix of the names of generated files
-d delta (the larger, the less likely to make distance-incresing moves) (default = 3)

The script creates three files:
1) chainA.pdb - chromosome-chain before condensation
2) chainB.pdb - chromosome-chain after condensation
3) chain - python-pickled chromosome-chain data after condensation

The script prints some self-explanatory information during the process.


position_the_chromosomes.py
Arguments:
-n number of steps for decondensation (default = 50000000)

The script creates multiple files:
1) final_chain_i.pdb - decondensed i-th chromosome-chain
2) positioned_chain_i.pdb - positioned condensed i-th chromosome-chain
3) lamin_laminfin.pdb - lamin
4) chain_i_step_k.pdb - i-th chromosome-chain after k steps of decondesation (after every 10 million steps)
5) final_chains - python-pickled chromosomes-chain data after decondensation

The script prints the progress (current step [out of n steps]).
