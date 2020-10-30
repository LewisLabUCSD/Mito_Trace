""" Lineage tracign simulation.

Running the first version of the lineage tracing simulation.
Files are simple_v02.yaml and wrapper_v02.yaml.
10000 iterations, 10000 num_cells.
Since this number is large, we are NOT going to save the full simulation
matrices, just the results
"""

from src.simulations import ParameterSweep

default_params_f = 'parameters/simple_v02.yaml'
sweep_params_f = 'parameters/sweep_v02.yaml'
sweep = ParameterSweep(default_params_f, sweep_params_f)

sweep.run_sweep()  # Runs the simulation
#sweep.save() # Do not save
sweep.plot_before_after_all_clones()
sweep.plot_sensitivity_and_dropout()  # plots results