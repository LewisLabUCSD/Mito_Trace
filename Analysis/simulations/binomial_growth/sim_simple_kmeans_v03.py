""" Lineage tracign simulation.

Running v03 version of the lineage tracing simulation.
Files are simple_v03.yaml and wrapper_v03.yaml.
1000 iterations, 10000 num_cells.

In this iteration, the initial population will not be limited to number
of cells but will be its own parameter. This will affect the growth.


Since this number is large, we are NOT going to save the full simulation
matrices, just the results
"""

from src.simulations.parametersweep import ParameterSweep

default_params_f = 'parameters/simple_v03.yaml'
sweep_params_f = 'parameters/sweep_v03.yaml'
sweep = ParameterSweep(default_params_f, sweep_params_f)

sweep.run_sweep()  # Runs the simulation
#sweep.save() # Do not save
sweep.plot_before_after_all_clones()
sweep.plot_sensitivity_and_dropout()  # plots results
sweep.plot_before_after_true_growth()
