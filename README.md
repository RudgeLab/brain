# brain
Reservoir computing with coupled genetic oscillators

Models:
danino - model of Danino oscillator.
danino2d - model of Danino oscillator which includes spatial diffusion.

Functions:
compute_period calculates period of oscillations.
compute_phase calculates phase shift between two oscillating signals.
prof_cos - function to generate sine wave.
prof_pulse - function to generate square wave.
compute_phases - compute phases in 2D space(to be confirmed).

Functions that run models:
SHIL_coupling tests a range of SHIL signal concentrations.
SHIL - to be confirmed.

run_coupling.m and running_model.m are example code to use SHIL_coupling and danino functions respectively.

analyse_all_mnist, prindle4, run_sim and mnist.nat are for reservoir computing. (to be confirmed).
