"""beamfea: a from-scratch 1D Euler-Bernoulli beam finite element solver.

Pipeline: config -> elements -> assembly -> solver -> postprocess / modal
-> plotting. All internal quantities are SI (m, N, Pa, kg, s).
"""

__version__ = "1.0.0"
