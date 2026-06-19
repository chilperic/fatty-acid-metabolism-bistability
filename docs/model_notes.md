# Model notes

The reusable implementation follows the four-state coarse-grained structure documented in the PhD thesis:

1. acetyl-CoA,
2. malonyl-CoA,
3. fatty-acid pool,
4. triglycerides.

The model uses constant influxes, proportional outfluxes, Michaelis-Menten terms, and non-competitive inhibitory interactions. The code is organized so that numerical simulation, steady-state approximation, and local stability diagnostics are separate operations.
