# Scientific and API migration to 0.2.0

The four-state equations and their pool balance are retained. The solver now checks input shapes, tolerances, finite derivatives and complete integration. A terminal state is called a steady state only after checking the infinity norm of its derivative. Stability classification applies to a specified equilibrium; a successful trajectory does not demonstrate bistability.

The default illustrative parameters are not calibrated experimental values. A publishable bistability claim still requires parameter provenance, multiple equilibria, continuation, stability analysis and robustness to tolerances. Cross-solver agreement and total-pool balance are tested.
