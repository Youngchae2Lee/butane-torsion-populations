# Prospective targeted follow-up, 2026-10-15

Before executing the next two jobs, freeze the existing cosine-order 3 and 2 models from commit 3090f6b. Compute new constrained optimizations at 15 and 45 degrees with the same route and atom-order convention as the seven-point set. Use its 180-degree electronic energy as reference. Do not refit either model or replace the original held-out metrics.

These angles were selected after inspecting the large 30-degree error. They are new quantum-chemical reference values for a targeted local diagnostic; they are not a randomly selected or representative global test set. Report each frozen-model residual and the two-point summary separately from the original three-point evaluation. Check normal termination, four convergence criteria, actual torsion within 0.01 degree, connectivity and source hashes. Do not launch further angles automatically.
