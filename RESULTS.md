# Main Simulation Results

The main experiment completed 1,200 independent datasets (12 conditions with 100 replicates each). Each dataset was evaluated by four methods, giving 4,800 method-level records. The pilot used another 120 datasets with seeds disjoint from the main experiment. Parameter settings are recorded in `frozen_config.json`.

## Mean misclustering rates

|Mean purchases|Heterogeneity|Preference|Count|Jaccard|Weighted|Bipartite SVD|
|---|---|---|---|---|---|---|
|6|high|broad_active|0.6287|0.6264|0.6287|0.6321|
|6|high|focused|0.4128|0.4231|0.4123|0.5970|
|6|low|broad_active|0.6054|0.6083|0.6056|0.6054|
|6|low|focused|0.5977|0.5994|0.5976|0.5975|
|6|medium|broad_active|0.6283|0.6291|0.6213|0.6260|
|6|medium|focused|0.5317|0.5392|0.5881|0.5884|
|12|high|broad_active|0.6294|0.6294|0.6297|0.6316|
|12|high|focused|0.1304|0.1315|0.4860|0.4775|
|12|low|broad_active|0.4853|0.4917|0.4855|0.4848|
|12|low|focused|0.3551|0.3633|0.3531|0.3546|
|12|medium|broad_active|0.6260|0.6284|0.6023|0.6028|
|12|medium|focused|0.2236|0.2265|0.3905|0.3701|

Assignment coverage was 100% for all methods. Confidence intervals are available in `outputs/main/summary.csv`; paired method differences are in `outputs/main/paired_differences.csv`.

These results apply to the specified sampling scheme with balanced product categories and equal within-category product popularity. They do not establish that any method is universally superior. Recovery was poor under the broad-active preference setting; these conditions remain in the study as part of the results.

## Validation

- All five unit tests passed, including checks of theoretical formulas against simulation.
- Main-experiment dataset counts, replicate counts, method pairing, purchase means, and both total-weight identities passed integrity checks.
- Regenerating a dataset from its recorded seed reproduced method metrics and theoretical values.
- Resuming the pilot did not duplicate completed records.
- Confidence intervals treat datasets as independent units. Paired differences are calculated by subtracting method errors within the same dataset before aggregation.

## Post hoc exploratory supplementary results

The original r=2 main results above are unchanged. After inspecting them and user-supplied exploratory descriptions, we recorded an additional design before fresh verification runs. The supplementary analyses use the original dependency environment, separate root seeds, and 100 replicates per cell. They comprise 1,200 r=3 sensitivity datasets and 1,800 fixed-information separation datasets.

At r=3, mean 12, medium heterogeneity, and broad active users, Count error is 62.46%, Jaccard 62.33%, Weighted 34.77%, and SVD 38.00%. The paired Count-minus-Weighted difference is 27.69 percentage points (95% Monte Carlo interval 26.13 to 29.25). Some focused r=3 cells become easy; the full grid is retained.

A later r=4 grid (the proposal's upper pilot value, root seed 23607, 1,200 further datasets) strengthens this pattern. In the same mean-12, medium-heterogeneity, broad-active condition, Count error is 61.05%, Jaccard 61.21%, Weighted 9.55%, and SVD 10.80%; the paired Count-minus-Weighted difference is 51.50 points (50.84 to 52.15). Count, Jaccard, and Weighted errors are below 5.4% in every focused r=4 cell. High-heterogeneity broad-active cells stay near the permutation reference because the remaining focused users buy only one or two products each.

The separation experiment fixes 1,000 informative users, then adds uniform background users. The following rows are illustrative; all 18 conditions are reported in `outputs/supplement/separation/` and the report appendix.

|Informative setting|Extra users|Purchases per extra user|Count error|Weighted error|SVD error|
|---|---|---|---|---|---|
|m=12, r=2|0|—|35.53%|35.53%|35.53%|
|m=12, r=2|50|48|59.05%|36.02%|36.33%|
|m=6, r=3|0|—|20.55%|20.55%|20.55%|
|m=6, r=3|50|24|56.83%|21.63%|22.32%|

The respective paired advantages of Weighted are 23.03 points (21.39 to 24.66) and 35.20 points (33.93 to 36.48). Informative matrices are identical across exposure levels within each replicate, while background users and purchases are nested. Total users and purchases increase, unlike the matched-total main design. Count and Weighted are exactly equivalent when all contributing users have the same purchase count; shared supplementary k-means seeds verify this control.

Twenty label permutations per dataset estimate a method-specific matching-adjusted reference. Mean reference errors span 62.86%–63.45% across supplementary cells; this is not a universal 63% chance-error constant. No independent-edge SBM threshold or universal preference crossover is claimed.

The combined conclusion is conditional: weighting can discard useful evidence from focused active buyers, and it can protect recovery against the tested preference-free active buyers when informative purchases are retained. These supplementary settings were chosen after observing results and must remain labeled exploratory.


## Post hoc SVD ablation

In focused, high-heterogeneity cells at mean 6, 80% of users buy one product. They create no product pairs, so Count and Weighted ignore them, but they still enter the bipartite SVD through the user and product degrees and the diagonal of M'M. `scripts/svd_ablation.py` regenerates the recorded datasets from their seeds, reproduces the saved SVD errors exactly, and reruns the unchanged SVD without single-purchase users (100 datasets per cell).

|Cell|Single-purchase users|Count|Weighted|SVD|SVD without single-purchase users|SVD difference (pp, 95% CI)|
|---|---|---|---|---|---|---|
|Mean 6, high, focused, r=2|800|41.28|41.23|59.70|41.27|18.43 (16.98 to 19.88)|
|Mean 12, high, focused, r=2|0|13.04|48.60|47.75|47.75|0.00 (0.00 to 0.00)|
|Mean 6, high, focused, r=3|800|7.62|7.62|39.34|7.62|31.72 (29.44 to 34.00)|

Removing single-purchase users makes SVD match Count in both mean-6 cells. The mean-12 cell has no single-purchase users and is unchanged; there SVD tracks Weighted because both downweight active users. This ablation was added after the supplementary analyses and is exploratory.
