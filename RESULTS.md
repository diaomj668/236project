# Simulation Results

The study uses two complementary designs. Making the most active users indiscriminate both adds noise and, when total purchases are fixed, removes informative purchases; no single design can hold both fixed. Design 1 (matched total) holds total purchases fixed and crosses preference strength r ∈ {2, 3, 4} with the 12-condition activity × preference × density grid, giving 36 conditions with 100 independent datasets each (3,600 datasets). Design 2 (fixed information) holds the informative users' purchases fixed and adds preference-free heavy buyers (18 conditions, 1,800 datasets). The r=2 grid was run first after a pilot of 120 datasets; the r=3 and r=4 grids followed with identical code and fresh seeds (root seeds 23603, 23604, and 23607). All rates below are mean misclustering percentages.

## Design 1 (matched total): r=2

|Mean|Activity|Preferences|Count|Jaccard|Weighted|SVD|
|---|---|---|---|---|---|---|
|6|Low|Focused|59.77|59.94|59.76|59.75|
|6|Low|Broad active|60.54|60.83|60.56|60.54|
|6|Medium|Focused|53.17|53.92|58.81|58.84|
|6|Medium|Broad active|62.83|62.91|62.13|62.60|
|6|High|Focused|41.28|42.31|41.23|59.70|
|6|High|Broad active|62.87|62.64|62.87|63.21|
|12|Low|Focused|35.51|36.33|35.31|35.46|
|12|Low|Broad active|48.53|49.17|48.55|48.48|
|12|Medium|Focused|22.36|22.65|39.05|37.01|
|12|Medium|Broad active|62.60|62.84|60.23|60.28|
|12|High|Focused|13.04|13.15|48.60|47.75|
|12|High|Broad active|62.94|62.94|62.97|63.16|

## Design 1 (matched total): r=3

|Mean|Activity|Preferences|Count|Jaccard|Weighted|SVD|
|---|---|---|---|---|---|---|
|6|Low|Focused|20.24|22.48|20.24|20.24|
|6|Low|Broad active|35.86|39.20|35.86|35.86|
|6|Medium|Focused|14.18|14.53|23.66|21.81|
|6|Medium|Broad active|62.85|62.71|59.81|60.45|
|6|High|Focused|7.62|7.76|7.62|39.34|
|6|High|Broad active|63.05|63.02|63.05|63.54|
|12|Low|Focused|1.60|1.62|1.60|1.60|
|12|Low|Broad active|4.51|4.58|4.51|4.51|
|12|Medium|Focused|2.15|2.16|2.12|1.83|
|12|Medium|Broad active|62.46|62.33|34.77|38.00|
|12|High|Focused|0.85|0.86|4.00|1.58|
|12|High|Broad active|62.92|62.87|62.17|62.95|

## Design 1 (matched total): r=4

|Mean|Activity|Preferences|Count|Jaccard|Weighted|SVD|
|---|---|---|---|---|---|---|
|6|Low|Focused|4.72|5.01|4.72|4.72|
|6|Low|Broad active|10.51|11.03|10.51|10.51|
|6|Medium|Focused|4.78|4.90|5.38|4.28|
|6|Medium|Broad active|62.66|62.64|53.95|56.47|
|6|High|Focused|2.07|2.12|2.07|10.43|
|6|High|Broad active|62.70|62.91|62.70|63.24|
|12|Low|Focused|0.15|0.14|0.15|0.15|
|12|Low|Broad active|0.49|0.50|0.49|0.49|
|12|Medium|Focused|0.28|0.29|0.15|0.13|
|12|Medium|Broad active|61.05|61.21|9.55|10.80|
|12|High|Focused|0.06|0.06|0.31|0.10|
|12|High|Broad active|62.86|62.95|60.74|62.10|

Assignment coverage was 100% for all methods. Confidence intervals are in `outputs/main/summary.csv` (r=2) and `outputs/supplement/sensitivity_r3/` and `outputs/supplement/sensitivity_r4/` (`summary.csv`, `paired_differences.csv`).

## Key comparisons in Design 1

- With focused active users, weighting is costly at r=2: at mean 12 and high heterogeneity, Count error is 13.04% and Weighted 48.60% (paired difference 35.57 points, 95% interval 34.14 to 36.99). Stronger preferences make focused conditions easy, so this cost shrinks; at r=4, Count, Jaccard, and Weighted stay below 5.4% in every focused condition.
- With broad active users at mean 12 and medium heterogeneity, the paired Count-minus-Weighted difference grows from 2.36 points at r=2 to 27.69 points at r=3 (26.13 to 29.25) and 51.50 points at r=4 (50.84 to 52.15). Weighted reaches 34.77% at r=3 and 9.55% at r=4 while Count stays near 61–63%.
- High-heterogeneity broad-active conditions stay near the random-label reference at every r, because the remaining focused users buy only one or two products each.

These results apply to the specified sampling scheme with balanced product categories and equal within-category product popularity. They do not establish that any method is universally superior.

## Design 2 (fixed information)

Design 2 was added after the matched-total results at r=2 showed that added background and lost informative purchases were confounded.

Design 2 fixes 1,000 informative users, then adds uniform background users. The following rows are illustrative; all 18 conditions are reported in `outputs/supplement/separation/` and the report appendix.

|Informative setting|Extra users|Purchases per extra user|Count error|Weighted error|SVD error|
|---|---|---|---|---|---|
|m=12, r=2|0|—|35.53%|35.53%|35.53%|
|m=12, r=2|50|48|59.05%|36.02%|36.33%|
|m=6, r=3|0|—|20.55%|20.55%|20.55%|
|m=6, r=3|50|24|56.83%|21.63%|22.32%|

The respective paired advantages of Weighted are 23.03 points (21.39 to 24.66) and 35.20 points (33.93 to 36.48). Informative matrices are identical across exposure levels within each replicate, while background users and purchases are nested. Total users and purchases increase, unlike Design 1. Count and Weighted are exactly equivalent when all contributing users have the same purchase count; shared k-means seeds verify this control.

Twenty label permutations per dataset estimate a method-specific matching-adjusted reference. Mean reference errors span 62.86%–63.45% across conditions; this is not a universal 63% chance-error constant. No independent-edge SBM threshold or universal preference crossover is claimed.

The combined conclusion is conditional: weighting can discard useful evidence from focused active buyers, and it can protect recovery against preference-free active buyers when informative purchases are retained.

## SVD ablation

In focused, high-heterogeneity cells at mean 6, 80% of users buy one product. They create no product pairs, so Count and Weighted ignore them, but they still enter the bipartite SVD through the user and product degrees and the diagonal of M'M. `scripts/svd_ablation.py` regenerates the recorded datasets from their seeds, reproduces the saved SVD errors exactly, and reruns the unchanged SVD without single-purchase users (100 datasets per cell).

|Cell|Single-purchase users|Count|Weighted|SVD|SVD without single-purchase users|SVD difference (pp, 95% CI)|
|---|---|---|---|---|---|---|
|Mean 6, high, focused, r=2|800|41.28|41.23|59.70|41.27|18.43 (16.98 to 19.88)|
|Mean 12, high, focused, r=2|0|13.04|48.60|47.75|47.75|0.00 (0.00 to 0.00)|
|Mean 6, high, focused, r=3|800|7.62|7.62|39.34|7.62|31.72 (29.44 to 34.00)|

Removing single-purchase users makes SVD match Count in both mean-6 cells. The mean-12 cell has no single-purchase users and is unchanged; there SVD tracks Weighted because both downweight active users.

## Validation

- All unit tests passed, including checks of theoretical formulas against simulation.
- Main-experiment dataset counts, replicate counts, method pairing, purchase means, and both total-weight identities passed integrity checks.
- Regenerating a dataset from its recorded seed reproduced method metrics and theoretical values.
- Resuming the pilot did not duplicate completed records.
- Confidence intervals treat datasets as independent units. Paired differences are calculated by subtracting method errors within the same dataset before aggregation.
