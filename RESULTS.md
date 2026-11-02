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
