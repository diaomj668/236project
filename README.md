# STAT 236: User Activity and Product Community Recovery

This project implements the main simulation design in the proposal. This directory contains the source code, experiment records, results, and reproducibility instructions.

## Project materials

- [Full report (PDF)](deliverables/report.pdf)
- [Editable report (Word)](deliverables/report.docx)
- [Results summary](RESULTS.md)
- [Original proposal](236%20proposal.pdf), preserved without content changes
- [Interactive webpage source](website/dist/index.html) and [portable website download](deliverables/website.zip)
- [Code and result download](deliverables/reproducibility.zip)
- [Example raw purchase matrix](data/example/purchase_matrix.csv), [user information](data/example/users.csv), and [true product categories](data/example/products.csv)

The repository is public with the owner's authorization. The demonstration webpage is published through GitHub Pages at **https://diaomj668.github.io/236project/**. It includes interactive results and downloadable report, code, and experiment tables.

## Quick start

Clone the repository and create a Python 3.11 environment:

```bash
git clone https://github.com/diaomj668/236project.git
cd 236project
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements-lock.txt
python -m unittest discover -s tests -v
python scripts/validate_project.py
```

The completed pilot and main results are included. Rebuild summaries and figures without running new simulations:

```bash
python run_experiments.py --stage main --analyze-only
```

To reproduce the main experiment using its recorded seeds in a new directory:

```bash
python run_experiments.py --stage main --frozen-config frozen_config.json \
  --seed 23603 --workers 4 --output outputs/reproduced_main
```

`requirements-lock.txt` records the verified dependency versions. `requirements.txt` specifies compatible version ranges. Statistical results depend on the recorded seeds and software environment; runtimes depend on the machine and numerical libraries. The resume command requires matching manifest metadata, including the recorded Python version, so a new output directory is preferable on a different machine.

## Raw simulated data

`data/example/` contains a reconstructed recorded main-experiment dataset with 1,000 users, 300 products, 12,000 purchases, focused preferences, and high activity heterogeneity. It is replicate 0 of that condition, with dataset seed 2344398831. It is not an additional independent replicate.

`purchase_matrix.csv` has user IDs in its first column and 300 binary product columns. `products.csv` supplies true category labels. `users.csv` records activity, preferred category, and preference ratio. `purchases.npz` provides the same arrays compactly, and `metadata.json` identifies the experimental condition. All identifiers and labels start at zero.

Regenerate this example, or export any recorded dataset:

```bash
python scripts/export_dataset.py
python scripts/export_dataset.py --dataset-seed 2344398831 --output data/reconstructed
```

The full 1,200 purchase matrices are not stored individually; their seeds and all analysis results are included.

## Files

- `model.py`: data generation, three projections, spectral clustering, bipartite SVD, optimal label matching, and theoretical calculations.
- `run_experiments.py`: experiment grids, independent random streams, parallel execution, checkpointing, summaries, paired comparisons, and pilot checks.
- `figures.py`: recovery curves, theoretical versus simulated contrasts, user influence, paired differences, and an illustrative network.
- `tests/test_model.py`: tests for matched means, purchase counts, weighting identities, isolated products, label permutations, theoretical moments, and recovery of a known partition.
- `frozen_config.json`: the initial parameter settings retained after the pilot passed its checks.
- `outputs/pilot/` and `outputs/main/`: pilot and main experiment records, including results for each dataset and figures.
- `RESULTS.md`: a summary of the completed main experiment and validation checks.

## Design and implementation

The main design uses n=1,000 users, p=300 products, K=3 categories, mean purchase counts of 6 and 12, three activity distributions (`low`, `medium`, and `high`), two preference settings (`focused` and `broad_active`), and r=2. This gives 12 experimental conditions.

Under medium heterogeneity, 80% of users purchase mean/2 products and 20% purchase 3 times the mean. Under high heterogeneity, 80% purchase mean/6 products and 20% purchase 13 times the mean divided by 3. Means are exactly matched across activity distributions. Product categories have equal sizes. Each user's preferred category is drawn independently and uniformly, giving approximately balanced preference groups. Activity is assigned independently of preferred category. Ties in selecting the most active 20% are broken randomly. Under `broad_active`, these users have r=1.

Products are sampled without replacement using NumPy `choice(replace=False, p=...)`, with weights r for products in the preferred category and 1 for other products. This corresponds to successive probability-proportional-to-size sampling without replacement, rather than independent Bernoulli sampling or Fisher's noncentral hypergeometric sampling.

All three projections retain their weighted edges without thresholding and set their diagonals to zero. Spectral clustering uses the eigenvectors corresponding to the K largest algebraic eigenvalues of D^(-1/2) W D^(-1/2), normalizes nonzero embedding rows, and applies k-means with `n_init=20`. The bipartite baseline removes zero-degree rows and columns, computes the SVD of Du^(-1/2) B Dp^(-1/2), retains the leading K right singular vectors, normalizes their nonzero rows, and uses the same k-means settings. The SVD baseline retains the implicit effect of the Gram matrix diagonal.

Products with zero strength or zero embedding rows receive label -1. Optimal label matching uses assigned products, but the misclustering rate divides by the total number of products, so unassigned products count as errors. Adjusted Rand index (ARI) is calculated among assigned products, and assignment coverage is reported separately. If fewer than K products have valid embeddings, all products are conservatively left unassigned.

All four methods use the same B within each dataset. Data generation and each method's k-means use independent random substreams. The default root seeds for pilot, calibration, and main experiments are 23601, 23602, and 23603, respectively. Dataset and k-means seeds are recorded. The manifest and dataset seed allow reconstruction without storing every B matrix.

The theory module uses a one-dimensional dynamic program for the number of selected products in the preferred category. It computes category second moments for the actual sequential sampling scheme. Conditional exchangeability of other products then gives within-category and between-category pair inclusion probabilities. The reported theoretical values are **unconditional expectations** over independently and uniformly drawn preferred categories. Random variation in preference-group sizes contributes to dataset variability. Count and Weighted expectations are exact up to numerical precision. Jaccard uses an explicitly labeled ratio-of-expectations approximation, not an exact expression for E[J]. Tests check the uniform-sampling formula at r=1 and compare general-r theoretical values against simulation.

## Outputs

- `results.csv`: one row per dataset and method, including misclustering error, ARI, coverage, within-category and between-category edge means, their contrast, strength coefficient of variation, runtime, theoretical values, and random seeds. Edge statistics are empty for the direct SVD baseline, which has no explicit projected network.
- `diagnostics.csv`: realized purchase counts, preference-group sizes, isolated products, weighting identities, Count-Jaccard edge-weight Spearman correlations over all pairs and over nonzero pairs, and clustering ARI.
- `summary.csv`: means, standard deviations, Monte Carlo standard errors, 95% Student-t confidence intervals for means, and valid replicate counts. A single-replicate smoke run does not estimate intervals.
- `paired_differences.csv`: all six paired method comparisons within each condition. Differences are method_a minus method_b; negative error differences favor method_a.
- `checkpoint.json`: the number of conditions satisfying the difficulty criterion, coverage, and whether the estimated serial generation-plus-algorithm time for 1,200 datasets is below 24 hours. The estimate excludes plotting, file writing, and process startup.
- `manifest.json`: parameters, root seed, package versions, sampling scheme, and interval definitions.
- `datasets/*.json`: atomic checkpoints for completed datasets, used to resume execution and rebuild result tables.
- `figures/`: PNG figures for the report.

Confidence intervals describe uncertainty in mean performance across independent generated datasets. They are not 95% prediction intervals for individual dataset errors. No multiple-comparison adjustment is applied. Raw edge contrasts do not establish spectral recovery guarantees, and projected edges are dependent.

## Additional commands

```bash
# Preference-strength sensitivity pilot
python run_experiments.py --stage pilot --r 1.25 2 4 --output outputs/pilot_r --workers 4

# Bounded calibration grid: purchase means 3, 6, 12, 24 and r=1.25, 2, 4
python run_experiments.py --stage calibrate --workers 4

# Rebuild summaries and figures from existing main-experiment CSV files
python run_experiments.py --stage main --analyze-only

# Check the complete pipeline with 12 datasets
python run_experiments.py --stage smoke --workers 2
```

At mean=3, the high-heterogeneity calibration setting cannot use mean/6 as an integer purchase count. It instead assigns one purchase to 80% of users and eleven purchases to 20%, preserving mean=3 exactly. The main experiment follows the proposal's table. Calibration does not automatically select the final design: freeze settings based on difficulty and feasibility, rather than the largest method advantage. A reviewed manifest can be supplied through `--frozen-config`. The optional unequal-product-popularity extension is not implemented; the main design retains equal sampling weights within each category.

Existing output directories require `--resume`, with an exactly matching configuration. Use a new `--output` directory when changing parameters. Changing `--workers` does not change datasets or k-means seeds. Algorithm runtimes include projection construction or matrix decomposition and clustering, but exclude diagnostics, theory, and generation. Shared machine load affects timing comparisons.

## Report and demonstration webpage

The complete English report is available in `deliverables/report.docx`, `deliverables/report.pdf`, and editable Markdown source `deliverables/report.md`. The report contains the original main analysis, clearly labeled exploratory supplementary analyses, references, and appendices. It includes the exact theoretical recurrence, all condition-level mean errors, paired confidence intervals, diagnostics, and limitations.

The responsive demonstration webpage is in `website/dist/index.html`. It works locally without external dependencies and includes interactive preference and density controls. To preview it, run `python -m http.server 8765 --directory website/dist` and visit `http://localhost:8765`. Report and reproducibility downloads are included under `website/dist/assets/`.

`scripts/build_report_assets.py` rebuilds the additional plots and compact website dataset from recorded results. `scripts/build_report.py` builds the Word document and requires `python-docx` in the document-authoring environment; it is not needed to run the statistical experiments. PDF rendering uses the bundled document renderer. Document-authoring dependencies are separate from the statistical requirements lockfile.

Before submission, all group members should review the report, confirm author names and the course requirement regarding overlap with existing research, and confirm the webpage is accessible to the intended class audience.

## GitHub automation and publication

The `Validate project` workflow runs unit tests and saved-result/asset checks on pushes and pull requests. The `Publish demonstration webpage` workflow is manual and does not run on push. Run it from the Actions tab after updating website assets to deploy `website/dist`. Both the repository and the demonstration webpage, including its report and code/results downloads, are public with the owner's authorization.

No course submission has been made. See [the submission checklist](SUBMISSION_CHECKLIST.md) for remaining group confirmations.

## Exploratory supplementary analyses

The original proposal and frozen r=2 main experiment are unchanged. After inspecting the main results and user-supplied exploratory descriptions, we specified two additional analyses in `outputs/supplement/design.json` before running their fresh verification simulations. They are post hoc exploratory analyses, not preregistered main analyses. All supplementary experiments use the same dependency environment as the original numerical work.

1. **Preference sensitivity:** the entire original 12-cell grid at r=3, with 100 new datasets per cell and root seed 23604.
2. **Separation experiment:** 1,000 informative users, with settings (m=12, r=2) and (m=6, r=3). Add A=50 or 200 uniform, preference-free users with D/m=1, 2, 4, or 8; include A=0 once per setting. Each of the 18 cells has 100 replicates. Root seed 23605 generates informative and background samples. Within a setting/replicate, informative purchases are identical across cells; background users and purchases are nested. Total users and purchases increase with background exposure, unlike the matched-total main design.

All four methods share a recorded k-means seed within each supplementary dataset, and across separation cells within a replicate. This helps verify the exact normalized-matrix equivalence of Count and Weighted when all contributing users have the same degree. The unchanged main study used its original independently recorded method seeds.

A method-specific permutation reference uses 20 random permutations of the balanced true product labels per dataset (seed namespace 23606), keeping fitted labels fixed and repeating optimal label matching. No universal 63% chance-error line or SBM detection threshold is assumed. Per-cell confidence intervals use independent replicate units; comparisons across exposure levels within the separation experiment are correlated. Paired method intervals are unadjusted for multiplicity.

```bash
python run_experiments.py --stage supplement --workers 4 --resume
python scripts/validate_supplement.py
```

On a machine with different Python/package metadata, use a new output folder to reproduce the same statistical design:

```bash
python run_experiments.py --stage supplement --workers 4 --output outputs/supplement_reproduced
```

A post hoc SVD ablation reruns the unchanged bipartite SVD on the recorded focused high-heterogeneity datasets after removing users with exactly one purchase. It reproduces the saved SVD errors before applying the change:

```bash
python scripts/svd_ablation.py --workers 4
```

Results are in `outputs/supplement/svd_ablation/` and `report-assets/svd_ablation.md`. Recovery figures in the report and webpage show the random-label reference as a dashed line.

Supplementary results, permutation references, paired intervals, and figures are under `outputs/supplement/`. `supplement.py` runs the analyses; `supplement_figures.py` rebuilds plots. The checks verify that the main-result and original-proposal SHA-256 hashes are unchanged.
