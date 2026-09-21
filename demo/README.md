# Offline synthetic demo

A small entry point into the repository's pathway-role statistics. **Every count is invented; there are no real genes, patient records or study-derived results in this example.**

## Run

From the repository root, with Python 3.10 or newer:

```bash
python3 demo/run_demo.py
python3 -m unittest discover -s tests -v
```

No third-party packages or downloads are required. The script always uses the bundled synthetic fixture and writes `results.csv` and `preview.svg` to `demo/output/`, which is ignored by Git. To choose another output directory:

```bash
python3 demo/run_demo.py --output-dir /tmp/sf3b4-synthetic-demo
```

## What it demonstrates

1. Validate non-negative integer counts and non-empty groups.
2. Call the **existing** `fisher_exact_one_sided_greater` function from `01_pan_cancer_1q/analyze_1q_pathway_signed_roles.py`.
3. Correct the three illustrative p-values with its existing `benjamini_hochberg` function.
4. Export a tidy CSV and a dependency-free SVG preview.

The tested 2 × 2 table is:

| Toy gene role | On 1q | Elsewhere |
| --- | --- | --- |
| Activator-like | `activator_1q` | `activator_other` |
| Repressor-like | `repressor_1q` | `repressor_other` |

The one-sided alternative asks whether the proportion on 1q is **greater in the activator-like group**. The bars show each group's fraction on 1q, not an effect-size estimate or confidence interval. Depletion is not significant in this test direction.

## Expected output

The fractions on 1q are 80% versus 20%, 50% versus 50%, and 20% versus 80%. The toy-enrichment p-value is approximately 0.0115, with BH-adjusted q approximately 0.0345 across the three scenarios.

Inspect the checked-in [reference CSV](example_output/results.csv) and [reference preview](example_output/preview.svg). Regenerate them with:

```bash
python3 demo/run_demo.py --output-dir demo/example_output
```

Output is deterministic within the same runtime; tests compare numerical reference values with floating-point tolerance. Each CSV row carries `data_type=synthetic_only_not_study_results`.

## Scope of validation

Tests compare the Fisher implementation against exact combinatorial calculations for 625 small tables, check BH correction, exercise signed-role parsing on a toy graph, reject malformed demo inputs, and verify offline output generation. This checks selected functions, not the full research workflow. It does not validate biological hypotheses, reproduce a manuscript figure, test real KEGG downloads or substitute for the controlled inputs needed by other analyses.
