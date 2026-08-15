# Validation

The normal package API uses only the private production implementation under
`R/` and `src/`. It does not load, export, or dispatch to V34.

`run_validation.R` provides the separate release-validation path. It checks
construction success, input immutability, known-true values, and optional
hidden V34 parity. It supports:

- `--scopes=corpus,known-true,anolis-empirical,anolis-simulated`
- `--population=all|ultrametric|nonultrametric`
- `--case-ids=id1,id2`
- `--include-label=REGEX` and `--exclude-label=REGEX`
- `--min-tips=N`, `--max-tips=N`, and `--exclude-tips=N1,N2`
- `--start=N` and `--limit=N`
- `--compare-v34=true|false`
- `--tolerance=1e-8`
- `--output=PATH`

Examples:

```sh
Rscript tools/validation/run_validation.R \
  --scopes=known-true \
  --output=tools/validation/output/known_true

Rscript tools/validation/run_validation.R \
  --scopes=corpus \
  --exclude-label='^generated_' \
  --max-tips=100 \
  --output=tools/validation/output/filtered_corpus

Rscript tools/validation/run_validation.R \
  --scopes=corpus,anolis-empirical,anolis-simulated \
  --compare-v34=true \
  --output=tools/validation/output/full_parity
```

The complete known-true inventory is
`reference/true/manifest_complete.csv`: 129 registered artifacts spanning PST,
PST+SS, PST+TT+SS, and original-metric scopes. The previous 91-row manifest is
retained beside it for provenance.

`historical/` contains compact proof summaries from the final pre-package V35
release run. Detailed logs and per-failure packets were intentionally omitted.
