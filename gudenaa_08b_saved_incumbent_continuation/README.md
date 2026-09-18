# 08b saved-incumbent continuation

This package continues the successful 08b constructive network repair from the
**actual final 08b incumbent**, rather than returning to the earlier
constructive start.

The ordinary 08b TSV files are retained as audit and mapping outputs. They do
not contain every continuous network-flow variable, so they are not by
themselves a complete solver checkpoint. The first helper below serializes the
full incumbent to an AMPL `.run` file. The second helper can then reset and
reload that checkpoint in a fresh AMPL session.

## Step 1 — save the current 08b incumbent

Do this while the completed 08b solution is still the current solution at the
`ampl:` prompt, before `reset;` or closing AMPL:

```ampl
include "gudenaa_08b_saved_incumbent_continuation/01_save_current_08b_incumbent.run";
```

It writes:

```text
Optimization/gudenaa_new_area_200m_global30_constructive_network_repair/
  gudenaa_new_area_200m_global30_constructive_network_repair_saved_incumbent.run
  gudenaa_new_area_200m_global30_constructive_network_repair_saved_incumbent_audit.tsv
```

Open or inspect the end of the generated incumbent file once. Its final line
must be:

```text
# CHECKPOINT_COMPLETE
```

The audit should report `network_shortfall_total_ha = 0` and the footprint and
quality values seen at the end of 08b.

## Step 2 — run the reset-safe continuation

Once Step 1 has completed, this command may be run immediately or later in a
fresh AMPL session:

```ampl
include "gudenaa_08b_saved_incumbent_continuation/02_continue_saved_08b_incumbent_8h4h.run";
```

The continuation deliberately performs `reset;`, then loads:

1. `Gudenaa_global30_constructive_network_repair_200m.mod`
2. the matching 08b `.dat`
3. the saved incumbent from Step 1

It then runs the lexicographic sequence:

| Stage | Objective | Time limit | Gap target |
|---|---|---:|---:|
| 2 | Minimize change from the notebook-07 footprint | 8 hours | 0.5% |
| 3 | Maximize biodiversity and connectivity quality at the Stage-2 footprint | 4 hours | 0.5% |

At the start of Stage 2, Gurobi should print a line like:

```text
Loaded user MIP start with objective ...
```

That objective should agree with
`saved_incumbent_audit.tsv::warm_footprint_change_ha` (approximately the final
08b footprint), **not** the original constructive value of 8772.19 ha.

The hard constraints remain unchanged: zero network shortfall, one global
certificate of at least 5,000 ha, one certificate of at least 1,000 ha in each
subarea, the global 30% nature budget, and the forestry/wet-soil transition
rules. Stage 3 retains the best Stage-2 footprint with only 0.01 ha of numerical
slack.

## Outputs

All continuation outputs use the suffix `continued_8h4h`, so the original 08b
files are not overwritten. Important files are:

- `*_continued_8h4h_start_audit.tsv`
- `*_continued_8h4h_stage_summary.tsv`
- `*_continued_8h4h_solution_summary.tsv`
- `*_continued_8h4h_network_summary.tsv`
- `*_continued_8h4h_selected_cells.tsv`
- `*_continued_8h4h_stage2_minimum_change_gurobi.log`
- `*_continued_8h4h_stage3_quality_gurobi.log`
- `*_continued_8h4h_final_incumbent.run`
- `*_continued_8h4h_completion_status.txt`

The final incumbent `.run` is another complete checkpoint and ends with the
same `# CHECKPOINT_COMPLETE` marker, so further continuation is possible
without relying on a live AMPL session.

## Important distinction

The `.mst` files written by Gurobi are useful solver-native records of the
integer solution at each stage. The AMPL `.run` checkpoints are the files used
for explicit reset-safe reloading here because they contain assignments in the
model's AMPL names and also preserve the network-flow values.
