# Patch-ledger findings used to design notebook 08a

The four supplied patch/network files are internally consistent. Their scopes
are different from the earlier `national_patches_vs_landuse_budget_audit.csv`,
which explains the initially surprising totals.

## Linked certificate-patch ledger

- Linked existing patches: **2,277**
- Candidate-to-patch links: **5,107**
- Candidate cells with at least one link: **2,742**
- Duplicate `(coarse_cell_id, patch_id)` links: **0**
- Full national area of linked patches: **6,498.340 ha**
- Area of those linked patches inside the planning polygon: **5,727.936 ha**
- Area of those linked patches outside the planning polygon: **770.404 ha**

The 770.404 ha outside the polygon is eligible only through a linked patch in
the **Global** certificate. It is not restoration area and is not part of any
30% subarea budget.

| Subarea | All national-patch nature inside (ha) | Linked certificate-patch nature inside (ha) | Existing but unlinked under contact rule (ha) | Land-use-raster current nature (ha) |
|---:|---:|---:|---:|---:|
| 1 | 1,841.748 | 1,800.196 | 41.551 | 1,841.800 |
| 2 | 1,025.788 | 962.550 | 63.238 | 1,024.990 |
| 3 | 1,969.483 | 1,860.373 | 109.110 | 1,969.420 |
| 4 | 1,357.144 | 1,104.817 | 252.327 | 1,149.440 |

“All national-patch nature” includes patches with no qualifying contact to a
candidate restoration unit. Those unlinked hectares are current nature, but
they cannot be claimed as part of an anchored network certificate under the
current contact definition. Notebook 08a now writes this distinction explicitly
to `patch_ledger_scope_audit.csv`.

## Source-ledger discrepancy

For Subareas 1–3, all-patch geometry and the current-nature land-use raster
agree closely. For Subarea 4 they differ by **207.704 ha (18.07%)**. Notebook
08a therefore keeps the two accounting systems separate:

- the four 30% budgets continue to use current nature from the land-use raster;
- network certificates use national patch geometry;
- the discrepancy is flagged, not silently corrected or substituted.

## Boundary handling

Eleven linked patches intersect more than one subarea. Notebook 08a stores one
area coefficient per `(patch_id, subarea_id)`. Thus each subarea certificate
counts only its own physical fragment, while the Global certificate counts the
full national patch once. The supplied fragment table reconciles to the linked
inside/outside audit to better than `1e-5` ha per patch.
