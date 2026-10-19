# Checks performed, 2026-10-01

Passed:

- Parsed all 17 R files with a real R interpreter through WebR, as well as an
  independent R syntax parser.
- Ran the included R comparison/status tests: unequal cell areas, expected
  area-weighted overlap, duplicate selection rejection, grid mismatch rejection,
  restored-area mismatch rejection, and accepted/rejected AMPL status classes.
- Loaded both AMPL models and solved small synthetic examples with AMPL/HiGHS
  1.15.1. Verified a whole-polygon 30% target with a subarea below 30%, the forest
  destination restriction, all five anchored certificates, and a positive
  optimal shortfall under an inadequate restoration budget.
- Exported all variable values, reloaded them into a fresh AMPL model, and
  recovered the same ecological objective.
- Evaluated the actual R `gp_stage()` run-file construction and executed its
  emitted AMPL code with Gurobi 13.0.3. Both a zero-shortfall ecological solve
  and a positive-shortfall feasibility solve passed. The latter exported partial
  network tables and a complete checkpoint, exercising the formerly manual
  handoff. The ecological generated run also passed with HiGHS.
- The generated runs check solver result codes, constraint and variable-bound
  residuals, and decision integrality before exporting an accepted incumbent.

Not executed:

- The complete spatial R pipeline against the actual Gudenå rasters, national
  patch data, shapefile and biodiversity surfaces. Those inputs were not attached.
- A full-scale Gurobi optimization, empirical runtime comparison, or exact
  10 m connectivity result for the real landscape.
- The six maps against the user's QGIS SLD and raster category table.
- A comparison with the independent binary model, which was not supplied.

The synthetic examples establish model/interface behavior, not full-data
feasibility or speed. Run the pipeline in the project environment to obtain the
real solver logs and independent spatial audit.
