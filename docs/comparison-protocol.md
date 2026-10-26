# Comparing CMAQ and AERMOD without comparing unlike quantities

## Two different kinds of verification

1. **Port verification:** Python versus the original implementation of the
   same routine/model on identical inputs. Close or exact agreement is a
   meaningful engineering target.
2. **Inter-model study:** CMAQ versus AERMOD on carefully harmonized scenarios.
   They solve different representations of atmospheric transport. Differences
   are scientific results to explain, not necessarily bugs or test failures.

The new comparator addresses the second task after data have been aligned.
It emits descriptive metrics, not a pass/fail threshold, equivalence label,
or model accuracy ranking.

## First case: passive-tracer source increment

Agree on the tracer, molecular weight, source geometry, release height,
temperature, exit velocity, time-varying emission rate, terrain, meteorology,
averaging windows, and common sampling support.

For the initial study, explicitly disable chemistry, dry deposition, wet
scavenging, and phase transfer for the tracer in both model setups. Review
surrogates: the CMAQ inert-tracer tutorial includes examples that still
undergo removal. Decide whether plume rise is disabled/fixed for a diagnostic
case or calculated by both models; do not conflate differing injection
heights with differing diffusion.

Use a source-isolated tracer or source-on minus source-off fields.
Document initial/boundary tracer values and subtraction. Do not compare a
CMAQ total field including other sources/background against an AERMOD
single-source increment. Signed increments are retained, not clipped.

## Required alignment

### Time

- Explicit UTC interval start/end; define averaging over that interval.
- Resolve AERMOD hour-ending conventions and local-standard-time offsets.
- Resolve CMAQ TFLAG and instantaneous versus averaged output semantics.
- Address spin-up, restart state, daylight-saving ambiguity, missing hours.
- Align the actual intervals, not merely a date/hour label.

### Horizontal and vertical support

AERMOD values at individual receptors do not directly equal a CMAQ
cell/layer mean. Options include:

- Average a documented receptor quadrature across the cell area and the
  layer thickness with appropriate weights.
- Use a separately justified common observation operator, explicitly labeling
  the point-versus-volume approximation and quantifying sensitivity.

For the preferred volume-average comparison, increase horizontal and vertical
sampling until the averaging error is characterized. Merely choosing the
nearest cell or averaging ground receptors is insufficient. Verify projections,
cell bounds, terrain reference, layer interfaces, and vertical weighting.
The initial CLI does not implement any of these operations.

### Species, units, and forcing

- Choose the same physical tracer, not a similarly named mechanism surrogate.
- Convert to a common mass concentration basis, including molecular weight,
  pressure/temperature, and dry/moist basis when applicable. Never assume a
  universal ppb-to-ug/m3 constant.
- Check whether a reported quantity is concentration, deposition flux, column
  abundance, aerosol mode sum, or a diagnostic.
- Harmonize meteorological forcing as far as scientifically possible;
  document remaining stationary/transient and profile differences.
- Match source emissions, total emitted mass, injection distribution, and
  active processes. Use sensitivity tests for grid numerical diffusion.

## Initial normalized JSON format

The files in [examples](../examples) are invented and clearly labeled.
No real model outputs are included.

Each dataset has exactly five top-level keys:

```json
{
  "schema_version": 1,
  "model": "CMAQ",
  "provenance": "Run/executable/input manifests and alignment procedure identifiers",
  "contract": {
    "tracer_id": "defined-physical-tracer",
    "source_id": "defined-source-configuration",
    "emissions_id": "emissions-manifest",
    "meteorology_case_id": "harmonized-forcing-study",
    "spatial_support_id": "projection-cell-bounds-and-averaging-procedure",
    "vertical_support_id": "layer-interfaces-and-weighting-procedure",
    "processes_id": "transport-only-no-removal",
    "units": "ug/m3",
    "quantity": "passive_tracer_source_increment"
  },
  "records": [
    {
      "start_utc": "2018-07-01T00:00:00Z",
      "end_utc": "2018-07-01T01:00:00Z",
      "support_id": "cell-001-layer-001",
      "value": 1.0
    }
  ]
}
```

The counterpart uses `model: "AERMOD"` and its own provenance. All contract
fields must match exactly. These IDs must refer to reviewed definitions and
manifests, not arbitrary labels chosen to bypass a mismatch. **Equal text
does not prove equal physical support or forcing.** The current CLI trusts
those user declarations and says so in its report.

Records pair by UTC start, UTC end, and support ID. Order may differ;
duplicate keys, overlapping intervals within a support, missing/extra
samples, nonfinite/non-numeric concentrations, zero/negative durations,
non-UTC/naive timestamps, unknown fields, and incompatible contracts fail.
No implicit intersection, interpolation, conversion, or missing-value fill
is performed. The in-memory normalized JSON tool is for modest aligned
datasets, not full multi-gigabyte CMAQ arrays.

The provenance string can point to executable/configuration/input hashes,
the extraction command, baseline subtraction and averaging procedure.
The comparison report also hashes each complete normalized input file.
Raw netCDF/AERMOD inputs are not independently verified by this CLI.

## Metrics

For each paired value define `d = CMAQ - AERMOD`, in ug/m3.

- Mean bias: arithmetic mean of d.
- MAE: arithmetic mean of abs(d).
- RMSE: square root of mean(d*d).
- Maximum absolute difference: maximum abs(d).

All pairs are equally weighted. For a volume-/area-/duration-weighted
domain statistic, perform a reviewed aggregation first or implement an
explicit weighted metric later; do not interpret the current mean as an
automatic mass or exposure budget. Different support sizes and interval
lengths can change interpretation even when pairing is correct.

Correlation, normalized bias, and ratio metrics are intentionally absent
from the first version; zero/near-zero tracer values and degenerate series
need explicit handling before such metrics are meaningful.

## Reproducible tool-only demonstration

```powershell
python -m cmaq_port compare '.\examples\synthetic-cmaq.json' '.\examples\synthetic-aermod.json'
```

Invented pairs are (3,1) and (1,3): errors +2 and -2. Expected bias is 0,
MAE=2, RMSE=2, maximum absolute difference=2. This verifies the tool, not
either atmospheric model.

## What a real study must report

Input and executable provenance; scientific configuration; sampling/
alignment choices; source mass balance; spin-up; paired coverage; actual
metrics and maps; grid/time/sampling sensitivity; limitations. No claim
of observed model agreement can be made until both runs and alignment
have actually been completed and checked.
