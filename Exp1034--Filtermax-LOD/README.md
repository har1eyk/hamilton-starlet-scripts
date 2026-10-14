# Exp1034 — FilterMax F5 DNA A260 screening LOD

## Read the DNA plate at A260

Use [read_dna_a260.py](read_dna_a260.py) for the **DNA dilution series in the UV-transparent
Corning 3635 plate**. It reads all 96 wells once at **260 nm** using the same acquisition
sequence and portrait geometry as the successful Orange G read. It runs independently of
the pipetting notebook and does not require an active preparation session.

From this experiment directory, using the `pylabrobot-env` Python environment:

```bash
export PYTHONPATH=/home/hamilton-robot/Documents/Hamilton-Starlet/pylabrobot
python read_dna_a260.py --dry-run
```

With no flags, the program also previews offline: no reader is constructed, no serial port
is opened, and no result directory is created. To acquire a plate, close SoftMax Pro and
other F5 connections, then run:

```bash
python read_dna_a260.py --run-reader
```

The reader defaults to serial **1191** on
`/dev/serial/by-id/usb-Prolific_Technology_Inc._USB-Serial_Controller_BFADb2A8J14-if00-port0`.
The supplied filter catalog resolves **260 nm, excitation slide 2, slot 1, 12 nm bandwidth**.
The script checks the connected reader, readiness and installed slide before moving the
tray; it never substitutes another wavelength. Use `--port`, `--expected-serial`, `--catalog`,
`--geometry` or `--output` to override those inputs. `--geometry` accepts a JSON object of
PLR `PlateGeometry` fields, as displayed by the preview.

Type `APPROVE DNA A260 F5 RUN` for that physical run. When the tray opens, load the plate in
**portrait orientation, A1 upper-right and nearest the instrument**, viewed from the front.
Remove the lid/seal, inspect for bubbles, seat the plate and clear the tray path before typing
`LOADED`. The reused reader profile has an absorbance height of **10.27 mm**; this is an optical
reader setting, independent of the STARlet's pipetting heights. It is not an independently
calibrated Corning profile. The driver closes the tray for measurement and ejects it afterward.

The confirmed map contains **72 DNA samples and 24 A31 matching-buffer blanks** at 50, 75, 100
and 200 µL. A30–A25 are labeled with their **planned** concentrations of 16, 8, 4, 2, 1 and
0.5 ng/µL. These labels and planned volumes are not measurements of the physical contents.
Sample IDs, source tubes and dispensing-replicate numbers retain the existing assignments.

Every acquisition creates a fresh `runs/dna-a260-f5-<timestamp>-<id>/` containing:

- `raw_a260.csv`: all 96 uncorrected readings, identities, timestamps and data-status flags.
- `blank_corrected_a260.csv`: raw readings and separate corrections using the six A31 blanks
  at the same planned volume; negative corrected readings are retained.
- `blank_summary.csv`: four blank means and **sample SDs**, blank well IDs and finite counts.
- `raw_backend_results.json`: the original reader results, preserving nonfinite values as strings.
- `well_manifest.csv` and `reader_metadata.json`: assignments, geometry, filter, instrument
  identity, temperature, hashes, timestamps and completion/failure status.

Missing/nonfinite readings are empty and flagged in CSV, never replaced with zero. If any of
the six blanks at one volume is unavailable, that group's mean, SD and corrections are left
empty and flagged; other volumes remain usable. Existing output directories are protected,
and interrupted measurements are never automatically repeated.

No path-length correction, reference subtraction, shaking or kinetic sequence is requested.
The program does not infer DNA concentration, purity or LOD from A260. Raw CSV and metadata
match the existing `analyze_a260.py` inputs, but `valid` and `exclusion_reason` are deliberately
left blank. Review a **copy** of the raw CSV and fill those fields before using that analyzer;
do not supply the already blank-corrected CSV, which would subtract background twice.

Offline tests use mocked readers only:

```bash
python -m unittest read_dna_a260_tests read_orange_g_tests
```

## Read the Orange G test plate on the FilterMax F5

Use [read_orange_g.py](read_orange_g.py) for the prepared **Orange G plate, with no DNA**.
It reuses all 96 assignments in `well_assignments.json`, including source-tube, planned-volume
and replicate IDs. It does not import the original DNA concentration labels or calculate DNA LOD.
This is one endpoint absorbance reading of every well, not a wavelength sweep or a grid scan
within each well.

The default is **450 nm absorbance**, using excitation slide **2**, slot **4**, bandwidth **8 nm**
in the saved `../helper-scripts/filtermax_slides.xml` catalog. The script checks the installed
slide before tray movement and does not substitute a fluorescence filter. Orange G's reported
visible maximum is [473–479 nm](https://www.sigmaaldrich.com/GB/en/product/aldrich/861286);
450 nm is an available off-peak test wavelength. Results should not be treated as equivalent
to the earlier OPTIMA measurement at 492 nm. Use `--wavelength` with another visible absorbance
filter only if it is described in the supplied catalog and physically installed.

From this experiment directory, check the plan without communicating with hardware:

```bash
export PYTHONPATH=/home/hamilton-robot/Documents/Hamilton-Starlet/pylabrobot
python read_orange_g.py --dry-run
```

For an actual read, close SoftMax Pro and other F5 connections, then run:

```bash
python read_orange_g.py --run-reader
```

The saved connection defaults are F5 serial **1191** and
`/dev/serial/by-id/usb-FTDI_US232R_FT4QLERM-if00-port0`. If the adapter has changed, supply
`--port /dev/serial/by-id/YOUR_F5_ADAPTER`; the script never guesses a different adapter.
`--expected-serial` selects a different known F5. Each physical run requires typing
`APPROVE ORANGE G F5 RUN` before opening its connection, followed by `LOADED` after loading.
The agent has only checked the script offline, without operating the reader.

The starting geometry is PLR's **Costar 96 clear portrait** preset, also used by the local
F5 helper. It is **not a separately calibrated Corning 3635 reader profile**. Review the
dimensions printed by `--dry-run`; use `--geometry /path/to/geometry.json` to supply a matching
`PlateGeometry` JSON object with the fields shown under `plate_geometry` in that output.
Do not use this preset if it does not match the physical plate. In portrait orientation, load
**A1 at the upper-right, nearest the instrument**, viewed from the front. Remove the lid/seal,
check for bubbles and seat the plate before typing `LOADED`. The current PLR driver closes
the tray for acquisition and ejects it on completion.

Each run creates a fresh `runs/orange-g-f5-.../` directory containing:

- `absorbance.csv`: 96 raw OD values with well IDs, sample IDs, source tubes and planned volumes.
- `raw_backend_results.json`: original PLR result, including temperature and measurement timestamp.
- `well_manifest.csv` and `reader_metadata.json`: assignments, settings, instrument identity,
  catalog and script hashes, completion state and any reported failure.

No pathlength correction or reference-wavelength subtraction is applied. Missing or nonfinite
values remain empty in CSV and are flagged, rather than converted to zero. Existing output
directories are never overwritten. Interrupted reads are not automatically repeated.

If **A31 contained dye-free matching buffer**, add `--blank-correct`. This retains the raw OD,
adds corrected OD and saves `blank_summary.csv`. Each well is corrected using only the six A31
blanks at its own planned volume (50, 75, 100 or 200 µL); negative results are retained. If any
blank is missing/nonfinite, that volume group's correction is withheld and flagged. These
are descriptive dye-test results; planned volumes are not measurements of remaining liquid.

Offline regression checks:

```bash
python -m unittest read_orange_g_tests
```

## Prepare the original experiment

Open **`Exp1034--Filtermax-LOD.ipynb`** and run its cells in order. It defaults to simulation,
using the same PyLabRobot LiquidHandler calls that the physical run uses. The operator confirmed the 300 µL A25 tip reach test. Complete plate dispensing accuracy
and mixing performance remain under evaluation.

The notebook preserves the supplied 96-well map and the single twofold dilution chain. The
standalone `exp1034.py` remains an offline planning entrypoint; `notebook_protocol.py` adds the
executable STARlet workflow. Keep the notebook and its Python/JSON sidecars together.

## Start directly at plate filling with prepared tubes

Open [Exp1034--Plate-only.ipynb](Exp1034--Plate-only.ipynb). Close the previous experiment's
kernel, enter the current total volumes in A31 and A30–A25, and run its cells in order. This
creates a new session/output directory, initializes tracking from those totals, and plans only
the 96 plate transfers with seven 300 µL tips. Preparation and refill are not replayed or marked
as having happened in the new session. Use an empty plate, fresh tips in A1–G1, and no attached tip.

Each DNA tube needs at least 1500 µL; A31 needs at least 2775 µL of clean matching buffer with
the default reserves. A25 may contain 1500 or 3000 µL: enter its actual current total, which
determines aspiration height. The prepared DNA concentrations must be 16, 8, 4, 2, 1 and
0.5 ng/µL in A30–A25. A32 is unused; its zero tracker in this mode is not a physical measurement.
Liquid-level detection remains off, so command completion does not establish physical delivery.

The plate-only notebook requests 4 mm lower tube aspirations, stopping at 3 mm above the inner
bottom to retain the checked insertion margin. It flags transfers limited by this clearance and
shortens surface following accordingly. Plate dispenses use a fixed 2.50 mm above the modeled
inner well bottom for every volume, with no surface-derived height or extra lowering. Volumes,
rates, dispense modes and tip grouping are unchanged. The review table shows the modeled
inner-bottom Z and dispense tip Z; offline checks do not verify physical clearance.

The equivalent API is `ExperimentSession(settings, new_output_dir,
prepared_source_volumes=volumes)`, followed by `await session.start()` and
`await session.execute_phase("plate")`. Do not reset the flags of a stopped or completed session.

## Fresh preparation after evaporation

Reset the deck to the initial quantities below, restart the notebook kernel, and run from the top.
The notebook creates a new `ExperimentSession` and unique run directory. Do not replay Friday's
session or use its plate-resume plan: the evaporated samples no longer match that checkpoint.
Friday's records remain preserved in `runs/physical-20261004T001957Z-fab79d55`.

Use A30 at 3000 µL total and measured 16 ng/µL, A31/A32 at 4000 µL buffer each, empty A29–A25
and plate, replenished racks, and no attached old tip. The fresh plan fills all 96 wells with seven
300 µL plate tips; it neither preserves the old A9 blank nor subtracts the old diagnostic aliquot.
A25 finishes at 1725 µL for this fresh run.

## Running a fresh experiment

1. Restart the Python kernel, then run the import and configuration cells. The notebook selects
   `/home/hamilton-robot/Documents/Hamilton-Starlet/pylabrobot` explicitly and prints its version
   and commit. It was checked against branch `branch-merge-for-Exp1034`, commit `1cfbf451b`.
2. Leave `simulate=True` for a complete offline run. Review the source-volume ledger, deck,
   tip budget, calculated heights and expanded manifest before the setup cell.
3. For a physical run, set `simulate=False` and enter the **measured** starting concentration
   and operator. `method_validation_reference` is optional: leave it blank for the first run.
   A blank reference is recorded as an exploratory run with no prior method validation supplied.
   If prior checks exist, enter a notebook reference, file path, or notes. The notebook checks
   the concentration, operator and modeled parameters before asking for `APPROVE EXP1034 RUN`,
   then initializes STARBackend. This is the repository's explicit approval for that physical run.
4. Run the preparation cell, acknowledge the clean-buffer refill, then run the plate-fill cell.
   A fresh output folder protects each run. Do not replay completed cells or partially failed steps.
5. Transfer the plate manually. The optional F5 cell is off by default and requires separate reader
   approval, an actual preparation record, an explicit port and verified reader profile. The
   analysis instructions follow in the notebook.

The notebook does not ask again for the carrier identities, track correction or plate dimensional
correspondence supplied by the user. The measured A30 concentration remains a required run input.
Physical accuracy, mixing and low-volume accessibility remain to be established during method
development; an exploratory run does not certify them.

## Corrected deck and resource provenance

All tracks are carrier left-edge tracks; module/site positions are zero-based. The complete model
assembles without bypassing collision checks or using the deprecated `rails` placement argument.

| Location | Resource and contents |
| --- | --- |
| Track 19 | `hamilton_tube_carrier_32_a00`, bare carrier, with eight McKesson 177-1505 12 × 75 mm glass tubes. A32, A31, A30, A29, A28, A27, A26, A25 map explicitly to sites 0–7, front to back. |
| Track 7 | `hamilton_mfx_carrier_L5_base`, with five `hamilton_mfx_plateholder_DWP_metal_tapped` modules. |
| Plate module 0 | Corning **3635 UV-transparent plate**, using `cor_96_wellplate_360uL_Fb` geometry, as confirmed dimensionally equivalent by the user. A1 is rear-left on the STAR deck. |
| Track 25, position 0 | `hamilton_96_tiprack_1000uL_filter`; 96 tips. |
| Track 25, position 1 | `hamilton_96_tiprack_300uL`; 96 nonfiltered tips, matching the confirmed catalog choice. |
| Track 25, position 2 | `hamilton_96_tiprack_50uL`; 96 nonfiltered tips, unused in this method. |

`TIP_CAR_480_A00` is implemented through its current equivalent `hamilton_tip_carrier_L5`.
`Cor_96_wellplate_360ul_Fb` is the deprecated alias of the plate factory used above.

The current PLR checkout contains the bare carrier definition. The McKesson file was located on
`add-12x75mm-test-tube`, commit **`d62db0a4234d3e87d3c1509a8e2a56b0d1ff1976`**. Its
`pylabrobot/resources/mckesson/tubes.py` is copied unchanged into `mckesson_tubes.py`, so the
notebook uses that exact resource without switching branches, merging, or editing the PLR checkout.
The explicit imports remain reviewable. The runner records the local definition's hash.

Tube mounting follows the supplied
`../Exp1031--Endotoxin-testing/endotoxin_assay_w_12x75mm_tubes.ipynb`: tube rim 25.6 mm above the
carrier, 1.0 mm bottom material thickness, and centers in the bare holders. The branch's measured
height-volume fit is `height = 3.622768533 + 0.01175032 × volume`, in mm and µL. Its intercept
at zero is a fit artifact, not a physical liquid surface or proof of accessible dead volume.

The resolved geometry, XYZ tube centers, cavity bottoms and rims are exported in
`resource_audit.json`. Aspiration heights account for the configured total loss allowance and
follow the source's modeled falling surface. Buffer-prefill dispense heights are 1 mm below the
destination's final surface; DNA transfers and plate fills remain above it. Mixing heights stay
submerged at the transient low volume. After physical setup,
actual channel Y reach is checked using the backend's discovered configuration before pipetting.

## Tip reach and budget

The operator-confirmed test reached 1 mm above A25's inner bottom with the 300 µL B1 tip,
aspirated 300 µL and retracted. This checks 73 mm insertion for the unchanged McKesson tube and
deck geometry. The planner retains a 2 mm margin within that checked depth, and enforces source
floors, bottom clearance and capacity. The deepest plate aspiration is approximately 70.03 mm.
This reach check does not establish volumetric accuracy, mixing or low-volume recovery.

With the notebook's **200 µL dead-volume floor + 25 µL total loss allowance per source**, a fresh run uses:

| Rack | Preparation tips | Plate tips | Total used | Remaining |
| --- | ---: | ---: | ---: | ---: |
| 1000 µL filtered | 5 | 0 | 5 | 91 |
| 300 µL nonfiltered | 0 | 7 | 7 | 89 |
| 50 µL nonfiltered | 0 | 0 | 0 | 96 |

One 300 µL tip fills all 24 blanks from A31. Each of A30–A25 uses its own 300 µL tip for all
12 wells at that concentration, including all four volumes. Tips change between source groups
and are discarded after their group. Each well receives one noncontact dispense of its complete
50, 75, 100 or 200 µL volume. The planner rejects settings beyond the checked 300 µL tip reach
instead of substituting 1000 µL tips. These counts apply to a fresh run with replenished racks.

The nominal preparation still uses only five tips: the first is retained from clean-buffer
prefilling through A30 mixing and the A30→A29 link; subsequent links get fresh tips. No DNA-exposed
tip re-enters blank buffer or a different concentration. Transfer and mixing aliquots are at most
750 µL by default. Air/prewet/over-aspiration commands are zero, and hidden liquid-class volume
correction is disabled. Both `preparation_aspirate_flow` and `preparation_dispense_flow`
default to **250 µL/s** for every preparation transfer, including buffer prefills and DNA transfers.
Mixing uses **400 µL/s in both directions**; plate transfers retain 100/120 µL/s aspiration/dispense.
The notebook explicitly sets these rates and displays them before setup. Rates are saved in
`planned_operations.json`. Setting a preparation rate to `None` explicitly opts into the installed
liquid-class flow default for that direction; other liquid-class parameters are not adopted.


## Preparation and reagent budget

A30 must contain **3000 µL TOTAL solution at measured 16 ng/µL**, including its DNA contribution.
This concentration is required by the supplied design, not an assumed measurement. The notebook
rejects a different concentration rather than silently labeling incorrect tubes. Stock preparation
from a concentrated source is outside this deck setup. Use identical DNA-free buffer for dilution
and blanks; A32 and A31 initially contain 4000 µL each, and A29–A25 are empty.

1. With one 1000 µL filtered tip, prefill each empty A29–A25 with 750 µL from A32 plus 750 µL
   from A31. Surface dispensing uses the McKesson height function at the destination's final total,
   minus `buffer_prefill_submersion=1.0` mm: 11.44 mm at 750 µL and 20.25 mm at 1500 µL,
   above the inner bottom. The tip contacts the liquid near the end of each addition. Reuse is
   permitted only while the contacted destination contains no DNA; the planner rejects otherwise.
2. Homogenize A30, then transfer A30→A29 as **two 750 µL aliquots**, both noncontact. Mix A29
   only after both aliquots have arrived, then discard the tip.
3. Repeat A29→A28→A27→A26→A25 with a fresh tip per link and destination mixing. Each mixture
   reaches 3000 µL before its onward transfer. Do not remove an extra aliquot from A25.
4. After the last DNA-exposed tip is discarded, restore A31 to **4000 µL TOTAL**. Nominal
   addition is **3750 µL**. Tracking changes only after the physical `REFILLED A31` acknowledgement.
   Earlier clean-buffer top-ups are inserted if configured buffer withdrawal floors require them.
5. Fill the exact manifest directly, with no in-well dilution: all blanks, then the 16, 8, 4, 2,
   1 and 0.5 ng/µL groups. Reuse one tip throughout each group and rotate volume order within it.
   All original well positions, volumes and replicate labels remain.

The default mix is 750 µL × 8 cycles, at 6 mm above the cavity bottom. It is configurable and
capacity/immersion checked; simulation does not validate physical homogeneity or bubble control.

| Tube | Initial µL | After complete chain µL | Plate withdrawal µL | Nominal final µL |
| --- | ---: | ---: | ---: | ---: |
| A32, buffer | 4000 | 250 | 0 | 250 |
| A31, buffer | 4000 | 250; refill +3750 → 4000 | 2550 | 1450 |
| A30, 16 ng/µL | 3000 | 1500 | 1275 | 225 |
| A29, 8 ng/µL | 0 | 1500 | 1275 | 225 |
| A28, 4 ng/µL | 0 | 1500 | 1275 | 225 |
| A27, 2 ng/µL | 0 | 1500 | 1275 | 225 |
| A26, 1 ng/µL | 0 | 1500 | 1275 | 225 |
| A25, 0.5 ng/µL | 0 | 3000 | 1275 | 1725 |

Dilution consumes 7500 µL buffer; blanks consume 2550 µL; total buffer delivered is **10,050 µL**.
The initial 8000 µL is short by at least 2050 µL. The selected refill supplies 11,750 µL buffer
and leaves 1700 µL in the buffer tubes. Total plate volume is **10,200 µL**.

The 225 µL DNA residual must cover inaccessible liquid plus all losses. The source floor is
checked on every aspiration, including during mixing; insufficiency fails before the run. The
loss allowance is conservative, not an assertion of measured loss. No unrequested ninth deck tube
or unaccounted starting DNA volume is introduced. Any extra actions require revised budgets.

## Execution records and failure handling

Each notebook run writes a fresh `runs/simulation-*` or `runs/physical-*` folder containing:

- `well_manifest.csv`: all 96 well/source/concentration/volume/DNA/replicate assignments.
- `dilution_ledger.csv`, `planned_operations.json`, `verification.json`: nominal volumes,
  capacities, mixed-tip budget, resolved heights and flow rates, and predicted source residuals.
- `resource_audit.json`: corrected deck, resource provenance and geometry.
- `execution.jsonl`: timestamps around every attempted/completed operation, aspiration and dispense.
- `preparation_record.json`: operator/method inputs, versions/hashes, refill acknowledgements,
  actual per-well dispensing timestamps, tracked residuals and completion/failure state. Physical
  sessions additionally query STARlet serial and firmware after approved setup.
- `simulation_backend.log`: Chatterbox call output for simulation only.

A failed aspiration or dispense latches the session. There is no automatic replay, tip discard,
tracking reset or recovery movement; inspect the actual deck and logs before planning recovery.
Completed phases cannot run twice in the same session. Recreating a session on a used physical
deck is not a recovery procedure. `session.stop()` closes the connection; it does not recover tips.

The refill prompt accepts `refilled` or `REFILLED A31`, ignoring case and extra spaces.
An unacknowledged refill pauses at that checkpoint without changing tracked volume; rerun the
same phase when ready. If A31 already contains 4000 µL total, acknowledge it without adding more.

For an existing live session stopped by `ValueError('Refill cancelled; no tracking update.')`,
use this cell without restarting the kernel, recreating `session`, or replaying preparation:

```python
import importlib
import notebook_protocol

importlib.reload(notebook_protocol)
await notebook_protocol.resume_refill(session)
print("A31 tracked total:", session.volumes()["A31"], "µL")
```

Enter `refilled` when A31 contains 4000 µL total of clean buffer. The helper verifies that only
the final A31 acknowledgement remains, preparation is complete, the plate is empty and tracking
matches preparation. It clears only that specific acknowledgement failure and keeps the existing
connection. It issues no robot commands. After it succeeds, continue with the plate-fill cell.

For recovery of a run whose physical contents still match its recorded state (not the evaporated
Friday samples), source-group tip reuse can be applied while retaining its kernel and connection:

```python
import importlib
import exp1034
import notebook_protocol

importlib.reload(exp1034)
importlib.reload(notebook_protocol)
await notebook_protocol.update_plate_tip_reuse(session)
```

The update preserves filled wells, source volumes, used-tip positions and the executed plan.
It archives the previous plan and replaces only the remaining plate operations. If a tip pickup
was cancelled, it checks live tip presence against PLR tracking before clearing that specific
failure. An attached empty tip from that pickup is scheduled for discard before the source groups
begin. Unknown tip state or an interrupted aspiration/dispense blocks this update. The update
itself issues no motion or pipetting commands. After it succeeds, run the plate-fill cell once.

## FilterMax handoff and analysis

Corning 3635 uses the 200 µL reference and recommended 75–200 µL range; 50 µL is an experimental
volume to validate. Do not substitute half-area 3679 for this design.
[Corning selection guide](https://www.corning.com/catalog/cls/documents/selection-guides/CLS-MP-014.pdf).

Transfer manually, inspect for bubbles, use consistent mixing and elapsed time, and remove an
ordinary lid/seal from the UV light path. No reader deck location or gripper command is assumed.
The optional notebook reader cell uses the installed `FilterMaxF5` integration and requires a
reader-specific verified 3635 profile, actual port and filter catalog. STARlet dimensional
correspondence does not determine the F5 optical read-height/orientation settings. Complete a
copy of `reader_profile.example.json`; do not mark guessed values verified.

The adapter verifies serial number/readiness and the installed excitation slide's 260 nm absorbance
filter before reading all 96 wells. No reference wavelength or pathlength correction is requested.
The local backend reports live validation at 450 nm; this does not establish A260 performance for
this experiment. Preserve the raw values, exact manifest, temperature, reader/firmware versions,
filter settings, timestamps and preparation record. The optional cell requires a separately typed
`APPROVE F5 RUN`; `LOADED` acknowledges manual plate loading before measurement.

Alternatively export unblanked A260 manually from the validated reader software. Preserve the
original export and settings. A reviewed CSV must contain `well,a260,valid,exclusion_reason,read_id,timestamp`
for each well, one reader repeat, and timezone-aware timestamps. Metadata must include
`wavelength:260`, `pathlength_correction_applied:false`, `reference_subtracted:false` and the
reader/preparation provenance. The analyzer deliberately rejects corrected data.

Copy `calibration_ranges.example.json`, select the observed linear range separately for each volume,
record the reason and set `reviewed:true`. Mark saturated/invalid readings false with a reason.
All six individual blanks are retained; at least three distinct valid concentrations must remain.

```bash
python analyze_a260.py --input reviewed_raw_a260.csv --ranges calibration_ranges.json \
  --reader-metadata runs/PHYSICAL-RUN/reader/reader_metadata.json --output runs/analysis-01
```

For each volume, the analyzer calculates blank mean and **sample SD (n−1), not SEM**, subtracts
that blank mean, fits the selected calibration with a free intercept, and reports:

`preliminary concentration LOD (ng/µL) = 3.3 × blank SD / slope (A260 per ng/µL)`

`preliminary mass LOD (ng) = concentration LOD × well volume (µL)`

All sample rep1–rep3 labels are **dispensing replicates from one dilution chain**, not independent
preparations. Repeated reads of the same wells are reader repeats. Equal concentrations have
different DNA masses at different volumes. This plate is screening, not an established LOD/LOQ;
confirmation requires fresh blanks and independently prepared low-concentration samples in later
runs. Raw-input copies/hashes, exclusions, selected ranges, blank values and fit statistics are saved.

## Offline verification

```bash
cd /home/hamilton-robot/Documents/hamilton-startlet-scripts/Exp1034--Filtermax-LOD
export PYTHONPATH=/home/hamilton-robot/Documents/Hamilton-Starlet/pylabrobot
python -m unittest -v exp1034_tests notebook_protocol_tests
python exp1034.py --a30-concentration 16 --output runs/nominal-review-01
```

The tests include exact-map invariants, mass/volume balance, refill cancellation, source insufficiency,
capacity, tip reach, live API signatures, full Chatterbox execution, refusal before unapproved
hardware construction, partial-step failure latching and the separately mocked reader/analysis path.
The notebook has also been executed from top to bottom in simulation. Simulation verifies software
behavior and modeled geometry; the candidate physical method remains subject to validation.

Scientific source: [Exp1034 DNA LOD in FilterMax F5](https://chatgpt.com/space/page_2505a363d96081918edc799c4d38a292).
The user's supplied revised well map and 3.3 × SD/slope analysis govern this experiment.
