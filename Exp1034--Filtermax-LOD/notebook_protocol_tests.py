"""No-I/O tests for the notebook's executable LiquidHandler workflow."""

import ast
import asyncio
import inspect
import json
import tempfile
import unittest
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import AsyncMock, Mock, patch

from pylabrobot.legacy.liquid_handling.backends import STARBackend
from pylabrobot.legacy.liquid_handling.errors import ChannelizedError

import notebook_protocol
from notebook_protocol import ExperimentSession, Settings
from exp1034_deck import tip_insertion_limit


def settings_from_notebook(path: Path) -> Settings:
  """Read literals and constant references in simulation mode without executing cells."""
  notebook = json.loads(path.read_text())
  for cell in notebook["cells"]:
    source = "".join(cell["source"])
    if cell["cell_type"] != "code" or "settings = Settings(" not in source:
      continue
    constants = {}
    for node in ast.parse(source).body:
      if (
        isinstance(node, ast.Assign)
        and len(node.targets) == 1
        and isinstance(node.targets[0], ast.Name)
        and isinstance(node.value, ast.Constant)
      ):
        constants[node.targets[0].id] = ast.literal_eval(node.value)
      if (
        isinstance(node, ast.Assign)
        and isinstance(node.value, ast.Call)
        and isinstance(node.value.func, ast.Name)
        and node.value.func.id == "Settings"
      ):
        values = {
          item.arg: constants[item.value.id]
          if isinstance(item.value, ast.Name)
          else ast.literal_eval(item.value)
          for item in node.value.keywords
        }
        values["simulate"] = True
        return Settings(**values)
  raise ValueError(f"No literal Settings constructor in {path.name}.")


async def failed_buffer_session(output: Path) -> ExperimentSession:
  """Simulate the recorded step-10 aspiration failure, including frontend tracker rollback."""
  session = ExperimentSession(
    Settings(
      channel=3,
      mix_cycles=1,
      preparation_aspiration_lowering=4,
      preparation_minimum_tip_height=1,
      plate_aspiration_lowering=4,
      plate_dispense_height=2.5,
    ),
    output,
  )
  await session.start()
  original_aspirate = session.lh.backend.aspirate

  async def fail_at_step_ten(*args, **kwargs):
    """Reject only the final A32 withdrawal before any simulated liquid movement."""
    if session.next_step == 9:
      raise ChannelizedError(errors={3: RuntimeError("fixture minimum-access rejection")})
    return await original_aspirate(*args, **kwargs)

  with patch.object(session.lh.backend, "aspirate", side_effect=fail_at_step_ten):
    try:
      await session.execute_phase("preparation")
    except ChannelizedError:
      return session
  raise AssertionError("Expected a step-10 aspiration failure")


class NotebookProtocolTests(unittest.IsolatedAsyncioTestCase):
  """Verify software gates, full simulated tracking and irreversible partial-step handling."""

  def test_settings_cell_reloads_cached_protocol(self) -> None:
    """Fresh-run settings must not reuse an older class cached in the notebook kernel."""
    notebook = json.loads((Path(__file__).parent / "Exp1034--Filtermax-LOD.ipynb").read_text())
    cell = next(c for c in notebook["cells"] if "settings = Settings(" in "".join(c["source"]))
    statements = []
    for node in ast.parse("".join(cell["source"])).body:
      statements.append(node)
      if isinstance(node, ast.Assign) and any(
        isinstance(target, ast.Name) and target.id == "settings" for target in node.targets
      ):
        break
    stale = Mock(side_effect=AssertionError("Cached Settings class used"))
    namespace = {"Settings": stale}
    with patch.object(notebook_protocol, "Settings", stale):
      exec(
        compile(ast.Module(body=statements, type_ignores=[]), "Settings cell", "exec"), namespace
      )
      self.assertIs(type(namespace["settings"]), notebook_protocol.Settings)
      self.assertIs(namespace["ExperimentSession"], notebook_protocol.ExperimentSession)
      self.assertEqual(namespace["settings"].preparation_minimum_tip_height, 4.3)
    stale.assert_not_called()

  async def test_number_of_tube_mixes_controls_each_dilution_stage(self) -> None:
    """The editable notebook count controls mixing after both aliquots, with valid bounds."""
    notebook = json.loads((Path(__file__).parent / "Exp1034--Filtermax-LOD.ipynb").read_text())
    settings_cell = next(
      cell for cell in notebook["cells"] if "settings = Settings(" in "".join(cell["source"])
    )
    source = "".join(settings_cell["source"])
    assignment = next(
      (line for line in source.splitlines() if line.startswith("number_of_tube_mixes = ")),
      None,
    )
    self.assertIsNotNone(assignment)
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      for count in (1, 4, 9):
        with self.subTest(count=count):
          settings_cell["source"] = source.replace(
            assignment, f"number_of_tube_mixes = {count}"
          ).splitlines(keepends=True)
          path = Path(directory) / "candidate.ipynb"
          path.write_text(json.dumps(notebook))
          settings = settings_from_notebook(path)
          self.assertEqual(settings.mix_cycles, count)
          session = ExperimentSession(settings, Path(directory) / f"run-{count}")
          for tube in ("A30", "A29", "A28", "A27", "A26", "A25"):
            mixes = [e for e in session.plan if e["action"] == "mix_cycle" and e["source"] == tube]
            self.assertEqual([e["cycle"] for e in mixes], list(range(1, count + 1)))
            if tube != "A30":
              index = session.plan.index(mixes[0])
              incoming = session.plan[index - 2 : index]
              self.assertEqual(
                [(e["action"], e["destination"], e["volume"]) for e in incoming],
                [("transfer", tube, 750), ("transfer", tube, 750)],
              )
          await session.start()
          with patch.object(session.lh, "aspirate", wraps=session.lh.aspirate) as aspirate:
            await session.execute_phase("preparation")
          self.assertEqual(aspirate.await_count, 20 + 6 * count)
          self.assertEqual(session.report["mixing"]["cycles"], count)
          await session.stop()
      for invalid in (0, -1, 10, 11, 2.5, 3.0, True, "3", None):
        with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "integer.*1.*9"):
          ExperimentSession(replace(settings, mix_cycles=invalid), Path(directory) / "invalid")

  async def test_full_notebook_uses_verified_plate_only_heights(self) -> None:
    """The full workflow copies verified plate heights without changing preparation."""
    root = Path(__file__).parent
    verified = settings_from_notebook(root / "Exp1034--Plate-only.ipynb")
    settings = settings_from_notebook(root / "Exp1034--Filtermax-LOD.ipynb")
    self.assertEqual(verified.plate_dispense_height, 2.5)
    self.assertEqual(verified.plate_aspiration_lowering, 4)
    self.assertEqual(settings.plate_dispense_height, verified.plate_dispense_height)
    self.assertEqual(settings.plate_aspiration_lowering, verified.plate_aspiration_lowering)
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      session = ExperimentSession(settings, Path(directory) / "full")
      baseline = ExperimentSession(
        replace(settings, plate_dispense_height=None, plate_aspiration_lowering=0),
        Path(directory) / "baseline",
      )
      self.assertEqual(
        [e for e in session.plan if e["phase"] != "plate"],
        [e for e in baseline.plan if e["phase"] != "plate"],
      )
      operations = [e for e in session.plan if e["phase"] == "plate" and e["action"] == "transfer"]
      totals = {}
      for event in operations:
        totals.setdefault(event["source"], event["source_before"])
      reference = ExperimentSession(
        verified, Path(directory) / "plate-only", prepared_source_volumes=totals
      )
      expected = [e for e in reference.plan if e["action"] == "transfer"]
      self.assertEqual(len(operations), 96)
      self.assertEqual({e["dispense_height"] for e in operations}, {2.5})
      self.assertEqual({e["volume"] for e in operations}, {50, 75, 100, 200})
      fields = (
        "source",
        "destination",
        "volume",
        "aspiration_height",
        "aspiration_end_height",
        "aspiration_minimum_height",
        "aspiration_height_limited",
        "surface_following_distance",
        "dispense_height",
        "dispense_height_mode",
        "dispense_tip_z",
        "tip_penetration_at_end",
        "dispense_jet",
        "dispense_blow_out",
      )
      for actual, original in zip(operations, expected):
        self.assertEqual(
          {key: actual[key] for key in fields}, {key: original[key] for key in fields}
        )
        self.assertGreaterEqual(
          actual["aspiration_end_height"], actual["aspiration_minimum_height"]
        )
      self.assertTrue(any(e["aspiration_height_limited"] for e in operations))
      await session.start()
      await session.execute_phase("preparation")
      await session.execute_phase("refill")
      with (
        patch.object(session.lh, "aspirate", wraps=session.lh.aspirate) as aspirate,
        patch.object(session.lh, "dispense", wraps=session.lh.dispense) as dispense,
      ):
        await session.execute_phase("plate")
      self.assertEqual(aspirate.await_count, 96)
      self.assertEqual(dispense.await_count, 96)
      for event, aspiration, dispensation in zip(
        expected, aspirate.call_args_list, dispense.call_args_list
      ):
        self.assertEqual(aspiration.kwargs["liquid_height"], [event["aspiration_height"]])
        self.assertEqual(
          aspiration.kwargs["surface_following_distance"], [event["surface_following_distance"]]
        )
        self.assertEqual(dispensation.kwargs["liquid_height"], [2.5])
        self.assertEqual(dispensation.kwargs["surface_following_distance"], [0.0])
        self.assertEqual(dispensation.kwargs["immersion_depth"], [0.0])
      self.assertEqual(session.volumes(), baseline.report["source_residuals"])
      self.assertEqual(session.completed_phases, ["preparation", "refill", "plate"])
      await session.stop()

  async def test_preparation_transfer_lowering_reaches_runtime(self) -> None:
    """Lower every buffer/DNA withdrawal, retaining clearance, mixing and plate geometry."""
    settings = Settings(
      preparation_aspiration_lowering=4,
      plate_aspiration_lowering=4,
      plate_dispense_height=2.5,
    )
    notebook_settings = settings_from_notebook(
      Path(__file__).parent / "Exp1034--Filtermax-LOD.ipynb"
    )
    self.assertEqual(notebook_settings.preparation_aspiration_lowering, 4)
    changed_fields = {
      "aspiration_height",
      "aspiration_end_height",
      "surface_following_distance",
      "tip_penetration_at_end",
      "aspiration_lowering_requested",
      "aspiration_start_lowering_applied",
      "aspiration_end_lowering_applied",
      "aspiration_height_limited",
    }
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      session = ExperimentSession(settings, Path(directory) / "lowered")
      baseline = ExperimentSession(
        replace(settings, preparation_aspiration_lowering=0), Path(directory) / "baseline"
      )
      transfers = []
      for event, original in zip(session.plan, baseline.plan):
        if event["phase"] != "preparation" or event["action"] != "transfer":
          self.assertEqual(event, original)
          continue
        transfers.append(event)
        self.assertIn(event["destination"], session.setup.tubes)
        self.assertEqual(event["volume"], 750)
        self.assertEqual(event["aspiration_lowering_requested"], 4)
        self.assertEqual(event["aspiration_minimum_height"], 4.3)
        for key in ("aspiration_height", "aspiration_end_height"):
          self.assertAlmostEqual(event[key], max(4.3, original[key] - 4))
        self.assertAlmostEqual(
          event["surface_following_distance"],
          event["aspiration_height"] - event["aspiration_end_height"],
        )
        self.assertGreaterEqual(event["surface_following_distance"], 0)
        self.assertLessEqual(
          event["tip_penetration_at_end"] + settings.tip_rim_margin,
          tip_insertion_limit(1000, session.setup.audit),
        )
        self.assertEqual(
          {k: v for k, v in event.items() if k not in changed_fields},
          {k: v for k, v in original.items() if k not in changed_fields},
        )
      self.assertEqual(len(transfers), 20)
      self.assertEqual(sum(e["buffer_prefill"] for e in transfers), 10)
      limited = [e for e in transfers if e["aspiration_height_limited"]]
      self.assertEqual(
        [(e["source"], e["destination"]) for e in limited], [("A32", "A25"), ("A31", "A25")]
      )
      self.assertTrue(all(e["aspiration_end_height"] == 4.3 for e in limited))
      a28 = [e for e in transfers if e["source"] == "A28"]
      self.assertEqual([e["destination"] for e in a28], ["A27", "A27"])
      for event, (start, end) in zip(
        a28, [(32.579970533, 23.767230533), (23.767230533, 14.954490533)]
      ):
        self.assertAlmostEqual(event["aspiration_height"], start)
        self.assertAlmostEqual(event["aspiration_end_height"], end)
      self.assertEqual(session.report, baseline.report)
      self.assertEqual(session.ledger, baseline.ledger)
      await session.start()
      operations = [
        e
        for e in session.plan
        if e["phase"] == "preparation" and e["action"] in ("transfer", "mix_cycle")
      ]
      with (
        patch.object(session.lh, "aspirate", wraps=session.lh.aspirate) as aspirate,
        patch.object(session.lh, "dispense", wraps=session.lh.dispense) as dispense,
      ):
        await session.execute_phase("preparation")
      self.assertEqual(aspirate.await_count, 68)
      self.assertEqual(dispense.await_count, 68)
      for event, aspiration, dispensation in zip(
        operations, aspirate.call_args_list, dispense.call_args_list
      ):
        self.assertEqual(aspiration.kwargs["liquid_height"], [event["aspiration_height"]])
        self.assertEqual(
          aspiration.kwargs["surface_following_distance"], [event["surface_following_distance"]]
        )
        self.assertEqual(aspiration.kwargs["immersion_depth"], [0])
        self.assertEqual(dispensation.kwargs["liquid_height"], [event["dispense_height"]])
        flow = 400 if event["action"] == "mix_cycle" else 250
        self.assertEqual(aspiration.kwargs["flow_rates"], [flow])
        self.assertEqual(dispensation.kwargs["flow_rates"], [flow])
      await session.execute_phase("refill")
      with patch.object(session.lh, "dispense", wraps=session.lh.dispense) as dispense:
        await session.execute_phase("plate")
      self.assertEqual(dispense.await_count, 96)
      self.assertTrue(all(c.kwargs["liquid_height"] == [2.5] for c in dispense.call_args_list))
      self.assertEqual(session.volumes(), baseline.report["source_residuals"])
      self.assertEqual(session.report["tip_counts"], {1000: 5, 300: 7})
      await session.stop()
      for value in (-1, float("nan"), float("inf")):
        with self.subTest(value=value), self.assertRaisesRegex(ValueError, "Invalid residual"):
          ExperimentSession(
            replace(settings, preparation_aspiration_lowering=value), Path(directory) / "invalid"
          )

  async def test_omitted_preparation_floor_matches_saved_notebook(self) -> None:
    """A Settings cell without the floor argument must still use the corrected path."""
    explicit = settings_from_notebook(Path(__file__).parent / "Exp1034--Filtermax-LOD.ipynb")
    arguments = asdict(explicit)
    del arguments["preparation_minimum_tip_height"]
    fallback = Settings(**arguments)
    self.assertEqual(fallback.preparation_minimum_tip_height, 4.3)
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      session = ExperimentSession(fallback, Path(directory) / "omitted")
      reference = ExperimentSession(explicit, Path(directory) / "explicit")
      self.assertEqual(session.plan, reference.plan)
      self.assertEqual(session.plan[9]["aspiration_end_height"], 4.3)
      self.assertEqual(session.plan[10]["aspiration_end_height"], 4.3)

  async def test_preparation_floor_avoids_rejected_buffer_endpoint(self) -> None:
    """Only the last two buffer endpoints change; DNA lowering and plate heights remain."""
    settings = settings_from_notebook(Path(__file__).parent / "Exp1034--Filtermax-LOD.ipynb")
    self.assertEqual(settings.preparation_minimum_tip_height, 4.3)
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      session = ExperimentSession(settings, Path(directory) / "adjusted")
      previous = ExperimentSession(
        replace(settings, preparation_minimum_tip_height=1), Path(directory) / "previous"
      )
      changed = []
      for event, original in zip(session.plan, previous.plan):
        if event["phase"] != "preparation" or event["action"] != "transfer":
          self.assertEqual(event, original)
          continue
        self.assertEqual(event["aspiration_height"], original["aspiration_height"])
        self.assertEqual(event["dispense_height"], original["dispense_height"])
        self.assertGreaterEqual(event["aspiration_end_height"], 4.3)
        if event["aspiration_end_height"] != original["aspiration_end_height"]:
          changed.append((event["source"], event["destination"]))
          self.assertEqual(original["aspiration_end_height"], 1)
          self.assertEqual(event["aspiration_end_height"], 4.3)
          bottom = session.setup.audit["tube_coordinates"][event["source"]]["cavity_bottom_z"]
          firmware_start = round((bottom + event["aspiration_height"]) * 10)
          firmware_follow = round(event["surface_following_distance"] * 10)
          self.assertEqual(firmware_start - firmware_follow, 1159)
        if not event["buffer_prefill"]:
          self.assertEqual(event["aspiration_start_lowering_applied"], 4)
          self.assertEqual(event["aspiration_end_lowering_applied"], 4)
      self.assertEqual(changed, [("A32", "A25"), ("A31", "A25")])
      self.assertEqual(session.report, previous.report)
      for invalid in (-1, 0, float("nan"), float("inf"), 10):
        with self.subTest(invalid=invalid), self.assertRaises(ValueError):
          ExperimentSession(
            replace(settings, preparation_minimum_tip_height=invalid),
            Path(directory) / "invalid",
          )

  async def test_prepared_sources_fill_empty_plate_without_preparation(self) -> None:
    """A fresh plate-only run honors actual A25 totals without inventing prior phases."""
    for a25_volume in (1500, 3000):
      with self.subTest(a25_volume=a25_volume), tempfile.TemporaryDirectory() as directory:
        volumes = {
          "A31": 4000,
          "A30": 1500,
          "A29": 1500,
          "A28": 1500,
          "A27": 1500,
          "A26": 1500,
          "A25": a25_volume,
        }
        with patch.object(
          STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")
        ):
          session = ExperimentSession(
            Settings(), Path(directory) / "plate-only", prepared_source_volumes=volumes
          )
          volumes["A25"] = 999  # Caller changes cannot alter the copied plan or starting trackers.
          self.assertEqual(session.report["tip_counts"], {300: 7})
          self.assertEqual({e["phase"] for e in session.plan}, {"plate"})
          a25_first = next(e for e in session.plan if e.get("source") == "A25")
          self.assertEqual(a25_first["source_before"], a25_volume)
          self.assertAlmostEqual(
            a25_first["aspiration_height"],
            session.setup.tubes["A25"].compute_height_from_volume(a25_volume - 25) - 2,
          )
          await session.start()
          self.assertEqual(session.completed_phases, [])
          self.assertEqual(session.volumes()["A25"], a25_volume)
          self.assertEqual(session.metadata["a30_total_volume"], 1500)
          self.assertEqual(session.metadata["entry_phase"], "plate")
          with self.assertRaisesRegex(ValueError, "next planned phase"):
            await session.execute_phase("preparation")
          await session.execute_phase("plate")
          self.assertEqual(session.completed_phases, ["plate"])
          self.assertEqual(session.volumes()["A25"], a25_volume - 1275)
          self.assertEqual(session.volumes()["A31"], 1450)
          self.assertEqual(len(session.metadata["per_well_dispensing_timestamps"]), 96)
          self.assertEqual(session.metadata["refill_acknowledgements"], [])
          self.assertNotIn("dilution_completed_at", session.metadata)
          for record in notebook_protocol.manifest():
            self.assertEqual(
              session.setup.plate.get_item(record.well).tracker.get_used_volume(),
              record.final_volume,
            )
          self.assertEqual(sum(s.has_tip() for s in session.setup.racks[1000].get_all_items()), 96)
          self.assertEqual(sum(s.has_tip() for s in session.setup.racks[300].get_all_items()), 89)
          with self.assertRaisesRegex(ValueError, "already completed"):
            await session.execute_phase("plate")
          await session.stop()
          with self.assertRaisesRegex(ValueError, "prepared_source_volumes"):
            await session.execute_phase("plate")

  async def test_lowered_plate_heights_reach_runtime_with_clearance(self) -> None:
    """Apply both requested offsets while stopping surface following at the reach floor."""
    volumes = {"A31": 4000, **{f"A{i}": 1500 for i in range(26, 31)}, "A25": 3000}
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      baseline = ExperimentSession(
        Settings(), Path(directory) / "baseline", prepared_source_volumes=volumes
      )
      session = ExperimentSession(
        Settings(plate_aspiration_lowering=4, plate_dispense_lowering=3),
        Path(directory) / "lowered",
        prepared_source_volumes=volumes,
      )
      operations = [e for e in session.plan if e["action"] == "transfer"]
      originals = [e for e in baseline.plan if e["action"] == "transfer"]
      self.assertEqual(len(operations), 96)
      for event, original in zip(operations, originals):
        self.assertAlmostEqual(event["aspiration_minimum_height"], 3)
        for key in ("aspiration_height", "aspiration_end_height"):
          self.assertAlmostEqual(event[key], max(3, original[key] - 4))
        self.assertAlmostEqual(
          event["surface_following_distance"],
          event["aspiration_height"] - event["aspiration_end_height"],
        )
        self.assertGreaterEqual(event["surface_following_distance"], 0)
        self.assertLessEqual(event["tip_penetration_at_end"] + 2, 73)
        self.assertAlmostEqual(event["dispense_height"], original["dispense_height"] - 3)
        self.assertAlmostEqual(event["dispense_height"], event["destination_final_surface"])
        self.assertGreaterEqual(event["dispense_height"], session.settings.minimum_tip_height)
        self.assertFalse(event["noncontact"])
        self.assertEqual(event["volume"], original["volume"])
        self.assertEqual(event["source_before"], original["source_before"])
        self.assertEqual(event["source_after"], original["source_after"])
        self.assertTrue(event["dispense_jet"])
        self.assertTrue(event["dispense_blow_out"])
      self.assertTrue(any(e["aspiration_height_limited"] for e in operations))
      self.assertTrue(any(e["surface_following_distance"] == 0 for e in operations))
      self.assertTrue(any(not e["aspiration_height_limited"] for e in operations))
      await session.start()
      with (
        patch.object(session.lh, "aspirate", wraps=session.lh.aspirate) as aspirate,
        patch.object(session.lh, "dispense", wraps=session.lh.dispense) as dispense,
      ):
        await session.execute_phase("plate")
      for event, aspiration, dispensation in zip(
        operations, aspirate.call_args_list, dispense.call_args_list
      ):
        self.assertEqual(aspiration.kwargs["liquid_height"], [event["aspiration_height"]])
        self.assertEqual(
          aspiration.kwargs["surface_following_distance"], [event["surface_following_distance"]]
        )
        self.assertEqual(dispensation.kwargs["liquid_height"], [event["dispense_height"]])
        self.assertEqual(aspiration.kwargs["flow_rates"], [100])
        self.assertEqual(dispensation.kwargs["flow_rates"], [120])
        self.assertEqual(dispensation.kwargs["jet"], [True])
        self.assertEqual(dispensation.kwargs["blow_out"], [True])
      self.assertEqual(aspirate.await_count, 96)
      self.assertEqual(dispense.await_count, 96)
      self.assertEqual(session.volumes(), baseline.report["source_residuals"])
      self.assertEqual(session.report["tip_counts"], {300: 7})
      for record in notebook_protocol.manifest():
        self.assertEqual(
          session.setup.plate.get_item(record.well).tracker.get_used_volume(), record.final_volume
        )
      await session.stop()

  async def test_fixed_plate_dispense_height_for_every_volume(self) -> None:
    """Every plate call uses 2.5 mm from the inner bottom, with no residual volume offset."""
    volumes = {"A31": 4000, **{f"A{i}": 1500 for i in range(26, 31)}, "A25": 3000}
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      settings = Settings(
        plate_aspiration_lowering=4, plate_dispense_lowering=3, plate_dispense_height=2.5
      )
      session = ExperimentSession(
        settings, Path(directory) / "fixed", prepared_source_volumes=volumes
      )
      original = ExperimentSession(
        replace(settings, plate_dispense_height=None),
        Path(directory) / "original",
        prepared_source_volumes=volumes,
      )
      operations = [e for e in session.plan if e["action"] == "transfer"]
      old_operations = [e for e in original.plan if e["action"] == "transfer"]
      self.assertEqual({e["volume"] for e in operations}, {50, 75, 100, 200})
      self.assertEqual({e["dispense_height"] for e in operations}, {2.5})
      for event, old_event in zip(operations, old_operations):
        self.assertEqual(event["aspiration_height"], old_event["aspiration_height"])
        self.assertEqual(event["aspiration_end_height"], old_event["aspiration_end_height"])
        self.assertEqual(event["dispense_lowering_applied"], 0)
        self.assertEqual(event["dispense_height_mode"], "fixed_above_inner_bottom")
        well = session.setup.plate.get_item(event["destination"])
        bottom = well.get_location_wrt(session.setup.deck).z + well.material_z_thickness
        self.assertEqual(event["destination_inner_bottom_z"], bottom)
        self.assertEqual(event["dispense_tip_z"], bottom + 2.5)
      await session.start()
      with patch.object(session.lh, "dispense", wraps=session.lh.dispense) as dispense:
        await session.execute_phase("plate")
      self.assertEqual(dispense.await_count, 96)
      for event, call in zip(operations, dispense.call_args_list):
        self.assertEqual(call.kwargs["liquid_height"], [2.5])
        self.assertEqual(call.kwargs["surface_following_distance"], [0.0])
        self.assertEqual(call.kwargs["immersion_depth"], [0.0])
        self.assertEqual(call.kwargs["vols"], [event["volume"]])
        self.assertEqual(call.kwargs["flow_rates"], [120])
      self.assertEqual(session.volumes(), original.report["source_residuals"])
      self.assertEqual(session.report["tip_counts"], {300: 7})
      self.assertEqual(session.completed_phases, ["plate"])
      for record in notebook_protocol.manifest():
        self.assertEqual(
          session.setup.plate.get_item(record.well).tracker.get_used_volume(), record.final_volume
        )
      await session.stop()
      for invalid in (0.5, -1, float("nan"), float("inf")):
        with self.subTest(height=invalid), self.assertRaisesRegex(ValueError, "bottom clearance"):
          ExperimentSession(
            replace(settings, plate_dispense_height=invalid), Path(directory) / "invalid"
          )
      with self.assertRaisesRegex(ValueError, "Dispense height exits"):
        ExperimentSession(replace(settings, plate_dispense_height=20), Path(directory) / "above")

  async def test_plate_height_offsets_leave_preparation_unchanged(self) -> None:
    """Plate-only corrections never move preparation transfers or mixing positions."""
    with tempfile.TemporaryDirectory() as directory:
      baseline = ExperimentSession(Settings(), Path(directory) / "baseline")
      adjusted = ExperimentSession(
        Settings(plate_aspiration_lowering=4, plate_dispense_lowering=3),
        Path(directory) / "adjusted",
      )
      self.assertEqual(
        [e for e in baseline.plan if e["phase"] != "plate"],
        [e for e in adjusted.plan if e["phase"] != "plate"],
      )
      fixed = ExperimentSession(Settings(plate_dispense_height=2.5), Path(directory) / "fixed")
      self.assertEqual(
        [e for e in baseline.plan if e["phase"] != "plate"],
        [e for e in fixed.plan if e["phase"] != "plate"],
      )
      for field in ("plate_aspiration_lowering", "plate_dispense_lowering"):
        for value in (-1, float("nan"), float("inf")):
          with (
            self.subTest(field=field, value=value),
            self.assertRaisesRegex(ValueError, "Invalid residual/loss/time setting"),
          ):
            ExperimentSession(replace(Settings(), **{field: value}), Path(directory) / "invalid")
      with self.assertRaisesRegex(ValueError, "bottom clearance"):
        ExperimentSession(Settings(plate_dispense_lowering=4), Path(directory) / "too-low")

  async def test_prepared_sources_validate_totals_before_hardware(self) -> None:
    """Missing, insufficient or invalid actual totals must not become assumed liquid levels."""
    volumes = {
      "A31": 4000,
      "A30": 1500,
      "A29": 1500,
      "A28": 1500,
      "A27": 1500,
      "A26": 1500,
      "A25": 1500,
    }
    candidates = [
      ({k: v for k, v in volumes.items() if k != "A25"}, "exactly A31"),
      ({**volumes, "A25": None}, "current total volume in A25"),
      ({**volumes, "A25": float("nan")}, "current total volume in A25"),
      ({**volumes, "A25": 6000}, "current total volume in A25"),
      ({**volumes, "A25": 1499}, "A25 needs at least 1500"),
      ({**volumes, "A31": 1500}, "A31 needs at least 2775"),
    ]
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      output = Path(directory) / "invalid"
      for candidate, message in candidates:
        with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
          ExperimentSession(Settings(), output, prepared_source_volumes=candidate)
        self.assertFalse(output.exists())

  async def test_plate_only_approval_and_existing_records(self) -> None:
    """Approval describes actual prepared totals, and old run directories remain protected."""
    volumes = {
      "A31": 4000,
      "A30": 1500,
      "A29": 1500,
      "A28": 1500,
      "A27": 1500,
      "A26": 1500,
      "A25": 1500,
    }
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      confirm = Mock(return_value="NO")
      session = ExperimentSession(
        Settings(simulate=False, measured_a30_concentration=16, operator="fixture"),
        Path(directory) / "plate-only",
        confirm=confirm,
        prepared_source_volumes=volumes,
      )
      with self.assertRaisesRegex(ValueError, "Not approved"):
        await session.start()
      prompt = confirm.call_args.args[0]
      self.assertIn("A25: 1500 uL at 0.5 ng/uL", prompt)
      self.assertIn("fresh 300 uL tips in A1-G1", prompt)
      self.assertIn("Only plate filling", prompt)
      self.assertFalse(session.output.exists())
      session.output.mkdir()
      previous = session.output / "previous-record.txt"
      previous.write_text("preserve")
      with self.assertRaisesRegex(ValueError, "fresh output directory"):
        await session.start()
      self.assertEqual(previous.read_text(), "preserve")

  async def test_full_liquid_handler_simulation(self) -> None:
    """All 96 real PLR well trackers and mixed tip inventories must match the plan."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      self.assertEqual(session.report["tip_counts"], {1000: 5, 300: 7})
      with patch.object(
        STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")
      ):
        await session.start()
      self.assertEqual(session.metadata["per_well_dispensing_timestamps"], {})
      self.assertNotIn("plate_resume", session.metadata)
      self.assertNotIn("a25_bottom_300ul_test", session.metadata)
      for phase in ("preparation", "refill", "plate"):
        await session.execute_phase(phase)
      self.assertEqual(session.volumes(), session.report["source_residuals"])
      self.assertEqual(session.volumes()["A25"], 1725)
      self.assertEqual(len(session.metadata["per_well_dispensing_timestamps"]), 96)
      for size, count in ((1000, 5), (300, 7)):
        rack = session.setup.racks[size]
        remaining = sum(spot.has_tip() for spot in rack.get_all_items())
        self.assertEqual(remaining, 96 - count)
      with self.assertRaisesRegex(ValueError, "already completed"):
        await session.execute_phase("plate")
      await session.stop()

  async def test_no_approval_no_backend(self) -> None:
    """Declined authorization prevents even STARBackend construction."""
    settings = Settings(
      simulate=False,
      measured_a30_concentration=16,
      operator="fixture",
      method_validation_reference="test-only",
    )
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(settings, Path(directory) / "run", confirm=lambda _: "NO")
      with patch(
        "notebook_protocol.STARBackend", side_effect=AssertionError("Hardware constructed")
      ):
        with self.assertRaisesRegex(ValueError, "Not approved"):
          await session.start()
      self.assertFalse(session.output.exists())

  async def test_missing_measurement_is_pre_io_error(self) -> None:
    """A planning concentration cannot be silently promoted to a measured stock."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(simulate=False), Path(directory) / "run")
      with patch(
        "notebook_protocol.STARBackend", side_effect=AssertionError("Hardware constructed")
      ):
        with self.assertRaisesRegex(ValueError, "measured_a30") as error:
          await session.start()
      self.assertIn("measured_a30_concentration=None", str(error.exception))
      self.assertIn("a30_concentration=16.0 is only the design target", str(error.exception))
      self.assertIn("rerun the full Settings cell to recreate session", str(error.exception))
      self.assertIn("operator is blank", str(error.exception))
      self.assertNotIn("method_validation_reference is blank", str(error.exception))
      self.assertIsNone(session.lh)
      self.assertFalse(session.output.exists())

  async def test_first_run_without_method_record_requires_approval(self) -> None:
    """Blank method notes permit first-run planning but never replace hardware approval."""
    for reference in ("", " \n "):
      with self.subTest(reference=reference), tempfile.TemporaryDirectory() as directory:
        confirm = Mock(return_value="NO")
        settings = Settings(
          simulate=False,
          measured_a30_concentration=16,
          operator="fixture",
          method_validation_reference=reference,
        )
        session = ExperimentSession(settings, Path(directory) / "run", confirm=confirm)
        with patch(
          "notebook_protocol.STARBackend", side_effect=AssertionError("Hardware constructed")
        ):
          with self.assertRaisesRegex(ValueError, "Not approved"):
            await session.start()
        confirm.assert_called_once()
        self.assertIn("no prior method validation supplied", confirm.call_args.args[0])
        self.assertIsNone(session.lh)
        self.assertFalse(session.output.exists())

  async def test_method_record_preserves_notes_without_inventing_validation(self) -> None:
    """Persist absent or supplied validation notes accurately in the run record."""
    for reference in ("", " \n ", "fixture notebook page 2"):
      with self.subTest(reference=reference), tempfile.TemporaryDirectory() as directory:
        session = ExperimentSession(
          Settings(method_validation_reference=reference), Path(directory) / "run"
        )
        with patch(
          "notebook_protocol.STARBackend", side_effect=AssertionError("Hardware constructed")
        ):
          await session.start()
        record = json.loads((session.output / "preparation_record.json").read_text())
        self.assertEqual(record["method_validation"]["reference"], reference.strip() or None)
        self.assertEqual(
          record["method_validation"]["status"],
          "reference_provided" if reference.strip() else "not_provided",
        )
        await session.stop()

  async def test_tip_reach_and_noncontact_heights(self) -> None:
    """Tip insertion stays checked; only clean-buffer prefills end below the final surface."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      for event in session.plan:
        if event["action"] in ("transfer", "mix_cycle"):
          shaft = tip_insertion_limit(event["tip_size"], session.setup.audit)
          self.assertLessEqual(
            event["tip_penetration_at_end"] + session.settings.tip_rim_margin, shaft
          )
        if event["action"] == "transfer":
          if event["buffer_prefill"]:
            self.assertAlmostEqual(event["dispense_height"], event["destination_final_surface"] - 1)
          else:
            self.assertGreater(event["dispense_height"], event["destination_final_surface"])
      dna_plate = [
        e
        for e in session.plan
        if e["action"] == "transfer" and e["phase"] == "plate" and e["source"] == "A30"
      ]
      self.assertTrue(all(e["tip_size"] == 300 for e in dna_plate))

  async def test_preparation_default_rates_and_contact_commands(self) -> None:
    """Run the complete preparation offline and check the actual frontend command arguments."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      operations = [
        e
        for e in session.plan
        if e["phase"] == "preparation" and e["action"] in ("transfer", "mix_cycle")
      ]
      with (
        patch.object(session.lh, "aspirate", wraps=session.lh.aspirate) as aspirate,
        patch.object(session.lh, "dispense", wraps=session.lh.dispense) as dispense,
      ):
        await session.execute_phase("preparation")
      self.assertEqual(aspirate.await_count, len(operations))
      self.assertEqual(dispense.await_count, len(operations))
      self.assertEqual(sum(e.get("buffer_prefill", False) for e in operations), 10)
      self.assertEqual(sum(e["action"] == "mix_cycle" for e in operations), 48)
      for event, aspiration, dispensation in zip(
        operations, aspirate.call_args_list, dispense.call_args_list
      ):
        mixing = event["action"] == "mix_cycle"
        prefill = event.get("buffer_prefill", False)
        self.assertEqual(aspiration.kwargs["flow_rates"], [400 if mixing else 250])
        self.assertEqual(dispensation.kwargs["flow_rates"], [400 if mixing else 250])
        self.assertEqual(dispensation.kwargs["jet"], [not (mixing or prefill)])
        self.assertEqual(dispensation.kwargs["blow_out"], [not mixing])
        if prefill:
          height = 11.435508533 if event["destination_after"] == 750 else 20.248248533
          self.assertAlmostEqual(dispensation.kwargs["liquid_height"][0], height)
          self.assertFalse(event["noncontact"])
        for call in (aspiration, dispensation):
          self.assertEqual(call.kwargs["hamilton_liquid_classes"], [None])
          self.assertEqual(call.kwargs["disable_volume_correction"], [True])
          self.assertEqual(call.kwargs["transport_air_volume"], [0])
          self.assertEqual(call.kwargs["blow_out_air_volume"], [0])
          self.assertEqual(call.kwargs["vols"], [750])
        self.assertEqual(aspiration.kwargs["pre_wetting_volume"], [0])
      self.assertEqual(
        session.volumes(),
        {"A32": 250, "A31": 250, **{f"A{i}": 1500 for i in range(26, 31)}, "A25": 3000},
      )
      for event in session.plan:
        if event["phase"] == "plate" and event["action"] == "transfer":
          self.assertEqual((event["aspirate_flow_rate"], event["dispense_flow_rate"]), (100, 120))
      await session.stop()

  async def test_preparation_flow_overrides_are_scoped_and_validated(self) -> None:
    """Optional numeric preparation rates must not alter mixing or plate dispensing."""
    with tempfile.TemporaryDirectory() as directory:
      settings = Settings(preparation_aspirate_flow=325, preparation_dispense_flow=350)
      session = ExperimentSession(settings, Path(directory) / "run")
      for event in session.plan:
        if event["action"] not in ("transfer", "mix_cycle"):
          continue
        expected = (
          (400, 400)
          if event["action"] == "mix_cycle"
          else (325, 350)
          if event["phase"] == "preparation"
          else (100, 120)
        )
        self.assertEqual((event["aspirate_flow_rate"], event["dispense_flow_rate"]), expected)
      for invalid in (0, -1, float("nan"), float("inf")):
        for settings in (
          Settings(preparation_aspirate_flow=invalid),
          Settings(preparation_dispense_flow=invalid),
        ):
          with (
            self.subTest(settings=settings),
            self.assertRaisesRegex(ValueError, "Preparation flow"),
          ):
            ExperimentSession(settings, Path(directory) / "invalid")
      with patch("notebook_protocol.get_star_liquid_class", return_value=None):
        with self.assertRaisesRegex(ValueError, "No preparation aspiration default"):
          ExperimentSession(
            Settings(preparation_aspirate_flow=None), Path(directory) / "missing_default"
          )
      with self.assertRaisesRegex(ValueError, "bottom clearance"):
        ExperimentSession(Settings(buffer_prefill_submersion=20), Path(directory) / "too_low")

  async def test_runtime_arguments_match_installed_star_api(self) -> None:
    """Check the real backend signatures, because Chatterbox accepts arbitrary kwargs."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      event = next(e for e in session.plan if e["action"] == "transfer")
      session.lh.aspirate = AsyncMock()
      session.lh.dispense = AsyncMock()
      await session.aspirate(event, mixing=False)
      await session.dispense(event, mixing=False)
      front_end = {"vols", "flow_rates", "liquid_height", "blow_out_air_volume"}
      for mock, api in (
        (session.lh.aspirate, STARBackend.aspirate),
        (session.lh.dispense, STARBackend.dispense),
      ):
        backend_args = {k: v for k, v in mock.call_args.kwargs.items() if k not in front_end}
        inspect.signature(api).bind(None, ops=[], **backend_args)
        self.assertEqual(backend_args["hamilton_liquid_classes"], [None])
        self.assertEqual(backend_args["transport_air_volume"], [0])
      await session.stop()

  async def test_preparation_resume_keeps_completed_steps_and_empty_tip(self) -> None:
    """Recover the live checkpoint without setup, tip changes or replaying filled tubes."""
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      session = await failed_buffer_session(Path(directory) / "run")
      original_plan = session.plan
      before = session.volumes()
      connection = session.lh
      tip = session.lh.head[3].get_tip()
      with (
        patch.object(session.lh, "aspirate", side_effect=AssertionError("Unexpected motion")),
        patch.object(session.lh, "dispense", side_effect=AssertionError("Unexpected motion")),
        patch.object(session.lh, "pick_up_tips", side_effect=AssertionError("Unexpected motion")),
        patch.object(session.lh, "discard_tips", side_effect=AssertionError("Unexpected motion")),
      ):
        await notebook_protocol.resume_preparation(session, no_liquid_aspirated=True)
      self.assertFalse(session.failed)
      self.assertIs(session.lh, connection)
      self.assertIs(session.lh.head[3].get_tip(), tip)
      self.assertEqual(tip.tracker.get_used_volume(), 0)
      self.assertEqual(session.volumes(), before)
      self.assertEqual(session.next_step, 9)
      self.assertEqual(session.settings.mix_cycles, 1)
      self.assertEqual(session.plan[:9], original_plan[:9])
      self.assertEqual(session.plan[9]["aspiration_end_height"], 4.3)
      self.assertEqual(session.plan[10]["aspiration_end_height"], 4.3)
      self.assertEqual(
        json.loads(
          (session.output / "planned_operations.before_preparation_resume.json").read_text()
        ),
        original_plan,
      )
      with self.assertRaises(ValueError):
        await notebook_protocol.resume_preparation(session, no_liquid_aspirated=True)
      with patch.object(session.lh, "aspirate", wraps=session.lh.aspirate) as aspirate:
        await session.execute_phase("preparation")
      self.assertEqual(aspirate.await_count, 18)
      self.assertEqual(aspirate.call_args_list[0].args[0], [session.setup.tubes["A32"]])
      self.assertEqual(aspirate.call_args_list[0].kwargs["vols"], [750])
      self.assertAlmostEqual(
        aspirate.call_args_list[0].kwargs["surface_following_distance"][0], 4.779330533
      )
      await session.execute_phase("refill")
      await session.execute_phase("plate")
      completed = [
        e["step"]
        for e in session.events
        if e.get("action") == "plan_step" and e.get("status") == "completed"
      ]
      self.assertEqual(completed, list(range(1, len(session.plan) + 1)))
      self.assertEqual(session.volumes(), session.report["source_residuals"])
      self.assertEqual(len(session.metadata["per_well_dispensing_timestamps"]), 96)
      await session.stop()

  async def test_preparation_resume_rejects_uncertain_or_changed_state(self) -> None:
    """A no-liquid acknowledgement cannot override partial aspiration or changed trackers."""
    with (
      tempfile.TemporaryDirectory() as directory,
      patch.object(STARBackend, "__init__", side_effect=AssertionError("Hardware constructed")),
    ):
      session = await failed_buffer_session(Path(directory) / "run")
      for value in (False, "yes", 1):
        with self.subTest(value=value), self.assertRaisesRegex(ValueError, "confirmation"):
          await notebook_protocol.resume_preparation(session, no_liquid_aspirated=value)
      session.setup.tubes["A32"].set_volume(999)
      with self.assertRaisesRegex(ValueError, "Source trackers changed"):
        await notebook_protocol.resume_preparation(session, no_liquid_aspirated=True)
      session.setup.tubes["A32"].set_volume(1000)
      session.lh.head[3].get_tip().tracker.set_volume(1)
      with self.assertRaisesRegex(ValueError, "original empty"):
        await notebook_protocol.resume_preparation(session, no_liquid_aspirated=True)
      session.lh.head[3].get_tip().tracker.set_volume(0)
      session.events[-2]["status"] = "completed"
      with self.assertRaisesRegex(ValueError, "no completed aspiration"):
        await notebook_protocol.resume_preparation(session, no_liquid_aspirated=True)
      session.events[-2]["status"] = "attempt"
      session.stopped = True
      with self.assertRaisesRegex(ValueError, "original failed live session"):
        await notebook_protocol.resume_preparation(session, no_liquid_aspirated=True)
      session.stopped = False
      self.assertTrue(session.failed)
      self.assertEqual(session.next_step, 9)
      self.assertEqual(session.plan[9]["aspiration_end_height"], 1)
      self.assertFalse(
        (session.output / "planned_operations.before_preparation_resume.json").exists()
      )
      await session.stop()

  async def test_partial_transfer_latches_failure(self) -> None:
    """A dispense failure after successful aspiration cannot replay the aspiration."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      session.lh.dispense = AsyncMock(side_effect=RuntimeError("fixture dispense failure"))
      with self.assertRaisesRegex(RuntimeError, "fixture dispense failure"):
        await session.execute_phase("preparation")
      self.assertTrue(session.failed)
      self.assertEqual(session.volumes()["A32"], 3250)
      with self.assertRaisesRegex(ValueError, "failed"):
        await session.execute_phase("preparation")
      with self.assertRaises(ValueError):
        await notebook_protocol.resume_refill(session)
      self.assertTrue(session.failed)
      await session.stop()

  async def test_refill_cancel_does_not_change_actual_tracker(self) -> None:
    """A physical refill acknowledgement precedes the change to the PLR tracker."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run", confirm=lambda _: "NO")
      await session.start()
      await session.execute_phase("preparation")
      session.settings = replace(session.settings, simulate=False)
      with self.assertRaisesRegex(ValueError, "Refill cancelled"):
        await session.execute_phase("refill")
      self.assertEqual(session.volumes()["A31"], 250)
      self.assertFalse(session.failed)
      session.confirm = Mock(return_value="refilled")
      await session.execute_phase("refill")
      self.assertEqual(session.volumes()["A31"], 4000)
      self.assertEqual(session.completed_phases, ["preparation", "refill"])
      self.assertEqual(len(session.metadata["refill_acknowledgements"]), 1)
      await session.stop()

  async def test_refill_acknowledgement_accepts_case_and_source(self) -> None:
    """Accept the named checkpoint's short confirmation but reject another tube."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      await session.execute_phase("preparation")
      session.settings = replace(session.settings, simulate=False)
      event = session.plan[session.next_step]
      for response in ("refilled A32", "yes", "", "NO"):
        with self.subTest(response=response):
          session.confirm = Mock(return_value=response)
          with self.assertRaisesRegex(ValueError, "Refill cancelled"):
            await session.execute_phase("refill")
          self.assertFalse(session.failed)
          self.assertEqual(session.volumes()["A31"], 250)
      for response in ("refilled", "REFILLED", "  rEfIlLeD   a31  ", "REFILLED A31"):
        with self.subTest(response=response):
          session.confirm = Mock(return_value=response)
          notebook_protocol.ExperimentSession.acknowledge_refill(session, event)
          self.assertEqual(session.volumes()["A31"], 4000)
          session.setup.tubes["A31"].set_volume(250)
      await session.stop()

  async def test_resume_existing_refill_keeps_connection_and_prepared_samples(self) -> None:
    """Recover the exact input-only failure without replaying preparation or robot commands."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      await session.execute_phase("preparation")
      session.settings = replace(session.settings, simulate=False)
      session.confirm = Mock(side_effect=ValueError("Refill cancelled; no tracking update."))
      with self.assertRaisesRegex(ValueError, "Refill cancelled"):
        await session.execute_phase("refill")
      self.assertTrue(session.failed)
      connection = session.lh
      before = session.volumes()
      next_step = session.next_step
      session.confirm = Mock(return_value="refilled")
      session.execute_event = AsyncMock(side_effect=AssertionError("Robot command replayed"))
      session.setup.tubes["A30"].set_volume(before["A30"] - 1)
      with self.assertRaisesRegex(ValueError, "Tracked sources differ"):
        await notebook_protocol.resume_refill(session)
      self.assertTrue(session.failed)
      session.confirm.assert_not_called()
      session.setup.tubes["A30"].set_volume(before["A30"])
      await notebook_protocol.resume_refill(session)
      session.execute_event.assert_not_called()
      self.assertIs(session.lh, connection)
      self.assertEqual(session.volumes(), {**before, "A31": 4000})
      self.assertEqual(session.next_step, next_step + 1)
      self.assertEqual(session.completed_phases, ["preparation", "refill"])
      self.assertFalse(session.failed)
      with self.assertRaises(ValueError):
        await notebook_protocol.resume_refill(session)
      self.assertEqual(len(session.metadata["refill_acknowledgements"]), 1)
      await session.stop()

  async def test_resume_refill_rejects_unrelated_failure(self) -> None:
    """A refill prompt recovery cannot clear a different failure or changed tracker state."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      await session.execute_phase("preparation")
      session.settings = replace(session.settings, simulate=False)
      session.confirm = Mock(side_effect=RuntimeError("unrelated failure"))
      with self.assertRaisesRegex(RuntimeError, "unrelated failure"):
        await session.execute_phase("refill")
      session.confirm = Mock(return_value="refilled")
      with self.assertRaisesRegex(ValueError, "only the refill acknowledgement"):
        await notebook_protocol.resume_refill(session)
      self.assertTrue(session.failed)
      session.confirm.assert_not_called()
      self.assertEqual(session.volumes()["A31"], 250)
      await session.stop()

  async def test_tip_reuse_update_preserves_partial_plate_and_reconciles_pickup(self) -> None:
    """Keep A9 filled once and continue after either physical outcome of a cancelled pickup."""
    for mounted in (False, True):
      with self.subTest(mounted=mounted), tempfile.TemporaryDirectory() as directory:
        session = ExperimentSession(Settings(), Path(directory) / "run")
        await session.start()
        await session.execute_phase("preparation")
        await session.execute_phase("refill")
        start = session.next_step
        transfer = dict(session.plan[start + 1])
        self.assertEqual(transfer["destination"], "A9")
        transfer["tip_size"] = 300
        legacy_steps = [
          {"action": "pick_up_tip", "tip_size": 300, "rack_spot": "A1"},
          transfer,
          {"action": "discard_tip", "tip_size": 300, "rack_spot": "A1"},
          {"action": "pick_up_tip", "tip_size": 300, "rack_spot": "B1"},
        ]
        for index, event in enumerate(legacy_steps, start + 1):
          event.update({"step": index, "phase": "plate"})
        session.plan = session.plan[:start] + legacy_steps
        for event in legacy_steps[:3]:
          await session.execute_event(event)
          session.next_step += 1
        original_pickup = session.lh.backend.pick_up_tips
        session.lh.backend.pick_up_tips = AsyncMock(side_effect=asyncio.CancelledError())
        session.lh.backend.request_tip_presence = AsyncMock(return_value=[mounted] + [False] * 7)
        with self.assertRaises(asyncio.CancelledError):
          await session.execute_phase("plate")
        self.assertTrue(session.failed)
        self.assertEqual(session.lh.head[0].has_tip, mounted)
        old_volumes = session.volumes()
        a9_timestamp = session.metadata["per_well_dispensing_timestamps"]["A9"]
        session.settings = replace(session.settings, simulate=False)
        if mounted:
          session.lh.backend.request_tip_presence.return_value = [False] * 8
          with self.assertRaisesRegex(ValueError, "Hardware tip presence differs"):
            await notebook_protocol.update_plate_tip_reuse(session)
          self.assertTrue(session.failed)
          self.assertEqual(session.volumes(), old_volumes)
          session.lh.backend.request_tip_presence.return_value = [True] + [False] * 7
        await notebook_protocol.update_plate_tip_reuse(session)
        session.lh.backend.pick_up_tips.assert_awaited_once()
        self.assertFalse(session.failed)
        self.assertEqual(session.volumes(), old_volumes)
        pending = session.plan[session.next_step :]
        self.assertFalse(any(e.get("destination") == "A9" for e in pending))
        self.assertEqual(sum(e["action"] == "pick_up_tip" for e in pending), 7)
        self.assertEqual(pending[0]["action"], "discard_tip" if mounted else "pick_up_tip")
        self.assertEqual(session.report["tip_counts"], {1000: 5, 300: 8 + int(mounted)})
        self.assertEqual([e["step"] for e in session.plan], list(range(1, len(session.plan) + 1)))
        session.lh.backend.pick_up_tips = original_pickup
        session.settings = replace(session.settings, simulate=True)
        await session.execute_phase("plate")
        self.assertEqual(session.volumes(), session.report["source_residuals"])
        self.assertEqual(session.metadata["per_well_dispensing_timestamps"]["A9"], a9_timestamp)
        self.assertEqual(len(session.metadata["per_well_dispensing_timestamps"]), 96)
        await session.stop()

  async def test_tip_reuse_update_does_not_clear_interrupted_dispense(self) -> None:
    """Never use regrouping to replay a transfer whose physical outcome is incomplete."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      await session.execute_phase("preparation")
      await session.execute_phase("refill")
      session.lh.dispense = AsyncMock(side_effect=asyncio.CancelledError())
      with self.assertRaises(asyncio.CancelledError):
        await session.execute_phase("plate")
      old_plan = list(session.plan)
      old_volumes = session.volumes()
      with self.assertRaisesRegex(ValueError, "not an aspiration or dispense"):
        await notebook_protocol.update_plate_tip_reuse(session)
      self.assertTrue(session.failed)
      self.assertEqual(session.plan, old_plan)
      self.assertEqual(session.volumes(), old_volumes)
      await session.stop()

  async def test_unverified_or_excessive_300ul_insertion_rejected(self) -> None:
    """The confirmed reach does not disable insertion limits or generalize to missing evidence."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      session.setup.audit["tip_reach_checks"]["300"]["maximum_penetration"] = 69
      with self.assertRaisesRegex(ValueError, "insertion limit"):
        notebook_protocol.geometry_for_plan(session.setup, session.plan, session.settings)
      session.setup.audit.pop("tip_reach_checks")
      with self.assertRaisesRegex(ValueError, "insertion limit"):
        notebook_protocol.geometry_for_plan(session.setup, session.plan, session.settings)

  async def test_diagnostic_resume_preserves_a9_and_accounts_for_loaded_tip(self) -> None:
    """Offline planning stays paused, then one approved resume fills exactly the remaining wells."""
    with tempfile.TemporaryDirectory() as directory:
      session = ExperimentSession(Settings(), Path(directory) / "run")
      await session.start()
      await session.execute_phase("preparation")
      await session.execute_phase("refill")
      start = session.next_step
      transfer = dict(session.plan[start + 1])
      legacy = [
        {"action": "pick_up_tip", "tip_size": 300, "rack_spot": "A1"},
        transfer,
        {"action": "discard_tip", "tip_size": 300, "rack_spot": "A1"},
        {"action": "pick_up_tip", "tip_size": 300, "rack_spot": "B1"},
      ]
      for index, event in enumerate(legacy, start + 1):
        event.update(step=index, phase="plate")
      session.plan = session.plan[:start] + legacy
      for event in legacy[:3]:
        await session.execute_event(event)
        session.next_step += 1
      original_pickup = session.lh.backend.pick_up_tips
      session.lh.backend.pick_up_tips = AsyncMock(side_effect=asyncio.CancelledError())
      session.lh.backend.request_tip_presence = AsyncMock(return_value=[True] + [False] * 7)
      with self.assertRaises(asyncio.CancelledError):
        await session.execute_phase("plate")
      session.lh.backend.pick_up_tips = original_pickup
      session.attached_tip = (300, "B1")
      await session.aspirate(
        {"source": "A25", "volume": 300, "aspiration_height": 1, "surface_following_distance": 0},
        mixing=False,
      )
      diagnostic = {
        "status": "completed",
        "stage": "aspirate_300ul_and_retract",
        "source": "A25",
        "tip_choice": "attached_B1",
        "volume": 300,
        "tip_volume_after": 300,
        "source_before": 3000,
        "source_after": 2700,
      }
      session.metadata["a25_bottom_300ul_test"] = diagnostic
      session.record({"action": "a25_bottom_300ul_test", **diagnostic})
      a9_timestamp = session.metadata["per_well_dispensing_timestamps"]["A9"]
      session.settings = replace(session.settings, simulate=False)
      session.lh.backend.request_tip_presence = AsyncMock(
        side_effect=AssertionError("Hardware I/O")
      )
      with patch.object(session.lh, "discard_tips", new_callable=AsyncMock) as discard:
        await notebook_protocol.update_plate_tip_reuse(session)
        discard.assert_not_awaited()
      session.lh.backend.request_tip_presence.assert_not_awaited()
      self.assertTrue(session.failed)
      self.assertEqual(session.metadata["plate_resume"]["status"], "prepared_not_authorized")
      self.assertEqual(session.report["source_residuals"]["A25"], 1425)
      self.assertEqual(session.report["tip_counts"], {1000: 5, 300: 9})
      checkpoint = json.loads((session.output / "plate_resume_checkpoint.json").read_text())
      pending = checkpoint["remaining_operations"]
      transfers = [e for e in pending if e["action"] == "transfer"]
      self.assertEqual(len(transfers), 95)
      self.assertEqual(sum(e["volume"] for e in transfers), 10000)
      self.assertTrue(all(e["tip_size"] == 300 for e in transfers))
      self.assertNotIn("A9", [e["destination"] for e in transfers])
      self.assertEqual(
        [e["rack_spot"] for e in pending if e["action"] == "pick_up_tip"],
        ["C1", "D1", "E1", "F1", "G1", "H1", "A2"],
      )
      with self.assertRaisesRegex(ValueError, "failed"):
        await session.execute_phase("plate")
      with self.assertRaisesRegex(ValueError, "already saved"):
        await notebook_protocol.update_plate_tip_reuse(session)
      session.settings = replace(session.settings, simulate=False)
      session.confirm = Mock(return_value="no")
      with patch.object(session.lh, "discard_tips", new_callable=AsyncMock) as discard:
        with self.assertRaisesRegex(ValueError, "Still paused"):
          await notebook_protocol.resume_plate_after_test(session)
        discard.assert_not_awaited()
      self.assertEqual(session.metadata["plate_resume"]["status"], "prepared_not_authorized")
      session.settings = replace(session.settings, simulate=True)
      session.setup.tubes["A25"].set_volume(2699)
      with self.assertRaisesRegex(ValueError, "Sources changed"):
        await notebook_protocol.resume_plate_after_test(session)
      session.setup.tubes["A25"].set_volume(2700)
      with patch.object(session.lh, "discard_tips", wraps=session.lh.discard_tips) as discard:
        await notebook_protocol.resume_plate_after_test(session)
      self.assertTrue(discard.call_args_list[0].kwargs["allow_nonzero_volume"])
      self.assertTrue(
        all(not call.kwargs["allow_nonzero_volume"] for call in discard.call_args_list[1:])
      )
      self.assertEqual(session.volumes(), session.report["source_residuals"])
      self.assertEqual(session.volumes()["A25"], 1425)
      self.assertEqual(len(session.metadata["per_well_dispensing_timestamps"]), 96)
      self.assertEqual(session.metadata["per_well_dispensing_timestamps"]["A9"], a9_timestamp)
      self.assertFalse(session.lh.head[0].has_tip)
      self.assertEqual(session.metadata["plate_resume"]["status"], "completed")
      with self.assertRaisesRegex(ValueError, "No unstarted resume"):
        await notebook_protocol.resume_plate_after_test(session)
      await session.stop()


if __name__ == "__main__":
  unittest.main()
