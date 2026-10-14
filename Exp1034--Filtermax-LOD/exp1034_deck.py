"""Corrected Exp1034 deck and source geometry, without opening instrument I/O."""

from dataclasses import asdict, dataclass

from pylabrobot.resources import (
  Coordinate,
  Plate,
  STARLetDeck,
  TipRack,
  Tube,
  cor_96_wellplate_360uL_Fb,
  hamilton_96_tiprack_50uL,
  hamilton_96_tiprack_300uL,
  hamilton_96_tiprack_1000uL_filter,
  hamilton_tip_carrier_L5,
)
from pylabrobot.resources.hamilton import (
  HamiltonSTARDeck,
  hamilton_mfx_carrier_L5_base,
  hamilton_mfx_plateholder_DWP_metal_tapped,
  hamilton_tube_carrier_32_a00,
)

from mckesson_tubes import McKesson_tube_5mL_Rb

TUBE_DEFINITION_COMMIT = "d62db0a42"
TUBE_DEFINITION_BRANCH = "add-12x75mm-test-tube"
TUBE_SITES = {"A32": 0, "A31": 1, "A30": 2, "A29": 3, "A28": 4, "A27": 5, "A26": 6, "A25": 7}


@dataclass
class DeckSetup:
  """Resources shared by planning, simulation and the approved hardware run."""

  deck: HamiltonSTARDeck
  plate: Plate
  tubes: dict[str, Tube]
  racks: dict[int, TipRack]
  audit: dict


def build_deck() -> DeckSetup:
  """Apply the user's corrected layout and Exp1031's explicit tube mounting."""
  deck = STARLetDeck()
  racks = {
    1000: hamilton_96_tiprack_1000uL_filter("tips_1000"),
    300: hamilton_96_tiprack_300uL("tips_300"),
    50: hamilton_96_tiprack_50uL("tips_50"),
  }
  tip_carrier = hamilton_tip_carrier_L5("TIP_CAR_480_A00")
  for index, rack in enumerate(racks.values()):
    tip_carrier[index].assign_child_resource(rack)
  deck.assign_child_resource(tip_carrier, track=25)
  carrier = hamilton_tube_carrier_32_a00("car_32_A00")
  tubes = {}
  for label, site in TUBE_SITES.items():
    tube = McKesson_tube_5mL_Rb(label)
    # The factory has no base thickness parameter. Use the measured value in Exp1031.
    tube._material_z_thickness = 1.0
    holder = carrier[site]
    bottom = carrier.get_size_z() + 25.6 - tube.get_size_z()
    holder.assign_child_resource(
      tube,
      location=Coordinate(
        (holder.get_size_x() - tube.get_size_x()) / 2,
        (holder.get_size_y() - tube.get_size_y()) / 2,
        bottom - holder.location.z,
      ),
    )
    tubes[label] = tube
  deck.assign_child_resource(carrier, track=19)
  modules = {i: hamilton_mfx_plateholder_DWP_metal_tapped(f"dwp_module_{i}") for i in range(5)}
  plate_carrier = hamilton_mfx_carrier_L5_base("plate_carrier", modules=modules)
  plate = cor_96_wellplate_360uL_Fb("Corning_3635_UV_plate")
  modules[0].assign_child_resource(plate)
  deck.assign_child_resource(plate_carrier, track=7)
  tips = {}
  for size, rack in racks.items():
    tip = rack.get_item("A1").get_tip()
    tips[str(size)] = {
      "nominal": tip.nominal_volume,
      "maximum": tip.maximal_volume,
      "filtered": tip.has_filter,
      "available": rack.num_items,
      "constructor_model": rack.model,
      "length": tip.get_size_z(),
      "fitting_depth": tip.fitting_depth,
      "exposed_length": tip.get_size_z() - max(tip.fitting_depth, tip.collar_height),
    }
  a1 = plate.get_item("A1").get_location_wrt(plate)
  h12 = plate.get_item("H12").get_location_wrt(plate)
  if not a1.x < h12.x or not a1.y > h12.y:
    raise ValueError("Unexpected plate orientation.")
  tube_coordinates = {}
  for label, tube in tubes.items():
    location = tube.get_location_wrt(deck)
    tube_coordinates[label] = {
      "site": TUBE_SITES[label],
      "bottom_origin": asdict(location),
      "center_x": location.x + tube.get_size_x() / 2,
      "center_y": location.y + tube.get_size_y() / 2,
      "cavity_bottom_z": location.z + tube.material_z_thickness,
      "rim_z": location.z + tube.get_size_z(),
    }
  audit = {
    "tips": tips,
    "tip_reach_checks": {
      "300": {
        "maximum_penetration": 73.0,
        "reference": "physical-20261004T001957Z-fab79d55/a25_bottom_300ul_test.py",
        "evidence": (
          "A25: 300 uL B1 tip reached 1 mm above the inner bottom, aspirated 300 uL "
          "and retracted; operator confirmed success. Applies to this tube and deck geometry."
        ),
      }
    },
    "deck_tracks": deck.num_tracks,
    "tip_carrier_track": 25,
    "tip_carrier_location": asdict(tip_carrier.location),
    "tube_carrier_track": 19,
    "plate_carrier_track": 7,
    "plate_carrier_model": plate_carrier.model,
    "plate_modules": {
      str(i): {"model": m.model, "location": asdict(m.location)} for i, m in modules.items()
    },
    "plate_reference": {
      "model": plate.model,
      "physical_catalog": "3635",
      "size": [plate.get_size_x(), plate.get_size_y(), plate.get_size_z()],
      "A1_lower_left": asdict(a1),
      "H12_lower_left": asdict(h12),
      "orientation": "landscape; A1 rear-left, H12 front-right on STAR deck",
      "physical_3635_correspondence_verified": True,
      "correspondence_evidence": "User confirmed 3635 has the same dimensions as the supplied resource.",
    },
    "tube_coordinates": tube_coordinates,
    "tube_geometry": {
      "outer_height": 75,
      "material_thickness": 1,
      "inner_height": 74,
      "rim_above_carrier": 25.6,
      "height_per_volume": 0.01175032,
      "height_intercept": 3.622768533,
      "definition_branch": TUBE_DEFINITION_BRANCH,
      "definition_commit": TUBE_DEFINITION_COMMIT,
      "mount_evidence": "Exp1031 endotoxin_assay_w_12x75mm_tubes.ipynb",
    },
    "placement_check": "PASS: no collision bypass or deprecated rails argument",
  }
  return DeckSetup(deck, plate, tubes, racks, audit)


def tube_height(volume: float, audit: dict) -> float:
  """Use the branch's measured linear fit; its nonzero intercept is not dead-volume validation."""
  geometry = audit["tube_geometry"]
  return geometry["height_intercept"] + geometry["height_per_volume"] * volume


def tip_insertion_limit(size: int, audit: dict) -> float:
  """Use an experiment-specific reach check when present, otherwise exposed shaft length."""
  checked = audit.get("tip_reach_checks", {}).get(str(size))
  return (
    checked["maximum_penetration"]
    if checked is not None
    else audit["tips"][str(size)]["exposed_length"]
  )


def choose_plate_tip(
  residual: float, audit: dict, immersion: float = 2, rim_margin: float = 2
) -> int:
  """Require the requested 300 uL plate tip to reach within the checked insertion limit."""
  height = tube_height(residual, audit) - immersion
  penetration = audit["tube_geometry"]["inner_height"] - height
  if penetration + rim_margin <= tip_insertion_limit(300, audit):
    return 300
  raise ValueError(f"300 uL tip exceeds checked reach at the modeled {residual:g} uL residual.")
