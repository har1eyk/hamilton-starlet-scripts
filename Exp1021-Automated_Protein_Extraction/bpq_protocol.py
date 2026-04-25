"""Reusable 96-well Bradford Protein Quantification protocol for PyLabRobot."""

from __future__ import annotations

from typing import Any, Sequence

CHANNELS = list(range(8))
PLATE_COLUMN_COUNT = 12
def _plate_col(plate: Any, col: int) -> Any:
    return plate[f"A{col}:H{col}"]


def _tip_col(rack: Any, col: int) -> Any:
    return rack[f"A{col}:H{col}"]


def _trough_channels(trough: Any) -> Any:
    return trough["A1"] * 8


def _well_height(plate: Any, volume: float) -> float:
    well = plate.get_item("A1")
    return max(1.0, well.compute_height_from_volume(max(0.0, volume)))


async def _pick_up_column_tips(lh: Any, rack: Any, col: int) -> None:
    await lh.pick_up_tips(_tip_col(rack, col), use_channels=CHANNELS)


async def _discard_column_tips(lh: Any) -> None:
    await lh.discard_tips(use_channels=CHANNELS)


async def _dispense_to_column(
    lh: Any,
    plate: Any,
    col: int,
    volume: float,
    *,
    reference_volume: float,
    flow_rate: float,
    blow_out: bool,
    **backend_kwargs: Any,
) -> None:
    liquid_height = max(1.0, _well_height(plate, reference_volume) - 1.0)
    await lh.dispense(
        _plate_col(plate, col),
        vols=[volume] * 8,
        use_channels=CHANNELS,
        liquid_height=[liquid_height] * 8,
        flow_rates=[flow_rate] * 8,
        blow_out=[1 if blow_out else 0] * 8,
        settling_time=[1] * 8,
        **backend_kwargs,
    )


async def _aspirate_plate_column(
    lh: Any,
    plate: Any,
    col: int,
    volume: float,
    *,
    liquid_height: float,
    flow_rate: float,
    **backend_kwargs: Any,
) -> None:
    await lh.aspirate(
        _plate_col(plate, col),
        vols=[volume] * 8,
        use_channels=CHANNELS,
        liquid_height=[liquid_height] * 8,
        flow_rates=[flow_rate] * 8,
        **backend_kwargs,
    )


async def _fill_sample_plates_with_reagent(
    lh: Any,
    plates: Sequence[tuple[int, Any]],
    bradford_reagent: Any,
    tiprack_1000: Any,
    *,
    reagent_volume: float,
    reagent_aspirate_volume: float,
) -> None:
    dispenses_per_refill = int(reagent_aspirate_volume // reagent_volume)
    if dispenses_per_refill < 1:
        raise ValueError("reagent_aspirate_volume must be at least one full reagent dispense.")

    total_plate_count = len(plates)
    # Keep reagent batching deterministic. The default STAR liquid classes silently alter
    # aspirate/dispense volumes and add transport air, which breaks exact multi-dispense math.
    reagent_backend_kwargs = {
        "hamilton_liquid_classes": [None] * 8,
        "disable_volume_correction": [True] * 8,
        "transport_air_volume": [0] * 8,
    }
    tips_picked_up = False
    try:
        await _pick_up_column_tips(lh, tiprack_1000, 1)
        tips_picked_up = True

        for plate_sequence_index, (plate_index, plate) in enumerate(plates, start=1):
            col = 1
            while col <= PLATE_COLUMN_COUNT:
                batch_start = col
                batch_end = min(PLATE_COLUMN_COUNT, col + dispenses_per_refill - 1)
                batch_dispense_count = batch_end - batch_start + 1
                batch_aspirate_volume = reagent_volume * batch_dispense_count

                await lh.aspirate(
                    _trough_channels(bradford_reagent),
                    vols=[batch_aspirate_volume] * 8,
                    use_channels=CHANNELS,
                    liquid_height=[2] * 8,
                    flow_rates=[250] * 8,
                    **reagent_backend_kwargs,
                )

                for batch_col in range(batch_start, batch_end + 1):
                    await _dispense_to_column(
                        lh,
                        plate,
                        batch_col,
                        reagent_volume,
                        reference_volume=reagent_volume,
                        flow_rate=200,
                        blow_out=False,
                        **reagent_backend_kwargs,
                    )

                print(
                    f"Bradford reagent: plate {plate_sequence_index}/{total_plate_count} "
                    f"(SamplePlate{plate_index}), "
                    f"columns {batch_start}-{batch_end}."
                )
                col = batch_end + 1
    finally:
        if tips_picked_up:
            await _discard_column_tips(lh)


async def _duplicate_elution_samples_to_plate(
    lh: Any,
    elution_plate: Any,
    destination_plate: Any,
    sample_tiprack_50: Any,
    source_cols: range,
    *,
    sample_aspirate_volume: float,
    sample_dispense_volume: float,
) -> None:
    for offset, source_col in enumerate(source_cols):
        destination_col_1 = offset * 2 + 1
        destination_col_2 = destination_col_1 + 1

        await _pick_up_column_tips(lh, sample_tiprack_50, source_col)
        await _aspirate_plate_column(
            lh,
            elution_plate,
            source_col,
            sample_aspirate_volume,
            liquid_height=1.0,
            flow_rate=35,
        )
        await _dispense_to_column(
            lh,
            destination_plate,
            destination_col_1,
            sample_dispense_volume,
            reference_volume=255,
            flow_rate=30,
            blow_out=False,
        )
        await _dispense_to_column(
            lh,
            destination_plate,
            destination_col_2,
            sample_dispense_volume,
            reference_volume=255,
            flow_rate=30,
            blow_out=False,
        )
        await _discard_column_tips(lh)
        print(
            f"Elution column {source_col} duplicated to assay columns "
            f"{destination_col_1} and {destination_col_2}."
        )


async def _triplicate_standards_to_plate(
    lh: Any,
    destination_plate: Any,
    bradford_stds: Any,
    standards_tiprack_50: Any,
    tip_col: int,
    *,
    standard_aspirate_volume: float,
    standard_dispense_volume: float,
) -> None:
    await _pick_up_column_tips(lh, standards_tiprack_50, tip_col)
    await _aspirate_plate_column(
        lh,
        bradford_stds,
        1,
        standard_aspirate_volume,
        liquid_height=1.0,
        flow_rate=35,
    )

    for destination_col in (10, 11, 12):
        await _dispense_to_column(
            lh,
            destination_plate,
            destination_col,
            standard_dispense_volume,
            reference_volume=255,
            flow_rate=30,
            blow_out=False,
        )

    await _discard_column_tips(lh)
    print(f"Bradford standards added to assay columns 10-12 using tip column {tip_col}.")


async def run_bpq_96w(
    lh: Any,
    elution_plate: Any,
    sample_plate_1: Any,
    sample_plate_2: Any,
    sample_plate_3: Any,
    bradford_reagent: Any,
    bradford_stds: Any,
    tiprack_1000: Any,
    sample_tiprack_50: Any,
    standards_tiprack_50: Any,
    *,
    add_reagent_sample_plate_1: bool = True,
    add_samples_plate_1: bool = True,
    add_stds_plate_1: bool = True,
    add_reagent_sample_plate_2: bool = True,
    add_samples_plate_2: bool = True,
    add_stds_plate_2: bool = True,
    add_reagent_sample_plate_3: bool = True,
    add_samples_plate_3: bool = True,
    add_stds_plate_3: bool = True,
    reagent_volume: float = 250,
    reagent_aspirate_volume: float = 800,
    sample_aspirate_volume: float = 12,
    sample_dispense_volume: float = 5,
    standard_aspirate_volume: float = 18,
    standard_dispense_volume: float = 5,
) -> None:
    """Run the Bradford Protein Quantification workflow.

    The notebook is responsible for constructing the deck, calling `await lh.setup(...)`,
    and passing the resources created in `BPQ_96w.ipynb`.
    """

    if sample_aspirate_volume < sample_dispense_volume * 2:
        raise ValueError("sample_aspirate_volume must cover both duplicate sample dispenses.")
    if standard_aspirate_volume < standard_dispense_volume * 3:
        raise ValueError("standard_aspirate_volume must cover all three standards dispenses.")
    if reagent_aspirate_volume <= 0 or reagent_volume <= 0:
        raise ValueError("Reagent volumes must be positive.")

    plate_configs = [
        {
            "plate_index": 1,
            "plate": sample_plate_1,
            "source_cols": range(1, 5),
            "standards_tip_col": 1,
            "add_reagent": add_reagent_sample_plate_1,
            "add_samples": add_samples_plate_1,
            "add_stds": add_stds_plate_1,
        },
        {
            "plate_index": 2,
            "plate": sample_plate_2,
            "source_cols": range(5, 9),
            "standards_tip_col": 2,
            "add_reagent": add_reagent_sample_plate_2,
            "add_samples": add_samples_plate_2,
            "add_stds": add_stds_plate_2,
        },
        {
            "plate_index": 3,
            "plate": sample_plate_3,
            "source_cols": range(9, 13),
            "standards_tip_col": 3,
            "add_reagent": add_reagent_sample_plate_3,
            "add_samples": add_samples_plate_3,
            "add_stds": add_stds_plate_3,
        },
    ]

    print("Starting BPQ 96-well protocol.")

    reagent_targets = [
        (config["plate_index"], config["plate"])
        for config in plate_configs
        if config["add_reagent"]
    ]
    if reagent_targets:
        await _fill_sample_plates_with_reagent(
            lh,
            reagent_targets,
            bradford_reagent,
            tiprack_1000,
            reagent_volume=reagent_volume,
            reagent_aspirate_volume=reagent_aspirate_volume,
        )
    else:
        print("Skipping Bradford reagent fill for all sample plates.")

    for config in plate_configs:
        plate_index = config["plate_index"]

        if not config["add_reagent"]:
            print(f"Skipping Bradford reagent addition for SamplePlate{plate_index}.")

        if config["add_samples"]:
            await _duplicate_elution_samples_to_plate(
                lh,
                elution_plate,
                config["plate"],
                sample_tiprack_50,
                config["source_cols"],
                sample_aspirate_volume=sample_aspirate_volume,
                sample_dispense_volume=sample_dispense_volume,
            )
            print(f"Completed elution sample transfer for SamplePlate{plate_index}.")
        else:
            print(f"Skipping elution sample transfer for SamplePlate{plate_index}.")

        if config["add_stds"]:
            await _triplicate_standards_to_plate(
                lh,
                config["plate"],
                bradford_stds,
                standards_tiprack_50,
                tip_col=config["standards_tip_col"],
                standard_aspirate_volume=standard_aspirate_volume,
                standard_dispense_volume=standard_dispense_volume,
            )
            print(f"Completed Bradford standards transfer for SamplePlate{plate_index}.")
        else:
            print(f"Skipping Bradford standards transfer for SamplePlate{plate_index}.")

    print("BPQ 96-well protocol complete.")
