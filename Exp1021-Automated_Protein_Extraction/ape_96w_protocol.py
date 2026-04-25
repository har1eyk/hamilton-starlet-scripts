"""Reusable 96-well Automated Protein Extraction protocol for PyLabRobot."""

from __future__ import annotations

import asyncio
import random
import time
from datetime import datetime
from typing import Any, Iterable, Sequence

from pylabrobot.liquid_handling.standard import Mix

CHANNELS = list(range(8))
COLUMN_COUNT = 12


def _plate_col(plate: Any, col: int) -> Any:
    return plate[f"A{col}:H{col}"]


def _trough_col(trough: Any) -> Any:
    return trough["A1"] * 8


def _column_numbers() -> Iterable[int]:
    return range(1, COLUMN_COUNT + 1)


def _tip_column(rack: Any, col: int) -> Any:
    return rack[f"A{col}:H{col}"]


def _well_height(plate: Any, volume: float) -> float:
    well = plate.get_item("A1")
    return max(0.0, well.compute_height_from_volume(max(0.0, volume)))


def _gripper_pair() -> tuple[int, int]:
    grip_one = random.randint(1, 6)
    return grip_one, grip_one + 1


def _format_timestamp(moment: datetime) -> str:
    return moment.strftime("%Y-%m-%d %I:%M:%S %p %Z")


def _human_readable_duration(seconds: float) -> str:
    if seconds < 1:
        return "less than 1 second"

    rounded_seconds = int(round(seconds))
    hours, remainder = divmod(rounded_seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    parts: list[str] = []

    if hours:
        parts.append(f"{hours} hour{'s' if hours != 1 else ''}")
    if minutes:
        parts.append(f"{minutes} minute{'s' if minutes != 1 else ''}")
    if secs or not parts:
        parts.append(f"{secs} second{'s' if secs != 1 else ''}")

    if len(parts) == 1:
        return parts[0]
    if len(parts) == 2:
        return f"{parts[0]} and {parts[1]}"
    return ", ".join(parts[:-1]) + f", and {parts[-1]}"


def _start_step(label: str) -> tuple[float, datetime]:
    started_at = datetime.now().astimezone()
    print(f"[START] {label} at {_format_timestamp(started_at)}.")
    return time.perf_counter(), started_at


def _record_step_timing(
    step_timings: list[tuple[str, float | None, str | None]],
    label: str,
    *,
    started_at: datetime | None = None,
    elapsed_seconds: float | None = None,
    detail: str | None = None,
    skipped: bool = False,
) -> None:
    if skipped:
        print(f"{label} was skipped.")
        step_timings.append((label, None, "Skipped."))
        return

    ended_at = datetime.now().astimezone()
    if started_at is not None:
        print(f"[END] {label} at {_format_timestamp(ended_at)}.")

    duration_text = _human_readable_duration(elapsed_seconds or 0)
    message = f"{label} took {duration_text} to complete."
    if detail:
        message += f" {detail}"
    print(message)
    step_timings.append((label, elapsed_seconds, detail))


def _pass_volumes(total_volume: float, *, max_single: float = 900, cleanup_volume: float = 50) -> list[float]:
    if total_volume <= 0:
        return []
    if total_volume < 500 or cleanup_volume <= 0:
        remaining = total_volume
        passes: list[float] = []
        while remaining > max_single:
            passes.append(max_single)
            remaining -= max_single
        passes.append(remaining)
        return passes

    passes = []
    remaining = total_volume
    while remaining > (max_single + cleanup_volume):
        passes.append(max_single)
        remaining -= max_single
    if remaining > cleanup_volume:
        passes.append(remaining - cleanup_volume)
    passes.append(cleanup_volume)
    return [vol for vol in passes if vol > 0]


async def _pick_up_column_tips(lh: Any, rack: Any, col: int) -> None:
    await lh.pick_up_tips(_tip_column(rack, col), use_channels=CHANNELS)


async def _return_column_tips(lh: Any, rack: Any, col: int) -> None:
    await lh.drop_tips(_tip_column(rack, col), use_channels=CHANNELS)


def _require_tip_racks(tip_racks: Sequence[Any]) -> tuple[Any, Any, Any]:
    if len(tip_racks) < 3:
        raise ValueError("run_ape_96w requires three 1000 uL tip racks.")
    return tip_racks[0], tip_racks[1], tip_racks[2]


def _ensure_core_pickup_tracking(lh: Any) -> None:
    if getattr(lh.backend, "num_arms", None) == 0:
        lh._resource_pickups = {0: None}


async def _move_plate_to_magnet(lh: Any, plate: Any, mag_plate: Any, *, grip_pair: tuple[int, int]) -> None:
    grip_one, grip_two = grip_pair
    if abs(grip_one - grip_two) != 1:
        raise ValueError("Co-Re gripper channels must be adjacent.")

    _ensure_core_pickup_tracking(lh)
    mag_plate.plate_z_offset = -5.5
    channel_1 = min(grip_one, grip_two)
    channel_2 = max(grip_one, grip_two)
    print(f"Using Co-Re gripper channels {channel_1} and {channel_2}.")

    await lh.move_plate(
        plate=plate,
        to=mag_plate,
        use_arm="core",
        core_front_channel=channel_1,
        core_grip_strength=50,
        pickup_distance_from_top=10,
        enable_recovery=True,
        return_core_gripper=False,
    )

    await lh.backend.core_check_resource_exists_at_location_center(
        location=plate.get_absolute_location(),
        resource=plate,
        gripper_y_margin=9,
        enable_recovery=True,
        audio_feedback=False,
    )
    await lh.backend.return_core_gripper_tools()


async def _move_plate_home(lh: Any, plate: Any, home: Any, *, grip_pair: tuple[int, int]) -> None:
    grip_one, grip_two = grip_pair
    if abs(grip_one - grip_two) != 1:
        raise ValueError("Co-Re gripper channels must be adjacent.")

    _ensure_core_pickup_tracking(lh)
    channel_1 = min(grip_one, grip_two)

    await lh.move_plate(
        plate=plate,
        to=home,
        use_arm="core",
        pickup_distance_from_top=10,
        core_front_channel=channel_1,
        core_grip_strength=50,
        enable_recovery=False,
        return_core_gripper=False,
    )
    await lh.backend.return_core_gripper_tools()


async def _high_volume_mix(
    lh: Any,
    plate: Any,
    plate_column: Any,
    *,
    mix_volume: float,
    total_volume: float,
    repetitions: int = 3,
) -> None:
    aspirate_height = max(1.0, _well_height(plate, total_volume) - 1)
    dispense_height = max(1.0, _well_height(plate, total_volume - mix_volume) - 1)
    for _ in range(repetitions):
        await lh.aspirate(
            plate_column,
            vols=[mix_volume] * 8,
            use_channels=CHANNELS,
            liquid_height=[aspirate_height] * 8,
            auto_surface_following_distance=True,
            swap_speed=[160] * 8,
            flow_rates=[400] * 8,
            minimum_height=[aspirate_height] * 8,
        )
        await lh.dispense(
            plate_column,
            vols=[mix_volume] * 8,
            use_channels=CHANNELS,
            liquid_height=[dispense_height] * 8,
            flow_rates=[500] * 8,
            auto_surface_following_distance=True,
            blow_out=[1] * 8,
            swap_speed=[160] * 8,
            settling_time=[1] * 8,
        )


async def _mix_column(
    lh: Any,
    plate: Any,
    plate_column: Any,
    *,
    mix_volume: float,
    dispense_height: float,
    repetitions: int,
    flow_rate: float,
) -> None:
    aspirate_height = max(1.0, _well_height(plate, mix_volume) - 1)
    await lh.aspirate(
        plate_column,
        vols=[mix_volume] * 8,
        use_channels=CHANNELS,
        liquid_height=[aspirate_height] * 8,
        auto_surface_following_distance=True,
        flow_rates=[flow_rate] * 8,
    )
    await lh.dispense(
        plate_column,
        vols=[mix_volume] * 8,
        use_channels=CHANNELS,
        liquid_height=[dispense_height] * 8,
        flow_rates=[flow_rate] * 8,
        mix=[Mix(volume=mix_volume, repetitions=repetitions, flow_rate=flow_rate)] * 8,
        blow_out=[1] * 8,
        settling_time=[1] * 8,
    )


async def _transfer_from_trough(
    lh: Any,
    trough: Any,
    destination_plate: Any,
    tip_rack: Any,
    *,
    total_volume: float,
    pass_volumes: Sequence[float],
    destination_height: float | Sequence[float],
    label: str,
    fixed_tip_column: int | None = None,
) -> None:
    if isinstance(destination_height, (int, float)):
        destination_heights = [float(destination_height)] * len(pass_volumes)
    else:
        destination_heights = [float(height) for height in destination_height]
        if len(destination_heights) != len(pass_volumes):
            raise ValueError("destination_height must match the number of pass_volumes.")

    picked_up_fixed_tips = False
    for col in _column_numbers():
        tip_col = col if fixed_tip_column is None else fixed_tip_column
        if fixed_tip_column is None or not picked_up_fixed_tips:
            await _pick_up_column_tips(lh, tip_rack, tip_col)
            picked_up_fixed_tips = True
        for transfer_volume, dispense_height in zip(pass_volumes, destination_heights):
            await lh.aspirate(
                _trough_col(trough),
                vols=[transfer_volume] * 8,
                use_channels=CHANNELS,
                liquid_height=[2] * 8,
            )
            await lh.dispense(
                _plate_col(destination_plate, col),
                vols=[transfer_volume] * 8,
                use_channels=CHANNELS,
                liquid_height=[dispense_height] * 8,
                blow_out=[1] * 8,
                settling_time=[1] * 8,
            )
        if fixed_tip_column is None:
            await _return_column_tips(lh, tip_rack, tip_col)
            picked_up_fixed_tips = False
        print(f"{label}: column {col}/12 ({total_volume:.0f} uL).")
    if fixed_tip_column is not None and picked_up_fixed_tips:
        await _return_column_tips(lh, tip_rack, fixed_tip_column)


async def _transfer_column_to_column(
    lh: Any,
    source_plate: Any,
    destination_plate: Any,
    tip_rack: Any,
    *,
    volume: float,
    source_height: float | None = None,
    destination_height: float | None = None,
    aspirate_flow_rate: float | None = None,
    dispense_flow_rate: float | None = None,
) -> None:
    for col in _column_numbers():
        await _pick_up_column_tips(lh, tip_rack, col)
        src_col = _plate_col(source_plate, col)
        dst_col = _plate_col(destination_plate, col)
        asp_height = source_height if source_height is not None else max(1.0, _well_height(source_plate, volume) - 1)
        disp_height = destination_height if destination_height is not None else max(1.0, _well_height(destination_plate, volume) - 1)
        aspirate_kwargs = {}
        dispense_kwargs = {}
        if aspirate_flow_rate is not None:
            aspirate_kwargs["flow_rates"] = [aspirate_flow_rate] * 8
        if dispense_flow_rate is not None:
            dispense_kwargs["flow_rates"] = [dispense_flow_rate] * 8

        await lh.aspirate(
            src_col,
            vols=[volume] * 8,
            use_channels=CHANNELS,
            liquid_height=[asp_height] * 8,
            auto_surface_following_distance=True,
            **aspirate_kwargs,
        )
        await lh.dispense(
            dst_col,
            vols=[volume] * 8,
            use_channels=CHANNELS,
            liquid_height=[disp_height] * 8,
            blow_out=[1] * 8,
            settling_time=[1] * 8,
            **dispense_kwargs,
        )
        await _return_column_tips(lh, tip_rack, col)
        print(f"Transferred column {col}/12.")


async def _mix_plate_columns(
    lh: Any,
    plate: Any,
    tip_rack: Any,
    *,
    rounds: int,
    mix_volume: float,
    total_volume: float | None = None,
    incubate_seconds: int,
    label: str,
    low_volume_repetitions: int = 2,
    low_volume_height: float = 6,
    low_volume_flow_rate: float = 400,
) -> None:
    for round_index in range(rounds):
        for col in _column_numbers():
            await _pick_up_column_tips(lh, tip_rack, col)
            plate_col = _plate_col(plate, col)
            if total_volume is None:
                await _mix_column(
                    lh,
                    plate,
                    plate_col,
                    mix_volume=mix_volume,
                    dispense_height=low_volume_height,
                    repetitions=low_volume_repetitions,
                    flow_rate=low_volume_flow_rate,
                )
            else:
                await _high_volume_mix(
                    lh,
                    plate,
                    plate_col,
                    mix_volume=mix_volume,
                    total_volume=total_volume,
                    repetitions=3,
                )
            await _return_column_tips(lh, tip_rack, col)
        print(f"{label}: round {round_index + 1}/{rounds}")
        if round_index < rounds - 1 and incubate_seconds > 0:
            await asyncio.sleep(incubate_seconds)


async def _remove_supernatant_to_plate(
    lh: Any,
    source_plate: Any,
    destination_plate: Any,
    tip_rack: Any,
    *,
    total_volume: float,
    destination_height: float,
    cleanup_volume: float = 50,
    discard_tips: bool = False,
) -> None:
    pass_volumes = _pass_volumes(total_volume, cleanup_volume=cleanup_volume)
    for col in _column_numbers():
        await _pick_up_column_tips(lh, tip_rack, col)
        source_col = _plate_col(source_plate, col)
        destination_col = _plate_col(destination_plate, col)
        remaining = total_volume
        for index, transfer_volume in enumerate(pass_volumes):
            is_cleanup_pass = index == len(pass_volumes) - 1 and cleanup_volume > 0
            aspirate_height = 0 if index == len(pass_volumes) - 1 and cleanup_volume > 0 else max(
                1.0, _well_height(source_plate, remaining) - 1
            )
            await lh.aspirate(
                source_col,
                vols=[transfer_volume] * 8,
                use_channels=CHANNELS,
                liquid_height=[aspirate_height] * 8,
                flow_rates=[50 if is_cleanup_pass else 100] * 8,
                auto_surface_following_distance=not is_cleanup_pass,
            )
            remaining = max(0.0, remaining - transfer_volume)
            await lh.dispense(
                destination_col,
                vols=[transfer_volume] * 8,
                use_channels=CHANNELS,
                liquid_height=[destination_height] * 8,
                blow_out=[1] * 8,
                settling_time=[1] * 8,
            )
        if discard_tips:
            await lh.discard_tips(use_channels=CHANNELS)
        else:
            await _return_column_tips(lh, tip_rack, col)
        print(f"Removed supernatant from column {col}/12.")


async def _remove_supernatant_to_trough(
    lh: Any,
    source_plate: Any,
    waste_trough: Any,
    tip_rack: Any,
    *,
    total_volume: float,
    cleanup_volume: float = 50,
    fixed_tip_column: int | None = None,
) -> None:
    pass_volumes = _pass_volumes(total_volume, cleanup_volume=cleanup_volume)
    picked_up_fixed_tips = False
    for col in _column_numbers():
        tip_col = col if fixed_tip_column is None else fixed_tip_column
        if fixed_tip_column is None or not picked_up_fixed_tips:
            await _pick_up_column_tips(lh, tip_rack, tip_col)
            picked_up_fixed_tips = fixed_tip_column is not None
        source_col = _plate_col(source_plate, col)
        remaining = total_volume
        for index, transfer_volume in enumerate(pass_volumes):
            liquid_height = 0 if index == len(pass_volumes) - 1 and cleanup_volume > 0 else max(
                1.0, _well_height(source_plate, remaining) - 1
            )
            await lh.aspirate(
                source_col,
                vols=[transfer_volume] * 8,
                use_channels=CHANNELS,
                liquid_height=[liquid_height] * 8,
                auto_surface_following_distance=index != len(pass_volumes) - 1,
            )
            remaining = max(0.0, remaining - transfer_volume)
            await lh.dispense(
                _trough_col(waste_trough),
                vols=[transfer_volume] * 8,
                use_channels=CHANNELS,
                liquid_height=[2] * 8,
                blow_out=[1] * 8,
                settling_time=[1] * 8,
            )
        if fixed_tip_column is None:
            await _return_column_tips(lh, tip_rack, tip_col)
        print(f"Cleared wash from column {col}/12.")
    if fixed_tip_column is not None and picked_up_fixed_tips:
        await _return_column_tips(lh, tip_rack, fixed_tip_column)


def _pause_for_manual_beads(bead_volume_note: float) -> None:
    print(
        "\n" + "=" * 72 +
        f"\nManual step required: add {bead_volume_note:.0f} uL magnetic beads to every well of the binding plate.\n"
        "The robot is paused and will not resume until you confirm the beads were added.\n" +
        "=" * 72 + "\n"
    )
    input("Press Enter after manually adding beads to all 96 wells... ")


async def run_ape_96w(
    lh: Any,
    binding_plate: Any,
    source_plate: Any,
    flowthrough_plate: Any,
    elution_plate: Any,
    mag_plate: Any,
    binding_trough: Any,
    wash1_trough: Any,
    wash2_trough: Any,
    waste_trough: Any,
    elution_trough: Any,
    tip_racks: Sequence[Any],
    *,
    fill_binding_plate: bool = True,
    pause_for_beads: bool = True,
    transfer_lysate: bool = True,
    mix_binding_step: bool = True,
    remove_binding_supernatant: bool = True,
    run_wash1: bool = True,
    run_wash2: bool = True,
    run_final_wash: bool = True,
    add_elution_buffer: bool = True,
    mix_elution_step: bool = True,
    recover_elution: bool = True,
    lysate_volume: float = 300,
    binding_buffer_volume: float = 1500,
    bead_volume_note: float = 50,
    wash_volume: float = 900,
    elution_volume: float = 120,
    binding_supernatant_cleanup_volume: float = 150,
    wash1_cycles: int = 1,
    binding_mix_rounds: int = 3,
    binding_incubation_seconds: int = 120,
    wash_mix_rounds: int = 3,
    magnet_settle_seconds: int = 20,
    wash_incubation_seconds: int = 120,
    elution_incubation_seconds: int = 120,
    final_wash_trough: Any | None = None,
) -> None:
    """Run the 96-well Automated Protein Extraction workflow.

    The notebook is responsible for creating `lh`, assigning resources to the deck, and calling
    `await lh.setup(skip_autoload=True)` before invoking this function.
    """

    fill_tip_rack, process_tip_rack, elution_tip_rack = _require_tip_racks(tip_racks)
    binding_plate_home = binding_plate.parent
    final_wash_source = wash2_trough if final_wash_trough is None else final_wash_trough
    step_timings: list[tuple[str, float | None, str | None]] = []
    protocol_start = time.perf_counter()
    protocol_started_at = datetime.now().astimezone()

    print(f"Starting APE 96-well protocol at {_format_timestamp(protocol_started_at)}.")

    if fill_binding_plate:
        step_start, step_started_at = _start_step("Fill binding plate")
        fill_passes = [binding_buffer_volume / 2, binding_buffer_volume / 2]
        first_binding_dispense_height = max(1.0, _well_height(binding_plate, fill_passes[0]) - 4.0)
        await _transfer_from_trough(
            lh,
            binding_trough,
            binding_plate,
            fill_tip_rack,
            total_volume=binding_buffer_volume,
            pass_volumes=fill_passes,
            destination_height=[first_binding_dispense_height, 24.0],
            label="Binding buffer added",
            fixed_tip_column=1,
        )
        _record_step_timing(
            step_timings,
            "Fill binding plate",
            started_at=step_started_at,
            elapsed_seconds=time.perf_counter() - step_start,
        )
    else:
        print("Skipping binding plate fill; assuming BindingPlate is already pre-filled.")
        _record_step_timing(step_timings, "Fill binding plate", skipped=True)

    if pause_for_beads and (fill_binding_plate or transfer_lysate or mix_binding_step):
        step_start, step_started_at = _start_step("Manual bead addition pause")
        _pause_for_manual_beads(bead_volume_note)
        _record_step_timing(
            step_timings,
            "Manual bead addition pause",
            started_at=step_started_at,
            elapsed_seconds=time.perf_counter() - step_start,
        )
    else:
        _record_step_timing(step_timings, "Manual bead addition pause", skipped=True)

    if transfer_lysate:
        step_start, step_started_at = _start_step("Transfer lysate to binding plate")
        await _transfer_column_to_column(
            lh,
            source_plate,
            binding_plate,
            process_tip_rack,
            volume=lysate_volume,
            destination_height=24,
            aspirate_flow_rate=50,
        )
        _record_step_timing(
            step_timings,
            "Transfer lysate to binding plate",
            started_at=step_started_at,
            elapsed_seconds=time.perf_counter() - step_start,
        )
    else:
        print("Skipping lysate transfer; assuming BindingPlate already contains lysate.")
        _record_step_timing(step_timings, "Transfer lysate to binding plate", skipped=True)

    if mix_binding_step:
        step_start, step_started_at = _start_step("Protein binding step")
        await _mix_plate_columns(
            lh,
            binding_plate,
            process_tip_rack,
            rounds=binding_mix_rounds,
            mix_volume=1000,
            total_volume=binding_buffer_volume + lysate_volume + bead_volume_note,
            incubate_seconds=binding_incubation_seconds,
            label="Binding mix",
        )
        binding_step_elapsed = time.perf_counter() - step_start
        binding_detail = (
            f"You performed {binding_mix_rounds} rounds because binding_mix_rounds={binding_mix_rounds}, "
            f"so each round took about {_human_readable_duration(binding_step_elapsed / binding_mix_rounds)}."
            if binding_mix_rounds > 0
            else "No binding rounds were run because binding_mix_rounds=0."
        )
        _record_step_timing(
            step_timings,
            "Protein binding step",
            started_at=step_started_at,
            elapsed_seconds=binding_step_elapsed,
            detail=binding_detail,
        )
    else:
        print("Skipping lysate/bead mixing.")
        _record_step_timing(step_timings, "Protein binding step", skipped=True)

    if remove_binding_supernatant:
        step_start, step_started_at = _start_step("Capture on magnet and remove binding supernatant")
        grip_pair = _gripper_pair()
        await _move_plate_to_magnet(lh, binding_plate, mag_plate, grip_pair=grip_pair)
        await asyncio.sleep(magnet_settle_seconds)
        await _remove_supernatant_to_plate(
            lh,
            binding_plate,
            flowthrough_plate,
            process_tip_rack,
            total_volume=binding_buffer_volume + lysate_volume + bead_volume_note,
            destination_height=22,
            cleanup_volume=binding_supernatant_cleanup_volume,
            discard_tips=True,
        )
        await _move_plate_home(lh, binding_plate, binding_plate_home, grip_pair=grip_pair)
        _record_step_timing(
            step_timings,
            "Capture on magnet and remove binding supernatant",
            started_at=step_started_at,
            elapsed_seconds=time.perf_counter() - step_start,
        )
    else:
        print("Skipping magnet capture and binding supernatant removal.")
        _record_step_timing(step_timings, "Capture on magnet and remove binding supernatant", skipped=True)

    wash_steps = [
        ("Wash 1", wash1_trough, run_wash1, wash1_cycles, wash_mix_rounds),
        ("Wash 2", wash2_trough, run_wash2, 1, wash_mix_rounds),
        ("Final wash", final_wash_source, run_final_wash, 1, wash_mix_rounds),
    ]
    wash_dispense_height = _well_height(binding_plate, wash_volume) + 4.0
    for wash_label, wash_source, enabled, stage_repetitions, stage_mix_rounds in wash_steps:
        step_label = f"{wash_label} step"
        if not enabled:
            print(f"Skipping {wash_label}.")
            _record_step_timing(step_timings, step_label, skipped=True)
            continue
        step_start, step_started_at = _start_step(step_label)
        for cycle_index in range(stage_repetitions):
            await _transfer_from_trough(
                lh,
                wash_source,
                binding_plate,
                fill_tip_rack,
                total_volume=wash_volume,
                pass_volumes=[wash_volume],
                destination_height=wash_dispense_height,
                label=f"{wash_label} added cycle {cycle_index + 1}",
                fixed_tip_column=1,
            )
            await _mix_plate_columns(
                lh,
                binding_plate,
                fill_tip_rack,
                rounds=stage_mix_rounds,
                mix_volume=min(wash_volume, 900),
                total_volume=None,
                incubate_seconds=wash_incubation_seconds, #no need to incubate here.
                label=f"{wash_label} mix cycle {cycle_index + 1}",
                low_volume_repetitions=2,
                low_volume_height=3, # wash 6
                low_volume_flow_rate=400,
            )
            grip_pair = _gripper_pair()
            await _move_plate_to_magnet(lh, binding_plate, mag_plate, grip_pair=grip_pair)
            await asyncio.sleep(magnet_settle_seconds)
            await _remove_supernatant_to_trough(
                lh,
                binding_plate,
                waste_trough,
                fill_tip_rack,
                total_volume=wash_volume,
                cleanup_volume=50,
                fixed_tip_column=1,
            )
            await _move_plate_home(lh, binding_plate, binding_plate_home, grip_pair=grip_pair)
        wash_step_elapsed = time.perf_counter() - step_start
        wash_detail = None
        if stage_repetitions > 1:
            wash_detail = (
                f"You performed {stage_repetitions} cycles, so each cycle took about "
                f"{_human_readable_duration(wash_step_elapsed / stage_repetitions)}."
            )
        _record_step_timing(
            step_timings,
            step_label,
            started_at=step_started_at,
            elapsed_seconds=wash_step_elapsed,
            detail=wash_detail,
        )

    if add_elution_buffer:
        step_start, step_started_at = _start_step("Add elution buffer")
        await _transfer_from_trough(
            lh,
            elution_trough,
            binding_plate,
            elution_tip_rack,
            total_volume=elution_volume,
            pass_volumes=[elution_volume],
            destination_height=5,
            label="Elution buffer added",
            fixed_tip_column=1,
        )
        _record_step_timing(
            step_timings,
            "Add elution buffer",
            started_at=step_started_at,
            elapsed_seconds=time.perf_counter() - step_start,
        )
    else:
        print("Skipping elution buffer addition.")
        _record_step_timing(step_timings, "Add elution buffer", skipped=True)

    if mix_elution_step:
        step_start, step_started_at = _start_step("Elution mix and incubation")
        await asyncio.sleep(elution_incubation_seconds)
        await _mix_plate_columns(
            lh,
            binding_plate,
            elution_tip_rack,
            rounds=1,
            mix_volume=min(100, elution_volume),
            total_volume=None,
            incubate_seconds=0,
            label="Elution mix",
            low_volume_repetitions=4,
            low_volume_height=2,
            low_volume_flow_rate=100,
        )
        await asyncio.sleep(elution_incubation_seconds)
        _record_step_timing(
            step_timings,
            "Elution mix and incubation",
            started_at=step_started_at,
            elapsed_seconds=time.perf_counter() - step_start,
        )
    else:
        print("Skipping elution mixing/incubation.")
        _record_step_timing(step_timings, "Elution mix and incubation", skipped=True)

    if recover_elution:
        step_start, step_started_at = _start_step("Recover elution")
        grip_pair = _gripper_pair()
        await _move_plate_to_magnet(lh, binding_plate, mag_plate, grip_pair=grip_pair)
        await asyncio.sleep(magnet_settle_seconds)
        await _remove_supernatant_to_plate(
            lh,
            binding_plate,
            elution_plate,
            elution_tip_rack,
            total_volume=elution_volume,
            destination_height=5,
            cleanup_volume=0,
        )
        await _move_plate_home(lh, binding_plate, binding_plate_home, grip_pair=grip_pair)
        _record_step_timing(
            step_timings,
            "Recover elution",
            started_at=step_started_at,
            elapsed_seconds=time.perf_counter() - step_start,
        )
    else:
        print("Skipping elution recovery.")
        _record_step_timing(step_timings, "Recover elution", skipped=True)

    total_elapsed = time.perf_counter() - protocol_start
    print("\nAPE 96-well timing summary:")
    for label, elapsed_seconds, detail in step_timings:
        if elapsed_seconds is None:
            print(f"- {label}: skipped")
            continue
        summary = f"- {label}: {_human_readable_duration(elapsed_seconds)}"
        if detail:
            summary += f" {detail}"
        print(summary)
    print(f"- Total protocol time: {_human_readable_duration(total_elapsed)}")
    print("APE 96-well protocol complete.")
