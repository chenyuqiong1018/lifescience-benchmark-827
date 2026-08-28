from opentrons import protocol_api

metadata = {
    "protocolName": "Audited 24-sample magnetic-bead cleanup",
    "author": "Codex benchmark execution",
    "description": "Frozen SOP translated to OT-2 Protocol API 2.16",
    "apiLevel": "2.16",
}

SAMPLE_WELLS = (
    "A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2",
    "A3", "B3", "C3", "D3", "A4", "B4", "C4", "D4",
    "A5", "B5", "C5", "D5", "A6", "B6", "C6", "D6",
)

REAGENTS = {
    "lysis": {"well": "A1", "initial": 2500, "dead": 300, "per_sample": 80, "uses": 1},
    "beads": {"well": "A2", "initial": 3500, "dead": 300, "per_sample": 120, "uses": 1},
    "wash": {"well": "A3", "initial": 9000, "dead": 300, "per_sample": 180, "uses": 2},
    "elution": {"well": "A4", "initial": 1500, "dead": 300, "per_sample": 40, "uses": 1},
}


def preflight() -> None:
    sample_count = len(SAMPLE_WELLS)
    if sample_count != 24 or len(set(SAMPLE_WELLS)) != sample_count:
        raise RuntimeError("Sample map is not the frozen set of 24 unique wells")
    for name, reagent in REAGENTS.items():
        consumed = reagent["per_sample"] * reagent["uses"] * sample_count
        if reagent["initial"] - consumed < reagent["dead"]:
            raise RuntimeError(f"{name} volume would fall below its dead volume")
    if sample_count * 6 > 192:
        raise RuntimeError("Tip policy requires more tips than the declared inventory")
    if (250 + 180 + 180) * sample_count > 195000:
        raise RuntimeError("Waste volume exceeds the declared reservoir capacity")
    if 50 + 80 + 120 > 2000:
        raise RuntimeError("Peak sample-well volume exceeds declared capacity")
    for volume in (30, 40, 80, 100, 120, 180, 250):
        if not 20 <= volume <= 300:
            raise RuntimeError(f"P300 range violation at {volume} uL")


def add_mix_with_fresh_tip(pipette, source, sample, add_ul, repetitions, mix_ul):
    pipette.pick_up_tip()
    pipette.aspirate(add_ul, source)
    pipette.dispense(add_ul, sample)
    pipette.mix(repetitions, mix_ul, sample)
    pipette.drop_tip()


def remove_with_fresh_tip(pipette, sample, waste, volume_ul):
    pipette.pick_up_tip()
    pipette.aspirate(volume_ul, sample)
    pipette.dispense(volume_ul, waste)
    pipette.drop_tip()


def wash_same_sample(protocol, pipette, wash_source, sample, waste, cycle):
    pipette.pick_up_tip()
    pipette.aspirate(180, wash_source)
    pipette.dispense(180, sample)
    protocol.delay(seconds=30, msg=f"Wash {cycle} dwell; same-sample tip remains attached")
    pipette.aspirate(180, sample)
    pipette.dispense(180, waste)
    pipette.drop_tip()


def run(protocol: protocol_api.ProtocolContext) -> None:
    preflight()

    reagent_labware = protocol.load_labware("nest_12_reservoir_15ml", "1", "reagents")
    tiprack_slot_4 = protocol.load_labware("opentrons_96_tiprack_300ul", "4", "tips_slot_4")
    magnetic_module = protocol.load_module("magnetic module gen2", "5")
    processing_plate = magnetic_module.load_labware(
        "nest_96_wellplate_2ml_deep", "processing_plate"
    )
    waste_labware = protocol.load_labware("nest_1_reservoir_195ml", "6", "liquid_waste")
    tiprack_slot_7 = protocol.load_labware("opentrons_96_tiprack_300ul", "7", "tips_slot_7")
    pipette = protocol.load_instrument(
        "p300_single_gen2", "right", tip_racks=[tiprack_slot_4, tiprack_slot_7]
    )

    source = {name: reagent_labware[data["well"]] for name, data in REAGENTS.items()}
    waste = waste_labware["A1"]
    samples = [processing_plate[well] for well in SAMPLE_WELLS]

    protocol.comment("Lysis additions and 5 x 100 uL mixes")
    for sample in samples:
        add_mix_with_fresh_tip(pipette, source["lysis"], sample, 80, 5, 100)
    protocol.delay(minutes=5, msg="Batch lysis incubation")

    protocol.comment("Bead-source resuspension, additions, and 10 x 180 uL sample mixes")
    for sample in samples:
        pipette.pick_up_tip()
        pipette.mix(10, 180, source["beads"])
        pipette.aspirate(120, source["beads"])
        pipette.dispense(120, sample)
        pipette.mix(10, 180, sample)
        pipette.drop_tip()
    protocol.delay(minutes=5, msg="Batch bead-binding incubation")

    magnetic_module.engage()
    protocol.delay(minutes=7, msg="Magnetic separation")

    protocol.comment("Supernatant removals")
    for sample in samples:
        remove_with_fresh_tip(pipette, sample, waste, 250)

    protocol.comment("Wash cycle 1 with one retained tip per sample")
    for sample in samples:
        wash_same_sample(protocol, pipette, source["wash"], sample, waste, 1)

    protocol.comment("Wash cycle 2 with one retained tip per sample")
    for sample in samples:
        wash_same_sample(protocol, pipette, source["wash"], sample, waste, 2)

    protocol.delay(minutes=2, msg="Air dry with magnet engaged")
    magnetic_module.disengage()

    protocol.comment("Elution additions and 10 x 30 uL mixes")
    for sample in samples:
        add_mix_with_fresh_tip(pipette, source["elution"], sample, 40, 10, 30)

