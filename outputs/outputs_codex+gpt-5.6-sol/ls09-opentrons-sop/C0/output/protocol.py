from opentrons import protocol_api

metadata = {
    "protocolName": "Frozen 24-sample magnetic-bead cleanup",
    "author": "Codex benchmark execution",
    "description": "Auditable implementation of the supplied local SOP",
    "apiLevel": "2.16",
}

requirements = {"robotType": "OT-2"}

SAMPLE_WELLS = [
    "A1", "B1", "C1", "D1", "A2", "B2", "C2", "D2",
    "A3", "B3", "C3", "D3", "A4", "B4", "C4", "D4",
    "A5", "B5", "C5", "D5", "A6", "B6", "C6", "D6",
]


def _preflight() -> None:
    reagent_budget = {
        "lysis_buffer": (2500, 300, 80 * 24),
        "bead_suspension": (3500, 300, 120 * 24),
        "wash_buffer": (9000, 300, 180 * 2 * 24),
        "elution_buffer": (1500, 300, 40 * 24),
    }
    if len(SAMPLE_WELLS) != 24 or len(set(SAMPLE_WELLS)) != 24:
        raise RuntimeError("Expected 24 unique frozen sample wells")
    for name, (initial, dead, required) in reagent_budget.items():
        if initial - required < dead:
            raise RuntimeError(f"Insufficient {name}: dead-volume boundary violated")
    if 6 * 24 > 192:
        raise RuntimeError("Insufficient tips")
    if max(50 + 80 + 120, 180, 40) > 2000:
        raise RuntimeError("Processing-well capacity exceeded")
    if 250 * 24 + 180 * 2 * 24 > 195000:
        raise RuntimeError("Liquid-waste capacity exceeded")
    for volume in (30, 40, 80, 100, 120, 180, 250):
        if not 20 <= volume <= 300:
            raise RuntimeError(f"P300 range violation: {volume} uL")


def run(protocol: protocol_api.ProtocolContext) -> None:
    _preflight()

    reagents = protocol.load_labware("nest_12_reservoir_15ml", "1", "reagents")
    tiprack_1 = protocol.load_labware("opentrons_96_tiprack_300ul", "4", "tips_1")
    magnetic_module = protocol.load_module("magnetic module gen2", "5")
    processing_plate = magnetic_module.load_labware(
        "nest_96_wellplate_2ml_deep", "processing_plate"
    )
    liquid_waste = protocol.load_labware("nest_1_reservoir_195ml", "6", "liquid_waste")
    tiprack_2 = protocol.load_labware("opentrons_96_tiprack_300ul", "7", "tips_2")
    pipette = protocol.load_instrument(
        "p300_single_gen2", "right", tip_racks=[tiprack_1, tiprack_2]
    )

    lysis = reagents["A1"]
    beads = reagents["A2"]
    wash = reagents["A3"]
    elution = reagents["A4"]
    waste = liquid_waste["A1"]
    samples = [processing_plate[name] for name in SAMPLE_WELLS]

    protocol.comment("Stage lysis: fresh tip per sample; add 80 uL and mix 5 x 100 uL")
    for sample in samples:
        pipette.pick_up_tip()
        pipette.aspirate(80, lysis.bottom(2))
        pipette.dispense(80, sample.bottom(5))
        pipette.mix(5, 100, sample.bottom(2))
        pipette.drop_tip()
    protocol.delay(minutes=5, msg="Lysis incubation")

    protocol.comment("Stage beads: resuspend source with each clean tip, add 120 uL, mix 10 x 180 uL")
    for sample in samples:
        pipette.pick_up_tip()
        pipette.mix(5, 200, beads.bottom(2))
        pipette.aspirate(120, beads.bottom(2))
        pipette.dispense(120, sample.bottom(5))
        pipette.mix(10, 180, sample.bottom(2))
        pipette.drop_tip()
    protocol.delay(minutes=5, msg="Bead-binding incubation")
    magnetic_module.engage(height_from_base=6.5)
    protocol.delay(minutes=7, msg="Magnetic separation")

    protocol.comment("Stage supernatant: fresh tip per sample; remove 250 uL to waste")
    original_aspirate_rate = pipette.flow_rate.aspirate
    pipette.flow_rate.aspirate = 30
    for sample in samples:
        pipette.pick_up_tip()
        pipette.aspirate(250, sample.bottom(0.5))
        pipette.dispense(250, waste.top(-5))
        pipette.drop_tip()

    protocol.comment("Stage wash 1: retain each sample's tip through its 30-second wait and removal")
    for sample in samples:
        pipette.pick_up_tip()
        pipette.aspirate(180, wash.bottom(2))
        pipette.dispense(180, sample.top(-5))
        protocol.delay(seconds=30, msg="Wash 1 dwell; tip remains attached")
        pipette.aspirate(180, sample.bottom(0.5))
        pipette.dispense(180, waste.top(-5))
        pipette.drop_tip()

    protocol.comment("Stage wash 2: retain each sample's tip through its 30-second wait and removal")
    for sample in samples:
        pipette.pick_up_tip()
        pipette.aspirate(180, wash.bottom(2))
        pipette.dispense(180, sample.top(-5))
        protocol.delay(seconds=30, msg="Wash 2 dwell; tip remains attached")
        pipette.aspirate(180, sample.bottom(0.5))
        pipette.dispense(180, waste.top(-5))
        pipette.drop_tip()
    pipette.flow_rate.aspirate = original_aspirate_rate

    protocol.delay(minutes=2, msg="Air dry while magnet remains engaged")
    magnetic_module.disengage()

    protocol.comment("Stage elution: fresh tip per sample; add 40 uL and mix 10 x 30 uL")
    for sample in samples:
        pipette.pick_up_tip()
        pipette.aspirate(40, elution.bottom(2))
        pipette.dispense(40, sample.bottom(3))
        pipette.mix(10, 30, sample.bottom(1))
        pipette.drop_tip()

