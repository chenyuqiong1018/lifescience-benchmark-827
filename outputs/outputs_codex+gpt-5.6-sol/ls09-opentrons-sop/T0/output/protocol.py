from opentrons import protocol_api

metadata = {
    "protocolName": "24-sample magnetic-bead cleanup (frozen SOP)",
    "description": "Six-tip-per-sample OT-2 translation with auditable preflight checks",
    "author": "Codex benchmark execution",
    "apiLevel": "2.16",
}

requirements = {"robotType": "OT-2"}

WELL_NAMES = tuple(
    f"{row}{column}" for column in range(1, 7) for row in ("A", "B", "C", "D")
)
PIPETTE_MIN_UL = 20
PIPETTE_MAX_UL = 300
TIP_INVENTORY = 192
PROCESSING_CAPACITY_UL = 2000
WASTE_CAPACITY_UL = 195000


def preflight() -> None:
    if len(WELL_NAMES) != 24 or len(set(WELL_NAMES)) != 24:
        raise RuntimeError("Frozen sample map must resolve to 24 unique wells")

    reagent_checks = (
        ("lysis_buffer", 2500, 300, 24 * 80),
        ("bead_suspension", 3500, 300, 24 * 120),
        ("wash_buffer", 9000, 300, 24 * 2 * 180),
        ("elution_buffer", 1500, 300, 24 * 40),
    )
    for name, initial, dead, consumption in reagent_checks:
        if initial - consumption < dead:
            raise RuntimeError(f"{name} would cross its declared dead volume")

    if 24 * 6 > TIP_INVENTORY:
        raise RuntimeError("The six-tip-per-sample policy exceeds inventory")
    if 24 * (250 + 2 * 180) > WASTE_CAPACITY_UL:
        raise RuntimeError("Liquid-waste capacity would be exceeded")
    if 50 + 80 + 120 > PROCESSING_CAPACITY_UL:
        raise RuntimeError("Processing-well capacity would be exceeded")
    for volume in (30, 40, 80, 100, 120, 180, 250):
        if volume < PIPETTE_MIN_UL or volume > PIPETTE_MAX_UL:
            raise RuntimeError(f"{volume} uL is outside the declared P300 range")


def add_and_mix(pipette, source, destination, add_volume, mix_repetitions, mix_volume):
    pipette.pick_up_tip()
    pipette.aspirate(add_volume, source.bottom(2))
    pipette.dispense(add_volume, destination.bottom(5))
    pipette.mix(mix_repetitions, mix_volume, destination.bottom(2))
    pipette.drop_tip()


def remove_to_waste(pipette, source, waste, volume):
    pipette.pick_up_tip()
    pipette.aspirate(volume, source.bottom(0.5))
    pipette.dispense(volume, waste.top(-5))
    pipette.drop_tip()


def wash_once(protocol, pipette, sample, wash_source, waste, label):
    pipette.pick_up_tip()
    pipette.aspirate(180, wash_source.bottom(2))
    pipette.dispense(180, sample.top(-5))
    protocol.delay(seconds=30, msg=f"{label}: tip retained for this sample")
    pipette.aspirate(180, sample.bottom(0.5))
    pipette.dispense(180, waste.top(-5))
    pipette.drop_tip()


def run(protocol: protocol_api.ProtocolContext) -> None:
    preflight()

    reagents = protocol.load_labware("nest_12_reservoir_15ml", "1", "reagents")
    tips_slot_4 = protocol.load_labware("opentrons_96_tiprack_300ul", "4", "tips_slot_4")
    magnet = protocol.load_module("magnetic module gen2", "5")
    plate = magnet.load_labware("nest_96_wellplate_2ml_deep", "processing_plate")
    waste_labware = protocol.load_labware("nest_1_reservoir_195ml", "6", "liquid_waste")
    tips_slot_7 = protocol.load_labware("opentrons_96_tiprack_300ul", "7", "tips_slot_7")
    pipette = protocol.load_instrument(
        "p300_single_gen2", "right", tip_racks=[tips_slot_4, tips_slot_7]
    )

    lysis = reagents["A1"]
    beads = reagents["A2"]
    wash = reagents["A3"]
    elution = reagents["A4"]
    waste = waste_labware["A1"]
    samples = [plate[name] for name in WELL_NAMES]

    protocol.comment("Lysis: 24 fresh tips")
    for sample in samples:
        add_and_mix(pipette, lysis, sample, 80, 5, 100)
    protocol.delay(minutes=5, msg="Lysis incubation")

    protocol.comment("Beads: resuspend source before each withdrawal; 24 fresh tips")
    for sample in samples:
        pipette.pick_up_tip()
        pipette.mix(10, 180, beads.bottom(2))
        pipette.aspirate(120, beads.bottom(2))
        pipette.dispense(120, sample.bottom(5))
        pipette.mix(10, 180, sample.bottom(2))
        pipette.drop_tip()
    protocol.delay(minutes=5, msg="Bead-binding incubation")

    magnet.engage()
    protocol.delay(minutes=7, msg="Magnetic separation")

    protocol.comment("Supernatant: 24 fresh tips")
    for sample in samples:
        remove_to_waste(pipette, sample, waste, 250)

    protocol.comment("Wash 1: one fresh tip per sample, retained through dwell and removal")
    for sample in samples:
        wash_once(protocol, pipette, sample, wash, waste, "Wash 1")

    protocol.comment("Wash 2: one fresh tip per sample, retained through dwell and removal")
    for sample in samples:
        wash_once(protocol, pipette, sample, wash, waste, "Wash 2")

    protocol.delay(minutes=2, msg="Air dry with magnet engaged")
    magnet.disengage()

    protocol.comment("Elution: 24 fresh tips")
    for sample in samples:
        add_and_mix(pipette, elution, sample, 40, 10, 30)

