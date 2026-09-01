from opentrons import protocol_api

metadata = {"apiLevel": "2.16", "protocolName": "24-sample magnetic-bead cleanup"}

def run(protocol: protocol_api.ProtocolContext):
    reagents = protocol.load_labware("nest_12_reservoir_15ml", "1")
    tips1 = protocol.load_labware("opentrons_96_tiprack_300ul", "4")
    mag = protocol.load_module("magnetic module gen2", "5")
    plate = mag.load_labware("nest_96_wellplate_2ml_deep")
    waste = protocol.load_labware("nest_1_reservoir_195ml", "6")
    tips2 = protocol.load_labware("opentrons_96_tiprack_300ul", "7")
    p300 = protocol.load_instrument("p300_single_gen2", "right", tip_racks=[tips1, tips2])
    wells = [plate[f"{r}{c}"] for c in range(1, 7) for r in "ABCD"]
    lysis, beads, wash, elution = (reagents["A1"], reagents["A2"], reagents["A3"], reagents["A4"])
    liquid_waste = waste["A1"]

    for well in wells:
        p300.pick_up_tip(); p300.transfer(80, lysis, well, new_tip="never"); p300.mix(5, 100, well); p300.drop_tip()
    protocol.delay(minutes=5)
    for well in wells:
        p300.pick_up_tip(); p300.mix(5, 180, beads); p300.transfer(120, beads, well, new_tip="never"); p300.mix(10, 180, well); p300.drop_tip()
    protocol.delay(minutes=5); mag.engage(); protocol.delay(minutes=7)
    for well in wells:
        p300.pick_up_tip(); p300.transfer(250, well.bottom(1), liquid_waste, new_tip="never"); p300.drop_tip()
    for well in wells:
        p300.pick_up_tip(); p300.transfer(180, wash, well, new_tip="never"); protocol.delay(seconds=30); p300.transfer(180, well.bottom(1), liquid_waste, new_tip="never"); p300.drop_tip()
    for well in wells:
        p300.pick_up_tip(); p300.transfer(180, wash, well, new_tip="never"); protocol.delay(seconds=30); p300.transfer(180, well.bottom(1), liquid_waste, new_tip="never"); p300.drop_tip()
    protocol.delay(minutes=2); mag.disengage()
    for well in wells:
        p300.pick_up_tip(); p300.transfer(40, elution, well, new_tip="never"); p300.mix(10, 30, well); p300.drop_tip()
