from opentrons import protocol_api

metadata = {
    "protocolName": "Combinatorial DNA Assembly",
    "description": """This protocol takes a master worklist generated from the assembly designer package and performs combinatorial DNA assembly reactions accordingly.\n\n
    Deck Layout:\n
    A2: Water trough with water in [A1],\n
    B2: mtp_dna_stocks,\n
    C2: mtp_mastermix,\n
    D2: mtp_source,\n
    reaction_plate: loaded directly in thermocycler""",
    "author": "SynCell Group - Dominic Kösters",
    "source": "OpentronsAI",
}

requirements = {"robotType": "Flex", "apiLevel": "2.27"}


def add_parameters(parameters: protocol_api.Parameters):
    parameters.add_csv_file(
        variable_name="master_worklist",
        display_name="Master Worklist CSV",
        description="Upload master_worklist.csv from ADesigner",
    )

    parameters.add_int(
        variable_name="gga_cycles",
        display_name="GGA Cycles",
        default=25,
        minimum=1,
        maximum=35,
        description="Number of Golden Gate Assembly cycles",
    )

    parameters.add_int(
        variable_name="gga_reaction_volume",
        display_name="GGA Volume (µL)",
        default=30,
        minimum=20,
        maximum=50,
        description="Total reaction volume per well",
    )

    parameters.add_str(
        variable_name="remove_dna_plate",
        display_name="Remove DNA Plate After Phase 1",
        default="no",
        choices=[
            {"display_name": "Yes - Pause to remove plate", "value": "yes"},
            {"display_name": "No - Continue without pause", "value": "no"},
        ],
        description="Pause after DNA dilution to store DNA stocks plate in freezer",
    )


def run(protocol: protocol_api.ProtocolContext):

    # ========================================
    # PARSE MASTER WORKLIST
    # ========================================
    def parse_worklist(csv_data):
        """Parse auto-generated master worklist with error handling."""
        rows = csv_data.parse_as_csv()

        if not rows or len(rows) < 2:
            raise ValueError("❌ CSV file is empty or contains only headers\n\nPlease verify your master worklist CSV has at least one data row.")

        header = rows[0]

        # Validate required columns exist
        required_columns = ["Well Nr. source plate", "Well Nr. destination plate", "Volumen", "Source_Worklist"]

        missing_columns = [col for col in required_columns if col not in header]
        if missing_columns:
            raise ValueError(f"❌ Missing required columns: {', '.join(missing_columns)}\n\nRequired columns: 'Well Nr. source plate', 'Well Nr. destination plate', 'Volumen', 'Source_Worklist'")

        # Get column indices
        idx = {col: header.index(col) for col in required_columns}
        idx["Part"] = header.index("Part") if "Part" in header else None

        # Group transfers by Source_Worklist
        worklists = {}
        for row_num, row in enumerate(rows[1:], start=2):
            try:
                wl_name = row[idx["Source_Worklist"]].strip()
                if wl_name not in worklists:
                    worklists[wl_name] = []

                worklists[wl_name].append({
                    "source": row[idx["Well Nr. source plate"]].strip(),
                    "destination": row[idx["Well Nr. destination plate"]].strip(),
                    "volume": float(row[idx["Volumen"]].strip()),
                    "part": row[idx["Part"]].strip() if idx["Part"] else "",
                })
            except (IndexError, ValueError, KeyError) as e:
                raise ValueError(f"❌ Error parsing row {row_num}: {e}\n\nPlease verify row {row_num} has valid data in all required columns.")

        if not worklists:
            raise ValueError("❌ No valid worklist data found in CSV\n\nPlease verify your CSV contains valid transfer data.")

        return worklists

    # ========================================
    # TIP AVAILABILITY CHECK
    # ========================================
    def check_tip_availability(worklists_dict):
        """Calculate and validate P50 tip requirements before protocol starts."""
        total_p50 = 0
        total_p1000 = 0

        def count_tips(worklist, new_tip_mode):
            """Count tips needed for a worklist based on tip strategy."""
            nonlocal total_p50, total_p1000

            if not worklist:
                return

            if new_tip_mode == "once":
                # Only 1 tip needed for entire worklist
                v = worklist[0]["volume"]
                if v <= 50:
                    total_p50 += 1
                else:
                    total_p1000 += 1
            else:  # new_tip == "always"
                # Count tips for each transfer
                for t in worklist:
                    if t["volume"] <= 50:
                        total_p50 += 1
                    else:
                        total_p1000 += 1

        # Count tips for each phase with correct tip strategy
        count_tips(worklists_dict.get("Worklist_DNA", []), "always")
        count_tips(worklists_dict.get("Worklist_W", []), "always")
        count_tips(worklists_dict.get("Worklist_PM_parts", []), "always")
        count_tips(worklists_dict.get("Worklist_GG_water", []), "once")
        count_tips(worklists_dict.get("Worklist_MM", []), "once")
        count_tips(worklists_dict.get("Worklist_PM", []), "always")
        count_tips(worklists_dict.get("Worklist_P", []), "always")
        count_tips(worklists_dict.get("Worklist_RBS", []), "always")
        count_tips(worklists_dict.get("Worklist_Tags", []), "always")
        count_tips(worklists_dict.get("Worklist_SP", []), "always")
        count_tips(worklists_dict.get("Worklist_GoI", []), "always")

        # Calculate available tips
        available_p50 = 4 * 96
        available_p1000 = 96

        # Display results
        protocol.comment("=" * 50)
        protocol.comment("🔍 TIP AVAILABILITY CHECK")
        protocol.comment("=" * 50)
        protocol.comment(f"P50 Tips:   Required={total_p50}, Available={available_p50}")
        protocol.comment(f"P1000 Tips: Required={total_p1000}, Available={available_p1000}")

        # Validate sufficiency
        if total_p50 <= available_p50 and total_p1000 <= available_p1000:
            protocol.comment("✅ Sufficient tips loaded")
            protocol.comment("=" * 50)
        else:
            protocol.comment("=" * 50)
            if total_p50 > available_p50:
                additional_racks = ((total_p50 - available_p50) // 96) + 1
                raise ValueError(f"❌ Insufficient P50 tips! Required: {total_p50}, Available: {available_p50}. Please load {additional_racks} additional P50 tip rack(s).")
            if total_p1000 > available_p1000:
                additional_racks = ((total_p1000 - available_p1000) // 96) + 1
                raise ValueError(f"❌ Insufficient P1000 tips! Required: {total_p1000}, Available: {available_p1000}. Please load {additional_racks} additional P1000 tip rack(s).")

    # ========================================
    # TIP RACK REPLACEMENT LOGIC
    # ========================================
    def replace_tip_racks():
        """Replace depleted tip racks with fresh ones using gripper."""
        try:
            protocol.comment("🔄 Replacing tip racks with gripper...")
            protocol.move_labware(tipracks_50[0], "B4", use_gripper=True)
            protocol.comment("  ✓ Moved B3 → B4")
            protocol.move_labware(tipracks_50_storage[0], "B3", use_gripper=True)
            protocol.comment("  ✓ Moved C4 → B3")
            protocol.move_labware(tipracks_50[0], "C4", use_gripper=True)
            protocol.comment("  ✓ Moved B4 → C4")
            protocol.move_labware(tipracks_50[1], "B4", use_gripper=True)
            protocol.comment("  ✓ Moved C3 → B4")
            protocol.move_labware(tipracks_50_storage[1], "C3", use_gripper=True)
            protocol.comment("  ✓ Moved D4 → C3")
            protocol.move_labware(tipracks_50[1], "D4", use_gripper=True)
            protocol.comment("  ✓ Moved B4 → D4")

            # Update tip rack references
            tipracks_50[0] = tipracks_50_storage[0]
            tipracks_50[1] = tipracks_50_storage[1]
            p50.reset_tipracks()
            protocol.comment("✅ Tip rack replacement complete")

        except Exception as e:
            protocol.comment(f"⚠️ Automated tip rack replacement failed: {e}")
            protocol.pause("⚠️ Tip rack replacement failed. Please manually:\n1. Remove depleted tip racks from B3 and C3\n2. Place fresh tip racks in B3 and C3\n3. Verify racks are properly seated\n4. Press Continue to resume protocol")
            tipracks_50[0] = tipracks_50_storage[0]
            tipracks_50[1] = tipracks_50_storage[1]
            p50.reset_tipracks()
            protocol.comment("✅ Manual tip rack replacement complete, resuming protocol")

    # ========================================
    # CALCULATE ACCUMULATED VOLUME PER WELL
    # ========================================
    def calculate_well_volumes(worklists_dict):
        """Calculate accumulated volume per destination well across all worklists."""
        well_volumes = {}
        for wl_name, worklist in worklists_dict.items():
            for transfer in worklist:
                dest_well = transfer["destination"]
                volume = transfer["volume"]
                if dest_well not in well_volumes:
                    well_volumes[dest_well] = 0
                well_volumes[dest_well] += volume
        return well_volumes

    # ========================================
    # TRANSFER LOGIC
    # ========================================
    def select_pipette(volume):
        return p50 if volume <= 50 else p1000

    def set_pipette_volume_mode(pip, volume):
        """Set appropriate volume mode for Flex pipettes based on transfer volume."""
        if pip == p50:
            # P50 has lowVolume mode for 1-5 µL (optimized for low volumes)
            if volume <= 5:
                pip.configure_for_volume(5)
            else:
                pip.configure_for_volume(50)

    def execute_worklist(worklist, source_plate, dest_plate, name="", mix=False, new_tip="always", well_volumes=None, mix_on_last_only=False):
        """Execute transfers from worklist with dynamic mix volume calculation and optimized flow rates."""
        if not worklist:
            protocol.comment(f"⚠️ {name}: empty - skipping")
            return

        protocol.comment(f"▶ {name}: {len(worklist)} transfers (tip={new_tip}, mix={mix})")

        # Validate single pipette for 'once' mode
        if new_tip == "once":
            pip = select_pipette(worklist[0]["volume"])
            if any(select_pipette(t["volume"]) != pip for t in worklist):
                raise ValueError(f"{name}: new_tip='once' requires same pipette for all transfers")
            pip.pick_up_tip()

        # Group transfers by destination well if mix_on_last_only is True
        if mix_on_last_only:
            dest_last_transfer = {}
            for i, t in enumerate(worklist):
                dest_last_transfer[t["destination"]] = i

        # Execute transfers
        for i, t in enumerate(worklist):
            pip = select_pipette(t["volume"])

            # Set appropriate volume mode for accurate low-volume transfers
            set_pipette_volume_mode(pip, t["volume"])

            # Check if we need to replace tip racks (only for P50)
            if pip == p50 and not pip.has_tip and len(pip.tip_racks[0].wells()) == 0 and len(pip.tip_racks[1].wells()) == 0:
                replace_tip_racks()

            # Determine source
            source = source_plate[t["source"]] if hasattr(source_plate, "wells") else source_plate
            dest = dest_plate[t["destination"]]

            # Transfer
            if new_tip == "always":
                pip.pick_up_tip()

            pip.aspirate(t["volume"], source)
            pip.dispense(t["volume"], dest)

            # Determine if we should mix
            should_mix = mix and (not mix_on_last_only or i == dest_last_transfer[t["destination"]])

            if should_mix:
                # Dynamic mix-volume-calculation based on accumulated volume in well
                if well_volumes and t["destination"] in well_volumes:
                    accumulated_volume = well_volumes[t["destination"]]
                    mix_vol = max(10, min(accumulated_volume * 0.7, pip.max_volume))
                else:
                    mix_vol = max(10, min(protocol.params.gga_reaction_volume * 0.7, pip.max_volume))

                # Save original flow rates
                original_aspirate_rate = pip.flow_rate.aspirate
                original_dispense_rate = pip.flow_rate.dispense
                original_blow_out_rate = pip.flow_rate.blow_out

                # Optimized flow rates for mixing (70% of standard flow rate)
                pip.flow_rate.aspirate = original_aspirate_rate * 0.7
                pip.flow_rate.dispense = original_dispense_rate * 0.7
                pip.flow_rate.blow_out = original_blow_out_rate * 0.1

                # Mix 5 times, 1 mm above bottom of plate for better mixing
                pip.mix(5, mix_vol, dest.bottom(1))

                # Restore original flow rates
                pip.flow_rate.aspirate = original_aspirate_rate
                pip.flow_rate.dispense = original_dispense_rate
                pip.flow_rate.blow_out = original_blow_out_rate

            if new_tip == "always":
                pip.drop_tip()

        if new_tip == "once":
            pip.drop_tip()

    # ========================================
    # LABWARE SETUP
    # ========================================
    reaction_plate_type = "biorad_96_wellplate_200ul_pcr"
    prep_plate_type = "nunc_96_wellplate_450ul"

    # Labware
    water_reservoir = protocol.load_labware("nest_12_reservoir_15ml", "A2")
    mtp_dna_stocks = protocol.load_labware(prep_plate_type, "B2")
    mtp_mastermix = protocol.load_labware(prep_plate_type, "C2")
    mtp_source = protocol.load_labware(prep_plate_type, "D2")

    # Modules - Load reaction plate directly in thermocycler
    tc_mod = protocol.load_module("thermocyclerModuleV2")
    mtp_reaction = tc_mod.load_labware(reaction_plate_type)

    # Pipettes - Active tip racks in B3 and C3
    tipracks_50 = [protocol.load_labware("opentrons_flex_96_tiprack_50ul", "B3"), protocol.load_labware("opentrons_flex_96_tiprack_50ul", "C3")]

    # Storage tip racks in C4 and D4
    tipracks_50_storage = [protocol.load_labware("opentrons_flex_96_tiprack_50ul", "C4"), protocol.load_labware("opentrons_flex_96_tiprack_50ul", "D4")]

    tiprack_1000 = protocol.load_labware("opentrons_flex_96_tiprack_1000ul", "D3")

    p50 = protocol.load_instrument("flex_1channel_50", "left", tip_racks=tipracks_50)
    p1000 = protocol.load_instrument("flex_1channel_1000", "right", tip_racks=[tiprack_1000])

    trash = protocol.load_trash_bin("A3")

    # ========================================
    # PARSE & VALIDATE WORKLISTS
    # ========================================
    protocol.comment("📋 Parsing master worklist...")

    wl = parse_worklist(protocol.params.master_worklist)

    # Extract worklists
    DNA = wl.get("Worklist_DNA", [])
    WATER = wl.get("Worklist_W", [])
    PM_PARTS = wl.get("Worklist_PM_parts", [])
    GG_WATER = wl.get("Worklist_GG_water", [])
    MM = wl.get("Worklist_MM", [])
    PM = wl.get("Worklist_PM", [])
    PROMOTER = wl.get("Worklist_P", [])
    RBS = wl.get("Worklist_RBS", [])
    TAGS = wl.get("Worklist_Tags", [])
    SP = wl.get("Worklist_SP", [])
    GOI = wl.get("Worklist_GoI", [])

    # Print found worklists
    protocol.comment("=" * 50)
    protocol.comment("📋 WORKLISTS FOUND:")
    protocol.comment("=" * 50)
    for wl_name in sorted(wl.keys()):
        protocol.comment(f"  ✓ {wl_name}: {len(wl[wl_name])} transfers")
    protocol.comment("=" * 50)

    # Summary
    protocol.comment(f"DNA:{len(DNA)} | Water:{len(WATER)} | PM_parts:{len(PM_PARTS)} | GG_water:{len(GG_WATER)} | MM:{len(MM)} | PM:{len(PM)} | Promoter:{len(PROMOTER)} | RBS:{len(RBS)} | Tags:{len(TAGS)} | SP:{len(SP)} | GoI:{len(GOI)}")

    # CHECK TIP AVAILABILITY BEFORE STARTING
    check_tip_availability(wl)

    # CALCULATE ACCUMULATED VOLUMES PER WELL
    phase1_volumes = {}
    for transfer in DNA + WATER:
        dest = transfer["destination"]
        if dest not in phase1_volumes:
            phase1_volumes[dest] = 0
        phase1_volumes[dest] += transfer["volume"]

    phase2_volumes = phase1_volumes.copy()
    for transfer in PM_PARTS:
        dest = transfer["destination"]
        if dest not in phase2_volumes:
            phase2_volumes[dest] = 0
        phase2_volumes[dest] += transfer["volume"]

    phase3_volumes = {}
    for transfer in GG_WATER + MM + PM + PROMOTER + RBS + TAGS + SP + GOI:
        dest = transfer["destination"]
        if dest not in phase3_volumes:
            phase3_volumes[dest] = 0
        phase3_volumes[dest] += transfer["volume"]

    # ========================================
    # PHASE 1: DNA DILUTION
    # ========================================
    protocol.comment("🧬 PHASE 1: DNA DILUTION")
    execute_worklist(WATER, water_reservoir["A1"], mtp_source, "Water", mix=False, new_tip="once")
    execute_worklist(DNA, mtp_dna_stocks, mtp_source, "DNA", mix=True, new_tip="always", well_volumes=phase1_volumes)

    # Optional pause to remove DNA stocks plate
    if protocol.params.remove_dna_plate == "yes":
        protocol.comment("=" * 50)
        protocol.comment("⏸️ PAUSING FOR DNA PLATE REMOVAL")
        protocol.comment("=" * 50)
        protocol.pause("Phase 1 complete. Remove mtp_dna_stocks plate from slot B2 and store in freezer. Press continue when ready to proceed.")
        protocol.comment("✅ Resuming protocol after DNA plate removal")

    # ========================================
    # PHASE 2: PLASMID MIX PREPARATION
    # ========================================
    protocol.comment("🧪 PHASE 2: PLASMID MIX")
    execute_worklist(PM_PARTS, mtp_source, mtp_source, "PM_parts", mix=True, new_tip="always", well_volumes=phase2_volumes, mix_on_last_only=True)

    # ========================================
    # PHASE 3: ASSEMBLY SETUP (OPTIMIZED)
    # ========================================
    protocol.comment("🔬 PHASE 3: ASSEMBLY SETUP")

    # Open thermocycler and pre-cool to first cycle temperature
    tc_mod.open_lid()
    tc_mod.set_block_temperature(16)
    protocol.comment("🌡️ Thermocycler pre-cooling during pipetting")

    # Pipette directly into thermocycler plate
    execute_worklist(GG_WATER, water_reservoir["A1"], mtp_reaction, "GG_water", mix=False, new_tip="once")
    execute_worklist(PM, mtp_source, mtp_reaction, "PlasmidMix", mix=False, new_tip="once")
    execute_worklist(PROMOTER, mtp_source, mtp_reaction, "Promoter", mix=False, new_tip="always")
    execute_worklist(RBS, mtp_source, mtp_reaction, "RBS", mix=False, new_tip="always")
    execute_worklist(SP, mtp_source, mtp_reaction, "SP", mix=False, new_tip="always")
    execute_worklist(GOI, mtp_source, mtp_reaction, "GoI", mix=True, new_tip="always")
    execute_worklist(TAGS, mtp_source, mtp_reaction, "Tags", mix=False, new_tip="always")
    execute_worklist(MM, mtp_mastermix, mtp_reaction, "Mastermix", mix=True, new_tip="always", well_volumes=phase3_volumes)

    # ========================================
    # THERMOCYCLER PROGRAM
    # ========================================
    protocol.comment("🌡️ THERMOCYCLER CYCLING")

    # Thermocycler control with error handling
    try:
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(95)
    except Exception as e:
        protocol.comment(f"⚠️ Thermocycler operation failed: {e}")
        protocol.pause("⚠️ Thermocycler error detected. Please:\n1. Check thermocycler lid is not obstructed\n2. Verify plate is properly seated in thermocycler\n3. Manually close lid if needed\n4. Press Continue to retry")
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(95)

    # Assembly cycling
    tc_mod.execute_profile(steps=[{"temperature": 37, "hold_time_seconds": 90}, {"temperature": 16, "hold_time_seconds": 180}], repetitions=protocol.params.gga_cycles, block_max_volume=protocol.params.gga_reaction_volume)

    # Post-cycling
    tc_mod.set_block_temperature(50, hold_time_minutes=5, block_max_volume=protocol.params.gga_reaction_volume)
    tc_mod.set_block_temperature(80, hold_time_minutes=10, block_max_volume=protocol.params.gga_reaction_volume)

    # Hold at 15°C and wait for operator
    tc_mod.set_block_temperature(15)
    protocol.pause("Plate is being held at 15°C in the thermocycler. Press continue when ready to finish the protocol.")

    # Deactivate and open after operator continues
    tc_mod.deactivate_lid()
    tc_mod.deactivate_block()
    tc_mod.open_lid()

    protocol.comment("✅ ASSEMBLY COMPLETED")