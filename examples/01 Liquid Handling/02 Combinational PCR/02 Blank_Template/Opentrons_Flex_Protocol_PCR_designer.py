from opentrons import protocol_api

metadata = {
    "protocolName": "Combinatorial DNA Assembly - PCR",
    "description": """This protocol takes a master worklist generated from the assembly designer package and performs combinatorial PCR reactions accordingly.\n\n
    Deck Layout:\n
    A2: Water trough with water in [A1],\n
    B2: mtp_primer_stocks,\n
    C2: mtp_source (for primer mixes and mastermix),\n
    D2: mtp_template_DNA,\n
    reaction_plate: loaded directly in thermocycler""",
    "author": "SynCell Group - Dominic Kösters",
    "source": "OpentronsAI",
}

requirements = {"robotType": "Flex", "apiLevel": "2.27"}


def add_parameters(parameters: protocol_api.Parameters):
    parameters.add_csv_file(
        variable_name="master_worklist",
        display_name="Master Worklist CSV",
        description="Upload master_worklist.csv from Assembly Designer",
    )

    parameters.add_int(
        variable_name="pcr_cycles",
        display_name="PCR Cycles",
        default=35,
        minimum=2,
        maximum=40,
        description="Number of PCR cycles",
    )

    parameters.add_int(
        variable_name="pcr_reaction_volume",
        display_name="PCR reaction volume",
        default=50,
        minimum=10,
        maximum=100,
        description="Total reaction volume per PCR reaction",
    )

    parameters.add_int(
        variable_name="pcr_annealing_temp",
        display_name="Annealing temperature",
        description="Annealing temperature for PCR",
        default=60,
        minimum=50,
        maximum=72,
        unit="°C",
    )

    parameters.add_int(
        variable_name="pcr_elongation_time",
        display_name="Elongation time",
        description="Elongation time for PCR",
        default=240,
        minimum=10,
        maximum=420,
        unit="s",
    )


def run(protocol: protocol_api.ProtocolContext):

    # ========================================
    # PARSE MASTER WORKLIST
    # ========================================
    def parse_worklist(csv_data):
        """Parse auto-generated master worklist with error handling."""
        rows = csv_data.parse_as_csv()

        if not rows or len(rows) < 2:
            raise ValueError(
                "❌ CSV file is empty or contains only headers\n\n"
                "Please verify your master worklist CSV has at least one data row."
            )

        header = rows[0]

        # Validate required columns exist
        required_columns = ["Well Nr. source plate", "Well Nr. destination plate", "Volumen", "Source_Worklist"]

        missing_columns = [col for col in required_columns if col not in header]
        if missing_columns:
            raise ValueError(
                f"❌ Missing required columns: {', '.join(missing_columns)}\n\n"
                f"Required columns: 'Well Nr. source plate', 'Well Nr. destination plate', 'Volumen', 'Source_Worklist'"
            )

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
                raise ValueError(
                    f"❌ Error parsing row {row_num}: {e}\n\n"
                    f"Please verify row {row_num} has valid data in all required columns."
                )

        if not worklists:
            raise ValueError(
                "❌ No valid worklist data found in CSV\n\n"
                "Please verify your CSV contains valid transfer data."
            )

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

        # Count tips for each worklist with correct tip strategy
        # Phase 1: Primer mix preparation in mtp_source
        count_tips(worklists_dict.get("Worklist_Primer1", []), "always")
        count_tips(worklists_dict.get("Worklist_Primer2", []), "always")

        # Phase 2: PCR reaction setup
        count_tips(worklists_dict.get("Worklist_Water", []), "once")
        count_tips(worklists_dict.get("Worklist_MM", []), "once")
        count_tips(worklists_dict.get("Worklist_PCR_Primer", []), "always")
        count_tips(worklists_dict.get("Worklist_Template", []), "always")

        # Calculate available tips (4 racks × 96 tips each for P50)
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
                raise ValueError(
                    f"❌ Insufficient P50 tips! Required: {total_p50}, Available: {available_p50}. "
                    f"Please load {additional_racks} additional P50 tip rack(s)."
                )

            if total_p1000 > available_p1000:
                additional_racks = ((total_p1000 - available_p1000) // 96) + 1
                raise ValueError(
                    f"❌ Insufficient P1000 tips! Required: {total_p1000}, Available: {available_p1000}. "
                    f"Please load {additional_racks} additional P1000 tip rack(s)."
                )

    # ========================================
    # TIP RACK REPLACEMENT LOGIC
    # ========================================
    def replace_tip_racks():
        """Replace depleted tip racks with fresh ones using gripper."""
        try:
            protocol.comment("🔄 Replacing tip racks with gripper...")

            # Step 1: Move B3 to B4 (temporary storage)
            protocol.move_labware(tipracks_50[0], "B4", use_gripper=True)
            protocol.comment("  ✓ Moved B3 → B4")

            # Step 2: Move C4 to B3 (fresh rack to active position)
            protocol.move_labware(tipracks_50_storage[0], "B3", use_gripper=True)
            protocol.comment("  ✓ Moved C4 → B3")

            # Step 3: Move B4 to C4 (used rack to storage)
            protocol.move_labware(tipracks_50[0], "C4", use_gripper=True)
            protocol.comment("  ✓ Moved B4 → C4")

            # Step 4: Move C3 to B4 (temporary storage)
            protocol.move_labware(tipracks_50[1], "B4", use_gripper=True)
            protocol.comment("  ✓ Moved C3 → B4")

            # Step 5: Move D4 to C3 (fresh rack to active position)
            protocol.move_labware(tipracks_50_storage[1], "C3", use_gripper=True)
            protocol.comment("  ✓ Moved D4 → C3")

            # Step 6: Move B4 to D4 (used rack to storage)
            protocol.move_labware(tipracks_50[1], "D4", use_gripper=True)
            protocol.comment("  ✓ Moved B4 → D4")

            # Update tip rack references
            tipracks_50[0] = tipracks_50_storage[0]
            tipracks_50[1] = tipracks_50_storage[1]

            # Reset pipette tip tracking
            p50.reset_tipracks()

            protocol.comment("✅ Tip rack replacement complete")

        except Exception as e:
            # Error recovery: pause for manual intervention
            protocol.comment(f"⚠️ Automated tip rack replacement failed: {e}")
            protocol.pause(
                "⚠️ Tip rack replacement failed. Please manually:\n"
                "1. Remove depleted tip racks from B3 and C3\n"
                "2. Place fresh tip racks in B3 and C3\n"
                "3. Verify racks are properly seated\n"
                "4. Press Continue to resume protocol"
            )

            # Update tip rack references to current racks in active positions
            tipracks_50[0] = tipracks_50_storage[0]
            tipracks_50[1] = tipracks_50_storage[1]

            # Reset pipette tip tracking
            p50.reset_tipracks()

            protocol.comment("✅ Manual tip rack replacement complete, resuming protocol")

    # ========================================
    # CALCULATE ACCUMULATED VOLUME PER WELL
    # ========================================
    def calculate_well_volumes(worklists_dict, worklist_names):
        """Calculate accumulated volume per destination well for specified worklists."""
        well_volumes = {}

        for wl_name in worklist_names:
            worklist = worklists_dict.get(wl_name, [])
            for transfer in worklist:
                dest_well = transfer["destination"]
                volume = transfer["volume"]

                if dest_well not in well_volumes:
                    well_volumes[dest_well] = 0
                well_volumes[dest_well] += volume

        return well_volumes

    # ========================================
    # TRANSFER LOGIC WITH LOW VOLUME SUPPORT
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
            # Create a dictionary to track last transfer to each destination
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
                    # Use 70% of accumulated volume in destination well
                    accumulated_volume = well_volumes[t["destination"]]
                    mix_vol = max(10, min(accumulated_volume * 0.7, pip.max_volume))
                else:
                    # Fallback: Use PCR reaction volume
                    mix_vol = max(10, min(protocol.params.pcr_reaction_volume * 0.7, pip.max_volume))

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
    reaction_plate_type = "opentrons_96_wellplate_200ul_pcr_full_skirt"
    prep_plate_type = "nunc_96_wellplate_450ul"

    # Labware
    water_reservoir = protocol.load_labware("nest_12_reservoir_15ml", "A2")
    mtp_primer_stocks = protocol.load_labware(prep_plate_type, "B2")
    mtp_source = protocol.load_labware(prep_plate_type, "C2")  # For primer mixes and mastermix
    mtp_template_DNA = protocol.load_labware(prep_plate_type, "D2")

    # Modules - Load reaction plate directly in thermocycler
    tc_mod = protocol.load_module("thermocyclerModuleV2")
    mtp_reaction = tc_mod.load_labware(reaction_plate_type)

    # Pipettes - Active tip racks in B3 and C3
    tipracks_50 = [
        protocol.load_labware("opentrons_flex_96_tiprack_50ul", "B3"),
        protocol.load_labware("opentrons_flex_96_tiprack_50ul", "C3"),
    ]

    # Storage tip racks in C4 and D4
    tipracks_50_storage = [
        protocol.load_labware("opentrons_flex_96_tiprack_50ul", "C4"),
        protocol.load_labware("opentrons_flex_96_tiprack_50ul", "D4"),
    ]

    tiprack_1000 = protocol.load_labware("opentrons_flex_96_tiprack_1000ul", "D3")

    p50 = protocol.load_instrument("flex_1channel_50", "left", tip_racks=tipracks_50)
    p1000 = protocol.load_instrument("flex_1channel_1000", "right", tip_racks=[tiprack_1000])

    trash = protocol.load_trash_bin("A3")

    # ========================================
    # PARSE & VALIDATE WORKLISTS
    # ========================================
    protocol.comment("📋 Parsing master worklist...")

    wl = parse_worklist(protocol.params.master_worklist)

    # Extract worklists based on new structure
    WATER = wl.get("Worklist_Water", [])
    MM = wl.get("Worklist_MM", [])
    PCR_PRIMER = wl.get("Worklist_PCR_Primer", [])
    PRIMER1 = wl.get("Worklist_Primer1", [])
    PRIMER2 = wl.get("Worklist_Primer2", [])
    TEMPLATE = wl.get("Worklist_Template", [])

    # Print found worklists
    protocol.comment("=" * 50)
    protocol.comment("📋 WORKLISTS FOUND:")
    protocol.comment("=" * 50)
    for wl_name in sorted(wl.keys()):
        protocol.comment(f"  ✓ {wl_name}: {len(wl[wl_name])} transfers")
    protocol.comment("=" * 50)

    # Summary
    protocol.comment(
        f"Water:{len(WATER)} | MM:{len(MM)} | PCR_Primer:{len(PCR_PRIMER)} | "
        f"Primer1:{len(PRIMER1)} | Primer2:{len(PRIMER2)} | Template:{len(TEMPLATE)}"
    )

    # CHECK TIP AVAILABILITY BEFORE STARTING
    check_tip_availability(wl)

    # CALCULATE ACCUMULATED VOLUMES PER WELL
    # Phase 1: Primer mix preparation (in mtp_source)
    phase1_volumes = calculate_well_volumes(wl, ["Worklist_Primer1", "Worklist_Primer2"])

    # Phase 2: PCR reaction setup (in mtp_reaction)
    phase2_volumes = calculate_well_volumes(wl, ["Worklist_Water", "Worklist_MM", "Worklist_PCR_Primer", "Worklist_Template"])

    # ========================================
    # PHASE 1: PRIMER MIX PREPARATION
    # ========================================
    protocol.comment("🧬 PHASE 1: PRIMER MIX PREPARATION")
    protocol.comment("Creating primer mixes in mtp_source plate...")

    # Transfer Primer1 from primer stocks to source plate
    execute_worklist(PRIMER1, mtp_primer_stocks, mtp_source, "Primer1 (Forward)", mix=False, new_tip="always")

    # Transfer Primer2 from primer stocks to source plate and mix
    execute_worklist(PRIMER2, mtp_primer_stocks, mtp_source, "Primer2 (Reverse)", mix=True, new_tip="always", well_volumes=phase1_volumes, mix_on_last_only=True)

    protocol.comment("✅ Primer mixes prepared in mtp_source")

    # ========================================
    # PHASE 2: PCR REACTION SETUP (OPTIMIZED)
    # ========================================
    protocol.comment("🧪 PHASE 2: PCR REACTION SETUP")

    # Open thermocycler and pre-cool
    tc_mod.open_lid()
    tc_mod.set_block_temperature(15)  # Pre-cool thermocycler block for pipetting
    protocol.comment("🌡️ Thermocycler pre-cooling during pipetting")

    # Pipette directly into thermocycler plate
    # Add water to reaction plate
    execute_worklist(WATER, water_reservoir["A1"], mtp_reaction, "Water", mix=False, new_tip="once")

    # Add mastermix from source plate (where it was pre-loaded)
    execute_worklist(MM, mtp_source, mtp_reaction, "Mastermix", mix=False, new_tip="once")

    # Add primer mixes from source plate (now containing mixed primers)
    execute_worklist(PCR_PRIMER, mtp_source, mtp_reaction, "Primer Mix", mix=False, new_tip="always")

    # Add template DNA and mix final reaction
    execute_worklist(TEMPLATE, mtp_template_DNA, mtp_reaction, "Template DNA", mix=True, new_tip="always", well_volumes=phase2_volumes)

    protocol.comment("✅ PCR reactions prepared")

    # ========================================
    # THERMOCYCLER PROGRAM
    # ========================================
    protocol.comment("🌡️ PCR CYCLING")

    # Thermocycler error handling
    try:
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(98)
    except Exception as e:
        protocol.comment(f"⚠️ Thermocycler operation failed: {e}")
        protocol.pause(
            "⚠️ Thermocycler error detected. Please:\n"
            "1. Check thermocycler lid is not obstructed\n"
            "2. Verify plate is properly seated in thermocycler\n"
            "3. Manually close lid if needed\n"
            "4. Press Continue to retry"
        )
        # Retry once
        tc_mod.close_lid()
        tc_mod.set_lid_temperature(98)

    # Initial denaturation
    tc_mod.set_block_temperature(temperature=98, hold_time_seconds=30, block_max_volume=protocol.params.pcr_reaction_volume)

    # PCR cycling
    pcr_profile = [
        {"temperature": 98, "hold_time_seconds": 30},
        {"temperature": protocol.params.pcr_annealing_temp, "hold_time_seconds": 10},
        {"temperature": 72, "hold_time_seconds": protocol.params.pcr_elongation_time},
    ]
    tc_mod.execute_profile(steps=pcr_profile, repetitions=protocol.params.pcr_cycles - 1, block_max_volume=protocol.params.pcr_reaction_volume)

    # Final extension
    tc_mod.set_block_temperature(72, hold_time_minutes=2, block_max_volume=protocol.params.pcr_reaction_volume)

    # Hold at 15°C and wait for operator
    tc_mod.set_block_temperature(15)
    protocol.pause("Plate is being held at 15°C in the thermocycler. Press continue when ready to finish the protocol.")

    # Deactivate and open after operator continues
    tc_mod.deactivate_lid()
    tc_mod.deactivate_block()
    tc_mod.open_lid()

    protocol.comment("✅ PCR PROTOCOL COMPLETED. Take out reaction plate from Thermocycler.")