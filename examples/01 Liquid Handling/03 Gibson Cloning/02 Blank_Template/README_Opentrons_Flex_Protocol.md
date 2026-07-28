# Opentrons Flex Protocol: Combinatorial DNA Assembly - Gibson

## Overview
This protocol automates combinatorial DNA assembly reactions using the Gibson Assembly method. It processes a master worklist generated from the Assembly Designer package and performs multi-phase liquid handling followed by thermocycler-controlled isothermal assembly reactions.

## Protocol Phases

### Phase 1: DNA Dilution
- Transfers DNA stocks from source plate to preparation plate
- Adds water for dilution
- Optional pause to store DNA stocks plate in freezer

### Phase 2: Plasmid Mix Preparation
- Combines DNA parts to create plasmid mixes
- Performs mixing after final transfer to each destination well

### Phase 3: Assembly Setup
- Adds water, mastermix, and DNA components directly into thermocycler plate
- Mixes final assembly reactions

### Phase 4: Thermocycler Program
- Executes Gibson assembly incubation at 50°C (default: 60 minutes)
- Holds at 15°C until operator retrieval

## Required Materials

### Hardware
- **Robot**: Opentrons Flex
- **Modules**: Thermocycler Module GEN2
- **Pipettes**:
  - Flex 1-Channel 50 µL (left mount)
  - Flex 1-Channel 1000 µL (right mount)

### Labware & Deck Layout
- Indicated when initializing the protocol on the Opentrons

### Consumables
- **Minimum 4 boxes** of Opentrons Flex 50 µL tips (protocol includes automated tip rack replacement)
- **1 box** of Opentrons Flex 1000 µL tips
- PCR-compatible 96-well plate (200 µL capacity)
  - Use **'biorad_96_wellplate_200ul_pcr'** (standard) or 'opentrons_96_wellplate_200ul_pcr_full_skirt' and adjust the definition in the python script if needed

## Protocol Parameters

### Required Upload
- **Master Worklist CSV**: Generated from Assembly Designer package
  - Must contain columns: `Well Nr. source plate`, `Well Nr. destination plate`, `Volumen`, `worklist_name`
  - Optional column: `Part` (for part identification)

### Customizable Parameters
1. **Gibson Incubation Time** (default: 60 minutes)
   - Range: 10-90 minutes
   - Controls duration of isothermal assembly at 50°C

2. **Gibson Reaction Volume** (default: 20 µL)
   - Range: 10-50 µL
   - Total volume per assembly reaction

3. **Remove DNA Plate After Phase 1** (default: No)
   - **Yes**: Protocol pauses after Phase 1 to allow DNA plate storage in freezer
   - **No**: Continues without pause

## Pre-Protocol Checklist

### Before Starting
1. ✅ Prepare master worklist CSV from Assembly Designer
2. ✅ Thaw DNA stocks, Gibson mastermix, and DNA parts on ice
3. ✅ Fill water reservoir (well A1) with nuclease-free water
4. ✅ Load DNA stocks into plate in slot B2
5. ✅ Load Gibson mastermix into plate in slot C2
6. ✅ Load diluted DNA parts into plate in slot D2
7. ✅ Place empty PCR plate directly in thermocycler
8. ✅ Load **4 boxes** of 50 µL tips (slots B3, C3, C4, D4)
9. ✅ Load 1 box of 1000 µL tips (slot D3)
10. ✅ Verify thermocycler lid is open and unobstructed

### Plate Preparation Guidelines
- **DNA Stocks Plate (B2)**: Arrange according to worklist `Well Nr. source plate` for `Worklist_DNA`
- **Mastermix Plate (C2)**: Load Gibson mastermix according to worklist requirements
- **Source Plate (D2)**: Pre-load any required DNA parts according to worklist

## Important Features

### Automated Tip Management
- Protocol automatically calculates required tips before starting
- **Validates tip availability** - will stop if insufficient tips loaded
- **Automated tip rack replacement** using gripper when P50 tips run low
- If automated replacement fails, protocol pauses for manual intervention

### Dynamic Mixing
- Mixing volumes automatically adjust based on accumulated liquid in each well
- Uses 70% of accumulated volume for optimal mixing
- Optimized flow rates (70% of standard) during mixing for better homogeneity

### Error Handling
- CSV validation with detailed error messages
- Thermocycler operation error recovery
- Graceful fallback to manual tip replacement if gripper fails

### Thermocycler Pre-cooling
- Block pre-cools to 15°C during Phase 3 pipetting
- Reduces total protocol time by ~5-10 minutes

## Expected Runtime
- **Phase 1 (DNA Dilution)**: ~10-15 minutes
- **Phase 2 (Plasmid Mix)**: ~5-10 minutes
- **Phase 3 (Assembly Setup)**: ~15-25 minutes
- **Phase 4 (Thermocycler)**: ~60-90 minutes (depends on incubation time)
- **Total**: ~1.5-2.5 hours

## Troubleshooting

### Common Issues

**"Insufficient tips" error**
- Solution: Load additional 50 µL tip racks as indicated in error message
- Protocol calculates exact requirements before starting

**"CSV parsing error"**
- Solution: Verify CSV has all required columns and valid data
- Check for empty rows or missing values
- Ensure worklist names match expected format (e.g., `Worklist_DNA`, `Worklist_MM`)

**Thermocycler lid won't close**
- Solution: Check for obstructions, verify plate is seated properly
- Protocol will pause and allow manual intervention

**Tip rack replacement fails**
- Solution: Protocol automatically pauses for manual replacement
- Remove depleted racks from B3/C3, place fresh racks, press Continue

### Protocol Pauses
The protocol may pause at:
1. **After Phase 1** (if enabled): Remove DNA stocks plate for freezer storage
2. **Thermocycler error**: Manual intervention required
3. **Tip replacement failure**: Manual tip rack replacement needed
4. **After thermocycler hold**: Retrieve completed reactions

## Support & Customization

### Modifying the Protocol
- **Incubation time**: Adjust via runtime parameters before starting (10-90 min range)
- **Reaction volume**: Adjust via runtime parameters (10-50 µL range)
- **Worklist structure**: Must maintain required CSV column format
- **Labware types**: Can be modified in protocol code (lines 280-283)


## Safety Notes
⚠️ **Important**:
- Always verify deck layout matches protocol requirements
- Ensure all labware is properly seated before starting
- Do not open robot door during pipetting operations
- Keep DNA stocks and Gibson mastermix on ice until loading
- Verify thermocycler plate is PCR-compatible and properly seated
- Gibson mastermix is temperature-sensitive - minimize time at room temperature

---
