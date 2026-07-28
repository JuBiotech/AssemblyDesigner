# Opentrons Flex Protocol: Combinatorial DNA Assembly - PCR

## Overview
This protocol automates combinatorial PCR reactions for DNA assembly. It processes a master worklist generated from the Assembly Designer package and performs multi-phase liquid handling followed by thermocycler-controlled PCR amplification.

## Protocol Phases

### Phase 1: Primer Mix Preparation
- Transfers forward primers (Primer1) from primer stocks to source plate
- Transfers reverse primers (Primer2) from primer stocks to source plate
- Mixes primer pairs in source plate wells

### Phase 2: PCR Reaction Setup
- Adds water to reaction plate in thermocycler
- Adds PCR mastermix to reaction wells
- Adds primer mixes from source plate
- Adds template DNA and mixes final reactions

### Phase 3: Thermocycler Program
- Initial denaturation at 98°C (2 minutes)
- PCR cycling with customizable parameters (default: 35 cycles)
- Final extension at 72°C (5 minutes)
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
1. **PCR Cycles** (default: 35)
   - Range: 1-40 cycles
   - Controls number of amplification cycles

2. **PCR Reaction Volume** (default: 50 µL)
   - Range: 10-100 µL
   - Total volume per PCR reaction

3. **Annealing Temperature** (default: 60°C)
   - Range: 50-72°C
   - Temperature for primer annealing step

4. **Elongation Time** (default: 240 seconds / 4 minutes)
   - Range: 10-420 seconds
   - Duration of extension step (adjust based on amplicon length)

## Pre-Protocol Checklist

### Before Starting
1. ✅ Prepare master worklist CSV from Assembly Designer
2. ✅ Thaw primer stocks, PCR mastermix, and template DNA on ice
3. ✅ Fill water reservoir (well A1) with nuclease-free water
4. ✅ Load primer stocks into plate in slot B2
5. ✅ Load PCR mastermix into plate in slot C2 (pre-loaded positions)
6. ✅ Load template DNA into plate in slot D2
7. ✅ Place empty PCR plate directly in thermocycler
8. ✅ Load **4 boxes** of 50 µL tips (slots B3, C3, C4, D4)
9. ✅ Load 1 box of 1000 µL tips (slot D3)
10. ✅ Verify thermocycler lid is open and unobstructed

### Plate Preparation Guidelines
- **Primer Stocks Plate (B2)**: Arrange forward and reverse primers according to worklist
- **Source Plate (C2)**: Pre-load PCR mastermix in designated wells (primer mixes will be created here)
- **Template DNA Plate (D2)**: Load template DNA according to worklist requirements

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
- Block pre-cools to 15°C during Phase 2 pipetting
- Reduces total protocol time by ~5-10 minutes

## Expected Runtime
- **Phase 1 (Primer Mix Preparation)**: ~5-10 minutes
- **Phase 2 (PCR Reaction Setup)**: ~15-25 minutes
- **Phase 3 (Thermocycler)**: ~2-3 hours (depends on cycle count and elongation time)
- **Total**: ~2.5-3.5 hours

## Troubleshooting

### Common Issues

**"Insufficient tips" error**
- Solution: Load additional 50 µL tip racks as indicated in error message
- Protocol calculates exact requirements before starting

**"CSV parsing error"**
- Solution: Verify CSV has all required columns and valid data
- Check for empty rows or missing values
- Ensure worklist names match expected format (e.g., `Worklist_Primer1`, `Worklist_Template`)

**Thermocycler lid won't close**
- Solution: Check for obstructions, verify plate is seated properly
- Protocol will pause and allow manual intervention

**Tip rack replacement fails**
- Solution: Protocol automatically pauses for manual replacement
- Remove depleted racks from B3/C3, place fresh racks, press Continue

### Protocol Pauses
The protocol may pause at:
1. **Thermocycler error**: Manual intervention required
2. **Tip replacement failure**: Manual tip rack replacement needed
3. **After thermocycler hold**: Retrieve completed reactions

## Support & Customization

### Modifying the Protocol
- **PCR cycles**: Adjust via runtime parameters before starting (1-40 cycles)
- **Reaction volume**: Adjust via runtime parameters (10-100 µL)
- **Annealing temperature**: Adjust via runtime parameters (50-72°C)
- **Elongation time**: Adjust via runtime parameters (10-420 seconds)
- **Worklist structure**: Must maintain required CSV column format
- **Labware types**: Can be modified in protocol code (lines 280-283)

## Safety Notes
⚠️ **Important**:
- Always verify deck layout matches protocol requirements
- Ensure all labware is properly seated before starting
- Do not open robot door during pipetting operations
- Keep primers, mastermix, and template DNA on ice until loading
- Verify thermocycler plate is PCR-compatible and properly seated
- PCR mastermix is temperature-sensitive - minimize time at room temperature
- Use nuclease-free water and maintain sterile technique

---
