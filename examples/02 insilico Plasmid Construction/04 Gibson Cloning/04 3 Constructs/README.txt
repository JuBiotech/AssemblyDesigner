This example builds two 3-TU constructs (C01, C02) that share an identical
TU1 and TU2, and only differ in TU3. It exists to demonstrate (and prevent a
regression of) two fixed bugs in the 3G/Gibson pipeline:

1. Building the *same* TU for two different constructs used to crash on
   Windows with "PermissionError: [WinError 32]" because the code tried to
   re-simulate and overwrite the same DNAcauldron report ZIP a second time.
   Fixed in assembly_designer/plasmidio/assembly.py: identical designs are
   now simulated once and reused.
2. AssemblyHistory.plot() always drew the backbone at the median row (which
   lands on top of the middle TU when the TU count is odd) and used a fixed
   figure size, so labels started overlapping from 3 TUs onward. Fixed in
   assembly_designer/plasmidio/history.py: the backbone is now always placed
   below every other row, and the figure auto-sizes with the number of rows.

Step 1:
- Add all relevant plasmid maps (Genbank or Snapgene) to the corresponding
  folders 01-06 (already populated here, reused from the other 3G/Gibson
  examples in this repo).

Step 2:
- Run "GGA_Parts_Parser.ipynb" if you change the parts in folders 01-06; it
  refreshes the "Info" sheet in the xlsx so you don't have to type part names.

Step 3:
- Plan each assembly in "gibson_designs_template.xlsx" (Constructs/TUs/PCR/
  Assembly sheets). Note that C01 and C02's TU1/TU2 rows are identical on
  purpose.

Step 4:
- Run "3G Assembly.ipynb" to assemble the two constructs and render the
  history plots into reports/History/.
