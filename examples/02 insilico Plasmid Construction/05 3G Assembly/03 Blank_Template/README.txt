Step 1:
- Add all relevant plasmid maps (Genbank or Snapgene) to the corresponding folders 01-06

Step 2:
- Run the "GGA_Parts_Parser.jpynb"
- It reads out the names of the added parts and prints them into the corresponding columns in the "Info" sheet to avoid typos when planning assemblies
- The "Info" sheet also contains information about the 3G principle and which primers are currently used for assemblies with 2-5 TUs

Step 3:
- Plan each assembly in the "gibson_designs_template.xlsx" file

Step 4:
- Run the "3G Assembly.jpynb" notebook to assemble the final 3G plasmids.
- Check that your plasmids were correctly assembled.
