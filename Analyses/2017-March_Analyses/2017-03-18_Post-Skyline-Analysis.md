Working with data post Skyline:

1. I have a .csv file; opened with Excel and saved as .xlsx
2. Inserted new columns, sorted based on "file name" and annotated each row with: sample #, site, treatment, trial 
3. Sorted based on Protein Name; since I erroneously layered the full fasta proteome on top of the stress proteome, I have duplicate "reads" - I can tell because there are entries with only a shortened name (e.g. cds.comp142816_c0_seq5), which correspond to my stress-proteome. I thus delete these entries.  
4. Used Text-to-Columns and the =CONCATENATE functions in Excel to isolate the Protein name without isoform, e.g. from cds.comp142816_c0_seq5|m.43494 to cds.comp142816_c0
5. Removed all entries with "N/A" area and retention time
6. Created a Pivot table with edited Protein Name as rows, Treatment and Site as columns, and Max Area as the cell value
7. From pivot table, pulled average Max Area (not including WB), then calculated the Eelgrass / Bare ratio. So, the ratio is the Average Max Area in Eelgrass / Average Max Area in Bare. 
8. Created another pivot table as before, this time with Average Area 
9. In both, identified the proteins that have an Eelgrass/Bare ratio >1.5, so proteins that were expressed 50% more in Eelgrass. 
10. Attained the annotations, and edited the protein names to match the shortened names in my data (via Excel functions). 
11. Loaded annotation file into Galaxy online; 
12. Created tab-delimited .txt file with Protein ID's that correspond to overexpressed proteins in eelgrass. Uploaded this into Galaxy.
12. Used Galaxy's "Text Manipulation -> Join" to combine annotation data
13. Pulled Uniprot Accession Number from the list of eelgrass overly expressed proteins, as identifed from the Eelgrass/Bare ration from the mean Average Area from the 4 Puget Sound sites.
14. Pasted Uniprot Accession numbers into DAVID 
15. Selected 
16. Generated a Annotation Summary Result, using "All Species" in my Gene List Manager 
17. 

