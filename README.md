# DB-Manager

This python script allow you to download a set of databases for proteomis automatically and generate the corresponding decoy databases in a corresponding structure: 
  - database (1)
    - 05-06-2015 (a version of the database depending of the date that was run the script)
    - latest     (a symbolink to the latest version of the database)
  - database (2)
    - 05-06-2015 (a version of the database depending of the date that was run the script)
    - latest     (a symbolink to the latest version of the database)

If the user provides the decoy option the (--decoy or -d) it will generate the decoy database: 
    - example of the command DecoyDatabase -in %s -out %s -reverse
    ** Note that the input and output file needs an structure like %s 
     
# Version and dependencies 
  - The current version is 0.0.1
  - The script works perfectly well with python 2.7 or higher and the native librarie.
  - The system should provide gzip command line option

# Databases

The program is now focus in human databases, it now supports:
  - [ ] Complete Human Uniprot Proteome: UniProt provides proteome sets of proteins thought to be expressed by organisms whose genomes have been completely sequenced.
  - [ ] Human SwissProt: The curated version of Human database in SwissProt
  - [ ] crap: The GPMDB contaminant database
  - [ ] max-crap: The MaxQuant contaminant database 
  - [ ] Human ENSEMBL: The ensembl human peptide database 
  - [ ] Human RefSeq: The Human Reference Sequence database from NCBI 

    
