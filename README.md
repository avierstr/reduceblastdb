## Reduceblastdb: a tool for reducing Blast databases to speed up local BLAST searches

reduceblastdb, a tool for reducing the BLAST databases based on taxonomy, sequence length and optional removing highly similar sequences from each species.

**As an initial indication of sizes:**

 - nt database: 445 GB 
 - Bacteria: 55 GB 
 - Bacteria max 2000 bp: 3.6 GB

(Only works on Linux and Mac because it uses multiprocessing. Mac with M1 processor can cause problems)

**Requirements:**
-   Python 3.8
-   edlib: Lightweight, super fast C/C++ library for sequence alignment using edit (Levenshtein) distance ([https://pypi.org/project/edlib/#description](https://pypi.org/project/edlib/#description))  
    (`python3 -m pip install edlib`)
-   biopyton (`sudo apt-get install python3-biopython`)
-   Blast executables (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (I am noticing differences in number of sequences that are processed depending on the version the Blast executables (2.12 and 2.15), so best to install the latest version!)
  
 **Options:**

There are 2 parts with options: **blastdb**  and **taxonomy**. **Blastdb** is to select parts of the blast databases to create a smaller local database. **Taxonomy** is to get information about the NCBI classification.

**Blastdb options:**

`-db, --database`: Which database ? nt, nr, LSU\_eukaryote\_rRNA, ...

`-o, --outfile`: Name of the outputfile saved in the outputfolder "reduced"

`-of, --outfolder`: Name of the output folder to save the databases. Custom created databases will be saved in the subfolder "reduced". If none is given, current working directory is used.

`-dl, --download`: Download the needed database from NCBI

`-min, --minlength`: Minimum sequence length

`-max, --maxlength`: Maximum sequence length

`-s, --select`: Select TaxIDs to make a sub-database, use "1" for all species use one or more TaxIDs or scientific names for one or more groups

`-u, --unknown`: Exclude non identified entries (unidentified, unknown,...)

`-g, --gene`: Select genes to make a sub-database, use one or more gene names that are used in the "sequence title" in GenBank

`-a, --append`: Append sequences to a previous created fasta file. (it adds the sequences to the "selected.fasta" file)

`-r, --reduce`: Reduce the number of sequences by removing highly similar sequences from each species

`-np, --nproc`: Number of processors to use in the "reduce" option. Default= 4

`-mbd, --makeblastdb`: Create the new Blast database from the selected or reduced file. It automatically selects the reduced file if present, if not it uses the selected file.


**Taxonomy options:**

the following commands have to be preceded by the word “tax”

`-of, --outfolder`: Name of the output folder to save the taxonmomy files. If none is given, current working directory is used.

`-tl, --taxidlist`: Create a TaxID list from NCBI with name, taxid and rank in alfabatic order (fast). List is in csv format. By default, everything below genus level is excluded to limit the file size. You can enter a taxonomic level to change the default or use "all" to use all levels. \[kingdom, phylum, class, order, family, genus, species, all\] You can limit the list by giving a TaxID or scientific name

`-tt, --taxtree`: Create a "Taxonomic tree" from NCBI with rank, name and taxid in alfabatic order per taxonomic level (very slow (12 hours for the complete database for all levels)). Information is saved in csv file. This is not a fylogenetic tree file. By default, everything below genus level is excluded to limit the file size. You can enter a taxonomic level to change the default or use "all" to use all levels. \[kingdom, phylum, class, order, family, genus, species, all\] You can limit the tree by giving a TaxID or scientific name

`-l, --lineage`: Create the full taxonomic lineage from a TaxID or scientific name

**How it works:**

**Blastdb:**

You need a blast database to start with. Downloading this takes some time and needs quite some harddisk space. NCBI has a lot of databases available ([https://ftp.ncbi.nlm.nih.gov/blast/db](https://ftp.ncbi.nlm.nih.gov/blast/db)). One of the biggest is the nucleotide (nt) that takes 450 GB of diskspace. (Taking a quick look at a available NCBI subdatabase like “ITS\_eukaryote\_sequences”, it seems that it does not contain all ITS sequences from eukaryotes, so I would suggest to start from the complete databases like nt or nr.)

Suppose you want to **download the** nucleotide (nt) **database**:

    python3 reduceblastdb.py blastdb -db nt -of nt_database -dl

This is interruption proof.  You can restart an interrupted download, it will continue where it was interrupted.  Once the database is downloaded, you can make a lot of **selections for sub-databases**. It needs the database and a name for outputfile:

Select all bacterial sequences:

    python3 reduceblastdb.py blastdb -db nt -of nt_database -o bacteria -s bacteria

Select all bacterial sequences with a maximum length of 2000 bp:

    python3 reduceblastdb.py blastdb -db nt -of nt_database -o bacteria -s bacteria -max 2000

Select all bacterial and fungal sequences with a maximum length of 2000 bp:

    python3 reduceblastdb.py blastdb -db nt -of nt\_database -o bacteria\_fungi -s bacteria fungi -max 2000

Select all bacterial 16S sequences:

    python3 reduceblastdb.py blastdb -db nt -of nt\_database -o bacteria\_16S -s bacteria -g 16S

Keep in mind that the “gene selection” is based on the information in the “sequence title” in Genbank. Genbank is not consistent in its descriptions in the title. '16S', '16S rRNA', '16S ribosomal RNA', 'small subunit ribosomal', '16S rDNA' are used gene descriptions in the titles for the bacterial 16S gene. Reduceblastdb.py has already a few genes covered in the script if you use one of the descriptions it will search for all variations. If you type 16S, it will search for the variations '16S', '16S rRNA', '16S ribosomal RNA', 'small subunit ribosomal', '16S rDNA' .

If they are not covered in the list below, you can enter all the variations for genes after the “-g” option.

['COI', 'CO1', 'COX1', 'COXI', 'cytochrome oxidase I', 'cytochrome oxidase 1', 'mitochondrion complete genome']

['COII', 'CO2', 'COX2', 'COXII', 'cytochrome oxidase 2', 'cytochrome oxidase II', 'mitochondrion complete genome']
['COIII', 'CO3', 'COX3', 'COXIII','cytochrome oxidase III', 'cytochrome oxidase 3', 'mitochondrion complete genome']

['12S', 'small subunit ribosomal']

['16S', '16S rRNA', '16S ribosomal RNA', 'small subunit ribosomal', '16S rDNA']

['18S', '18S rRNA', '18S ribosomal', 'small subunit ribosomal', 'SSU rRNA']

['23S', 'large subunit ribosomal']

['26S', 'large subunit ribosomal']

['28S', 'LSU', 'large subunit ribosomal', 'large subunit rRNA','large subunit of rRNA']

['ITS', 'internal transcribed spacer']

['NADH1', 'nad1', 'NADH dehydrogenase subunit 1', 'nd1', 'NADH dehydrogenase subunit I', 'NDI', 'nicotinamide adenine dinucleotide dehydrogenase subunit 1']

['NADH2','nad2', 'NADH dehydrogenase subunit 2', 'ND2', 'NADH-ubiquinone oxidoreductase chain 2']

['psbA', 'photosystem II']

['psaA', 'photosystem I']

['rbcL', 'ribulose 1,5 bisphosphate carboxylase/oxygenase large']

['tufA', 'elongation factor Tu']

['rpb2', 'RNA polymerase II', 'RNA polymerase second', 'rpbII']

['ef1', 'TEF', 'elongation factor 1', 'elongation factor alpha', 'EF-1']


**Unknown sequences** can be excluded. This is open for discussion. For now it will exclude all entries where the title begins with \['unidentified', 'unknown', 'unspecified', 'untyped', 'ungrouped', 'undetermined', 'undescribed', 'uncultured', 'uncultivated', 'unclassified'\].

    python3 reduceblastdb.py blastdb -db nt -of nt_database -o bacteria -s bacteria -u
    
But there are a lot of other entries that are also not informative for some researchers and are desired for others.  The following are NOT excluded when the “-u” option is used. \[soil, sludge, sediment, salal, rumen, root, rhizosphere, rape rhizosphere, rainbow trout, psychrophilic, prokaryote enrichment, primary endosymbiont, nitrogen fixing, mycorrhizal, mucus, mixed, methylotrophic, methanotrophic, methanotroph, methanogenic, metal-contaminated, mesophilic, mercury-resistant, marine, maize, low G+C, leaf litter, iron-reducing, intestinal, ice core, humic, halophilic, haloarchaeon, haloalkaliphilic, groundwater, grassland, glacial, gamma proteobacterium, fungal, fuel, freshwater, fossil, forest, foliar, filamentous, fern, extreme, eukaryote, eubacterial, ericoid, epsilon, epacrid, environmental, enrichment, endosymbiont, endophytic, endocytic, ectomycorrhizal, earthworm, drinking, denitrifying, delta, cyanobacterium enrichment, crenarchaeote enrichment, candidate division, beta proteobacterium, barley, bacterium enrichment, archaeon enrichment, anammox bacterium enrichment, alpha proteobacterium enrichment, activated sludge, actinobacterium enrichment)\]

There is a **append** “-a” option available. Suppose you want a combined database with bacterial 16S sequences and fungal ITS sequences. If you use the “-g 16S” option for bacteria and fungi, it will search for “small subunit ribosomal” witch is also used for the 18S gene in fungi which you don’t want. Then it is better to perform 2 separate selections with the append option.

    python3 reduceblastdb.py blastdb -db nt -of nt\_database -o bacteria\_fungi -s bacteria -g 16S

    python3 reduceblastdb.py blastdb -db nt -of nt\_database -o bacteria\_fungi -s fungi -g ITS -a

The **reduce** option “-r” will read the selected sequences from above (saved in xxx\_selected.fasta) and create a lot of xxx.spec files in 100 tmp folders. Each file contains the sequences of one species. Next, the program will read each file and compare the sequences in the file with each other. If some sequences are >= 97% similar, it will keep only the longest sequence. That way there will be an extra reduction in number of sequences. The result is saved in xxx\_reduced.fasta. Depending on the amount of selected sequences, this can take a long time (for example more than 2 weeks on 8 cores for all bacterial sequences). It is “interruption proof”: you can interrupt the reduction process and restart it later, it will continue where it was interupted. The script is saving 2 files (xxx\_removeset.pick and xxx\_uniqueset_pick) which contain the accessions that can be removed and which are unique. If you download a newer database later, the script will know which it can remove immediately when you perform a new reduction on the new dataset. This will save a lot of time.

    python3 reduceblastdb.py blastdb -db nt -of nt\_database -o bacteria\_fungi -r -np 8

The last step is to **create a new Blast database** from the selected or reduced dataset. This is done with the “-mbd” option. The script will look for a xxx\_reduced.fasta file and create a database. If it does not find the xxx\_reduced.fasta, it will use the xxx_selected.fasta file.

    python3 reduceblastdb.py blastdb -db nt -of nt\_database -o bacteria\_fungi -mbd
    
All the above commands can be performed in one go:

    python3 reduceblastdb.py blastdb -db nt -dl -o bacteria -s bacteria -max 2000 -r -np 8 -mbd

**Taxonomy:**

Taxonomy contains some extra options which are not directly connected to the reduce database functions. Genbank uses a unique Taxonomic Identifier (TaxID) for each entry (for example ‘species’) or node (for example ‘family’) ([https://www.ncbi.nlm.nih.gov/books/NBK53758](https://www.ncbi.nlm.nih.gov/books/NBK53758)).

Suppose you want to know the **lineage** from the root to a certain species or taxonomic rank:

    python3 reduceblastdb.py tax -l cordulegaster

    Rank Scientific name TaxID
    
    No Rank cellular organisms 131567
    Superkingdom Eukaryota 2759
    Clade Opisthokonta 33154
    Kingdom Metazoa 33208
    Clade Eumetazoa 6072
    Clade Bilateria 33213
    Clade Protostomia 33317
    Clade Ecdysozoa 1206794
    Clade Panarthropoda 88770
    Phylum Arthropoda 6656
    Clade Mandibulata 197563
    Clade Pancrustacea 197562
    Subphylum Hexapoda 6960
    Class Insecta 50557
    Clade Dicondylia 85512
    Subclass Pterygota 7496
    Infraclass Palaeoptera 33339
    Order Odonata 6961
    Suborder Epiprocta 2510002
    Infraorder Anisoptera 6962
    Superfamily Cavilabiata 70899
    Family Cordulegastridae 107811
    Genus Cordulegaster 126172

If you want to produce an **alphabetic list** with everything from the order Odonata:

    python3 reduceblastdb.py tax -tl odonata

Partial screenshot:

![treelist](https://github.com/avierstr/reduceblastdb/blob/main/reduceblastdb_2.jpg)

If you want the same with a more **tree-like** look for the order Odonata:

    python3 reduceblastdb.py tax -tt odonata

Partial screenshot:

![tree](https://github.com/avierstr/reduceblastdb/blob/main/reduceblastdb_1.jpg)

**Release notes:**

 2024/07/08:
 - initial release.

