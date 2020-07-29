# Teals_Strain-Level_Metagenome_Pipeline
<b>Purpose:</b>Using a mixture of my Universal Reference Database, single- and co-assembly, novel software, and dereplication, this pipeline gives you detailed community functions and strain(ish)-level MAGs.

The fundamental process is always the same for metagenome pipelines: read QC, assembly, read coverage, bin.<br> 
The devil is in the details...

## A New Perspetive On Metagenomes
Having already spent 3 years fixing metatranscriptomics by:
1. curating the taxonomy and phylogeny of all 1.8 million named species so they all conform with taxonomic standards, have no synonyms (eg Bacillus is a stick-insect and a bacterial genus), and divide without any flying buttresses (google it) into a tree of life with 8 levels (kingdom, phylum, class, order, family, genus, species, strain).
   - https://github.com/TealFurnholm/Universal-Taxonomy-Database
2. creating a universal reference database (URDB) containing all genes, gene functions, gene phylogeny/common ancestor, and contamination annotations (~15% of public genomes are significantly contaminated or mis-annotated). In the process creating an additional 110K new strain-level entries in the taxonomy database. (backstory: NCBI stupidly decided it was too hard to give each sequenced strain its own number, because taking the last number and adding +1 is such a challenge ::hard eyeroll::).
   - https://github.com/TealFurnholm/Meta-NGS_Reference_Database
3. having created this universal database, I know that >75% of all sequenced genes have orthologs w/>95% nucleotide identity. If you are going to do any meta-omics correctly, <b>YOU HAVE TO MULTIMAP YOUR READS</b>! We had to create new alignment software that can multi-map hundreds of millions of reads to 191 million distinct (547M total) genes without indexing (the reference database is 220GB = too big to all fit into RAM) and not have it take forever. This one I had to get help on.
   - https://github.com/kshedden/muscato.<p></p>

So I had a really different perspective shifting into metagenomes, having seen and fixed features of all the genes and taxonomy and functions of life (animals, plants, bacteria, viruses...).

<em><b>Aside: I was doing this work for my postdoc. My boss was forced to quit just as we were getting this all together and shortly afterwards I was diagnosed with advanced metastatic cancer - so you probably won't see this in publication. Just email me if you want a copy of the database/software.</b></em>

## Metagenomics is Borked.
After our lab dissolved, I started new contract work on human scalp samples, which had very reduced species/genus diversity but VERY high strain heterogenity. I'm not a glutton for punishment, I tried a variety of pipelines and software and at the end ... 6 bins. Yep. Billions of reads and only 6 bins. Having come from working on psoriasis, I KNEW the skin would have WAY more than 6 species. 
### 1. ALL assemblers are wrong.
As per usual you trim your reads, remove human contamination (I don't care if you are doing environmental samples, do this anyway), and remove low entropy reads that mess with assembly. 
- Gene orthology means conserved means assemblers are likely to split at the most important thing in the metagenome... genes.
- It is coding laziness or mental laziness (I'm a programmer I know) which is why they chose to split at branches. This is NOT how you handle conserved<->novel transitions. You have to trace the reads so you get the longest unambiguous stretch. 
- Since no software like this exists, the best you can do is combine your best assemblers (MetaSpades and MegaHit) and merge/overlap the output. You get different genes which give you very different bin outcomes BTW, so you REALLY should <ins>use multiple assemblers!</ins>

### 2. Binning loses you 90% of your data
Yep. Right out of the gate, their contig length cut-off loses you millions of genes. Remember, the average bacterial gene length is 800nt. Average, as in middle of the bell curve, and the contig length cut-offs are between 1000-2000 for the different binners. I've already harped on 
### 3. Bin quality checkers are wrong.
There is 
