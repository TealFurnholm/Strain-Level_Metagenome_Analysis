# What you NEED to know about microbiomics! (seriously, read this)
Between my years working on host/human and environmental metatranscriptomics, metabolomics, and metagenomics, I have been thoroughly disappointed by existing software and pipelines for microbiome analysis. In metagenomics, focus is on bins rather than all contig, BUT, >80% of contigs/data are discarded when binning (contigs too short because of poor assembly practices/software). Due to orthology, both 30-70% of reads will multimap (bad for read depth calculation) and contigs will multi-bin (bins less complete), but pipeline authors and programmers have failed to take this into account. *Collectively this results in ridiculously few decent quality bins.* Most are way more contaminated/less complete than they think because they use "single-copy" genes for QC. <br> <br>
Also, in my experience I can tell you that the gene-level (ie. strains and functional pathways) is where it matters. Think about the strain heterogeneity between samples, the gene orthology & HGT, strain specific microbe-microbe/host interactions, microbial metabolic feeding and quorum regulation. You may have two people with the same disease, where the microbes are interfering with two different steps in the same pathway, or worse, different pathways (eg. disrupt anti-inflammatory regulator vs increase pro-inflammatory lipid). Right now scientists average samples & compare conditions in whatever R package. But if only 10% strains/50% genes are conserved between samples, and if you compare gene-wise (which you have to because of read normalization) or species-wise (also because of normalization), you are throwing away a huge chunk of the Big Picture/Answer. I don't have the final answer, which probably involves constructing full metabolic and immunological pathways with each player within each sample, and only then seeing if the samples are comparable by conserved mechanism or species. Still, I can get you quite a bit farther with bins, functions and communities ... <br>
 
#### Unfortunately, I also have metastatic cancer and probably won't get to publish these things (long story, bosses quitting, contract work, bad doctors getting me killed .. )
SO. If you want to do the microbiomics field any good: 
- 5% contamination is NOT OK (1 in 20 genes is contamination? really? No, bad!) aim for < 0.1% at the bare minimum
- don't trust any software blindly, NO ONE checks the software/code logic before publication, even if the programmer wanted (and we DO!). Try every alternative software at least once, and play with the settings - they will give you substantially (thousands of genes) better or worse results.
- take the time to name and describe your public sequence deposits with the consideration for all of science, now and future (at least make sure no other species has that name!)
- culture, co-culture, and look under the microscope = DONT NEGLECT THE WET LAB. In silico is not a silver bullet. 
  - Check out the Cytation 5 machine by BioTek, use non-cytotoxic live fluorescent stains to actually SEE how microbes are interacting with host/other microbes!

## Anyhoo, I've developed databases, software, and pipelines to fix the many issues and gaps in the analysis of microbiomes I've found so far, including:
1. a Taxonomy-Phylogeny Database for all 2.2 million species/strains of bacteria, viruses, plasmids, eukaryotes, archaea...etc, which is curated to be consistent, logical, and produces perfect untangled phylogenetic trees
2. my own Universal (as in every sequenced species of anything) Gene Database with 
  - all known functions for each of the 529 million genes (191 distinct sequences - there's lots of orthology)
  - the phylogenetic lineage of the last common ancestor and each individual species
  - after aligning 100nt pieces of every gene to every other gene and comparing their %GC and phylogeny, identifying any genome contaminants/HGT
3. various in house software for removing low-entropy reads, correctly normalizing reads for creating quantitative-comparative community phylogenies and doing functional analyses, and (probably what you'll care most about) creating and QC-ing metagenome bins.

### Outputs: 
1. detailed community genes, functions, and phylogeny
2. differential analysis and visualization
3. many more complete & strain(ish)-level MAGs

# So check out the wiki to get started! https://github.com/TealFurnholm/Teals_Strain-Level_Metagenome_Pipeline/wiki
