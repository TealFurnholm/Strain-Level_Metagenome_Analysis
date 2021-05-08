### OVERVIEW
#### This pipeline is for paired short-read metagenome analysis.<br>
It uses a "universal database" and novel in-house scripts to create strain-level bins, as well as doing comprehensive linked community+functional primary and secondary analyses. This pipeline analyses each sample separately, and if 

### Requirements
 ** You should be doing this in a unix/linux shell environment, not Windows, unless you want to install a VM. We use a container called "comics" on our system, so you will need to change "comics" to whatever your specific system call is for the various software. We also have a 96 core system, of which I usually use 20 or 40 threads (excep the perl scripts only use 1) and a few hundred GB of RAM.
 * If you don't already have them on your system, install 
   - trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
   - bbtools: https://sourceforge.net/projects/bbmap/
   - diamond: https://github.com/bbuchfink/diamond
   - perl: https://www.perl.org/get.html
   - megahit: https://github.com/voutcn/megahit
   - spades: https://github.com/ablab/spades
   - quast: http://quast.sourceforge.net/quast
   - metabat: https://bitbucket.org/berkeleylab/metabat/src/master/
 * You'll also need to grab all the scripts and template files: https://github.com/TealFurnholm/Teals_Strain-Level_Metagenome_Pipeline/archive/master.zip
 * Get the human genome (even if you have environmental metagenome!!!)
>     wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.p13.genome.fa.gz

### samp_list.txt
The pipeline is set up to run many samples listed in "samp_list.txt" through each step of the pipeline. 
The "samp_list.txt" is just the unique sample prefixes before the fwd/rev: **something_123_fwd.fastq should be something_123**
#### Example of a samp_list.txt
control_1<br>
control_2<br>
control_3<br>
test_1<br>
test_2<br>
test_3<br>

## Begin: https://github.com/TealFurnholm/Teals_Strain-Level_Metagenome_Pipeline/wiki
