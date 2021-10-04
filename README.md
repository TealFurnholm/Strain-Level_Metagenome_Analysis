### OVERVIEW
* The tail-end of the pipeline - binning and secondary analysis - are being worked on as of 5/17/21. This message will disappear once it is complete.
#### This pipeline is for paired short-read metagenome analysis.<br>
It uses a "universal database" and novel in-house scripts to create strain-level bins, as well as doing comprehensive linked community+functional primary and secondary gene analyses. This pipeline analyses each sample separately, but there is a co-assembly option, and the binning process can be run with multiple samples.

### Requirements
 * You should be doing this in a unix/linux shell environment
 * I usually use a max of 40 threads (excep the perl scripts only use 1)
 * You'll want/need a few hundred GB of RAM
 * If you don't already have them on your system, install perl and make sure it can run from any folder: https://www.perl.org/get.html

### Set up your samp_list.txt
The pipeline is set up to run many samples listed in "samp_list.txt" through each step of the pipeline. 
The "samp_list.txt" is just the unique sample prefixes before the fwd/rev: **something_123_fwd.fastq** should be listed as **something_123**
#### Example of a samp_list.txt
control_1<br>
control_2<br>
control_3<br>
test_1<br>
test_2<br>
test_3<br>

### Start [Here](https://github.com/TealFurnholm/Teals_Strain-Level_Metagenome_Pipeline/wiki)
