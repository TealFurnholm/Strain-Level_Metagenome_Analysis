use warnings;
use Statistics::Basic qw(:all nofill);
$time = localtime;
$time = uc($time);
$time =~ /^[A-Z]+\s+([A-Z]+)\s+\S+\s+\S+\s+(\d\d\d\d)/;
$month = $1; $year = $2;
$version=$1."_".$2;

$argv = join(" ", @ARGV);
if($argv !~ /\w/){
        print "You must specify the following: \n";
        print "\t-d /path/to/reference_databases/files/ \n";
        print "\t-s the_sample_prefix\n";
        print "\t-i /path/to/your/alignment/m8/files/ (if current directory, you can ignore this flag)\n";
        print "\tEx. for Samp_1234_PROTS.m8 in this directory -d ~/Ref_DBs/ -i ../ALIGNMENTS/ -s Samp_1234 \n";
        print "\tThe script uses the taxonomy, function names, and UniRef100 databases. See my GitHub for details\n";
        print "You can optionally specify the:\n";
        print "\t-c \t Minimum \% query coverage (default -c=70)\n";
        print "\t-m \t Minimum \% identity alignment (default -m=40)\n";
        print "\t-x \t Minimum number of genes to keep cellular organisms (default -x=3). Note: does not include viruses or plasmids\n";
        print "\t-k \t Top score margin, between 0 and 1. \n";
        print "\t\t Scores are \%identity x \%coverage (max possible score is 100x100=10000).\n";
        print "\t\t The -k margin keeps hits where: hit score >= top_score x [margin].\n";
        print "\t\t ex: 10000 x 0.75 = keep all hits with score >= 7500\n";
        print "\t\t Top score hits only = 1. Ignore hit scores = 0. (default = 0.75)\n\n";
        print "These options allow that the reference database may not have your sample's strains, just various relatives with mixed gene orthology\n";
        print "Lower stringency helps stop false negatives. Higher helps increase speed and reduce noise/false positives.\n";
        print "Many Excess/spurious hits and species are screened out later in the script.\n\n\n";
        die;
}


if($argv =~ /\-d\s+(\S+)/){
        $dir=$1;
        if($dir !~ /\/$/){$dir.='/';}
        $dir='/geomicro/data22/teals_pipeline/BOSS/'; ###REMOVE AFTER DEBUG
        opendir(DIR, $dir) or die "Could not open $dir\n";
        @FILES = grep(/TAXONOMY\_DB.*\.txt/i, readdir DIR);
        $intax=$dir.$FILES[0]; @FILES = ();

        opendir(DIR, $dir) or die "Could not open $dir\n";
        @FILES = grep(/UNIREF100_INFO.*\.txt/i, readdir DIR);
        $ininfo=$dir.$FILES[0]; @FILES = ();

        opendir(DIR, $dir) or die "Could not open $dir\n";
        @FILES = grep(/Function_Names.*\.txt/i, readdir DIR);
        $infn=$dir.$FILES[0]; @FILES = ();

        if(!-s $intax){ print "Cant find $intax in $dir. Please make sure $dir contains taxonomy database: https://github.com/TealFurnholm/Universal-Taxonomy-Database \n"; die;}
        if(!-s $ininfo){print "Cant find $ininfo in $dir. Please make sure $dir contains uniref100 info: https://github.com/TealFurnholm/Universal_Microbiomics_Alignment_Database \n"; die;}
        if(!-s $infn){  print "Cant find $infn in $dir. Please make sure $dir contains function names file: https://github.com/TealFurnholm/Universal_Biological_Compounds_Database/wiki \n"; die;}
}


if($argv =~ /\-s\s+(\S+)/){$samp=$1;}
else{
        print "A sample prefix [sample]_READSvsGENES.rpkm, [sample]_GENES.fna, [sample]_GENES.m8 ...etc is required\n";
        print "Specify sample ex: -s Toxin_Day2\n"; die;
}
if($argv =~ /\-i\s+(\S+)/){$sdir=$1;} else{$sdir='./';}
$samp=$sdir.$samp;
$ingrc  = $samp."_READSvsGENES.rpkm";           if(!-s $ingrc){ print "missing file $ingrc.\n"; die; }
$inccov = $samp."_READSvsCONTIGS.rpkm";         if(!-s $inccov){ print "missing file $inccov.\n"; die; }
#CHECK FOR GENES OR PROTEIN ALIGNMENTS
   if(-s $samp."_PROTS.m8"){$ingvu  = $samp."_PROTS.m8";}
elsif(-s $samp."_GENES.m8"){$ingvu  = $samp."_GENES.m8";}
else{ print "missing samp_PROTS.m8 or samp_GENES.m8.\n"; die;}
   if(-s $samp."_PROTEINS.faa"){$ingen  = $samp."_PROTEINS.faa";}
elsif(-s $samp."_GENES.fna"){   $ingen  = $samp."_GENES.fna";}
else{ print "missing samp_PROTEINS.faa or samp_GENES.fna \n"; die;}
   if(-s $samp."_MCDD.fa"){ $incon = $samp."_MCDD.fa";}
else{ print "missing contig fasta file samp_MCDD.fa \n"; die;}
if($argv =~/\-k\s+([\d\.]+)/){$top=$1;}else{$top = 0.75;}
if($argv =~/\-x\s+(\d+)/){$mingen=$1;} else{$mingen = 3;}
if($argv =~/\-m\s+(\d+)/){$minid=$1;}  else{$minid = 40;}
if($argv =~/\-c\s+(\d+)/){$mincov=$1}  else{$mincov = 70;}


#SET OUTPUT FILES
$samp =~ s/.*[\/\\]//;
$outgene  = $samp."_CONTIG_GENE_INFO_".$version."Q.txt";
$outfuncs = $samp."_COMMUNITY_FUNCTIONS_".$version."Q.txt";
$outcyto  = $samp."_COMMUNITY_PHYLOGENY_".$version."Q.cyto";
$outtrans = $samp."_TRANSPORTER_COMP_".$version."Q.txt";
$outmrx   = $samp."_INTERACTION_MATRIX_".$version."Q.txt";
open(DEBUG, ">", $samp."_CONAN_".$version."Q.debug")||die;

print "sample $samp taxondb $intax funcnm $infn dbinfo $ininfo samppfx $samp topscore $top mingen $mingen minid $minid mincov $mincov\n";



######## INPUT DATA ###########
################################
my ($congens, $genpos, $genprt, $genpgc, $genelength) = GetConGenPos();
        %CON_GENES=%$congens;           $kc1 = keys %CON_GENES;         if($kc1 < 1){print "no $kc1\n"; die;} #{$contig}{$gene}=1;
        %GENE_POS= %$genpos;            $kc2 = keys %GENE_POS;          if($kc2 < 1){print "no $kc2\n"; die;} #{$gene}=$start;
        %GENE_PART=%$genprt;            $kc3 = keys %GENE_PART;         if($kc3 < 1){print "no $kc3\n"; die;} #{$gene}=$part;
        %GENE_PGC= %$genpgc;            $kc4 = keys %GENE_PGC;          if($kc4 < 1){print "no $kc4\n"; die;} #{$gene}=$prgc;
        %GENE_LEN= %$genelength;        $kc5 = keys %GENE_LEN;          if($kc5 < 1){print "no $kc5\n"; die;} #{$gene}=$len;


my ($generc, $generpkm, $totreads) = GetGeneReadCounts();
        %GENE_RD=%$generc;              $kc6 = keys %GENE_RD;           if($kc6 < 1){print "no $kc6\n"; die;} #{$gene}= depth
        %GENE_KM=%$generpkm;            $kc7 = keys %GENE_KM;           if($kc7 < 1){print "no $kc7\n"; die;} #{$gene}= rpkm

my ($cpgc, $conlen) = GetContigPGC();
        %CON_PGC=%$cpgc;                $kc8 = keys %CON_PGC;           if($kc8 < 1){print "no $kc8\n"; die;} #{$contig} = %GC
        %CON_LEN=%$conlen;              $kc9 = keys %CON_LEN;           if($kc9 < 1){print "no $kc9\n"; die;} #{$contig} = $len

my ($ccov, $crpkm) = GetContigCov();
        %CON_RD=%$ccov;                 $kc10 = keys %CON_RD;           if($kc10 < 1){print "no $kc10\n"; die;} #{$contig} = depth
        %CON_KM=%$crpkm;                $kc11 = keys %CON_KM;           if($kc11 < 1){print "no $kc11\n"; die;} #{$contig} = rpkm

my $funcnm = GetFuncNames();
        %FUNC_NM=%$funcnm;              $kc12 = keys %FUNC_NM;          if($kc12 < 1){print "no $kc12\n"; die;} #{$id} = $name

my ($phytid, $lintids) = GetTaxonInfo();
        %PHY=%$phytid;                  $kc13 = keys %PHY;              if($kc13 < 1){print "no $kc13\n"; die;} #{$tid} = $lin
        %LIN_TID=%$lintids;             $kc14 = keys %LIN_TID;          if($kc14 < 1){print "no $kc14\n"; die;} #{$lev} = $tid

my ($gentophit, $ghits, $gmxs) = InDiamond();
        %GENE_HITS=%$gentophit;         $kc15 = keys %GENE_HITS;        if($kc15 < 1){print "no $kc15\n"; die;} #{$gene}{$hit} = $sco
        %HIT_GENES=%$ghits;             $kc16 = keys %HIT_GENES;        if($kc16 < 1){print "no $kc16\n"; die;} #{$hit}{$gene} = $sco
        %GENE_SCO=%$gmxs;               $kc17 = keys %GENE_SCO;         if($kc17 < 1){print "no $kc17\n"; die;} #{$gene} = $sco
PrintKC();

my ($lingen, $genlin, $genfun, $gennm, $tlgen, $htlns) = GetRefGeneInfo();
        %LIN_GENES=%$lingen;            $kc18 = keys %LIN_GENES;        if($kc18 < 1){print "no $kc18\n"; die;} #{$lin}{$gene}=$sco
        %GENE_LINS=%$genlin;            $kc19 = keys %GENE_LINS;        if($kc19 < 1){print "no $kc19\n"; die;} #{$gene}{$lin}=$sco
        %GENE_FUNCS=%$genfun;           $kc20 = keys %GENE_FUNCS;       if($kc20 < 1){print "no $kc20\n"; die;} #{$gene}{$i}{$id}++; i=7..24
        %GENE_NAME=%$gennm;             $kc21 = keys %GENE_NAME;        if($kc21 < 1){print "no $kc21\n"; die;} #{$gene}=$name
        %TOT_LIN_GENES=%$tlgen;         $kc22 = keys %TOT_LIN_GENES;    if($kc22 < 1){print "no $kc22\n"; die;} #{$lin}++
        %HIT_LINS=%$htlns;              $kc23 = keys %HIT_LINS;         if($kc23 < 1){print "no $kc23\n"; die;} #{$hit}{$lin}

my $linscos = CleanLineages();
        %LIN_SCO=%$linscos;             $kc24 = keys %LIN_SCO;          if($kc24 < 1){print "no $kc23\n"; die;}

my ($cnlca, $gnlca, $cnlin) = GetContigLCA();
        %CON_LCA=%$cnlca;               $kc25 = keys %CON_LCA;          if($kc25 < 1){print "no $kc24\n"; die;} #{$contig}=$conlca
        %GEN_LCA=%$gnlca;               $kc26 = keys %GEN_LCA;          if($kc26 < 1){print "no $kc25\n"; die;} #{$gene}  =$genelca
        %CON_LINS=%$cnlin;              $kc27 = keys %CON_LINS;         if($kc27 < 1){print "no $kc26\n"; die;} #{$contig}{$lin}=1

PrintKC();



##############################################
############ OUTPUT CONTIG INFO ##############
##############################################

print "output gene-contig info\n";
open(OUTGENE,">",$outgene)||die;
print OUTGENE "gene\tcontig\tname\tpartial\tsco\t";
print OUTGENE "start\tgenlen\tconlen\tcongencnt\t";
print OUTGENE "gene_rpkm\tcontig_rpkm\t";
print OUTGENE "gene_depth\tcontig_depth\t";
print OUTGENE "gene_pgc\tcontig_pgc\t";
print OUTGENE "gen_str\tgentids\tgene_lca\t";
print OUTGENE "con_str\tcontids\tcon_lca\t";
print OUTGENE "SIGALPEP\tTMS\tDNA\tMETAL\tTCDB\tLOCATION\tCOG\tPFAM\tTIGR\tGO\tIPR\tEC\tKEGG\tRHEA\tBIOCYC\tREACTANTS\tPRODUCTS\tTRANS_CPD\n";

foreach my $contig (sort(keys %CON_PGC)){

        $conlca = $CON_LCA{$contig};
        $constr  = keys %{$CON_LINS{$contig}};
        %CTIDS=();
        foreach my $lin (keys %{$CON_LINS{$contig}}){ $CTIDS{$LIN_TID{$lin}}=1;}
        @TLIST=();
        foreach my $tid (sort(keys %CTIDS)){ push(@TLIST,$tid); }
        $contids = join(";",@TLIST);
        if($contids !~ /\d/){$contids=0;}
        $congcnt = keys %{$CON_GENES{$contig}};
        $conlen = $CON_LEN{$contig};
        $conpgc = $CON_PGC{$contig};
        $concov = $CON_RD{$contig};
        $conrpkm = $CON_KM{$contig};

        if(!exists($CON_GENES{$contig})){
                @FUNC_LIST=(); for my $i (0..17){ push(@FUNC_LIST, "0"); }
                $funcs=join("\t",@FUNC_LIST);
                print OUTGENE "0\t$contig\tNO_GENES\t0\t0\t";
                print OUTGENE "0\t0\t$conlen\t0\t";
                print OUTGENE "0\t$conrpkm\t";
                print OUTGENE "0\t$concov\t";
                print OUTGENE "0\t$conpgc\t";
                print OUTGENE "0\t0\t0\t";
                print OUTGENE "0\t0\t$conlca\t";
                print OUTGENE "$funcs\n";

                #GET CYTO STUFF
                $conlca=~/^(.)/;
                $king = $1;
                @TAXON=split(";",$conlca);
                for($i=0; $i<=$#TAXON; $i++){
                        $type = $king."_".$i;
                        if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
                        else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
                        $CYTO_CCON{$phyla}++;
                        $CYTO_CGEN{$phyla}+=0;
                        $CYTO_CLEN{$phyla}+=$conlen;
                        $CYTO_RPKM{$phyla}+=$conrpkm;
                        $CYTO_CPGC{$phyla}+=$conpgc;
                        $CYTO_CCOV{$phyla}+=$concov;
                }
                delete($CON_PGC{$contig});
                delete($CON_LCA{$contig});
                delete($CON_LINS{$contig});
                delete($CON_LEN{$contig});
                delete($CON_RD{$contig});
                delete($CON_KM{$contig});
                delete($CON_GENES{$contig});
                next;
        }

        foreach my $gene (keys %{$CON_GENES{$contig}}){
                $genlca = "NOVEL";
                $genstr = 0;
                $genname = "NOVEL_GENE";
                $gensco  = 0;
                $gentids = 0;
                if(exists($GENE_SCO{$gene})){
                        #score means it has a hit in the URDB
                        $genlca = $GEN_LCA{$gene};
                        $genstr = keys %{$GENE_LINS{$gene}};
                        %GTIDS=();
                        foreach my $lin (keys %{$GENE_LINS{$gene}}){ $GTIDS{$LIN_TID{$lin}}=1;}
                        @TLIST=(); foreach my $tid (sort(keys %GTIDS)){ push(@TLIST,$tid); }
                        $gentids = join(";",@TLIST);
                        $genname = $GENE_NAME{$gene};
                        $gensco  = $GENE_SCO{$gene};
                }
                $start  = $GENE_POS{$gene};
                $genlen = $GENE_LEN{$gene};
                $genpart = $GENE_PART{$gene};
                $gencov  = $GENE_RD{$gene};
                $genrpkm = $GENE_KM{$gene};
                $genpgc  = $GENE_PGC{$gene};

                %FIDS=(); %RXNS=(); %CRX=(); %CPX=(); %CTX=(); %TCDB=();
                @FUNC_LIST=();
                for my $i (7..24){
                        $j=$i-7;
                        foreach my $id (sort{$GENE_FUNCS{$gene}{$i}{$b}<=>$GENE_FUNCS{$gene}{$i}{$a}} keys %{$GENE_FUNCS{$gene}{$i}}){
                                if($id!~/\w/){next;}
                                $FUNC_LIST[$j].=$id.";";
                                if($i < 22){$FIDS{$id}=1;}
                                if($i=~/19|20|21/){$RXNS{$id}=1;}
                                if($i==11){$TCDB{$id}=1;} #TCDB ONLY USES CHEBI CPDS
                                if($i==22 && $id =~ /CHEBI/){$CRX{$id}=1;}
                                if($i==23 && $id =~ /CHEBI/){$CPX{$id}=1;}
                                if($i==24 && $id =~ /CHEBI/){$CTX{$id}=1;}
                }       }

                #GET COMPOUND~REACTION~TRANSPORTER
                %LINS=();
                foreach my $lin (keys %{$GENE_LINS{$gene}}){if($lin =~/^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)\;.+\;.+\;.+\;.+\;.+\;/ && $lin !~/UNCLASSIFIED[^\;]+GENUS/){$LINS{$lin}=1;}}
                foreach my $id (keys %CTX){ #TCDB CPD
                        foreach my $tcdb (keys %TCDB){
                                $CPD_RPKM{$id}{$gene}{$tcdb}=$genrpkm;
                                foreach my $lin (keys %LINS){$CPD_LIN_TCDB{$id}{$lin}{$tcdb}++;}
                }       }
                foreach my $id (keys %CRX){ #REACTANT CPD
                        foreach my $rxn (keys %RXNS){
                                $CPD_RPKM{$id}{$gene}{$rxn}=$genrpkm;
                                foreach my $lin (keys %LINS){$CPD_LIN_RNT{$id}{$lin}{$rxn}++;}
                }       }
                foreach my $id (keys %CPX){ #PRODUCT CPD
                        foreach my $rxn (keys %RXNS){
                                $CPD_RPKM{$id}{$gene}{$rxn}=$genrpkm;
                                foreach my $lin (keys %LINS){$CPD_LIN_PDT{$id}{$lin}{$rxn}++;}
                }       }

                #GET SUMMARY FUNCTION OUTPUT
                foreach my $id (keys %FIDS){
                        $FUNC_GENES{$id}++;
                        $FUNC_SCOS{$id}+=$gensco;
                        $FUNC_LCA{$id}{$genlca}++;
                        foreach my $lin (keys %{$GENE_LINS{$gene}}){
                                $FUNC_LINS{$id}{$lin}=1;
                                if(!exists($FUNC_CNT{$id}) || $FUNC_CNT{$id} < $gensco){
                                        $FUNC_CNT{$id}=$gensco; $FUNC_TOP{$id}=$lin;
                                }
                        }
                        $FUNC_RPKM{$id}+=$genrpkm;
                }

                #GET GENE FUNCTION OUTPUT
                for my $i (0..17){if($FUNC_LIST[$i]!~/\w/){$FUNC_LIST[$i]=0;}else{$FUNC_LIST[$i]=~s/^\;+|\;+$//g;}}
                $funcs=join("\t",@FUNC_LIST);


                # PRINT OUTPUTS
                print OUTGENE "$gene\t$contig\t$genname\t$genpart\t$gensco\t";
                print OUTGENE "$start\t$genlen\t$conlen\t$congcnt\t";
                print OUTGENE "$genrpkm\t$conrpkm\t";
                print OUTGENE "$gencov\t$concov\t";
                print OUTGENE "$genpgc\t$conpgc\t";
                print OUTGENE "$genstr\t$gentids\t$genlca\t";
                print OUTGENE "$constr\t$contids\t$conlca\t";
                print OUTGENE "$funcs\n";

                delete($GENE_POS{$gene});
                delete($GENE_LEN{$gene});
                delete($GENE_PART{$gene});
                delete($GENE_RD{$gene});
                delete($GENE_KM{$gene});
                delete($GENE_PGC{$gene});
                delete($GENE_SCO{$gene});
                delete($GENE_NAME{$gene});
                delete($GEN_LCA{$gene});
                delete($GENE_LINS{$gene});
                delete($GENE_HITS{$gene});
        }

        #GET CYTO STUFF
        $conlca=~/^(.)/; $king=$1;
        @TAXON=split(";",$conlca);
        for($i=0; $i<=$#TAXON; $i++){
                $type = $king."_".$i;
                if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
                else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
                $CYTO_CCON{$phyla}++;
                $CYTO_CGEN{$phyla}+=$congcnt;
                $CYTO_CLEN{$phyla}+=$conlen;
                $CYTO_RPKM{$phyla}+=$conrpkm;
                $CYTO_CPGC{$phyla}+=$conpgc;
                $CYTO_CCOV{$phyla}+=$concov;
        }

        delete($CON_PGC{$contig});
        delete($CON_LCA{$contig});
        delete($CON_LINS{$contig});
        delete($CON_LEN{$contig});
        delete($CON_RD{$contig});
        delete($CON_KM{$contig});
        delete($CON_GENES{$contig});
}


qx{wget -N https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz};
qx{gunzip -f names.tsv.gz};
open(INCPDNM, "names.tsv")||die;
while(<INCPDNM>){
        if($_!~/\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t", $_);
        $cpd="CHEBI:".$stuff[1];
        $name=$stuff[4];
        $name=~/([A-Z]+)/;
        $len=length($1);
        $name=CleanNames($name);
        if(exists($CPD_NAMES{$cpd})){
                $CPD_NAMES{$cpd}=~/([A-Z]+)/;
                $old=length($1);
                if($len>$old){$CPD_NAMES{$cpd}=$name;}
        }
        else{$CPD_NAMES{$cpd}=$name;}
}


print "Creating Interaction Cytos\n";
foreach my $cpd (keys %CPD_LIN_TCDB){ #loop through all transported compounds
        %T_LINS=();
        %R_LINS=();
        %P_LINS=();
        $clin='';
        @TEST=(); @GOOD=();
        foreach my $lin (keys %{$CPD_LIN_TCDB{$cpd}}){push(@TEST,$lin);}
        while($TEST[0]=~/\w/){
                $clin=shift(@TEST);
                if(!grep(/\Q$clin\E/, @TEST) && !grep(/\Q$clin\E/, @GOOD)){ push(@GOOD,$clin); $T_LINS{$clin}=1;}
        }
        @TEST=(); @GOOD=();
        foreach my $lin (keys %{$CPD_LIN_RNT{$cpd}}){push(@TEST,$lin);}
        while($TEST[0]=~/\w/){
                $clin=shift(@TEST);
                if(!grep(/\Q$clin\E/, @TEST) && !grep(/\Q$clin\E/, @GOOD)){ push(@GOOD,$clin); $R_LINS{$clin}=1;}
        }
        @TEST=(); @GOOD=();
        foreach my $lin (keys %{$CPD_LIN_PDT{$cpd}}){push(@TEST,$lin);}
        while($TEST[0]=~/\w/){
                $clin=shift(@TEST);
                if(!grep(/\Q$clin\E/, @TEST) && !grep(/\Q$clin\E/, @GOOD)){ push(@GOOD,$clin); $P_LINS{$clin}=1;}
        }

        #GET TOTAL INTERACTIONS - at least one must have the compound transporter
        foreach my $pl (keys %P_LINS){
                $c=0;
                if(exists($T_LINS{$pl})){$c=1;}
                foreach my $rl (keys %R_LINS){
                        $d=0;
                        if(exists($T_LINS{$rl})){$d=1;}
                        if($c==1 || $d==1){ $LIN_INTS{$pl}{$rl}{$cpd}=1;}
                }
        }

        #PRODUCT RXNS/TCDB
        %PTC_LINS=(); %PC_LINS=();
        foreach my $lin (keys %P_LINS){
                if(exists($T_LINS{$lin})){print "cpd $cpd PTC lin $lin\n";      $PTC_LINS{$lin}=1; delete($T_LINS{$lin});}
                else{$PC_LINS{$lin}=1; print "cpd $cpd PC lin $lin\n";}
                delete($P_LINS{$lin});
        }
        undef(%P_LINS);

        %RTC_LINS=(); %RC_LINS=();
        foreach my $lin (keys %R_LINS){
                if(exists($T_LINS{$lin})){$RTC_LINS{$lin}=1; delete($T_LINS{$lin});}
                else{$RC_LINS{$lin}=1;}
                delete($R_LINS{$lin});
        }
        undef(%R_LINS);

        %TC_LINS=();
        foreach my $lin (keys %T_LINS){$TC_LINS{$lin}=1; delete($T_LINS{$lin});}
        undef(%T_LINS);


        #CREATE separate trees of product and reactant lineages attached to compound
        $cpdnm='';
        $cpdnm=$cpd."_".$CPD_NAMES{$cpd};
        foreach my $lin (keys %TC_LINS){
                $t_rpkm=0;
                $sco=0.1;
                foreach my $tcdb (keys %{$CPD_LIN_TCDB{$cpd}{$lin}}){foreach $gene (keys %{$CPD_RPKM{$cpd}}){$t_rpkm+=$CPD_RPKM{$cpd}{$gene}{$tcdb}; }}
                $lin=~/^(.)/; $king=$1;
                @TAXON=split(";",$lin);
                for($i=0; $i<=$#TAXON; $i++){
                        $rn="T_".$cpd;
                        $source=$rn."_".$TAXON[$i-1];
                        $target=$rn."_".$TAXON[$i];
                        if($i == 0){
                                $type = $king."_-2";
                                $phyla = "$type\tROOT\t$cpdnm";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_TRN{$phyla}++;
                                $INT_CAT{$phyla}="T";
                                $INT_RPKM{$phyla}+=$p_rpkm;
                        }
                        if($i == 0){
                                $type = $king."_-1";
                                $phyla = "$type\t$cpdnm\t$rn";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_TRN{$phyla}++;
                                $INT_CAT{$phyla}="T";
                                $INT_RPKM{$phyla}+=$p_rpkm;
                        }
                        $type = $king."_".$i;
                        if($i == 0){$phyla = "$type\t$rn\t$target";}
                        else{$phyla = "$type\t$source\t$target";}
                        $INT_SCO{$phyla}+=$sco;
                        $INT_CNT{$phyla}++;
                        $INT_TRN{$phyla}++;
                        $INT_CAT{$phyla}="T";
                        $INT_RPKM{$phyla}+=$p_rpkm;
        }       }

        foreach my $lin (keys %PTC_LINS){
                $p_rpkm=0;
                $sco=1;
                foreach my $tcdb (keys %{$CPD_LIN_TCDB{$cpd}{$lin}}){foreach $gene (keys %{$CPD_RPKM{$cpd}}){$p_rpkm+=$CPD_RPKM{$cpd}{$gene}{$tcdb}; }}
                foreach my $rxn  (keys %{$CPD_LIN_PDT{$cpd}{$lin}}){ foreach $gene (keys %{$CPD_RPKM{$cpd}}){$p_rpkm+=$CPD_RPKM{$cpd}{$gene}{$rxn}; }}
                $lin=~/^(.)/; $king=$1;
                @TAXON=split(";",$lin);
                for($i=0; $i<=$#TAXON; $i++){
                        $rn="P_".$cpd;
                        $source=$rn."_".$TAXON[$i-1];
                        $target=$rn."_".$TAXON[$i];
                        if($i == 0){
                                $type = $king."_-2";
                                $phyla = "$type\tROOT\t$cpdnm";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_TRN{$phyla}++;
                                $INT_CAT{$phyla}="P";
                                $INT_RPKM{$phyla}+=$p_rpkm;
                        }
                        if($i == 0){
                                $type = $king."_-1";
                                $phyla = "$type\t$cpdnm\t$rn";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_TRN{$phyla}++;
                                $INT_CAT{$phyla}="P";
                                $INT_RPKM{$phyla}+=$p_rpkm;
                        }
                        $type = $king."_".$i;
                        if($i == 0){$phyla = "$type\t$rn\t$target";}
                        else{$phyla = "$type\t$source\t$target";}
                        $INT_SCO{$phyla}+=$sco;
                        $INT_CNT{$phyla}++;
                        $INT_TRN{$phyla}++;
                        $INT_CAT{$phyla}="P";
                        $INT_RPKM{$phyla}+=$p_rpkm;
        }       }
        foreach my $lin (keys %PC_LINS){
                $p_rpkm=0;
                $sco=0.5;
                foreach my $rxn  (keys %{$CPD_LIN_PDT{$cpd}{$lin}}){ foreach $gene (keys %{$CPD_RPKM{$cpd}}){$p_rpkm+=$CPD_RPKM{$cpd}{$gene}{$rxn}; }}
                $lin=~/^(.)/; $king=$1;
                @TAXON=split(";",$lin);
                for($i=0; $i<=$#TAXON; $i++){
                        $rn="P_".$cpd;
                        $source=$rn."_".$TAXON[$i-1];
                        $target=$rn."_".$TAXON[$i];
                        if($i == 0){
                                $type = $king."_-2";
                                $phyla = "$type\tROOT\t$cpdnm";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_CAT{$phyla}="P";
                                $INT_RPKM{$phyla}+=$p_rpkm;
                        }
                        if($i == 0){
                                $type = $king."_-1";
                                $phyla = "$type\t$cpdnm\t$rn";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_CAT{$phyla}="P";
                                $INT_RPKM{$phyla}+=$p_rpkm;
                        }
                        $type = $king."_".$i;
                        if($i == 0){$phyla = "$type\t$rn\t$target";}
                        else{$phyla = "$type\t$source\t$target";}
                        $INT_SCO{$phyla}+=$sco;
                        $INT_CNT{$phyla}++;
                        $INT_CAT{$phyla}="P";
                        $INT_RPKM{$phyla}+=$p_rpkm;
        }       }

        foreach my $lin (keys %RTC_LINS){
                $r_rpkm=0;
                $sco=1;
                foreach my $tcdb (keys %{$CPD_LIN_TCDB{$cpd}{$lin}}){foreach $gene (keys %{$CPD_RPKM{$cpd}}){$r_rpkm+=$CPD_RPKM{$cpd}{$gene}{$tcdb}; }}
                foreach my $rxn  (keys %{$CPD_LIN_RNT{$cpd}{$lin}}){ foreach $gene (keys %{$CPD_RPKM{$cpd}}){$r_rpkm+=$CPD_RPKM{$cpd}{$gene}{$rxn}; }}
                $lin=~/^(.)/; $king=$1;
                @TAXON=split(";",$lin);
                for($i=0; $i<=$#TAXON; $i++){
                        $rn="R_".$cpd;
                        $source=$rn."_".$TAXON[$i-1];
                        $target=$rn."_".$TAXON[$i];
                        if($i == 0){
                                $type = $king."_-2";
                                $phyla = "$type\tROOT\t$cpdnm";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_TRN{$phyla}++;
                                $INT_CAT{$phyla}="R";
                                $INT_RPKM{$phyla}+=$r_rpkm;
                        }
                        if($i == 0){
                                $type = $king."_-1";
                                $phyla = "$type\t$cpdnm\t$rn";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_TRN{$phyla}++;
                                $INT_CAT{$phyla}="R";
                                $INT_RPKM{$phyla}+=$r_rpkm;
                        }
                        $type = $king."_".$i;
                        if($i == 0){$phyla = "$type\t$rn\t$target";}
                        else{$phyla = "$type\t$source\t$target";}
                        $INT_SCO{$phyla}+=$sco;
                        $INT_CNT{$phyla}++;
                        $INT_TRN{$phyla}++;
                        $INT_CAT{$phyla}="R";
                        $INT_RPKM{$phyla}+=$r_rpkm;
        }       }
        foreach my $lin (keys %RC_LINS){
                $r_rpkm=0;
                $sco=0.5;
                foreach my $rxn  (keys %{$CPD_LIN_RNT{$cpd}{$lin}}){ foreach $gene (keys %{$CPD_RPKM{$cpd}}){$r_rpkm+=$CPD_RPKM{$cpd}{$gene}{$rxn}; }}
                $lin=~/^(.)/; $king=$1;
                @TAXON=split(";",$lin);
                for($i=0; $i<=$#TAXON; $i++){
                        $rn="R_".$cpd;
                        $source=$rn."_".$TAXON[$i-1];
                        $target=$rn."_".$TAXON[$i];
                        if($i == 0){
                                $type = $king."_-2";
                                $phyla = "$type\tROOT\t$cpdnm";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_CAT{$phyla}="R";
                                $INT_RPKM{$phyla}+=$r_rpkm;
                        }
                        if($i == 0){
                                $type = $king."_-1";
                                $phyla = "$type\t$cpdnm\t$rn";
                                $INT_SCO{$phyla}+=$sco;
                                $INT_CNT{$phyla}++;
                                $INT_CAT{$phyla}="R";
                                $INT_RPKM{$phyla}+=$r_rpkm;
                        }
                        $type = $king."_".$i;
                        if($i == 0){$phyla = "$type\t$rn\t$target";}
                        else{$phyla = "$type\t$source\t$target";}
                        $INT_SCO{$phyla}+=$sco;
                        $INT_CNT{$phyla}++;
                        $INT_CAT{$phyla}="R";
                        $INT_RPKM{$phyla}+=$r_rpkm;
}       }       }


open(OUTRANS, ">", $outtrans)||die "unable to open $outtrans: $!\n";
print OUTRANS "type\tSource\tTarget\tCategory\tLin_Cnt\tLin-w-Tr_Cnt\tType_sum_scor\tsum_gene_rpkm\n";
foreach my $phyla (sort(keys %INT_SCO)){
        $sco=$INT_SCO{$phyla};
        $cat=$INT_CAT{$phyla};
        $cnt=$INT_CNT{$phyla};
        $trn=$INT_TRN{$phyla};
        $rpkm=$INT_RPKM{$phyla};
        print OUTRANS "$phyla\t$cat\t$cnt\t$trn\t$sco\t$rpkm\n";
}


open(OUTMRX, ">", $outmrx)||die;
print OUTMRX "PRODUCT_LIN\tREACTANT_LIN\tTotal_Ints\tP-R_Ints\tPR_Cpd_Names\tPR_Cpd_Ids\tR-P_Ints\tRP_Cpd_Names\tRP_Cpd_Ids\n";
foreach my $pl (sort(keys %LIN_INTS)){
        foreach my $rl (sort(keys %{$LIN_INTS{$pl}})){
                @FNAMES=(); @FCPDS=(); @RNAMES=(); @RCPDS=();
                foreach my $cpd (sort(keys %{$LIN_INTS{$pl}{$rl}})){push(@FNAMES,$CPD_NAMES{$cpd}); push(@FCPDS,$cpd);}
                $fnames=join(";",@FNAMES);
                $fcpds=join(";",@FCPDS);
                $fcc=@FCPDS;
                foreach my $cpd (sort(keys %{$LIN_INTS{$rl}{$pl}})){push(@RNAMES,$CPD_NAMES{$cpd}); push(@RCPDS,$cpd);}
                $rnames=join(";",@RNAMES);
                $rcpds=join(";",@RCPDS);
                $rcc=@RCPDS;
                $tot=$fcc+$rcc;
                print OUTMRX "$pl\t$rl\t$tot\t$fcc\t$fnames\t$fcpds\t$rcc\t$rnames\t$rcpds\n";
}       }



print "output cytoscape file\n";
open(OUTCYTO, ">", $outcyto)||die;
print OUTCYTO "type\tsource\ttarget\ttot_contigs\ttot_genes\tsum_contig_length\tsum_rpkm\tavg_read_depth\tavg_%_GC\n";
foreach my $phyla (sort(keys %CYTO_CCON)){
        $csumcon=$CYTO_CCON{$phyla};
        $csumtgc=$CYTO_CGEN{$phyla};
        $csumlen=$CYTO_CLEN{$phyla};
        $csumrkm=$CYTO_RPKM{$phyla}; $csumrkm=~s/(\.\d\d).*/$1/;
        $cavgcov=$CYTO_CCOV{$phyla}/$csumcon; $cavgcov=~s/(\.\d\d).*/$1/;
        $cavgpgc=$CYTO_CPGC{$phyla}/$csumcon; $cavgpgc=~s/(\.\d\d).*/$1/;
        print OUTCYTO "$phyla\t$csumcon\t$csumtgc\t$csumlen\t$csumrkm\t$cavgcov\t$cavgpgc\n";
}
undef(%CYTO_CCON);
undef(%CYTO_CGEN);
undef(%CYTO_CLEN);
undef(%CYTO_RPKM);
undef(%CYTO_CCOV);
undef(%CYTO_CPGC);





print "output function summary\n";
open(OUTFUNC,">",$outfuncs)||die;
print OUTFUNC "fid\tname\tgene_count\tsum_rpkm\tstrains\tsum_score\ttop_score\ttop_lin\tlca\n";
foreach my $id (sort(keys %FUNC_GENES)){
        $name = $FUNC_NM{$id};
        $fgc  = $FUNC_GENES{$id};
        $fsco = $FUNC_SCOS{$id};
        $frkm = $FUNC_RPKM{$id};
        $flc  = keys %{$FUNC_LINS{$id}};
        $toplin = $FUNC_TOP{$id};
        $topsco = $FUNC_CNT{$id};
        @PHYL=();
        foreach my $lca (keys %{$FUNC_LCA{$id}}){ push(@PHYL, $lca); }
        $flca=MakeLCA(@PHYL);
        print OUTFUNC "$id\t$name\t$fgc\t$frkm\t$flc\t$fsco\t$topsco\t$toplin\t$flca\n";
}






##########################################
############# SUBROUTINES ################
##########################################
sub GetConGenPos{
        $/=">";
        %hash1=();
        %hash2=();
        %hash3=();
        %hash4=();
        %hash5=();
        print "input $ingen\n";
        open(INGEN,$ingen)||die;
        while(<INGEN>){
                if($_ !~ /\w/){next;}
                $_=uc($_);
                @stuff=split("\n",$_);
                $header = shift(@stuff);
                $header =~ /^((\S+)\_\d+)[\s\#]+(\d+)[\s\#]+(\d+)[\s\#]+.*PARTIAL\=(\d+)/i;
                $gene   =$1;
                $contig =$2;
                $start  =$3;
                $len    =$4-$3;
                $part   =$5;
                if($gene!~/\w/ || $contig!~/\w/ || $start!~/\w/ || $len!~/\w/ || $part !~/\w/){print "gene $gene contig $contig start $start len $len part $part\nmissing something $header\n";die;}


                $seq=join("",@stuff);
                $seq=~s/[^A-Z]//g;
                my $count = $seq =~ tr/GC/gc/;
                $len=length($seq);
                $prgc=$count*100/$len;
                $prgc =~ s/(\.\d\d).*/$1/;

                $hash1{$contig}{$gene}=1;
                $hash2{$gene}=$start;
                $hash3{$gene}=$part;
                $hash4{$gene}=$prgc;
                $hash5{$gene}=$len;
        }
        close(INGEN);
        return(\%hash1, \%hash2, \%hash3, \%hash4, \%hash5);
}


sub GetGeneReadCounts{
        $/="\n";
        %hash1=();
        %hash2=();
        print "input $ingrc\n";
        open(INGRC,$ingrc)||die;
        while(<INGRC>){
                if($_ !~ /\w/){next;}
                if($_ =~ /^\#/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                @stuff=split("\t",$_,-1);
                $depth=$stuff[3];
                $depth =~ s/(\.0*\d\d).*/$1/;
                $rpkm=$stuff[5];
                $rpkm =~ s/(\.0*\d\d).*/$1/;
                $stuff[0] =~ s/\s.*//;
                $hash1{$stuff[0]}=$depth;
                $hash2{$stuff[0]}=$rpkm;
        }
        close(INGRC);
        return(\%hash1, \%hash2);
}


sub GetContigPGC{
        $/=">";
        %hash1=();
        %hash2=();
        print "input $incon\n";
        open(INCON, $incon)||die;
        while(<INCON>){
                if($_ !~ /\w/){next;}
                $_=uc($_);
                @stuff=split("\n",$_);
                $contig=shift(@stuff);
                $seq=join("",@stuff);
                $seq=~s/[^A-Z]//g;
                $gc = $seq =~ tr/GC/gc/;
                $len = length($seq);
                $pgc = $gc*100/$len;
                $pgc =~ s/(\.\d\d).*/$1/;
                $hash1{$contig}=$pgc;
                $hash2{$contig}=$len;
        }
        close(INCON);
        return(\%hash1, \%hash2);
}


sub GetContigCov{
        $/="\n";
        %hash1=();
        %hash2=();
        print "input $inccov\n";
        open(INCC, $inccov)||die;
        while(<INCC>){
                if($_ !~ /\d/){next;}
                $_=uc($_);
                $_ =~ s/[\r\n]//;
                @stuff = split("\t", $_, -1);
                $contig=$stuff[0];
                if(!exists($CON_PGC{$contig})){next;}
                $depth=$stuff[3];
                $depth =~ s/(\.0*\d\d).*/$1/;
                $hash1{$contig}=$depth;
                $rpkm=$stuff[5];
                $rpkm =~ s/(\.0*\d\d).*/$1/;
                $hash2{$contig}=$rpkm;
        }
        close(INCC);
        return(\%hash1, \%hash2);
}


sub GetFuncNames{
        $/="\n";
        %hash1=();
        print "input $infn\n";
        open(INFN, $infn)||die;
        while(<INFN>){
                if($_ !~ /\w/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                (my $id, my $name)=split("\t",$_,-1);
                $id=~s/\s//g;
                $hash1{$id}=$name;
        }
        close(INFN);
        return(\%hash1);
}


sub GetTaxonInfo{
        $/="\n";
        %hash1=();
        %hash2=();
        print "input $intax\n";
        open(INTAX, $intax)||die;
        while(<INTAX>){
                if($_ !~ /\w/){next;}
                $_ =~ s/[\r\n]//;
                $_=uc($_);
                @stuff = split("\t", $_, -1);
                $tid =  shift(@stuff);
                $lin = join(";",@stuff);
                if(!exists($hash2{$lin})){$hash2{$lin}=$tid;}
                $hash1{$tid}=$lin;
        }
        close(INTAX);
        return(\%hash1, \%hash2);
}


sub InDiamond{
        $/="\n";
        %hash1=();
        %hash2=();
        %hash3=();
        print "input $ingvu\n";
        open(INGVU, $ingvu)||die;
        while(<INGVU>){
                if($_ !~ /\w/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                @stuff=split("\t",$_,-1);
                $gene = $stuff[0];
                $hit  = $stuff[2];
                $pid  = $stuff[9];
                if($pid < $minid){next;}
                $cov  = $stuff[11];
                if($cov < $mincov){next;}
                $sco  = $pid*$cov;

                if(!exists($hash3{$gene})){$hash3{$gene}=$sco;}
                if($sco > $hash3{$gene}){
                        $hash3{$gene}=$sco;  #set new top score hit for the gene
                        foreach my $h (keys %{$hash1{$gene}}){
                                #loop through prior hits for the gene, remove any not within the new top 90%
                                if($hash1{$gene}{$h} < $hash3{$gene}*$top){ delete($hash1{$gene}{$h}); delete($hash2{$h}{$gene}); }
                        }
                }
                if($sco >= $hash3{$gene}*$top){
                        $hash1{$gene}{$hit}=$sco;
                        $hash2{$hit}{$gene}=$sco;
                }
        }
        close(INGVU);
        return(\%hash1, \%hash2, \%hash3);
}


sub GetRefGeneInfo{
        $/="\n";
        %hash1=();
        %hash2=();
        %hash3=();
        %hash4=();
        $on=0;
        if($ininfo=~/\.gz$/){open(INFO, "gunzip -c $ininfo |")||die; }
        else{open(INFO,$ininfo)||die;}
        while(<INFO>){
                if($_ !~ /\w/){next;}
                $_=uc($_);
                $_=~s/[\r\n]+//;
                @stuff=split("\t",$_,-1);
                $hit = $stuff[0];
                $name = $stuff[1];

                #functions
                $allfuncs='';
                if(exists($HIT_GENES{$hit})){
                        @IDS=();
                        $stuff[12]=~s/\|/\;/g;
                        for my $i (7..24){
                                if($stuff[$i]!~/\w/){next;}
                                if($i==8){ $IDS[0]=$stuff[$i]; }
                                else{@IDS=split(";", $stuff[$i]);}
                                foreach my $id (@IDS){
                                        if($i==10){$id="METAL-BINDING:".$id;}
                                        if($i==12){$id="LOCATION:".$id;}
                                        $allfuncs.=$id.";";
                                        foreach my $gene (keys %{$HIT_GENES{$hit}}){ $hash3{$gene}{$i}{$id}++; }
                }       }       }

                #taxa
                @tids = split(";",$stuff[5]); %LINS=();
                foreach my $tid (@tids){
                        #use lineage because duplicate lins with different tids
                        $lin=$PHY{$tid};
                        if($lin !~ /\w/){print DEBUG "missing $tid\n"; next;}
                        if(exists($LINS{$lin})){next;}
                        $LINS{$lin}=1;
                        $hash5{$lin}++;
                        if(!exists($HIT_GENES{$hit})){next;}
                        $hash6{$hit}{$lin}=1;
                        foreach my $gene (keys %{$HIT_GENES{$hit}}){
                                $sco=$HIT_GENES{$hit}{$gene};
                                $hash1{$lin}{$gene}=$sco;       #total contig genes for tid
                                $hash2{$gene}{$lin}=$sco;       #for looping thru gene tids
                                $hash4{$gene}=$name;
                }       }

                if(exists($HIT_GENES{$hit})){delete($HIT_GENES{$hit});}

                if($on%100000==0){
                        $hr = keys %HIT_GENES;
                        $time = localtime;
                        print "on $on GetRefGeneInfo remaining $hr hit $hit time $time\nfuncs $allfuncs\n";
                } $on++;
                if(keys %HIT_GENES < 1){last;}
        }
        undef(%HIT_GENES);
        return(\%hash1, \%hash2, \%hash3, \%hash4, \%hash5, \%hash6);
}


sub CleanLineages{
        print "summarize lineage gene scores\n";
        %hash1=();
        %AG=();
        %TOSS=();

        #### REMOVE LOW QUALITY LINEAGES ####
        $slkc = keys %LIN_GENES;
        %AG  =%GENE_LINS;
        %TOSS=%LIN_GENES;
        foreach my $lin (keys %LIN_GENES){
                $uniq=0;
                foreach my $gene (keys %{$LIN_GENES{$lin}}){
                        $hash1{$lin}+=$LIN_GENES{$lin}{$gene}; #score
                        $gkc = keys %{$GENE_LINS{$gene}};
                        #has only 1 matching organism and is very high score hit and is >100aa (longer than single read)
                        if($gkc==1 && $LIN_GENES{$lin}{$gene} >= 9900 && $GENE_LEN{$gene}>=100){ $uniq=1; }
                }
                if($uniq>0){
                        if(keys %{$LIN_GENES{$lin}} >= $mingen || $lin =~ /^MONA/ ){
                                delete($TOSS{$lin});
                                foreach my $gene (keys %{$LIN_GENES{$lin}}){delete($AG{$gene});}
                        }
                }
        }
        $akc = keys %AG;
        $tokc = keys %TOSS;

        #now toss only has orgs that don't have a unique gene with a high score
        #%AG only contains genes from organisms without unique genes
        foreach my $lin ( sort{$hash1{$b}<=>$hash1{$a}} keys %TOSS ){
                $hasag=0;
                if(keys %AG > 0){
                        #delete %AG in order of best scoring lineage
                        foreach my $gene (keys %{$LIN_GENES{$lin}}){
                                if(exists($AG{$gene})){ delete($TOSS{$lin}); delete($AG{$gene}); $hasag=1;}
                        }
                }
                if($hasag==0){ #organism has no unique genes and scores lower for all other genes
                        foreach my $gene (keys %{$LIN_GENES{$lin}}){ #remove lineage from gene
                                if(keys %{$GENE_LINS{$gene}} > 1){ delete($GENE_LINS{$gene}{$lin}); }
                        }
                        delete($LIN_GENES{$lin});
                        delete($TOT_LIN_GENES{$lin});
                        delete($hash1{$lin});
                }
        }

        #now screen out genes bad hits
        foreach my $gene (keys %GENE_HITS){
                foreach my $hit (keys %{$GENE_HITS{$gene}}){
                        $goodhit=0;
                        foreach my $lin (keys %{$HIT_LINS{$hit}}){
                                if(exists($GENE_LINS{$gene}{$lin})){$goodhit=1;}
                        }
                        if($goodhit==0){delete($HIT_LINS{$hit}); delete($GENE_HITS{$gene}{$hit});}
                }
        }
        undef(%HIT_LINS);
        $elkc = keys %LIN_GENES;
        $hkc  = keys %hash1;
        print "startlins $slkc nouniqlins $tokc nouniqlins genes $akc endlin_genes $elkc linscor $hkc\n";
        return(\%hash1);
}



sub GetContigLCA{
        print "Clean Up Gene Hits and Get Contig LCA\n";
        %hash1=(); #$hash1{$contig}=$conlca;
        %hash2=(); #$hash2{$gene}=$genelca;
        %hash3=(); #$hash3{$contig}{$lin}=1;

        #clean up gene hits using whole contig, get contig lca
        foreach my $contig (sort(keys %CON_PGC)){
                $gc = keys %{$CON_GENES{$contig}};

                # CONTIG HAS NO GENES
                if(keys %{$CON_GENES{$contig}} < 1){ $nogencon++; $hash1{$contig}="NOVEL"; }

                # ADD NON-NOVEL GENES TO GLIST
                %GLIST=();
                foreach my $gene (keys %{$CON_GENES{$contig}}){
                        if(keys %{$GENE_HITS{$gene}} > 0){$GLIST{$gene}=1;}
                        else{$hash2{$gene}="NOVEL";}
                }

                # CONTIG HAS ONLY NOVEL GENES
                if(keys %GLIST < 1){ $novelonly++; $hash1{$contig}="NOVEL"; next;}

                # MAKE GENE LCA
                @CONLINS=(); @CMONA=();
                foreach my $gene (keys %GLIST){
                        @PHYL=(); @PMONA=();
                        foreach my $lin (keys %{$GENE_LINS{$gene}}){
                                if(!exists($LIN_GENES{$lin})){next;}
                                if($lin =~ /^MONA/){ push(@PMONA, $lin); push(@CMONA, $lin); }
                                               else{ push(@PHYL,$lin); push(@CONLINS, $lin); }
                                $hash3{$contig}{$lin}=1;
                        }

                        # if gene has MONA then it is MONA gene
                        if($PMONA[0] =~ /\w/){$genelca=MakeLCA(@PMONA);}
                        else{$genelca=MakeLCA(@PHYL);}
                        $hash2{$gene}=$genelca;
                }

                # MAKE CONTIG LCA
                #if every gene on contig is viral gene - then viral contig
                   if($CONLINS[0]!~/\w/){$conlca = MakeLCA(@CMONA); $monacontig++;}
                elsif($CONLINS[0]=~/\w/){$conlca = MakeLCA(@CONLINS);}
                else{$conlca="NOVEL";}
                $hash1{$contig}=$conlca;
        }
        print "contigs no genes: $nogencon\n";
        print "contigs only novel genes: $novelonly\n";
        print "contigs from MONA: $monacontig\n";
        return(\%hash1, \%hash2, \%hash3);
}


sub MakeLCA{
        my @ARR;
        %seen=();
        @array1 = @_;
        @array1 = grep { !$seen{$_}++ } @array1;

        #get the kingdoms, JIC lineage is NCA
        %LET=();
        foreach my $lin (@array1){
                if($lin !~ /^(BACTERIA|ARCHAEA|MONA|EUKARYOTA)/i){next;}
                $lin =~ /^(\w)/; $LET{$1}=1;
        }
        @LET=();
        foreach my $let (sort(keys %LET)){push(@LET,$let);}
        $let=join("",@LET);

        $LCA='';
        $len1 = @array1;
        if($len1 == 1){$LCA = $array1[0]; }
        elsif($len1 > 1){
                $first = $array1[0];
                @levels = split(";", $first);
                for($i=0; $i<=$#levels; $i++){
                        $alevel=$levels[$i];
                        @matched = grep(/\Q$alevel\E/i, @array1);
                        $len2 = @matched;
                        if($len2 == $len1){push(@ARR, $alevel);}
                        else{last;}
                }
                $len3 = @ARR;
                if($len3 > 1){$LCA = join(";", @ARR);}
                elsif($len3==1){$LCA = $ARR[0];}
                else{$LCA = "NCA"; }
        }
        else{$LCA = "NCA"; }

        #add kingdoms to NCA
        if($LCA eq "NCA"){ $LCA.="-".$let; }
        return($LCA);
}

sub CleanNames{
        $name = $_[0];

        #remove pointless ambiguators
        $name =~ s/(CANDIDATUS|CANDIDATUAS|CANDIDATE|VOUCHERED|UNDESCRIBED|UNSCREENED|UNKNOWN|UNCULTIVATED|UNCULTURED)\s*/UNCHARACTERIZED\_/g;
        $name =~ s/(UNIDENTIFIED|UNCLASSIFIED|CONTAMINATION|SCREENED|UNASSIGNED|PUTATIVE|HYPOTHETICAL)\s/UNCHARACTERIZED\_/g;
        $name =~ s/(UNCHARACTERIZED\_)+/UNCHARACTERIZED\_/g;
        $name =~ s/\-*LIKE\s*/\_/g;

        #remove junk punctuation/standardize
        $name =~ s/\s+/_/g;
        $name =~ s/[^\w]+/_/g;
        $name =~ s/\_+/\_/g;
        $name =~ s/(^\_+|\_+$)//g;

        return($name);
}



sub PrintKC{
print "1 $kc1\n";
print "2 $kc2\n";
print "3 $kc3\n";
print "4 $kc4\n";
print "5 $kc5\n";
print "6 $kc6\n";
print "7 $kc7\n";
print "8 $kc8\n";
print "9 $kc9\n";
print "10 $kc10\n";
print "11 $kc11\n";
print "12 $kc12\n";
print "13 $kc13\n";
print "14 $kc14\n";
print "15 $kc15\n";
print "16 $kc16\n";
print "17 $kc17\n";
print "18 $kc18\n";
print "19 $kc19\n";
print "20 $kc20\n";
print "21 $kc21\n";
print "22 $kc22\n";
print "23 $kc23\n";
print "24 $kc24\n";
print "25 $kc25\n";
print "26 $kc26\n";
print "27 $kc27\n";
#print "28 $kc28\n";
#print "29 $kc29\n";
#print "30 $kc30\n";
#print "31 $kc31\n";
#print "32 $kc32\n";
#print "33 $kc33\n";
#print "34 $kc34\n";
}
