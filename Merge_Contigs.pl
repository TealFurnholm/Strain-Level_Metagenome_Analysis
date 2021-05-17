#BEFORE YOU BEGIN
#USE bbtools dedupe to remove containments and duplicates and orient the sequences into overlaps
   #   comics dedupe in=multi_assembly_contigs.fasta out=dedup.fasta tuc mid=99 e=2 minscaf=200 fo mo=200 c pc=t cc=t fcc=t fmj=t absorbcontainment=t numaffixmaps=3 overwrite=t
#USE bbtools again on the dedupe output to generate the .dot file - use more  permissive edit to allow for indels, will screen later
   #   comics dedupe in=dedup.fasta out=dedup2.fasta tuc e=10 rnc=t ngn=t fo c pc=f absorbcontainment=f mo=200 numaffixmaps=3 dot=graph.dot overwrite=t
#the .dot file outputs based on input contigs and will not keep contig names if you cluster (bug in dedupe)



#use File::Copy;
#use warnings;
$userpid = 0.995;
$percid = 1-$userpid;
$sample = $ARGV[0];
$sample =~ s/\s+//;

$incon  = $sample."_dedup4.fa";
$indot  = $sample."_graph4.dot";
$prod   = $sample."_prod_X.fa";
$output = $sample."_MERGED_CONTIGS.fa";
$log    = $sample."_MERGED_CONTIGS.log";
$debug  = $sample."_MERGED_CONTIGS.debug";
open(INCON, $incon)||die;
open(INDOT, $indot)||die;
open(OUTPUT, ">", $output)||die;
open(LOG, ">", $log)||die;
open(DEBUG, ">", $debug)||die;


#FIRST: GET OVERLAP INFO AND SEQ IDS
$on=0;
$max=0;
$totseqlen=0;
$/=">";
while(<INCON>){
        if($_ !~ /\w/){next;}
        $on++;
        $_ = uc($_);
        @stuff = split("\n",$_);
        $header = shift(@stuff);
        $header =~ /(CLUSTER\_(\d+)\,CONTIG\_(\d+))/;
        $contig = $1;
        $cluster = $2;
        $connum = $3;
        $seq = join("", @stuff);
        $seq =~ s/[^A-Z]//g;
        $ID2SEQ{$cluster}{$connum}=$seq;
        $CON2SEQ{$contig}=$seq;
        $ALLSEQS{$seq}=$cluster;
        $NAME{$seq}= "CLUSTER-".$cluster."_".$connum;
        $totseqlen+=length($seq);
        if(length($seq)>$max){$max=length($seq);}
}
$sc = keys %NAME;
$ckc = keys %ID2SEQ;
$avglen = $totseqlen/$on;
print    "contigs:\t$on\ndistinct:\t$sc\nlongest:\t$max\naverage:\t$avglen\ntotallen:\t$totseqlen\nclusters:\t$ckc\n";
print LOG "contigs:\t$on\ndistinct:\t$sc\nlongest:\t$max\naverage:\t$avglen\ntotallen:\t$totseqlen\nclusters:\t$ckc\n";



#GET OVERLAPS
$/="\n";
open(INDOT, $indot)||die;
while(<INDOT>){
        if($_ !~ /\-\>/){next;}
        $_ = uc($_);
        $line=$_;
        $_ =~ /.*(CLUSTER\_(\d+)\S+)\s.*(CLUSTER\S+)\s.*\=\"(.+)\"/i;
        $con1 = $1;
        $cluster = $2;
        $con2 = $3;
        @stuff = split(",",$4);
        $seq1 = $CON2SEQ{$con1}; if($seq1 !~ /^[A-Z]+$/){ print DEBUG "INDOT: $con1 no matching seq in fasta line $line \n"; die; }
        $seq2 = $CON2SEQ{$con2}; if($seq2 !~ /^[A-Z]+$/){ print DEBUG "INDOT: $con2 no matching seq in fasta line $line \n"; die; }
        $s1len = length($seq1)-1;
        $s2len = length($seq2)-1;
        #$type = $stuff[0];
        $olen = $stuff[1];
        $edit = $stuff[2];
        $maxbad = $olen*$percid;
        if($edit > $maxbad){ $toobad++; next; }

        if($stuff[4] < $stuff[6] && $stuff[5] < $stuff[7]){
                $CLUSTERS{$cluster}{$seq2}{$seq1}=$stuff[6]."|".$stuff[7]."|".$s1len."|".$stuff[4]."|".$stuff[5]."|".$s2len."|".$edit."|".$olen;
                $KEEP{$cluster}{$seq1}++; $KEEP{$cluster}{$seq2}+=0;
                $ML{$seq1}=1; $ML{$seq2}=1;
        }
        elsif($stuff[4] > $stuff[6] && $stuff[5] > $stuff[7]){
                $CLUSTERS{$cluster}{$seq1}{$seq2}=$stuff[4]."|".$stuff[5]."|".$s1len."|".$stuff[6]."|".$stuff[7]."|".$s2len."|".$edit."|".$olen;
                $KEEP{$cluster}{$seq2}++; $KEEP{$cluster}{$seq1}+=0;
                $ML{$seq1}=1; $ML{$seq2}=1;
        }
        else{ $notoverlap++; next; }
}
$ckc = keys %CLUSTERS;
$skc = keys %ML;

print LOG "toobad $toobad nooverlap $notoverlap olclusters $ckc olseqs $skc\n";
print    "toobad $toobad nooverlap $notoverlap olclusters $ckc olseqs $skc\n";
undef(%CON2SEQ);



#OUTPUT SEQS WITHOUT ANY OVERLAPS
$cntout=0;
open(OUTPUT, ">>", $output)||die;
open(LOG, ">>", $log)||die;
foreach my $seq (sort(keys %ALLSEQS)){
        if(!exists($ML{$seq})){
                $cntout++;
                print LOG "type 1 $NAME{$seq}\n";
                print OUTPUT ">$NAME{$seq}\n$seq\n";
                delete($ALLSEQS{$seq});
                delete($ID2SEQ{$ALLSEQS{$seq}});
                delete($NAME{$seq});
        }
}
$akc = keys %ALLSEQS;
undef(%ML);
print LOG "output $cntout unoverlaped contigs, there are $akc contigs remaining to be merged\n";
print "output $cntout unoverlaped contigs, there are $akc contigs remaining to be merged\n";



#THIRD: FOR EACH CLUSTER - MERGEOVERLAP
print "looping through clusters\n";
open(OUTPUT, ">>", $output)||die;
open(LOG, ">>", $log)||die;
$on=0;
$totclust = keys %CLUSTERS;
$starttime = localtime;
foreach my $cluster (sort{$CLUSTERS{$a} <=> $CLUSTERS{$b}} keys %CLUSTERS){

        $time = localtime;
        if($on%100000==0){print "cluster $cluster on $on of $totclust start $starttime time $time perfect $perfol prodigal $prodol\n";}
        $on++;

        #### ~ SCREEN %KEEP FOR HEADS ~ ####
        $min=100000000; %HEADS=();
        foreach my $seq (sort{$KEEP{$cluster}{$a} <=> $KEEP{$cluster}{$b}} keys %{$KEEP{$cluster}}){
                if($KEEP{$cluster}{$seq} < $min ){ $min = $KEEP{$cluster}{$seq}; }
                if($KEEP{$cluster}{$seq}==$min){$HEADS{$seq}=1;}
                else{last;}
        }
        $sthedkc = keys %HEADS;

       #### ~  MERGE ALL OVERLAPS ~ ####
        %DONE=();
        while(keys %HEADS > 0){
                foreach my $seqx (keys %HEADS){
                        my @SEQ1 = split("",$seqx);
                        $NAME{$seqx} =~ /(\d+)$/; #grab end contig number
                        $id = $1;
                        $seq1 = $ID2SEQ{$cluster}{$id};
                        $f1=0;
                        foreach my $seq2 (keys %{$CLUSTERS{$cluster}{$seq1}}){
                                my @SEQ2 = split("",$seq2);
                                $NAME{$seq2} =~ /(\d+)$/; #grab contig number
                                $id2 = $1;

                                #GET OVERLAP - 726|1419|1419|0|693|739|1|694
                                $regex1 = $CLUSTERS{$cluster}{$seq1}{$seq2};
                                $regex1 =~ /^(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|(\d+)/;
                                $for = $4; $fiv = $5; $six = $6; $olr = $8;
                                $tre = $#SEQ1; $two = $#SEQ1; $one = $tre-$olr+1;

                                #CHECK BEST CONSERVED AND MERGE
                                if($one > 0){$header = join('', @SEQ1[0..$one-1]);} else{$header = '';}
                                $mid1 = join('', @SEQ1[$one..$two]);
                                $mid2 = join('', @SEQ2[$for..$fiv]);
                                if($fiv < $six){$trailer = join('', @SEQ2[$fiv+1..$six]);} else{$trailer = '';}
                                $ns1 = $header.$mid1.$trailer;
                                $ns2 = $header.$mid2.$trailer;


                                #IF OL IS PERFECT, ELSE KEEP BOTH
                                if($mid1 eq $mid2){$perfol++; $NAME{$ns1} = $NAME{$seqx}."|".$id2; $HEADS{$ns1}=1; if(keys %{$CLUSTERS{$cluster}{$seq2}} < 1){$DONE{$ns1}=$NAME{$ns1};}}
                                else{   $prodol++;
                                        $NAME{$ns1} = $NAME{$seqx}."|a".$id2; $HEADS{$ns1}=1;
                                        $NAME{$ns2} = $NAME{$seqx}."|b".$id2; $HEADS{$ns2}=1;
                                        if(keys %{$CLUSTERS{$cluster}{$seq2}} < 1){$DONE{$ns1}=$NAME{$ns1}; $DONE{$ns2}=$NAME{$ns2};}
                                }
                                $f1=1;
                                delete($ALLSEQS{$seq1}); #makes sure we know which contigs really got merged
                                delete($ALLSEQS{$seq2});
                        }
                        if($f1==0){ $DONE{$seqx}=$NAME{$seqx};} #were no more tails to merge
                        delete($HEADS{$seqx}); #EITHER HAD OVERLAP SO IS NOW A NEW HEAD OR NO OVERLAP SO IS DONE = DELETE
                }
        }

        #OUTPUT MERGED SEQS
        $dkc = keys %DONE;
        foreach my $seq (keys %DONE){
                $name = $NAME{$seq};
                print OUTPUT ">$name\n$seq\n";
                print LOG "cluster $cluster type:2 start:$skc heads:$sthedkc merged:$dkc final $name\n";
        }
        delete($KEEP{$cluster});
        delete($CLUSTERS{$cluster});
        delete($ID2SEQ{$cluster});
}

$akc = keys %ALLSEQS;
print "total $akc unprocessed contigs remaining after merge\nperfect-ol $perfol\nbranches $prodol\n";
print LOG "total $akc unprocessed contigs remaining after merge\nperfect-ol $perfol\nbranches $prodol\n";
foreach my $seq (keys %ALLSEQS){
        print OUTPUT ">$NAME{$seq}\n$seq\n";
        print LOG "$NAME{$seq} type:4 residual unmerged for debug\n";
}





##### ~~ SUBROUTINES ~~ #####
#############################

sub do_prodigal{
        my $nsX = $_[0];
        open(PROD, ">", $prod)||die;
        print PROD ">nsX\n$nsX\n";
        $compl=0; $conf=0; $nlen =0;
        $Q = qx(comics prodigal -c -q -p meta -f gff -i $prod);
        @Q = split("\n", $Q);
        foreach my $l (@Q){
                if($l =~ /^\#/){next;}
                $l =~ /^(\S+)\s+\S+\s+\S+\s+(\d+)\s+(\d+).*partial\=(\d+).*conf\=(\d+)/i;
                $len = $3-$2;
                if($4 !~ /1/){$compl++;}
                $conf+=$5;
                $nlen += $len;
        }
        return($conf, $compl, $nlen);
}






