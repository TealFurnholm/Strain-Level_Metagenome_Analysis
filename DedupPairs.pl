#use warnings;

$input1 = $ARGV[0];
$input2 = $ARGV[1];
$output1 = $ARGV[2];
$output2 = $ARGV[3];
$log = $output1;
$log =~ s/\_*(fwd|for|1\_).*?$/_dedup.log/i;
open(LOG, ">", $log)||die;


if($ARGV[0] !~ /\w/){
        print "Usage: perl DedupPairs.pl in_fwd.fastq.gz in_rev.fastq.gz out_fwd.fastq.gz out_rev.fastq.gz\n";
        print "Note: files may be .fastq or .fastq.gz [compressed]\n";
        die;
}

if($input1 =~ /\.gz$/i){
        open(INPUT1, "gunzip -c $input1 |") or die "gunzip $input1: $!";
        open(INPUT2, "gunzip -c $input2 |") or die "gunzip $input2: $!";
}
else{   open(INPUT1, $input1)||die; open(INPUT2, $input2)||die; }

if($output1 =~ /\.gz$/i){
        open(OUTPUT1,'>:gzip', $output1)||die;
        open(OUTPUT2,'>:gzip', $output2)||die;}
else{   open(OUTPUT1, ">", $output1)||die;
        open(OUTPUT2, ">", $output2)||die;}


$on=0;
$start = localtime;
$badline=0;
@RAND=("A","T","G","C");
while( my $line = <INPUT1> . <INPUT1> . <INPUT1> . <INPUT1> . <INPUT2> . <INPUT2> . <INPUT2> . <INPUT2> ) {
        if($badline > 3){print "on $on File parsing issue, check file format\n"; die;}

        $totreads++;
        $line = uc($line);
        @stuff = split("\n", $line);
        if($stuff[1] !~ /^[A-Z]+$/ || $stuff[5] !~ /^[A-Z]+$/){print "non-nucleotide detected on $on at $stuff[1] \n or \n $stuff[5]\n"; $badline++; next;}
        $stuff[4] =~ /^(\S+)(\s|\_).*?/;
        $read=$1;
        if($stuff[0] !~ /$read/ || $read !~ /\w/){ print "paired name issues $stuff[0] ne $stuff[4] read $read badline $badline\n"; $badline++; next;}
        if(exists($SEQS{$stuff[1]}{$stuff[5]})){$dup++; next;}

        $SEQS{$stuff[1]}{$stuff[5]}=$read;
        $SCOS{$stuff[1]}=$stuff[3];
        $SCOS{$stuff[5]}=$stuff[7];
        $bothgood++;
        if($on%1000000==0){$time = localtime;
                print "on $on time $time total pairs $totreads unique pairs $bothgood duplicates $dup\n";
        }
        $on++;
}

$fkc = keys %SEQS;
foreach my $seq1 (keys %SEQS){
        if($seq1 !~ /^[A-Z]+$/){$bll++; next;}
        $rkc += keys %{$SEQ{$seq1}};
        foreach my $seq2 (keys %{$SEQS{$seq1}}){
                if($seq2 !~ /^[A-Z]+$/){$blr++; next;}
                $totout++;
                $nm1=$SEQS{$seq1}{$seq2}."_1";
                $nm2=$SEQS{$seq1}{$seq2}."_2";
                print OUTPUT1 "$nm1\n$seq1\n+\n$SCOS{$seq1}\n";
                print OUTPUT2 "$nm2\n$seq2\n+\n$SCOS{$seq2}\n";
        }
}
$end = localtime;

print "$input1 on $on start $start end $end total pairs $totreads unique pairs $bothgood duplicates $dup fkc $fkc rkc $rkc totout $totout bll $bll blr $blr\n";
print LOG "$input1 on $on start $start end $end total pairs $totreads unique pairs $bothgood duplicates $dup fkc $fkc rkc $rkc totout $totout bll $bll blr $blr\n";

