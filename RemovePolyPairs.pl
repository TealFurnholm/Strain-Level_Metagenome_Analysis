use warnings;

$input1 = $ARGV[0]; 
$input2 = $ARGV[1];
$minlen = $ARGV[2];
$output1 = $ARGV[3];
$output2 = $ARGV[4];

if($input1 !~ /(\_1\_|fwd)/i){ print "Please specify forward read file with either _1_ or _fwd_ in the file name\n"; }
if($input2 !~ /(\_2\_|rev)/i){ print "Please specify reverse read file with either _2_ or _rev_ in the file name\n"; }
if($output1 !~ /\w/ || $output2 !~ /\w/ || $output1 eq $output2){ print "Please specify forward and reverse output file names\n"; }

if($input1 =~ /\.gz$/i){
	open(INPUT1, "gunzip -c $input1 |") or die "gunzip $input1: $!";
        open(INPUT2, "gunzip -c $input2 |") or die "gunzip $input2: $!";
}
else{	open(INPUT1, $input1)||die; open(INPUT2, $input2)||die; }

if($output1 =~ /\.gz$/i){
	open(OUTPUT1,'>>:gzip', $output1)||die;
	open(OUTPUT2,'>>:gzip', $output2)||die;
}
else{	open(OUTPUT1, ">>", $output1)||die; open(OUTPUT2, ">>", $output2)||die; }


$on=0;
$badline=0;
$start = localtime;
while( my $line = <INPUT1> . <INPUT1> . <INPUT1> . <INPUT1> . <INPUT2> . <INPUT2> . <INPUT2> . <INPUT2> ) {

	$totreads++;
   	$line = uc($line);
	@stuff = split("\n", $line);
	if($stuff[1] !~ /^[A-Z]+$/ || $stuff[5] !~ /^[A-Z]+$/){print "non-nucleotide detected on $on at $stuff[1] \n or \n $stuff[5]\n"; $badline++; next;}
	$stuff[4]=~/^(\S+)\s/; $read=$1;
	if($stuff[0] !~ /$read/){ print "paired fastq files not matched pairs $stuff[0] ne $stuff[4]\n"; $badline++; next;}
	if($badline > 3){print "File parsing issue, make sure all lines in file match fastq format and pairs are matched\n"; die;}


	#REMOVE AMBIGS
	$stuff[1] =~ s/[^ATGC]/N/g;
	if($stuff[1] =~ /^N|N$/){
		@seq1 = split("", $stuff[1]);
		@sco1 = split("", $stuff[3]);
		while($seq1[0] eq "N"){ shift(@seq1); shift(@sco1); }
		@seq1 = reverse(@seq1);
		while($seq1[0] eq "N"){ shift(@seq1); shift(@sco1); }
		@seq1 = reverse(@seq1);
		$stuff[1] = join("", @seq1);
		$stuff[3] = join("", @sco1);
	}
        $stuff[5] =~ s/[^ATGC]/N/g;
        if($stuff[5] =~ /^N|N$/){
                @seq2 = split("", $stuff[5]);
                @sco2 = split("", $stuff[7]);
                while($seq2[0] eq "N"){ shift(@seq2); shift(@sco2); }
                @seq2 = reverse(@seq2);
                while($seq2[0] eq "N"){ shift(@seq2); shift(@sco2); }
                @seq2 = reverse(@seq2);
                $stuff[5]=join("", @seq2);
		$stuff[7] = join("", @sco2);
        }


	#REMOVE LOW ENTROPY READS
	%DI = ();
	@SEQ = split("", $stuff[1]);
	for my $x (0..$#SEQ-1){
	        $y = $x+1;
	        $di = join("", @SEQ[$x..$y]);
	        $DI{$di}=1;
	}
	$kc1 = keys %DI;
        %DI2 = ();
        @SEQ = split("", $stuff[5]);
        for my $x (0..$#SEQ-1){
                $y = $x+1;
                $di = join("", @SEQ[$x..$y]);
                $DI2{$di}=1;
        }
        $kc2 = keys %DI2;

	
	#REMOVE POLY AND SHORT REPEATS
	$seq = $stuff[1];
	$sco = $stuff[3];
	while($seq =~ /^.?A{20,}|T{20,}|G{20,}|C{20,}|(.{2,6})\1{5,}/){
        if($seq =~ /^(A{20,}|T{20,}|G{20,}|C{20,})/){ $beg = @+[0]; }
        elsif( $seq =~ /^(.{2,6})\1{5,}/){ $beg = @+[0]; }
        elsif($seq =~ /^.?(A{20,}|T{20,}|G{20,}|C{20,})/){ $beg = @+[0]; }
        elsif( $seq =~ /^.?(.{2,6})\1{5,}/){ $beg = @+[0]; }
        else{last;}
        $beg = @+[0];
        @seq = split("", $seq);
        @seq = @seq[$beg..$#seq];
        @sco = split("", $sco);
        @sco = @sco[$beg..$#sco];
        $seq = join("", @seq);
        $sco = join("", @sco);
	}

	while($seq =~ /A{20,}|T{20,}|G{20,}|C{20,}|(.{2,6})\1{5,}.?$/){
        if($seq =~ /(A{20,}|T{20,}|G{20,}|C{20,})$/){ $beg = @+[0]; }
        elsif( $seq =~ /(.{2,6})\1{5,}$/){ $beg = @+[0]; }
        elsif($seq =~ /(A{20,}|T{20,}|G{20,}|C{20,}).?$/){ $beg = @+[0]; }
        elsif( $seq =~ /(.{2,6})\1{5,}.?$/){ $beg = @+[0]; }
        else{last;}
        $end = @-[0];
        $end--;
        @seq = split("", $seq);
        @seq = @seq[0..$end];
        @sco = split("", $sco);
        @sco = @sco[0..$end];
        $seq = join("", @seq);
        $sco = join("", @sco);
	}
	$stuff[1] = $seq;
	$stuff[3] = $sco;


	#REMOVE POLY REVERSE
	$seq = $stuff[5];
	$sco = $stuff[7];
        while($seq =~ /^.?A{20,}|T{20,}|G{20,}|C{20,}|(.{2,6})\1{5,}/){
        if($seq =~ /^(A{20,}|T{20,}|G{20,}|C{20,})/){ $beg = @+[0]; }
        elsif( $seq =~ /^(.{2,6})\1{5,}/){ $beg = @+[0]; }
        elsif($seq =~ /^.?(A{20,}|T{20,}|G{20,}|C{20,})/){ $beg = @+[0]; }
        elsif( $seq =~ /^.?(.{2,6})\1{5,}/){ $beg = @+[0]; }
        else{last;}
        $beg = @+[0];
        @seq = split("", $seq);
        @seq = @seq[$beg..$#seq];
        @sco = split("", $sco);
        @sco = @sco[$beg..$#sco];
        $seq = join("", @seq);
        $sco = join("", @sco);
        }


        while($seq =~ /A{20,}|T{20,}|G{20,}|C{20,}|(.{2,6})\1{5,}.?$/){
        if($seq =~ /(A{20,}|T{20,}|G{20,}|C{20,})$/){ $beg = @+[0]; }
        elsif( $seq =~ /(.{2,6})\1{5,}$/){ $beg = @+[0]; }
        elsif($seq =~ /(A{20,}|T{20,}|G{20,}|C{20,}).?$/){ $beg = @+[0]; }
        elsif( $seq =~ /(.{2,6})\1{5,}.?$/){ $beg = @+[0]; }
        else{last;}
        $end = @-[0];
        $end--;
        @seq = split("", $seq);
        @seq = @seq[0..$end];
        @sco = split("", $sco);
        @sco = @sco[0..$end];
        $seq = join("", @seq);
        $sco = join("", @sco);
        }
        $stuff[5] = $seq;
        $stuff[7] = $sco;


	#OUTPUT READS
	my $NFs = $stuff[1] =~ tr/Nn/NN/; 
	my $NRs = $stuff[5] =~ tr/Nn/NN/;
	   if($kc1 < 7 && $kc2 < 7){													$lowent++; next;} #both low entropy
	elsif($NFs > length($stuff[1])*0.05 && $NRs > length($stuff[5])*0.05){			$ambig++; next;} #both too many ambigs
	elsif(length($stuff[1]) < $minlen && length($stuff[5]) < $minlen){				$short++; next;} #both too short
	elsif($kc1 < 7 || $NFs > length($stuff[1])*0.05 || length($stuff[1]) < $minlen){#forward read bad
	   if($kc2 < 7 || $NRs > length($stuff[1])*0.05 || length($stuff[5]) < $minlen){$mixed++; next;} #both bad
	   $revgood++;
	}
	elsif($kc2 < 7 || $NRs > length($stuff[1])*0.05 || length($stuff[5]) < $minlen){#reverse read bad
	   if($kc1 < 7 || $NFs > length($stuff[1])*0.05 || length($stuff[1]) < $minlen){$mixed++; next;} #both bad
	   $forgood++;
	}
	else{$bothgood++;}

	$stuff[0]=~s/\s.*//;
	$SEQS{$stuff[1]}{$stuff[5]}=$stuff[0];
	$SCOS{$stuff[1]}=$stuff[3];
	$SCOS{$stuff[5]}=$stuff[7];

	if($on%1000000==0){$time = localtime; print "on $on time $time tot $totreads bothgood $bothgood forgood $forgood revgood $revgood bothlowent $lowent bothambig $ambig bothshort $short mixed $mixed\n";} 
	$on++;
}
$end = localtime;
print "on $on start $start end $end minlen $minlen\n";
print "tot $totreads bothgood $bothgood forgood $forgood revgood $revgood bothlowent $lowent bothambig $ambig bothshort $short mixed $mixed\n";


foreach my $seq1 (keys %SEQS){
	foreach my $seq2 (keys %{$SEQS{$seq1}}){
		$nm1=$SEQS{$seq1}{$seq2}."_1";
		$nm2=$SEQS{$seq1}{$seq2}."_2";
		print OUTPUT1 "$nm1\n$seq1\n+\n$SCOS{$seq1}\n";
		print OUTPUT2 "$nm2\n$seq2\n+\n$SCOS{$seq2}\n";
	}
}
