#use warnings;

#INPUT GENE READ COUNTS
print "input read counts\n";
$inrc="MERGED_CONTIGS_COR_GENES.rpkm";
open(INRC,$inrc)||die;
while(<INRC>){
        if($_ !~ /\w/){next;}
        if($_ =~ /^\#/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_,-1);
        $rc=$stuff[4];
        $rpkm=$stuff[5];
        $totreads+=$rc;
        $stuff[0] =~ s/\s.*//;
        $GENE_RC{$stuff[0]}=$rc;
        $GENE_KM{$stuff[0]}=$rpkm;
}

#IN TAXONOMY
#TID     GENES   AVG_GLEN        AVG_%_GC        K 
print "INPUT TAXONOMY\n";
#$intax = 'TAXONOMY_DB_2020.txt';
$intax = "/geomicro/data2/tealfurn/MUSCATO/TAXONOMY_DB_2020_mod2021_withCOUNTS.txt";
open(INTAX, $intax)||die;
while(<INTAX>){
    	if($_ !~ /\w/){next;}
    	@stuff = split("\t", $_,-1);
    	$stuff[$#stuff]=~ s/[\n\r]+//;
    	$tid = shift(@stuff);
	$gc=shift(@stuff);
	$glen=shift(@stuff);
	$pgc=shift(@stuff);
	$lin=join(";",@stuff);
    	$PHY{$tid}=$lin;
	$LIN2TID{$lin}=$tid;
	if($gc > 0 && $pgc > 0){
		if(!exists($LIN2PGC{$lin}) || $gc > $LIN2TGC{$lin}){$LIN2PGC{$lin}=$pgc; $LIN2TGC{$lin}=$gc;}
	}
}
#have to separate from tax input because of duplicate taxon ids
foreach my $lin (keys %LIN2PGC){
	@LEVS=split(";",$lin);
        for my $i (0..$#LEVS){
         	$j = $#LEVS-$i;
                $nlin = join(";",@LEVS[0..$j]);
                $PHY_PGC{$nlin}+=$LIN2PGC{$lin};
		$PHY_CNT{$nlin}++;
		if($nlin eq "BACTERIA;SPIROCHAETES;SPIROCHAETIA;LEPTOSPIRALES;LEPTOSPIRACEAE;LEPTOSPIRA;LEPTOSPIRA_MEYERI"){ 
			$avg=$PHY_PGC{$nlin}/$PHY_CNT{$nlin};
			print "foundleptospira pgc $PHY_PGC{$nlin} avg $avg lin $lin\n"; 
		}
        }
}


#qseqid 	qlen 	sseqid 		slen 	qstart 	qend 	sstart 	send 	evalue 		pident 	mismatch qcovhsp scovhsp
#CLUSTER-0_0_1  2562    CLID:84788615   475     703     1947    60      474     4.5e-112        49.9    208      48.6    87.4
print "input diamond hits\n";
$ingvu = "SAMPLE_53600_GENESvsURDB.m8"; 
open(INGVU, $ingvu)||die;
while(<INGVU>){
	if($_ !~ /\w/){next;}
	$_=uc($_);
	$_=~s/[\r\n]+//;
	@stuff=split("\t",$_,-1);
	$gene = $stuff[0];
	$hit  = $stuff[2];
	$eval = $stuff[8];
	$pid  = $stuff[9];
	$cov  = $stuff[11];
	$sco  = $pid*$cov;
	$TG{$gene}=1;
	#remove very low quality hits - spot checked: ok
	if($sco < 1000 && $eval > 0.0000001){ next; }
	if(!exists($GENEMAXSCO{$gene})){$GENEMAXSCO{$gene}=$sco;}
	if($sco > $GENEMAXSCO{$gene}){ $GENEMAXSCO{$gene}=$sco; }
	if($sco >= $GENEMAXSCO{$gene}*0.5){ $GENETOPHIT{$gene}{$hit}=$sco; $HITS{$hit}{$gene}=$sco;}
}
$kc = keys %GENEMAXSCO;
$tc = keys %TG;
print "totcongen $tc good $kc\n";
undef(%TG);


print "input URDB info\n";
$outmis = "missingt_tids.txt";
open(OUTMIS, ">", $outmis)||die;
$ininfo = "ALL_GENES_INFO.txt";
open(INFO, $ininfo)||die;
while(<INFO>){
	if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        @stuff=split("\t",$_,-1);
	if(!exists($HITS{$stuff[0]})){next;}
	$hit = $stuff[0];
	$name = $stuff[1];
	$HITNAME{$hit}=$name;

	#taxa
	@tids = split(";",$stuff[3]);
	foreach my $tid (@tids){ #using lineage because can be duplicate lins with different tids
		if($PHY{$tid}!~/\w/){ print OUTMIS "$tid\n"; next;}
		$lin=$PHY{$tid};
		$HIT2LINS{$hit}{$lin}=1;
		foreach my $gene (keys %{$HITS{$hit}}){
			$LIN2GENES{$lin}{$gene}=$HITS{$hit}{$gene}; #total contig genes for tid
			$GENE2LINS{$gene}{$lin}=1; #for looping thru gene tids
		}
	}

	#functions
	for my $f (9..14,18){ $stuff[$f]=~s/\s+//g; }
	if($stuff[9] =~/\w/){@FID=split(";",$stuff[9] ); foreach my $id (@FID){ foreach my $gene (keys %{$HITS{$hit}}){ $GENEKEGG{$gene}{$id}++; $GENEFUNCS{$gene}{$id}=1; $FUNCGENES{$id}{$gene}=1; }}}
	if($stuff[11]=~/\w/){@FID=split(";",$stuff[10]); foreach my $id (@FID){ foreach my $gene (keys %{$HITS{$hit}}){ $GENECOG{$gene}{$id}++;  $GENEFUNCS{$gene}{$id}=1; $FUNCGENES{$id}{$gene}=1; }}}
	if($stuff[10]=~/\w/){@FID=split(";",$stuff[11]); foreach my $id (@FID){ foreach my $gene (keys %{$HITS{$hit}}){ $GENEPFAM{$gene}{$id}++; $GENEFUNCS{$gene}{$id}=1; $FUNCGENES{$id}{$gene}=1; }}}
	if($stuff[12]=~/\w/){@FID=split(";",$stuff[12]); foreach my $id (@FID){ foreach my $gene (keys %{$HITS{$hit}}){ $GENEGO{$gene}{$id}++;   $GENEFUNCS{$gene}{$id}=1; $FUNCGENES{$id}{$gene}=1; }}}
	if($stuff[13]=~/\w/){@FID=split(";",$stuff[13]); foreach my $id (@FID){ foreach my $gene (keys %{$HITS{$hit}}){ $GENEIPR{$gene}{$id}++;  $GENEFUNCS{$gene}{$id}=1; $FUNCGENES{$id}{$gene}=1; }}}
	if($stuff[14]=~/\w/){@FID=split(";",$stuff[14]); foreach my $id (@FID){ foreach my $gene (keys %{$HITS{$hit}}){ $GENEEC{$gene}{$id}++;   $GENEFUNCS{$gene}{$id}=1; $FUNCGENES{$id}{$gene}=1; }}}
	if($stuff[18]=~/\w/){@FID=split(";",$stuff[18]); foreach my $id (@FID){ foreach my $gene (keys %{$HITS{$hit}}){ $GENERXN{$gene}{$id}++;  $GENEFUNCS{$gene}{$id}=1; $FUNCGENES{$id}{$gene}=1; }}}
}


#get the sum score (%ID x coverage) of the found genes for each tid
print "summarize lineage gene scores\n";
foreach my $lin (keys %LIN2GENES){
	foreach my $gene (keys %{$LIN2GENES{$lin}}){
		$LINSUMSCO{$lin}+=$LIN2GENES{$lin}{$gene};
	}
}


print "input contig genes\n";
$/=">";
$incg = "MERGED_CONTIGS_COR_GENES.fna";
open(INCG, $incg)||die;
while(<INCG>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        @stuff = split("\n",$_,-1);
        $header = shift(@stuff);
        $header =~ /((CLUSTER\-\d+\_\d+)\_\d+)\D+(\d+)\D+(\d+)\D.*PARTIAL\=(\d+)/;
        $gene=$1;
        $contig=$2;
        $start=$3;
        $end=$4;
        $part=$5;
        $seq = join("",@stuff);
        $seq =~ s/[^A-Z]+//g;
	my $count = $seq =~ tr/GC/gc/;
	$len=length($seq);
	$prgc=$count*100/$len;
        $CONGENS{$contig}{$gene}=1;
        $CONGENPART{$gene}=$part;
        $CONGENLEN{$gene}=$len;
        $CONGENPOS{$gene}=$start."\t".$end;
	$CONGENPGC{$gene}=$prgc;
}


print "input contigs and refine gene tids and hits\n";
$incon = "MERGED_CONTIGS_COR.fasta";
$on=0;
open(INCON, $incon)||die;
while(<INCON>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        @stuff = split("\n",$_,-1);
        $contig = shift(@stuff);
	$TOTCONTIGS{$contig}=1;
	if(keys %{$CONGENS{$contig}} < 1){ next; } # add print output no genes
	$contig=~/(CLUSTER\-\d+)/;
        $cluster=$1;
	$seq = join("",@stuff);
        $seq =~ s/[^A-Z]+//g;
	$conlen = length($seq);
	$CONLENS{$contig}=$conlen;

	#FOR GENES WITH URDB HITS: ADD TO GLIST AND ADD SCORES OF GENETIDS TO TSCO
	%TSCO=(); %GLIST=();
	foreach my $gene (keys %{$CONGENS{$contig}}){
		if(exists($GENETOPHIT{$gene})){ $GLIST{$gene}=1;
			foreach my $lin (keys %{$GENE2LINS{$gene}}){ $TSCO{$lin}=$LINSUMSCO{$lin};}
		}
	}

	#FOR GENES WITH URDB HITS: REMOVE EXTRA TIDS
	#IF CONTIG HAS MULTIPLE GENES: TRY TO FIND MINIMUM TIDS THAT EXPLAIN GENES
	$maxsco=0; $minsco=0; %GCNT = %GLIST; 
	if(keys %GLIST > 0){
		#IF ONLY ONE GENE: GET TOP TID AND REMOVE TIDS < 90% OF TOP SCORE
		if(keys %GLIST == 1){ 
			foreach my $lin (sort{$TSCO{$b}<=>$TSCO{$a}} keys %TSCO){
                                if($TSCO{$lin} > $maxsco){$maxsco = $TSCO{$lin};}
                                if($TSCO{$lin} < $maxsco*0.9){delete($TSCO{$lin});}
                        }
                }
		else{ 	#REDUCE TID LIST TO MINIMUM NEEDED TO EXPLAIN GENES
			foreach my $lin (sort{$TSCO{$b}<=>$TSCO{$a}} keys %TSCO){
				if(keys %GCNT < 1){ 
					if($TSCO{$lin}>=$minsco){next;}
					else{delete($TSCO{$lin});} #remove lineages below minimum
				}
				foreach my $gene (keys %GCNT){ if(exists($GENE2LINS{$gene}{$lin})){ delete($GCNT{$gene}); $minsco=$TSCO{$lin};}}
				if($TSCO{$lin} > $maxsco){$maxsco=$TSCO{$lin};}
		    	}
		}
		
		#NOW LINS OF GENES NARROWED DOWN, NARROW DOWN GENE HITS AND LINS
		foreach my $gene (keys %GLIST){
			foreach my $hit (keys %{$GENETOPHIT{$gene}}){
				$good=0;
				foreach my $lin (keys %{$HIT2LINS{$hit}}){ 
					if(exists($TSCO{$lin})){$good++;}
					else{delete($GENE2LINS{$gene}{$lin});}
				}
				if($good==0){delete($GENETOPHIT{$gene}{$hit});}
			}
			$CLUST2GENES{$cluster}{$gene}=1;
		}
	}
	$on++;
}


#NARROW DOWN AND GET CLUSTER LCAS
print "refine cluster lineages\n";
foreach my $cluster (sort(keys %CLUST2GENES)){
	$tgc = keys %{$CLUST2GENES{$cluster}};
	%LINS=(); @TOTLINS=(); %COUNTS=(); %CLEVS=();
	foreach my $gene (keys %{$CLUST2GENES{$cluster}}){
		foreach my $lin (keys %{$GENE2LINS{$gene}}){
			$LINS{$lin}++; push(@TOTLINS,$lin);
			@LEVS=split(";",$lin);
			for my $i (0..$#LEVS){
                       		$j = $#LEVS-$i;
				$nlin = join(";",@LEVS[0..$j]);
				$COUNTS{$nlin}{$gene}=1;
				$CLEVS{$nlin}=$j;
				if($tgc>=5){ $kc=keys %{$COUNTS{$nlin}}; print "cluster $cluster gene $gene i $i j $j kc $kc tgc $tgc nlin $nlin\n";}
			}
		}
	}

	#find lins with at least half the genes
	@GOODLEVS=(); $bl=0;
	foreach my $nlin (sort{$CLEVS{$b}<=>$CLEVS{$a}} keys %CLEVS){
		$kc=keys %{$COUNTS{$nlin}};
		if($bl>$CLEVS{$nlin}){last;}
		if($kc > $tgc/2){push(@GOODLEVS,$nlin); $bl=$CLEVS{$nlin};}
	}
	if($GOODLEVS[0]!~/\w/){ $clusterlca = MakeLCA(@TOTLINS); }
	else{ $clusterlca = MakeLCA(@GOODLEVS); }
	$CLUSTLCA{$cluster}=$clusterlca;
	if($tgc>=5){$tl=@TOTLINS; $gl=@GOODLEVS; print "cluster $cluster tl $tl @TOTLINS\ncluster $cluster gl $gl @GOODLEVS\ncluster $cluster lca $clusterlca\n";}
}
undef(%CLUST2GENES);


print "output gene summary info\n";
$outgene="CONTIG_GENE_INFO.txt";
open(OUTGENE,">",$outgene)||die;
print OUTGENE "gene\tcontig\tstart\tend\tpartial\tgene_pgc\tlev_avg_pgc\tgenlen\tconlen\trpm\trpkm\tname\ttophit\tsco\tkid\tcid\tpid\tgid\tiid\teid\trid\tgtids\tctids\tgenelca\tconlca\tlusterlca\n";
foreach my $contig (keys %TOTCONTIGS){

	#GET CLUSTER LCA
	$contig=~/(CLUSTER\-\d+)/;
        $cluster=$1;
	$clusterlca="NOVEL";	
	if($CLUSTLCA{$cluster}){$clusterlca=$CLUSTLCA{$cluster};}
        $clusterlca =~ /^(.)/; $king = $1;
        @TAXON=split(";",$clusterlca);
        for($i=0; $i<=$#TAXON; $i++){
                $type = $king."_".$i;
                if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
                else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
                $TCYTO{$phyla}=1;
                $CYTO_CLUSTCONS{$phyla}+=keys %{$CONGENS{$contig}};
        }

	#GET CONTIG LCA BY CLEAN GENE LINS
        %LINS=(); @TOTLINS=(); %COUNTS=(); %CLEVS=();
	foreach my $gene (keys %{$CONGENS{$contig}}){
		if(!exists($GENE2LINS{$gene})){next;}
		foreach my $lin (keys %{$GENE2LINS{$gene}}){
                        $LINS{$lin}++; push(@TOTLINS,$lin);
                        @LEVS=split(";",$lin);
                        for my $i (0..$#LEVS){
                                $j = $#LEVS-$i;
                                $nlin = join(";",@LEVS[0..$j]);
                                $COUNTS{$nlin}{$gene}=1;
                                $CLEVS{$nlin}=$j;
                        }
                }
        }
        #find lins with at least half the genes
        @GOODLEVS=(); $bl=0;
        foreach my $nlin (sort{$CLEVS{$b}<=>$CLEVS{$a}} keys %CLEVS){
                $kc=keys %{$COUNTS{$nlin}};
                if($bl>$CLEVS{$nlin}){last;}
                if($kc > $tgc/2){push(@GOODLEVS,$nlin); $bl=$CLEVS{$nlin};}
        }
	  if($GOODLEVS[0] =~ /\w/){$conlca = MakeLCA(@GOODLEVS); }
	elsif($TOTLINS[0] =~ /\w/){$conlca = MakeLCA(@TOTLINS); }
	else{$conlca="NOVEL";}
	$conlen = $CONLENS{$contig};
        $conlca =~ /^(.)/; $king = $1;
        @TAXON=split(";",$conlca);
        for($i=0; $i<=$#TAXON; $i++){
                $type = $king."_".$i;
                if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
                else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
                $TCYTO{$phyla}=1;
                $CYTO_CCON{$phyla}++;
                $CYTO_CGEN{$phyla}+=keys %{$CONGENS{$contig}};
                $CYTO_CLEN{$phyla}+=$conlen;
        }
	@CTIDS=(); foreach my $lin (@TOTLINS){ push(@CTIDS,$LIN2TID{$lin}); }
	$ctids=join(";",@CTIDS);

	#PRINT AND NEXT IF NO GENES ON CONTIG
	if(!exists($CONGENS{$contig})){print OUTGENE "0\t$contig\t0\t0\t0\t0\t0\t$conlen\t0\t0\t\t\t\t\t\t\t\t\t\t\t\t$ctids\t\t$conlca\t$clusterlca\n"; next;}


	#OUTPUT GENE INFO
	foreach my $gene (keys %{$CONGENS{$contig}}){
		$genelca="NOVEL"; $gtids='';
		$kid=''; $cid=''; $pid=''; $gid=''; $iid=''; $eid=''; $rid=''; 
		$name="NOVEL_GENE";
		$sco=0;
		$tophit='';
        	$partial = $CONGENPART{$gene};
        	$glen 	 = $CONGENLEN{$gene};
        	$pos 	 = $CONGENPOS{$gene};
		$rpm 	 = $GENE_RC{$gene}*1000000/$totreads;
		$rpkm 	 = $GENE_KM{$gene};		
		$prgc	 = $CONGENPGC{$gene};
                $avgpgc  = 0;


		if(exists($GENEMAXSCO{$gene})){
			#MAKE GENE LCA
		 	@GPHYL=(); @GTIDS=();
			foreach my $lin (keys %{$GENE2LINS{$gene}}){
				push(@GPHYL,$lin); push(@GTIDS, $LIN2TID{$lin});
				#GET TID SCORE FOR ALL GENE FUNCTIONS LATER
				foreach my $id (keys %{$GENEFUNCS{$gene}}){ $FUNCLINSCO{$id}{$lin}+=$GENEMAXSCO{$gene}; $FUNCLINRPKM{$id}{$lin}+=$rpkm; }
			}
			$genelca = MakeLCA(@GPHYL);
			$gtids=join(";",@GTIDS);
			
			#get gene funcs
				foreach my $id (sort {$GENEKEGG{$gene}{$b}<=>$GENEKEGG{$gene}{$a}} keys %{$GENEKEGG{$gene}}){ $kid=$id; last;}
                		foreach my $id (sort {$GENECOG{$gene}{$b} <=>$GENECOG{$gene}{$a}}  keys %{$GENECOG{$gene}}){  $cid=$id; last;}
                	$c=0; 	foreach my $id (sort {$GENEPFAM{$gene}{$b}<=>$GENEPFAM{$gene}{$a}} keys %{$GENEPFAM{$gene}}){ if($c<5){$c++; $pid.=$id.";";} last;} $pid=~s/\;+$//;
                	$c=0; 	foreach my $id (sort {$GENEGO{$gene}{$b}  <=>$GENEGO{$gene}{$a}}   keys %{$GENEGO{$gene}}){   if($c<5){$c++; $gid.=$id.";";} last;} $gid=~s/\;+$//;
                	$c=0; 	foreach my $id (sort {$GENEIPR{$gene}{$b} <=>$GENEIPR{$gene}{$a}}  keys %{$GENEIPR{$gene}}){  if($c<5){$c++; $iid.=$id.";";} last;} $iid=~s/\;+$//;
                	$c=0; 	foreach my $id (sort {$GENEEC{$gene}{$b}  <=>$GENEEC{$gene}{$a}}   keys %{$GENEEC{$gene}}){   if($c<5){$c++; $eid.=$id.";";} last;} $eid=~s/\;+$//;
                	$c=0; 	foreach my $id (sort {$GENERXN{$gene}{$b} <=>$GENERXN{$gene}{$a}}  keys %{$GENERXN{$gene}}){  if($c<5){$c++; $rid.=$id.";";} last;} $rid=~s/\;+$//;

			#get top gene hit, name, and score
			foreach my $hit (sort{$GENETOPHIT{$gene}{$b}<=>$GENETOPHIT{$gene}{$a}} keys %{$GENETOPHIT{$gene}}){ #loop through gene hits 
				foreach my $lin (sort {$LINSUMSCO{$b} <=> $LINSUMSCO{$a}} keys %{$GENE2LINS{$gene}}){ #loop through good contig tids 
					if(exists($HIT2LINS{$hit}{$lin})){ $tophit=$hit; last; } #if that gene top hit has a top tid - set it as top and end
				}
				if($tophit=~/\w/){last;}
			}
			$name=$HITNAME{$tophit};
			$sco=$GENEMAXSCO{$gene};
			if($PHY_CNT{$genelca} >= 1){$avgpgc  = $PHY_PGC{$genelca}/$PHY_CNT{$genelca};}
		}
		
		#PRINT OUT GENE INFO
		print OUTGENE "$gene\t$contig\t$pos\t$partial\t$prgc\t$avgpgc\t$glen\t$conlen\t$rpm\t$rpkm\t$name\t$tophit\t$sco\t$kid\t$cid\t$pid\t$gid\t$iid\t$eid\t$rid\t$gtids\t$ctids\t$genelca\t$conlca\t$clusterlca\n";

		#GET CYTO
       		$genelca =~ /^(.)/; $king = $1;
		@TAXON=split(";",$genelca);
        	for($i=0; $i<=$#TAXON; $i++){
                	$type = $king."_".$i;
                	if($i == 0){$phyla = "$type\tROOT\t$TAXON[$i]";}
                	else{$phyla = "$type\t$TAXON[$i-1]\t$TAXON[$i]";}
			$TCYTO{$phyla}=1;
			$CYTO_GLCA{$phyla}++;
			$CYTO_GLEN{$phyla}+=$glen;
			$CYTO_GSCO{$phyla}+=$sco;
			$CYTO_GRPM{$phyla}+=$rpm;
			$CYTO_GRKM{$phyla}+=$rpkm;
		}
	}
}

$/="\n";
$infn = "Function_Names.txt";
open(INFN, $infn)||die;
while(<INFN>){
        if($_ !~ /\w/){next;}
        $_=uc($_);
        $_=~s/[\r\n]+//;
        (my $id, my $name)=split("\t",$_,-1);
	$id=~s/\s//g;
	$id=~s/^EC\://;
	$FUNC2NAME{$id}=$name;
}



print "output functions\n";
$outfuncs = "SAMPLE_53600_funcs.txt";
open(FUNC,">",$outfuncs)||die;
print FUNC "FID\tName\tGC\tLINS\tsumsco\tsumrpkm\tLCA\n";
foreach my $id (sort(keys %FUNCGENES)){
	$gc = keys %{$FUNCGENES{$id}};
	$tc = keys %{$FUNCLINSCO{$id}};
	@LINS=();
	$sumsco=0; $sumrpkm=0;
	foreach my $lin (keys %{$FUNCLINSCO{$id}}){
		$sumsco+=$FUNCLINSCO{$id}{$lin};
		$sumrpkm+=$FUNCLINRPKM{$id}{$lin};
		push(@LINS,$lin);
	}
	$funclca = MakeLCA(@LINS);
	print FUNC "$id\t$FUNC2NAME{$id}\t$gc\t$tc\t$sumsco\t$sumrpkm\t$funclca\n";
}

print "output cytoscape file\n";
$outcyto = "SAMPLE_53600.cyto";
open(OUTCYTO, ">", $outcyto)||die;
print OUTCYTO "type\tsource\ttarget\tgsumgc\tgsumlen\tgsumsco\tgsumrpm\tgsumrpkm\tcsumcon\tcsumtgc\tcsumlen\tclusterlca\n";
foreach my $phyla (sort(keys %TCYTO)){
	$gsumgc=$CYTO_GLCA{$phyla};
	$gsumlen=$CYTO_GLEN{$phyla};
	$gsumsco=$CYTO_GSCO{$phyla};
	$gsumrpm=$CYTO_GRPM{$phyla};
	$gsumrpkm=$CYTO_GRKM{$phyla};
	$csumcon=$CYTO_CCON{$phyla};
	$csumtgc=$CYTO_CGEN{$phyla};
	$csumlen=$CYTO_CLEN{$phyla};
	$csumclu=$CYTO_CLUSTCONS{$phyla};
	print OUTCYTO "$phyla\t$gsumgc\t$gsumlen\t$gsumsco\t$gsumrpm\t$gsumrpkm\t$csumcon\t$csumtgc\t$csumlen\t$csumclu\n";
}




sub MakeLCA{
        my @TAXON;
        %seen=();
        @PHYL = @_;
        @PHYL = grep { !$seen{$_}++ } @PHYL;
	if(grep {/^MONA/} @PHYL){
		@TMP=(); 
		foreach my $lin (@PHYL){ 
			if($lin =~ /^MONA/){ push(@TMP,$lin); } 
		}
		 @PHYL=@TMP;
	}
        $len1 = @PHYL;
        if($len1 == 1){$LCA = $PHYL[0]; }
        elsif($len1 >1){
                $first = $PHYL[0];
                @levels = split(";", $first);
                for($i=0; $i<=$#levels; $i++){
                        $alevel=$levels[$i];
                        @matched = grep(/\Q$alevel\E/i, @PHYL);
                        $len2 = @matched;
                        if($len2 == $len1){push(@TAXON, $alevel);}
                        else{last;}
                }
                $len3 = @TAXON;
                if($len3 > 1){$LCA = join(";", @TAXON);}
                elsif($len3==1){$LCA = $TAXON[0];}
                else{$LCA = "NCA"; }
        }
        else{$LCA = "NCA"; }
        return($LCA);
}	

