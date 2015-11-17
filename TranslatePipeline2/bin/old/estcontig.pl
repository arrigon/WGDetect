 #use warnings;
 #my $out = <STDIN>;
  #parse genewise output
 #   open (WISE, $out);

     my @seqs_nuc=();
    
    LOOP1 : while(<STDIN>){
	
	my $seq='';
	
	if (/^>\S+\[([-\d]+)\:([-\d]+)\]\.sp$/){ #header for the first CDS chunk on either the forward or reverse strands. The terminal 'sp' indicates that it is a CDS but not the protein 
	    
	    # records coordinates
	    my $start=$1<$2?$1:$2;
	    my $end=$1<$2?$2:$1;
	 
	    
	    while(<STDIN>){
		if (/\/\//){ # // means end of annotation for the strand: put the data in a hash, store the hash in an array  and return to the LOOP1 flag to parse the annotation of the reverse strand or terminate parsing
		    my %hash=('seq'=>$seq,
			      'start'=>$start,
			      'end'=>$end
			      );
		    push @seqs_nuc,\%hash;
		    next LOOP1;
		}
		elsif(/^>\S+\[(\d+)\:(\d+)\]\.sp$/){ #header for a new chunk of CDS
		    $end=$2>$1?$2:$1; #update end coordinate
		    next;
		}
		else { # we are in the sequence: records
		    chomp;
		   s/\s//g;
		    $seq.=$_;
		 
		}
	    }
	}
    }
    
  #  close WISE;
    
        # select the longest CDS among the 2 annotations (fow. vs. rev.)
    
    if (length($seqs_nuc[0]{seq}) > length($seqs_nuc[1]{seq})){
	$seq = $seqs_nuc[0]{seq};
	$start = $seqs_nuc[0]{start};
	$end = $seqs_nuc[0]{end};
	print $seq;
    }else{
	$seq = $seqs_nuc[1]{seq};
	$start = $seqs_nuc[1]{start};
	$end = $seqs_nuc[1]{end};
	print $seq;
    }
  
 
