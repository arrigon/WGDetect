# singleLinkageClustering_KMD.pl
# KM Dlugosch ~ Aug 2010
#
# Assumes an input file which is a list of tab-separated pairs of numerical IDs, with or without 'x's on either side of the IDs (will be removed in output)
#
# run: perl singleLinkageClustering.pl <file_of_tab_separated_hit_pairs>

my $blasthits = $ARGV[0];
my $INIFAS = $ARGV[1];
my $outfolder = $ARGV[2];

###### Convert accs into numbers
################################
open (INFILE, $blasthits) or die "Couldn't open file $blasthits: $!\n";
open (OUT, ">$outfolder/.accs.converted.txt");
print "Converting accs into numbers...\r";

## Hash converting accs to numbers
my %accstokey = ();
my $numkey = 0;

## Hash keeping track of relative orientations
my %pairstoorient = ();


## Parse filtered blast report (assumes three columns: ref query orient)
while (my $line = <INFILE>) {
  chomp $line;
  my @array = split (/\t/, $line);

  # convert accs to keys
  if(!defined($accstokey{$array[0]})){
    $accstokey{$array[0]} = $numkey;
    $numkey++;
    }
  if(!defined($accstokey{$array[1]})){
    $accstokey{$array[1]} = $numkey;
    $numkey++;
    }
  print OUT "$accstokey{$array[0]}\t$accstokey{$array[1]}\n";

  # feed %orient
  $orient{$accs[0]}{$accs[1]} = $accs[2];

  }
close(OUT);
print "Converting accs into numbers... done\n";


###### Produce gene families (script KMD)
#########################################
print "Producing gene families...\n";
my (@families, @sizes);
my $tmp = "$outfolder/.accs.converted.txt";
unless (-e "$tmp") { die "\n\tError: File <$tmp> not found.\n\n"; }
&slClustering($tmp);
print "Producing gene families... done\n";

## Print clusters
open FAMS, ">$outfolder/$blasthits.clusters"; 
print FAMS "Cluster\tSize\tMembers\n";
for (my $i = 0; $i < @families; $i++) { print FAMS "$i\t$sizes[$i]\t$families[$i]\n"; }
close FAMS;
unlink("$outfolder/.accs.converted.txt");

###### Reconverting numbers into accessions
print "Converting numbers into accs...\r";
open (OUT, ">$outfolder/$blasthits.families");

## revert the accstokey hash,
%keytoaccs = reverse %accstokey;
for (my $i = 0; $i < @families; $i++) {
  $tmp = $families[$i];
  @members = split / /, $tmp;

  foreach $dude (@members){
    print OUT "$i\t$keytoaccs{$dude}\n";
    }
  }
close(OUT);
print "Converting numbers into accs... done\n";


###### Produce corresponding fasta files
########################################

## open FASTA file containing all sequences corresponding to blast hits
print "Producing gene family fasta files...\n";

## import fasta file
open(FILE, "$INIFAS");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## load it into hash. Note that %fasta contains your sequences. Can be reused elsewhere.
my %fasta;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  
  # clean sequence header as blast would do
  my $header = $tmp[0];
  $header =~ /([\w|\_]*)\s*/;
  $header = $1;

  $fasta{$header} = $seq; 
  }

## produce Gene family FASTA files
chdir("$outfolder/families/");


# loop over all gene families
for (my $i = 0; $i < @families; $i++) {

  # get members of that family
  $tmp = $families[$i];
  @members = split / /, $tmp;
  
  # prepare output file
  my $outfile = "Fam.$i.fas";
  open OUT, ">$outfile";
  
  # print output
  foreach $dude (@members){
    $seqid = $keytoaccs{$dude};
    $seq = $fasta{$seqid};
    print OUT ">$seqid\n$seq\n";   
    }
  close(OUT);
  }
print "Producing gene family fasta files... done\n";





################################################
#################################### SUBROUTINES
sub slClustering {
	my $file = $_[0];
	my ($lowID, %groups, @pointers);
	my ($highVal, $highKey) = (0,0);

	open HITS, "<$file";
	while(<HITS>) { 
		chomp $_; 
		$_ =~ s/x//g;  
		#read the information from each pair into a hash %groups where the key is the lowest ID # of the pair
		my @tabs = split /\t/, $_;
		# pop @tabs; TODO this was for a case where I had an alignment score after the pair of IDs
		@tabs = sort {$a <=> $b} @tabs;
		if ($groups{$tabs[0]}) { $groups{$tabs[0]} .= " $tabs[1]"; } # if the key already exists, append the new member to the value
		else { $groups{$tabs[0]} = $tabs[1]; }
		if ($tabs[1] > $highVal) { $highVal = $tabs[1]; }
		if ($tabs[0] > $highKey) { $highKey = $tabs[0]; }
	}
	close HITS;

	#create an array of pointers indicating cluster membership for each ID and use it to unite linked groups, iterate until unchanged
	@pointers = (0 .. $highVal);
	my ($newPointerRef);
	my $changed = 'Y';
	while ($changed eq 'Y') {  #keep doing this until the clusters are unchanged
		$newPointerRef = regroup(\@pointers, \%groups, $highVal, $highKey);  #send to the subroutine 'regroup', which will return an array reference $newPointerRef
		my $mismatches = 0;
		for (my $j = 0; $j < @pointers; $j++) { if ($pointers[$j] != $$newPointerRef[$j]) { $mismatches++; } }
		@pointers = @$newPointerRef;
		if ($mismatches == 0) { $changed = 'N'; }
	}
	print "Done.\n";
	(%groups, @$newPointerRef) = ();

	#use the pointers to create a new hash of clusters
	$groups{0} = 0;
	for (my $i = 1; $i < @pointers; $i++) { 
		if ($pointers[$i] == 0) { $groups{$pointers[$i]} .= " $i"; }
		elsif ($groups{$pointers[$i]}) {$groups{$pointers[$i]} .= " $i"; } 
		else { $groups{$pointers[$i]} = $i; } 
	}
	(@pointers) = ();

	# get @sizes and @families with space separated families for each line
	(@families, @sizes) = ();
	foreach my $fam (keys %groups) {
		my @temp = split / /, $groups{$fam};
		if (scalar@temp > 1) {
			push @sizes, scalar(@temp);
			push @families, $groups{$fam};
# 			print "Family is called <$fam> and includes <$groups{$fam}>\n";
		}
	}
	%groups = ();
} # end sub slClustering



sub regroup {
	my ($pointerRef, $groupRef, $highVal, $highKey) = @_;
	my (@newPointers);
	print "Iterating...\n";

	#update array of pointers indicating cluster membership for each ID
	@newPointers = @$pointerRef;
	for (my $i = 0; $i <= $highKey; $i++) {
		if ($$groupRef{$i}) {
			my @entries = split / /, $$groupRef{$i};
			# find the lowest ID among the hash and all of their entries in the pointer list and use it as the cluster identifier
			my $lowID = $i;
			if ($newPointers[$i] < $lowID) { $lowID = $newPointers[$i]; }
			foreach my $ID (@entries) { if ($newPointers[$ID] < $lowID) { $lowID = $newPointers[$ID]; } }
			foreach my $ID (@entries) { $newPointers[$ID] = $lowID; }
			$newPointers[$i] = $lowID;
		}
	}
	return (\@newPointers);
} # end sub regroup

