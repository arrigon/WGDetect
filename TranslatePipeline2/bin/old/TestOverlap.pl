#!/bin/perl

## Test single comparison
$test = ovlp(5, 12, 11, 15);
print("$test\n");


## Test multiple comparisons
my %coords;
$coords{"A"}{"start"} = 0;
$coords{"B"}{"start"} = 5;
$coords{"C"}{"start"} = 11;
$coords{"D"}{"start"} = 10;
$coords{"A"}{"stop"} = 4;
$coords{"B"}{"stop"} = 11;
$coords{"C"}{"stop"} = 19;
$coords{"D"}{"stop"} = 19;

## Test of overlap among proteins
my $ref = \%coords;
my %RES = %{ovlppw($ref)};

## Cleaning overlapping proteins (keeps longest ones)
foreach $pws (keys(%RES)){
  $rmve = $RES{$pws}{"rmve"};
  if($coords{$rmve}){
    delete($coords{$rmve});
    }
  }

print "DONE";

## sub for testing CDS overlaps. 
# Input: hash containing start / end coordinates of all proteins found in focal sequence
#
# produces ref to hash
# where keys = pws comparison
# and contains keep / rmve status of each tested protein
#
###
sub ovlppw {
  my %coords = %{$_[0]};
  my @prots = keys(%coords);
  my %RES;
  my $pws = 0;
  foreach my $idx1 (0..$#prots){
    foreach my $idx2 (0..$#prots){
      if($idx2 > $idx1){
	my $prt1 = @prots[$idx1];
	my $prt2 = @prots[$idx2];

	my $s1 = $coords{$prt1}{"start"};
	my $e1 = $coords{$prt1}{"stop"};
	my $s2 = $coords{$prt2}{"start"};
	my $e2 = $coords{$prt2}{"stop"};
	
	my $test = ovlp($s1, $e1, $s2, $e2);
	print "test $prt1 vs $prt2: ovlp = $test\n";

	if($test == 1){
	  if(abs($e1 - $s1) >= abs($e2 - $s2)){
	    $RES{$pws}{"keep"} = $prt1;
	    $RES{$pws}{"rmve"} = $prt2;
	    } else {
	    $RES{$pws}{"keep"} = $prt2;
	    $RES{$pws}{"rmve"} = $prt1;
	    }
	  $pws++;
	  }
	}
      }
    }
  my $adrs = \%RES;
  return $adrs;
  }

# overlap test give start / end positions of two proteins
sub ovlp {
  $s1 = $_[0];
  $e1 = $_[1]; 
  $s2 = $_[2]; 
  $e2 = $_[3];
  if($s2 < $e1 && $s1 < $e2){
    return 1;
    } else {
    return 0;
    }
  }