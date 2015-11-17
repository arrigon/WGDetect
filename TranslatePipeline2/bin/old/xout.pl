## Translates sequence into amino acids or leaves them as is
## Usage perl dna xout.pl < STDIN
## Usage perl prot xout.pl < STDIN  
# STDOUT is requested sequence
use warnings;
use Bio::SeqIO;
use Bio::Seq; 

$type = $ARGV[0];
$header = $ARGV[1];
$out = $ARGV[2];

## Sequence converter: produce sequence string
my $seq= <STDIN>;
my $seqobjnuc = ();
$seqobjnuc = Bio::PrimarySeq->new (-seq => $seq);

my $seqobjprot = $seqobjnuc -> translate;

if ($seqobjprot->seq=~/X/){
  my (@triplet)=();
  my @prot=split//,$seqobjprot->seq;

  #split nuc seq in codon (triplet), and omit stop codon
  my $tmpseq=$seqobjnuc->seq;
  while($tmpseq=~/(\w{3})/g){
    my $codon=$1;
    next if $codon=~/N\w\w/i || $codon=~/\wN\w/i || $codon=~/\w\wN/i;
    push @triplet,$codon;
    }

  #update Bio::Seq object 
  $seqobjnuc->seq(join('',@triplet));
  $seqobjprot->seq($seqobjnuc->translate->seq);
  }

if($type eq "dna"){
  $outseq = $seqobjnuc->seq();
  } else {
  $outseq = $seqobjprot->seq();
  } 

open(OUT, ">>$out");
print OUT ">$header\n$outseq\n";
close(OUT)