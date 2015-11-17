## removes any sequences coming after stop codon.

use warnings;
use Bio::SeqIO;
use Bio::Seq; 

my $seq= <STDIN>;
my $seqobjnuc = ();
$seqobjnuc = Bio::PrimarySeq->new (-seq => $seq );

my $seqobjprot = $seqobjnuc -> translate;

if ($seqobjprot->seq=~/\D*/){
  my (@triplet)=();
  my @prot=split//,$seqobjprot->seq;

  #split nuc seq in codon (triplet), and omit stop codon
  my $tmpseq=$seqobjnuc->seq;
  while($tmpseq=~/(\w{3})/g){
    my $codon=$1;
    next if $codon=~/TAG/i || $codon=~/TGA/i || $codon=~/TAA/i;
    push @triplet,$codon;
    }

  #update Bio::Seq object 
  $seqobjnuc->seq(join('',@triplet));
  $seqobjprot->seq($seqobjnuc->translate->seq);
  }
$seqout = $seqobjnuc->seq();
print "$seqout\n";