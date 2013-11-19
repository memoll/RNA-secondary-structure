#!/usr/bin/perl
use strict; use warnings;
use fasta_hash;

#############################################################################
## Read in UTR sequences, calculate bp frequences for each, then output a  ##
## randomly generated seqence per UTR in FASTA format.                     ##
## To run:                                                                 ##
##        ./randseq.pl < orig.fasta > new.fasta                            ##
#############################################################################

#subroutine: find bp frequencies
sub acugProbs{
    my $seq = $_[0];
    my $n = length($seq);
    my @seq = split //,$seq; ## make array for easy looping
    my %freqs = ();
    foreach my $c (@seq) {
        $freqs{$c}++;
    }
    ## transform to frequencies
    while ( my ($key, $value) = each(%freqs) ) {
        $freqs{$key} = $value/$n;
    }
    # return as reference.
    return \%freqs;
}


#####################################################################
### Subroutine for generating random sequence. Pass in reference  ###
### to nucleotide frequencies, and from these generating a random ###
### sequence of a specified length.                               ###
#####################################################################
sub genSeq{
    my ($freqs_ref, $n_ref) =  @_;
    ## dereference freqency hash
    my %freqs = %$freqs_ref;
    my $randseq = '';
    my ($prob, $tmp);
    ## loop that generates a random base each time from a uniform
    ## distribution and appends to a string.
    for (my $i = 0; $i < ${$n_ref}; $i++){
        $prob = rand;
        $tmp = 0;
        foreach my $key (keys %freqs) {
            if ($prob > $tmp && $prob < $tmp + $freqs{$key}) {
                $randseq .= $key;
            }
            $tmp += $freqs{$key};
        }
    }
    ## return reference to sequence
    return (\$randseq);
}

my %sequence_data;
my %fasta;

# Open fasta file from system
my $fasta_file="data/Splice_to_start.txt";
my $fh;

open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
while (fasta_hash::read_fasta_sequence($fh, \%sequence_data, 0)) {
    $fasta{$sequence_data{header}} = $sequence_data{seq};
}

## loop through fasta, generate new seqeunces and print
my ($freqs, $n, $randseq);
for my $c (keys %fasta) {
    my $freqs = acugProbs($fasta{$c});
    my $n = length($fasta{$c});
    my $randseq = genSeq($freqs, \$n);
    print ">$c\n";
    print "${$randseq}\n\n";
}
