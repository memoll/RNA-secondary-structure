#!/usr/bin/perl
use strict;
use warnings;
use rnafold;
use fasta_hash;

#### DATA ####
### Read in data from system and save reference to variable.
my $fasta_file="./data/Splice_to_Start_Rand.txt";
my $fasta_file_fold="./data/5pUTR_37_rand.txt";
my $utr5_28="./data/5pUTR_28_rand.txt";
#my $fasta_file="./data/Splice_to_start.txt";
#my $fasta_file_fold="./data/splice_to_start_dotbracket_37.txt";
#my $utr5_28="./data/5pUTR_28.txt";

## 5'UTR fasta with no fold information
my $fh;
open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";

## 5'UTR fasta with with secondary structures at 37C
my $fh2;
open($fh2, $fasta_file_fold) or die "can't open $fasta_file_fold: $!\n";

## 5'UTR fasta with with secondary structures at 28C
my $fh3;
open($fh3, $utr5_28) or die "can't open $utr5_28: $!\n";

##############

## create appropriate hash files for each fasta file.
my %sequence_data;
my %fasta; ## no folded structures
while (fasta_hash::read_fasta_sequence($fh, \%sequence_data, 0)) {
    $fasta{$sequence_data{header}} = $sequence_data{seq};
}

my %fasta_struct; ## folted at 37C
while (fasta_hash::read_fasta_sequence($fh2, \%sequence_data, 1)) {
    $fasta_struct{$sequence_data{header}} = $sequence_data{struct};
}

my %fasta_utr5_28; # folded at 28C
while (fasta_hash::read_fasta_sequence($fh3, \%sequence_data, 1)) {
    $fasta_utr5_28{$sequence_data{header}} = $sequence_data{struct};
}

## Construct a new folded RNA object
my @genes = keys %fasta;
my $mfe37; my $mfe28; my $mdist; my $mdist_norm; my $hpins37; my $hpins28;
my $len;
print "Accession\tmfe37\tmfe28\tmdist\tmdistnorm\thpins37\thpins28\tlength";
print "\n";
while(my $acc = shift(@genes)) {
    my $RNA = rnafold::new($acc, 'names'=>['UTR5', 'UTR3'], 'temps'=>[37, 28]);
    $RNA->read_seq(\%fasta, 'UTR5');
    $RNA->read_dotbracket(\%fasta_struct, 'UTR5', 37);
    $RNA->read_dotbracket(\%fasta_utr5_28, 'UTR5', 28);

    $mfe37 = $RNA->mfe('UTR5', 37);
    $mfe28 = $RNA->mfe('UTR5', 28);
    $mdist = $RNA->mountain_dist('UTR5', 37, 28, 0);
    $mdist_norm = $RNA->mountain_dist('UTR5', 37, 28, 1);
    $hpins37 = $RNA->num_hairpins('UTR5', 37);
    $hpins28 = $RNA->num_hairpins('UTR5', 28);
    $len = $RNA->length('UTR5');

    print "$acc $mfe37 $mfe28 $mdist $mdist_norm $hpins37 $hpins28 $len";
    print "\n";
}
