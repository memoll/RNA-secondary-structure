#!/usr/bin/perl
use strict; use warnings;
use List::MoreUtils qw(each_array);
use POSIX;
package rnafold;

## package needs to store accession (name) of gene, keys for folding params,
## temperatures folded at, seqence, dotbrackets and MFE.
sub new {
    my $self = shift @_;
    my %params = @_;

    return bless{
        'ACCESSION' => $self,
        'names' => $params{'names'},
        'temps' => $params{'temps'},
        'sequence' => {},
        'dotbracket' => {},
        'MFE' => {}
    }
}

## returns the accession when called on object
sub accession { return $_[0]->{'ACCESSION'}; }

## return the sequence if called
sub sequence {
    if(defined( $_[0]->{$_[1]}->{'sequence'} )) {
        return $_[0]->{$_[1]}->{'sequence'};
    }
    else { kill "No sequence stored!"; }
}

## return the dotbracket if called
sub dotbracket {
    if(defined( $_[0]->{$_[1]}->{$_[2]}->{'dotbracket'} )) {
        return $_[0]->{$_[1]}->{$_[2]}->{'dotbracket'};
    }
    else { kill "No dotbracket structure stored here!"; }
}

## return the minimum free energy if called
sub mfe {
    if(defined( $_[0]->{$_[1]}->{$_[2]}->{'MFE'} )) {
        return $_[0]->{$_[1]}->{$_[2]}->{'MFE'};
    }
    else { kill "No MFE found!"; }
}

sub new_temp { unshift(@{$_[0]->{'temps'}},$_[1]); }

## need functions for gettin sequence, dot-bracket, expression,
## look up sequence from hash reference with accession as key
sub read_seq {
    (my $hashref, my $name) = ($_[1], $_[2]);
    my $seq = $hashref->{$_[0]->{'ACCESSION'} };
    if (defined($seq)) {
        $_[0]->{$name}->{'sequence'} = $seq;
        return $_[0]->{$name}->{'sequence'};
    }
    ## return 0 if there is no sequence for this name.
    else {return 0;}
}

##############################################################################
## Submodule for reading in the dotbracket file from an appropriate file.   ##
## Note that at the end of the line containing the dotbracket plot, there is##
## also the minimum free energy enclosed in parenthesis.                    ##
##############################################################################
sub read_dotbracket {
    (my $hashref, my $name, my $temp) = ($_[1], $_[2], $_[3]);
    my $db = $hashref->{$_[0]->{'ACCESSION'} };
    if (defined($db)) {
        my @vals = split(" ", $db, 2); # splits dotbracket plot from MFE
        $vals[1] =~ tr/\(\)/ /; ## remove the parenthesis around the MFE
        $_[0]->{$name}->{$temp}->{'dotbracket'} = $vals[0];
        $_[0]->{$name}->{$temp}->{'MFE'} = $vals[1];
    }
    else { kill "File isn't proper format."}
}

## methods: mountain distance, query motif, hausdorff distance, tree distance?
# (needs RNAforester system function), count hairpins, symmetric set diff,
# printing functions, length,

##############################################################################
## Returns length of structure.                                             ##
##############################################################################
sub length {
    return length($_[0]->{$_[1]}->{'sequence'});
}

##############################################################################
## subroutine will scale through secondary structure and +1 for (, -1 for ),##
## and + 0 for .                                                            ##
##############################################################################
sub mountain_dist {
    ## mountain function for a structure
    sub mountain {
        my $structure = $_[0];
        my @mount;
        my @struct = split //,$structure;
        my $i = 0;
        while(my $ele = shift(@struct)) {
            if ($ele eq "\.") { push @mount,$i; }
            elsif ($ele eq "\(") { push @mount,$i++; }
            else { push @mount,$i--; }
        }
        return @mount;
    }
    my($name, $temp1, $temp2, $normalized) = ($_[1], $_[2], $_[3], $_[4]);

    ## check to make sure comparing on equal length structures.
    my $struct1 = $_[0]->{$_[1]}->{$temp1}->{'dotbracket'};
    my $struct2 = $_[0]->{$_[1]}->{$temp2}->{'dotbracket'};
    if(CORE::length($struct1) != CORE::length($struct2)) {
        die "Can only compare structures of same length\n";
    }

    my @mount1 = mountain($struct1);
    my @mount2 = mountain($struct2);

    ## mountain distance: at each position, take difference of mountain
    ## function at that position between two structures. Use modulus of
    ## of difference for manhattan distance (can also use euclidean).
    ## Sum across all positions.

    ## Find distance. This loops through each element of both arrays and
    ## takes the modulus of the difference.
    my $dist = 0;
    my $it = List::MoreUtils::each_array( @mount1, @mount2 );
    while (my ($x, $y) = $it->() ) {
        $dist += abs($x-$y)**2;
    }

    ## if distance requested not to be normalized by length, just return this.
    if($normalized != 1) { return $dist; }
    ## normalize
    else {
        my $n = CORE::length($struct1);
        ## Need to normalize by the maximal distance for a string of length n.
        ## This is the distance between a completely open string and a string
        ## that is a hairpin loop with only three unpaired bases.

        my $s_open = '.'x$n; ## open string

        ## s_max: hairpin loop of length n with three unpaired bases in center
        my $s_max = ('(')x(POSIX::floor(($n-3)/2))."...";
        $s_max .= (")")x(POSIX::floor(($n-3)/2));
        if (CORE::length($s_max) < $n) { $s_max .= "."; }

        my @open_mount = mountain($s_open);
        my @max_mount  = mountain($s_max);
        ## as before, find distance.
        my $maxdist = 0;
        my $itx = List::MoreUtils::each_array( @open_mount, @max_mount );
        while (my ($a, $b) = $itx->() ) {
            $maxdist += abs($a-$b)**2;
        }
        return sqrt($dist)/sqrt($maxdist);
    }
}

## Symmetric set difference. Very crude measure, will finish if time
## time permitted.
sub sym_diff{
    my($name, $temp1, $temp2) = ($_[1], $_[2], $_[3]);
}

##############################################################################
## This subroutine reads reads in a name and temperature, will query the    ##
## appropriate dotbracket string and use a regex to find number of hairpins.##
## Hair pins are defined as two parenthesis enclosing between 3 and 5 .'s   ##
##############################################################################
sub num_hairpins{
    my($name, $temp) = ($_[1], $_[2]);
    my $count = 0;
    $count++ while ($_[0]{$name}{$temp}{'dotbracket'} =~ m/\(\.{1,5}\)/g );
    return $count;
}

1;
