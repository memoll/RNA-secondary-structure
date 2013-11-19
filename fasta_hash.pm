#!/usr/bin/perl
use strict;
use warnings;
package fasta_hash;

###############################################################################
## This module reads in a fasta file and will parse it into a hash.          ##
## Super convenient. Accessions are keys and sequences or structures are     ##
## values.                                                                   ##
## Code hacked to suit the needs of RNA dotbracket files from this URL:      ##
## http://code.izzid.com/2011/10/31/How-to-read-a-fasta-file-in-perl.html    ##
##                                                                           ##
## Arguements: fh = pointer to fasta file.                                   ##
##             seq_info = where to save reference                            ##
##             folded = flag on whether normal fasta file or with dotbrackets##
###############################################################################

sub read_fasta_sequence {
    my ($fh, $seq_info, $folded) = @_;

    $seq_info->{seq} = undef; # clear out previous sequence

    if ($folded == 1) { $seq_info->{struct} = undef; }
    # put the header into place
    $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

    my $file_not_empty = 0;
    while (<$fh>) {
        $file_not_empty = 1;
        next if /^\s*$/;  # skip blank lines
        chomp;

        if (/^>/) { # fasta header line
            my $h = $_;
            $h =~ s/^>//;
            if ($seq_info->{header}) {
                $seq_info->{next_header} = $h;
                return $seq_info;
            }
            else { # first time through only
                $seq_info->{header} = $h;
            }
        }
        else {
            ## if user passes a fasta with dot bracket data file
            if($folded == 1) {
                ## if current line begins with A,C,G,U, save to sequence
                if (/^(A|C|G|U)/) { $seq_info->{seq} .= $_; }
                ## or if it begins with (,),., save to struct
                elsif (/^(\(|\)|\.)/) { $seq_info->{struct} .= $_; }
            }
            else { s/\s+//; $seq_info->{seq} .= $_; }
        }
    }

    if ($file_not_empty) {
        return $seq_info;
    }
    else {
        # clean everything up
        $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

        return;
    }
}

1;
