#!/usr/bin/perl

# duplex_finder.pl --- find miRNA* sequence for a given miRNA and its precursor

# Copyright (C) 2014  KANG <kanglmfATgmailDOTcom>

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Commentary:

# Usage:

# duplex_finder.pl [--pre|-p precursor] [--mir|-m miRNA]

# precursor: miRNA precursor sequence
# miRNA:     miRNA sequence

# NOTE:

# By default, this program will search for miRNA* that has 2-nt 3'
# overhang, but some miRNAs in miRBase do not have this feature, such
# as zma-miR159b and zma-miR166c, both of which have only 1-nt
# overhang. So, change the constant 'OVERHANG' if necessary.


BEGIN {
    # check RNAfold
    # my $rnafold = '/usr/bin/RNAfold';
    die "RNAfold not found\n" unless (qx(which RNAfold) =~ /RNAfold/);
}

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use constant OVERHANG => 2;

# ==========================================================================

my $VERSION = '0.2';

my ($pre, $mir) = check_argv();
my ($mir_start_pos, $mir_end_pos) = get_mir_pos($pre, $mir);
my $part = get_structure_parts($pre, $mir, $mir_start_pos, $mir_end_pos);
my ($star_seq, $star_start_pos, $star_end_pos) = get_folding_info($part);

print "#star_seq\tstar_start\tstar_end\tstar_len\n";
print "$star_seq\t$star_start_pos\t$star_end_pos\t" . length($star_seq) . "\n";

# ==========================================================================

sub check_argv {
    my ($pre, $mir, $help);
    GetOptions(
               "pre|p=s" => \$pre,
                "mir|m=s" => \$mir,
                "help|h!" => \$help
              );
    usage() if $help;
    usage() unless $pre and $mir;
    $pre =~ /[ATCGU]+/i or die "Invalid precursor seq.\n";
    $mir =~ /[ATCGU]+/i or die "Invalid miRNA seq.\n";
    return $pre, $mir;
}

sub get_mir_pos {
    my ($pre, $mir) = @_;
    my ($mir_start_pos, $mir_end_pos);
    if ($pre =~ /$mir/i) {
        $mir_start_pos = length($`) + 1; # 1-based position on precursor
        $mir_end_pos = length($`) + length($&);
    } else {
        die "miRNA NOT matched on precursor.\n";
    }
    return $mir_start_pos, $mir_end_pos;
}

sub get_structure_parts {
    my ($pre, $mir, $mir_start_pos, $mir_end_pos) = @_;
    my $noPS = (qx(RNAfold --version) =~ /RNAfold 2.\d/) ? "--noPS" : "-noPS";
    my $structure = (split /\s/, qx(echo $pre|RNAfold $noPS|sed 1d))[0];
    $structure =~ y/(/L/;       #| to make REGEX more clear
    $structure =~ y/)/R/;       #| rather than '\(' or '\)'
    # die if there is no loop or >2 loops
    die "Invalid precursor structure.\n"
        if ($structure !~ /L.*?R/ or $structure =~ /L.*?R.*?L.*?R/);
    my $mir_struct = substr($structure, $mir_start_pos - 1,
                            $mir_end_pos - $mir_start_pos + 1
                           );
    my ($struct5, $struct3, $pre5, $pre3, $loop_len, $mir_arm);
    if ($structure =~ /^(.+L)(.*?)(R.+)$/) {
        # this REGEX must succeed (for data extraction only)
        $struct5 = $1;
        $loop_len = length $2;
        $struct3 = $3;
        $pre5 = substr($pre, 0, length($1));
        $pre3 = substr($pre, length($1) + length($2));
        $mir_arm = ($mir_end_pos <= length($pre5)) ? '5p' :
            ($mir_start_pos > (length($pre) - length($pre3))) ? '3p' : undef;
        die "miRNA not in 5p or 3p arm.\n" unless (defined $mir_arm);
    } else {
        die "Precursor structure REGEX failed.\n";
    }
    my $part = {
                'struct5'       => $struct5,
                'struct3'       => $struct3,
                'pre5'          => $pre5,
                'pre3'          => $pre3,
                'loop_len'      => $loop_len,
                'mir_arm'       => $mir_arm,
                'mir_struct'    => $mir_struct
               };
    # return $struct5, $struct3, $pre5, $pre3, $loop_len, $mir_arm, $mir_struct;
    return $part;
}

sub get_folding_info {
    my $part = shift;
    my ($struct5, $struct3, $pre5, $pre3, $loop_len, $mir_arm, $mir_struct)
        = @$part{qw/struct5 struct3 pre5 pre3 loop_len mir_arm mir_struct/};

    # common variables
    my $right_base_count;
    my ($mir_5p_unmatch, $mir_mid_len, $mir_3p_unmatch);

    my $base_counter1;
    my $base_counter1_match;
    my $star_counter1;
    my $star_end_pos;

    my $base_counter2;
    my $base_counter2_match;
    my $star_counter2;
    my $star_start_pos;

    my $star_seq;


    $struct5 = reverse $struct5;
    $pre5 = reverse $pre5;
    if ($mir_arm eq '5p') {
        $right_base_count = length($struct5) - $mir_end_pos;
        my $mir_struct_rev = reverse $mir_struct;
        if ($mir_struct_rev =~ /^(\.*)(L.+L)(\.*)$/) {
            $mir_3p_unmatch = length $1;
            $mir_mid_len    = length $2;
            $mir_5p_unmatch = length $3;
        } else {
            die "\$mir_struct_rev not matched.\n";
        }
        #---------------------------------------------------------------
        # calculate star_start_pos (in 3p arm)
        $base_counter1 = $right_base_count + $mir_3p_unmatch;
        # $base_counter1_match includes the right-most match base within miRNA
        $base_counter1_match = (substr($struct5, 0, $base_counter1) =~ y/L//);
        $star_counter1 = ($struct3 =~ /^((R\.*){$base_counter1_match})/) ? length($1) + 1 : undef;
        die "Invalid \$star_counter1 (miRNA: 5p)\n" unless (defined $star_counter1);
        my $star_start_pos_in_pre3 = $star_counter1 - $mir_3p_unmatch + OVERHANG; # 1-based position
        $star_start_pos = length($pre5) + $loop_len + $star_start_pos_in_pre3;
        #---------------------------------------------------------------
        # calculate star_end_pos (in 3p arm)
        $base_counter2 = $base_counter1 + $mir_mid_len;
        # $base_counter2_match includes the left-most match base within miRNA
        $base_counter2_match = (substr($struct5, 0, $base_counter2 - 1) =~ y/L//);
        $star_counter2 = ($struct3 =~ /^((R\.*){$base_counter2_match})/) ? length($1) + 1 : undef;
        die "Invalid \$star_counter2 (miRNA: 5p).\n" unless (defined $star_counter2);
        my $star_end_pos_in_pre3 = $star_counter2 + $mir_5p_unmatch + OVERHANG; # 1-based position
        $star_end_pos = length($pre5) + $loop_len + $star_end_pos_in_pre3;
        #---------------------------------------------------------------
        $star_seq = substr($pre3, $star_start_pos_in_pre3 - 1,
                           $star_end_pos_in_pre3 - $star_start_pos_in_pre3 + 1);
    } elsif ($mir_arm eq '3p') {
        my $mir_start_pos_in_pre3 = $mir_start_pos - length($pre5) - $loop_len;
        #my $mir_end_pos_in_pre3 = $mir_start_pos_in_pre3 + length($mir) - 1;
        my $right_base_count = $mir_start_pos_in_pre3 - 1;
        if ($mir_struct =~ /^(\.*)(R.+R)(\.*)$/) {
            $mir_5p_unmatch = length $1;
            $mir_mid_len    = length $2;
            $mir_3p_unmatch = length $3;
        } else {
            die "miRNA structure not matched.\n";
        }
        #---------------------------------------------------------------
        # calculate star_end_pos (in 5p arm)
        $base_counter1 = $right_base_count + $mir_5p_unmatch;
        $base_counter1_match = (substr($struct3, 0, $base_counter1) =~ y/R//);
        $star_counter1 = ($struct5 =~ /^((L\.*){$base_counter1_match})/) ? length($1) + 1 : undef;
        die "Invalid \$star_counter1 (miRNA: 3p).\n" unless (defined $star_counter1);
        my $star_end_pos_in_pre5 = $star_counter1 - $mir_5p_unmatch - OVERHANG;
        $star_end_pos = length($pre5) - $star_end_pos_in_pre5 + 1; # NOTE the trailing base
        #---------------------------------------------------------------
        # calculate star_start_pos (in 5p arm)
        $base_counter2 = $base_counter1 + $mir_mid_len;
        $base_counter2_match = (substr($struct3, 0, $base_counter2 - 1) =~ y/R//);
        $star_counter2 = ($struct5 =~ /^((L\.*){$base_counter2_match})/) ? length($1) + 1 : undef;
        die "Invalid \$star_counter2 (miRNA: 3p).\n" unless (defined $star_counter2);
        my $star_start_pos_in_pre5 = $star_counter2 + $mir_3p_unmatch - OVERHANG;
        $star_start_pos = length($pre5) - $star_start_pos_in_pre5 + 1;
        $pre5 = reverse $pre5;
        #---------------------------------------------------------------
        $star_seq = substr($pre5, $star_start_pos - 1,
                           $star_end_pos - $star_start_pos + 1);
    }
    return $star_seq, $star_start_pos, $star_end_pos;
}

sub usage {
    my $prog = basename $0;
    my $usage = <<EOF;
$prog v$VERSION
Find miRNA* sequence for a given precursor and an miRNA with 2-nt 3' overhang
Usage:
$prog [--pre|-p precursor] [--mir|-m miRNA]
EOF
    print STDERR $usage;
    exit 1;
}
