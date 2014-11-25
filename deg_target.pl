#!/usr/bin/perl

# deg_target.pl --- target prediction and degradome-based target identification
# for miRNAs

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

# deg_target.pl [--srna|-s FILE] [--mrna|-m FILE] [--fasta|-a FILE] [args...]

# Please invoke this script without argument for more usage detail

# NOTE:

# If any miRNA is shorter than 21 nt, then another miRNA fasta input for
# fasta36 alignment should be generated to make those miRNAs >= 21 nt.

# To prepare fasta input for miRNAs shorter than 21 nt:
# cat $srna_fasta | perl -wpe '/^[^>]/ && (length() < 22) && s/$/"N" x (22 - length())/e' >$srna_fasta.N
# then use $srna_fasta.N for fasta36 alignment, and the original $fa for this
# script.

# To do fasta36 alignment:
# fasta36 -C 15 -E 30000 -f -20 -g -20 -i -m 0 -n -q -r +15/-10 -U -W 4 $srna_fasta $mrna_fasta >$fasta_result
# where '-E 30000' is tested to be sufficient for not omitting putative targets
# in maize.

# To do bowtie alignment:
# bowtie-build -f -r $mrna_fasta $bowtie_prefix
# bowtie -f -v 0 -k 2 -m 5 -p 16 --norc --best --fullref --al ${srna_fasta}.vmode0 $bowtie_prefix $srna_fasta ${srna_fasta}.hit

# To run this script:
# perl deg_target.pl -f deg -s $srna_fasta -m $mrna_fasta -b $bowtie_result -l 27 -a $fasta36_result -p 0.1


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Spec;
# use Data::Dumper;

my $VERSION = '0.2';
my $argv = join ' ', @ARGV;
my $pwd = $ENV{'PWD'};
my ($func, $srna_db, $mrna_db, $fasta_result, $cutoff, $predict_out,
    $bowtie, $deg_length, $p_value, $deg_out, $target_out, $plot_dir, $help);

GetOptions(
    "func|f=s"          => \$func,          # default: 'deg'
    "srna|s=s"          => \$srna_db,       # required
    "mrna|m=s"          => \$mrna_db,       # required
    "fasta|a=s"         => \$fasta_result,  # required
    "cutoff|c=s"        => \$cutoff,        # default: 4.5
    "predictout|r=s"    => \$predict_out,   # default: 'out_prediction'

    "bowtie|b=s"        => \$bowtie,        # required for deg
    "deglen|l=i"        => \$deg_length,    # required for deg
    "pvalue|p=f"        => \$p_value,       # default: 0.05
    "degout|d=s"        => \$deg_out,       # default: 'out_degradome_hit'
    "targetout|t=s"     => \$target_out,    # default: 'out_targets'
    "plotdir|o=s"       => \$plot_dir,      # default: 'out_plot'

    "help|h!"           => \$help           # optional
);

my $bp;
$bp->{'AU'} = [0,   '|'];
$bp->{'UA'} = [0,   '|'];
$bp->{'GC'} = [0,   '|'];
$bp->{'CG'} = [0,   '|'];
$bp->{'GU'} = [0.5, 'o'];
$bp->{'UG'} = [0.5, 'o'];
$bp->{'AC'} = [1,   ' '];
$bp->{'CA'} = [1,   ' '];
$bp->{'AG'} = [1,   ' '];
$bp->{'GA'} = [1,   ' '];
$bp->{'UC'} = [1,   ' '];
$bp->{'CU'} = [1,   ' '];
$bp->{'A-'} = [1,   ' '];
$bp->{'U-'} = [1,   ' '];
$bp->{'G-'} = [1,   ' '];
$bp->{'C-'} = [1,   ' '];
$bp->{'-A'} = [1,   ' '];
$bp->{'-U'} = [1,   ' '];
$bp->{'-G'} = [1,   ' '];
$bp->{'-C'} = [1,   ' '];
$bp->{'AA'} = [1,   ' '];
$bp->{'UU'} = [1,   ' '];
$bp->{'CC'} = [1,   ' '];
$bp->{'GG'} = [1,   ' '];

#=======================================================================
# BEGIN main program

check_argv();
message('START');
check_fasta_format($srna_db, $mrna_db);
my $srna_fasta_ref = read_srna($srna_db);
my $prediction = PREDICT_TARGET($srna_db, $mrna_db);
if ($func eq 'predict') {
    message('END');
    exit 0;
}

# BEGIN degradome-based target identification
my ($mrna_seq_info, $mrna_valid_bases) = get_mrna_info($mrna_db);
my ($cleavage_info, $category_proportion) = parse_bowtie_output($bowtie);
my $deg_targets = GET_DEG_TARGET($prediction, $cleavage_info);
message('END');

# END main program
#=======================================================================

sub check_argv {
    usage() if ($help);
    $func = 'deg' unless (defined $func);
    if ($func eq 'deg') {
        unless (defined $srna_db) {
            print STDERR "sRNA FASTA file not specified.\n";
            usage();
        }
        unless (defined $mrna_db) {
            print STDERR "mRNA FASTA file not specified.\n";
            usage();
        }
        unless (defined $fasta_result) {
            print STDERR "FASTA alignment result not specified.\n";
            usage();
        }
        $cutoff = 4.5 unless (defined $cutoff);
        if ($cutoff < 0 or $cutoff > 4.5) {
            print STDERR "Prediction score cutoff should be >=0 and <=4.5.\n";
            usage();
        }
        $predict_out = 'out_prediction' unless (defined $predict_out);
        #---------------------------------------------------------------
        unless (defined $bowtie) {
            print STDERR "Bowtie alignment file not specified.\n";
            usage();
        }
        unless (defined $deg_length) {
            print STDERR "Degradome read length not specified.\n";
            usage();
        }
        $p_value = 0.05 unless (defined $p_value);
        if ($p_value <= 0 or $p_value > 1) {
            print STDERR "P-value should be >0 and <=1.\n";
            usage();
        }
        $deg_out    = 'out_degradome_hit' unless (defined $deg_out);
        $target_out = 'out_targets' unless (defined $target_out);
        $plot_dir   = 'out_plot' unless (defined $plot_dir);

        # print command line summary
        print STDERR <<EOF;
# $0 $argv

------------------------------------------------------------------------
sRNA fasta input:          $srna_db
mRNA fasta input:          $mrna_db
FASTA alignment input:     $fasta_result
Alignment score cutoff:    $cutoff
Prediction output:         $predict_out
Bowtie mapping input:      $bowtie
Degradome read length:     $deg_length
P-value cutoff:            $p_value
Degradome hit output:      $deg_out
Final targets output:      $target_out
T-plot output dir:         $plot_dir
------------------------------------------------------------------------

EOF
    } elsif ($func eq 'predict') {
        unless (defined $srna_db) {
            print STDERR "sRNA FASTA file not specified.\n";
            usage();
        }
        unless (defined $mrna_db) {
            print STDERR "mRNA FASTA file not specified.\n";
            usage();
        }
        unless (defined $fasta_result) {
            print STDERR "FASTA alignment result not specified.\n";
            usage();
        }
        $cutoff = 4.5 unless (defined $cutoff);
        if ($cutoff < 0 or $cutoff > 4.5) {
            print STDERR "Prediction score cutoff should >=0 and <=4.5.\n";
            usage();
        }
        $predict_out = 'out_prediction' unless (defined $predict_out);
        # the following vars should not be specified in 'prediction' mode
        if (defined $bowtie) {
            print STDERR "Warning: \$bowtie is ignored in prediction mode.\n";
        }
        if (defined $deg_length) {
            print STDERR "Warning: \$deg_length is ignored in prediction mode.\n";
        }
        if (defined $p_value) {
            print STDERR "Warning: \$p_value is ignored in prediction mode.\n";
        }
        if (defined $deg_out) {
            print STDERR "Warning: \$deg_out is ignored in prediction mode.\n";
        }
        if (defined $target_out) {
            print STDERR "Warning: \$target_out is ignored in prediction mode.\n";
        }
        if (defined $plot_dir) {
            print STDERR "Warning: \$plot_dir is ignored in prediction mode.\n";
        }
        print STDERR <<EOF;
# $0 $argv

------------------------------------------------------------------------
sRNA fasta input:          $srna_db
mRNA fasta input:          $mrna_db
FASTA alignment input:     $fasta_result
Alignment score cutoff:    $cutoff
Prediction output:         $predict_out
------------------------------------------------------------------------

EOF
    } else {
        print STDERR "Unrecognized mode specification.\n";
        usage();
    }
}

sub message {
    my $stage = shift;
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) =
        localtime(time);
    my $time = sprintf("%04d-%02d-%02d %02d:%02d:%02d",
        $year + 1900, $mon + 1, $mday, $hour, $min, $sec);
    print STDERR "# Job $stage at:\t$time\n";
}

# STOP if '>fasta_id' line in the DB file contains '\s' or '|', since FASTA3x
# will output partial seq id (in 'm8C' format) in this condition

sub check_fasta_format {  #(@)
    my @data = @_;
    my $num_of_seq = 0;     # just for running statistics
    for my $file (@data) {
        open (F, $file) or die "Can not open fasta file $file!\n";
        while (<F>) {
            next unless (/^>/);
            $num_of_seq++;
            chomp;
            die "Bad FASTA file format in $file line $.\n" if (/\s|\|/);
        }
        close F;
        print STDERR "Number of sequences in $file:\t$num_of_seq\n";
        $num_of_seq = 0;
    }
}

sub read_srna {   #($)
    my $file = shift;
    my $id_seq;
    open (F, $file) or die $!;
    while (<F>) {
        chomp (my $id = $_);
        $id =~ s/^>//;  # assuming the format of id line is fine
        chomp (my $seq = <F>);
        $id_seq->{$id} = $seq;
    }
    close F;
    return $id_seq;
}

sub PREDICT_TARGET {  #($$)
    my ($srna_db, $mrna_db) = @_;
    my $prediction;     # returned ref
    my $num_of_targets = 0;
    open (PREDICT, '>', $predict_out) or die "Can not open $predict_out!\n";
    print PREDICT "#sRNA_ID\tmRNA_ID\tscore\tmRNA_start\tmRNA_end\t\@cleavage_site\n";
    #-------------------------------------------------------------------
    my $prev_mrna_id;           # one mRNA may match >1 positions
    #my $fasta_cmd = "$fasta_bin -C 15 -E 30000 -f -20 -g -20 -i -m 0 -n -q -r +15/-10 -U -W 4 $srna_db $mrna_db >$fasta_result";
    open (FASTA, $fasta_result) or die "Can not open FASTA result!\n";
    while (<FASTA>) {
        if (/^>>([^\-]{2}.+?)\s/o) {  # line 1, new mRNA ID
            my $mrna_id     = $1;
            $prev_mrna_id   = $mrna_id;
            my ($srna_id, $srna_match_seq, $mrna_match_seq, $mrna_coord_start, $mrna_coord_end)
                = read_nine_lines($srna_fasta_ref, \*FASTA); # all vars are refs
            # all are undefs if the 1st returned val is undef
            next unless (defined $srna_id);
            my $info = align_info($srna_match_seq, $mrna_match_seq,
                                $mrna_coord_start, $mrna_coord_end);
            next unless (defined $info);
            $$srna_match_seq =~ y/AUCG/UAGC/;   # seq should be converted back after align_info()
            # push VALUE of var, or, all vars are the same ref address (SCALAR(...))
            push @{$prediction->{$$srna_id}},
            [$$srna_match_seq, $mrna_id, $$mrna_match_seq, $$mrna_coord_start, $$mrna_coord_end, $info];
            print PREDICT "$$srna_id\t$mrna_id\t$info->{'total_score'}\t$$mrna_coord_start\t$$mrna_coord_end\t\@$info->{'cleavage_site'}\n\n";
            print PREDICT "sRNA 3' $$srna_match_seq 5'\n";
            print PREDICT "        $info->{'align_char'}\n";
            print PREDICT "mRNA 5' $$mrna_match_seq 3'\n\n";
            $num_of_targets++;
            undef $srna_id; undef $srna_match_seq; undef $mrna_match_seq;
            undef $mrna_coord_start; undef $mrna_coord_end; undef $info;
            #$prediction->{$srna_id}->{'mRNA_ID'}           = $mrna_id;
            #$prediction->{$srna_id}->{'mRNA_match_seq'}    = $$mRNA_match_seq;
            #$prediction->{$srna_id}->{'mRNA_coord_start'}  = $$mrna_coord_start;
            #$prediction->{$srna_id}->{'mRNA_coord_end'}    = $$mrna_coord_end;
            #$prediction->{$srna_id}->{'total_score'}       = $info->{'total_score'};
            #$prediction->{$srna_id}->{'GU_score'}          = $info->{'GU_score'};
            #$prediction->{$srna_id}->{'core_score'}        = $info->{'core_score'};
            #$prediction->{$srna_id}->{'align_char'}        = $info->{'align_char'};
            #$prediction->{$srna_id}->{'cleavage_site'}     = $info->{'cleavage_site'};
        } elsif (/^>--$/o) {          # line 1, previous mRNA ID
            my $mrna_id = $prev_mrna_id;
            my ($srna_id, $srna_match_seq, $mrna_match_seq, $mrna_coord_start, $mrna_coord_end)
                = read_nine_lines($srna_fasta_ref, \*FASTA);
            next unless (defined $srna_id);
            my $info = align_info($srna_match_seq, $mrna_match_seq,
                                $mrna_coord_start, $mrna_coord_end);
            next unless (defined $info);
            $$srna_match_seq =~ y/AUCG/UAGC/; # seq should be converted back after align_info()
            push @{$prediction->{$$srna_id}},
            [$$srna_match_seq, $mrna_id, $$mrna_match_seq, $$mrna_coord_start, $$mrna_coord_end, $info];
            print PREDICT "$$srna_id\t$mrna_id\t$info->{'total_score'}\t$$mrna_coord_start\t$$mrna_coord_end\t\@$info->{'cleavage_site'}\n\n";
            print PREDICT "sRNA 3' $$srna_match_seq 5'\n";
            print PREDICT "        $info->{'align_char'}\n";
            print PREDICT "mRNA 5' $$mrna_match_seq 3'\n\n";
            $num_of_targets++;
            undef $srna_id; undef $srna_match_seq; undef $mrna_match_seq;
            undef $mrna_coord_start; undef $mrna_coord_end; undef $info;
        }
    }
    close FASTA;
    close PREDICT;
    print STDERR "Total predicted targets for all sRNAs:\t $num_of_targets\n";
    return $prediction;
}

sub read_nine_lines {
    my $srna_fasta_ref = shift;
    my $fh = shift;
    readline $fh;
    my $coord_line = readline $fh;  # line 3
    my ($srna_coord_end, $srna_coord_start, $mrna_coord_start1, $mrna_coord_end1);
    if ($coord_line =~ /\((\d+)-(\d+):(\d+)-(\d+)\)$/o) {
        $srna_coord_end     = $1;
        $srna_coord_start   = $2;
        $mrna_coord_start1  = $3;
        $mrna_coord_end1    = $4;
    } else {
        die "sRNA/mRNA coordinate line $. REGEX error in FASTA result.\n"
    }
    readline $fh;
    readline $fh;
    my $srna_line = readline $fh;   # line 6
    my ($srna_id, $srna_line_space, $srna_match_seq);
    if ($srna_line =~ /^(\S+)(\s+)(\S+)/o) {
        $srna_id            = $1;
        $srna_line_space    = $2;
        $srna_match_seq     = $3;
        return undef if (($srna_match_seq =~ y/-//) > 1);
    } else {
        die "sRNA line $. error in FASTA result.\n";
    }
    my $srna_match_length = length $srna_match_seq;
    # NOTE: for miRNAs shorter then 20 nt (or just 20 nt), fasta36 will *NOT*
    # align efficiently (and possibly not properly as well), therefore, extra
    # 'N's should be added to the end of miRNA to make sure that each miRNA
    # sequence is >= 21 nt.
    (my $srna_match_seq_real = $srna_match_seq) =~ s/N//gi;
    my $srna_trailing_N_count = $srna_match_length - length $srna_match_seq_real;
    my $srna_match_length_real = length $srna_match_seq_real;
    my $mrna_start_pos = length($srna_id) + length($srna_line_space) + $srna_trailing_N_count;
    readline $fh;
    my $mrna_line = readline $fh;   # line 8
    my $mrna_match_seq = substr($mrna_line, $mrna_start_pos, $srna_match_length_real);
    return undef if ($mrna_match_seq =~ /\s/);  # some mRNA match string may contain space!
    my $srna_mrna_string = $srna_match_seq_real . $mrna_match_seq;
    return undef if (($srna_mrna_string =~ y/-//) > 1); # consider both seq
    # $mrna_coord_start calculation should be based on $srna_length
    $srna_id =~ s/-+$//;
    die "Invalid sRNA ID in FASTA result line $..\n" unless (exists $srna_fasta_ref->{$srna_id});
    my $srna_length = length $srna_fasta_ref->{$srna_id};
    return undef if ($srna_match_length_real < $srna_length);    # some alignment may include only partial sRNA seq
    # NOTE: srna_length is the real length of miRNA, which does not include
    # the trailing 'N's.
    my $mrna_coord_start = $mrna_coord_start1 - ($srna_length - $srna_coord_end);
    my $mrna_coord_end = $mrna_coord_end1 + ($srna_coord_start - 1);
    readline $fh;
    readline $fh;
    return \$srna_id, \$srna_match_seq_real, \$mrna_match_seq, \$mrna_coord_start, \$mrna_coord_end;
}

sub align_info {  #($$$$)
    my ($srna_match_seq_rc, $mrna_match_seq) = @_[0, 1]; # seq refs
    (my $srna_seq_r = $$srna_match_seq_rc) =~ y/AUCG/UAGC/;
    my $mrna_seq = $$mrna_match_seq;
    my ($mrna_coord_start, $mrna_coord_end) = @_[2, 3]; # coord refs
    my ($mrna_start, $mrna_end) = ($$mrna_coord_start, $$mrna_coord_end);
    #-------------------------------------------------------------------
    my ($normal_score, $core_score, $gu_score, $align_char);
    my @srna_bases      = reverse(split //, $srna_seq_r);   # for convenience\
    my @mrna_bases      = reverse(split //, $mrna_seq);     # of looping below
    my $adjacent        = 0;
    my $is_base_srna    = 0;    # calculation of $cleavage_site is related to
    my $is_base_mrna    = 0;    # both valid base (AUCG) of sRNA and mRNA.
    my $BASE_POS        = 0;    # the real base position of sRNA
    my $cleavage_site;
    my $last_pair;
    for my $pos (1 .. scalar(@srna_bases)) {
        my $srna_base    = $srna_bases[$pos - 1];
        my $mrna_base    = $mrna_bases[$pos - 1];
        my $pair         = $bp->{$srna_base . $mrna_base};
        return undef unless (defined $pair);
        # NOTE: mRNA seq may exist non-AUCG! (eg. GHKMNRY)
        $BASE_POS += ($srna_base eq '-') ? 0 : 1; # incrementing before calculation
        $adjacent++  if (($BASE_POS > 12) and ($pair->[0] != 0)
                        and defined($last_pair) and ($last_pair->[0] != 0));
        $normal_score += (($BASE_POS >= 2) and ($BASE_POS <= 12)) ? 0 : $pair->[0];
        $core_score += (($BASE_POS >= 2) and ($BASE_POS <= 12)) ? ($pair->[0] * 2) : 0;
        $gu_score += (($pair eq 'GU') or ($pair eq 'UG')) ? 0.5 : 0;
        $align_char .= $pair->[1];
        $last_pair = $pair;
        if ($is_base_srna == 10) {      # pos = 11 now
            return undef if ($pair->[0] != 0 or $last_pair->[0] != 0);
            $cleavage_site = $mrna_end - $is_base_mrna + 1;
            # $cleavage_site is 1-based position of mRNA pairing with
            # the 10th base of sRNA. That is, RISC will cleave before
            # this position, which is used later to identify
            # 'real' targets by combining the degradome mapping data.
        }
        $is_base_srna += ($srna_base eq '-') ? 0 : 1; # increment after 'if'.
        $is_base_mrna += ($mrna_base eq '-') ? 0 : 1; # increment after 'if'.
    }
    $align_char = reverse $align_char;
    my $total_score = $normal_score + $core_score;
    return undef if (   $total_score    > $cutoff
                    or  $gu_score       > 4
                    or  $core_score     > 2
                    or  $adjacent       > 1);
    my $info;
    $info->{'total_score'}      = $total_score;
    $info->{'GU_score'}         = $gu_score;
    $info->{'core_score'}       = $core_score;
    $info->{'align_char'}       = $align_char;
    $info->{'cleavage_site'}    = $cleavage_site;
    return $info;
}

sub get_mrna_info { #($)
    my $mrna_db = shift;
    my $mrna_seq_info;      # returned ref
    my ($mrna_id, $mrna_len);
    my $mrna_valid_bases;   # returned ref
    my $valid_base_count;
    open (MRNA, $mrna_db) or die "Can not open $mrna_db\n";
    while (<MRNA>) {
        chomp;
        if (/^>/) {
            $mrna_len = 0;
            ($mrna_id = $_) =~ s/^>//;
        } else {
            $mrna_len += length $_;
            $valid_base_count += y/ATCGatcg//;  # may differ from length
        }
        $mrna_seq_info->{$mrna_id} = $mrna_len;
    }
    close MRNA;
    $mrna_valid_bases = $valid_base_count - ((scalar (keys %$mrna_seq_info)) * $deg_length);
    return $mrna_seq_info, $mrna_valid_bases;
}

sub parse_bowtie_output { #($)
    my $bowtie = shift;
    my $cleavage;       # cleavage sites per mRNA
    my $cleavage_info;  # updated cleavage sites info per mRNA, returned
    my $category;
    my $site;           # count of each cleavage site
    my $hit_per_deg;    # number of hit(s) on mRNA per degradome read
    my $category_proportion;
    open (BOWTIE, $bowtie) or die "Can not open bowtie result!\n";
    while (<BOWTIE>) {
        my @col = (split)[0 .. 4];
        my ($deg_id, $deg_count, $strand, $mrna_id, $map_start);
        if (($col[1] eq '+') or ($col[1] eq '-')) { # DO NOT use regex here
            ($deg_id, $strand, $mrna_id, $map_start) = @col[0 .. 3];
            $deg_count = 1;
        } elsif (($col[2] eq '+') or ($col[2] eq '-')) {
            ($deg_id, $deg_count, $strand, $mrna_id, $map_start) = @col;
        } else {
            die "Error: invalid bowtie result format in $bowtie.\n";
        }
        next unless ($strand eq '+');
        # $map_start in bowtie is 0-based, now convert it to 1-based position
        $map_start = $map_start + 1; # real start site (hereafter in this sub)
        push @{$cleavage->{$mrna_id}}, [$deg_id, $map_start];
        $site->{$mrna_id}->{$map_start} += $deg_count;
        $hit_per_deg->{$deg_id} += $deg_count;
    }
    close BOWTIE;
    # after gathering info above, sort by $mrna_id is never needed
    $category->{'cat0'} = 0;
    $category->{'cat1'} = 0;
    $category->{'cat2'} = 0;
    $category->{'cat3'} = 0;
    $category->{'cat4'} = 0;
    open (DEG, ">", $deg_out) or die "Can not open $deg_out!\n";
    for my $mrna_id (keys %$cleavage) {
        my $mrna_len = $mrna_seq_info->{$mrna_id};
        print DEG ">$mrna_id\t$mrna_len\n";
        my $deg_info = $cleavage->{$mrna_id};  # ref of AoA
        my $map_hits;
        for my $each_info (@$deg_info) { # some $map_start may be equal
            my $deg_id      = $each_info->[0];
            my $map_start   = $each_info->[1];
            my $deg_count   = $hit_per_deg->{$deg_id};
            $map_hits->{$map_start}->{'raw'}  = $site->{$mrna_id}->{$map_start};
            $map_hits->{$map_start}->{'norm'} += (1 / $deg_count);
            $map_hits->{$map_start}->{'uniq'} += ($deg_count == 1) ? 1 : 0;
        }
        #---------------------------------------------------------------
        # FIXME: category is defined by raw_read, not norm_read
        my @norm_hits = sort {$b <=> $a} map {$_->{'norm'}} values %$map_hits;
        my $norm_max = shift @norm_hits;
        # @norm_hits now exclude the 1st max hit and thus can be used to check
        # if there are other hits that are equal to the max hit, to determine
        # the category (is cat1?) of each hit by the grep() below.
        my $items = @norm_hits;
        my $norm_median =
            ($items == 0) ? undef :
            ($items == 1) ? (($norm_max + $norm_hits[0]) / 2) :
            (($items > 1) and ($items % 2 == 0)) ?
            (($norm_hits[$items/2] + $norm_hits[$items/2-1]) / 2) :
            ($norm_hits[int(@norm_hits / 2)]);
        my $is_cat1 = (defined $norm_median) ?
            grep {$_ == $norm_max} @norm_hits : undef; # determine cat1
        #---------------------------------------------------------------
        for my $each_info (@$deg_info) {
            my $deg_id      = $each_info->[0];
            my $map_start   = $each_info->[1];
            my $hit_raw     = $map_hits->{$map_start}->{'raw'};
            my $hit_norm    = $map_hits->{$map_start}->{'norm'};
            my $hit_uniq    = $map_hits->{$map_start}->{'uniq'};
            my $each_line   = "$deg_id\t$map_start\t$hit_raw\t$hit_norm\t$hit_uniq";
            my $mrna_category = undef;
            if ($hit_raw == 1) {    # then this hit must be of Cat_4
                print DEG "$each_line\tCat_4\n";
                $category->{'cat4'}++;
                $mrna_category = 'cat4';
            } elsif ($hit_norm == $norm_max) {    # hit >= 1 and hereafter
                # first check if defined
                if ((! defined $is_cat1) or $is_cat1 == 0) {
                    print DEG "$each_line\tCat_0\n";
                    $category->{'cat0'}++;
                    $mrna_category = 'cat0';
                } else {
                    print DEG "$each_line\tCat_1\n";
                    $category->{'cat1'}++;
                    $mrna_category = 'cat1';
                }
            } elsif ($hit_norm >= $norm_median) {
                print DEG "$each_line\tCat_2\n";
                $category->{'cat2'}++;
                $mrna_category = 'cat2';
            } else {
                print DEG "$each_line\tCat_3\n";
                $category->{'cat3'}++;
                $mrna_category = 'cat3';
            }
            push @{$cleavage_info->{$mrna_id}},
                [$deg_id, $map_start, $hit_raw, $hit_norm, $mrna_category];
        }
    }
    close DEG;
    $category_proportion->{'cat0'} = $category->{'cat0'} / $mrna_valid_bases;
    $category_proportion->{'cat1'} = $category->{'cat1'} / $mrna_valid_bases;
    $category_proportion->{'cat2'} = $category->{'cat2'} / $mrna_valid_bases;
    $category_proportion->{'cat3'} = $category->{'cat3'} / $mrna_valid_bases;
    $category_proportion->{'cat4'} = $category->{'cat4'} / $mrna_valid_bases;
    return $cleavage_info, $category_proportion;
}
#-----------------------------------------------------------------------
# TODO: indexing by degradome map info or target prediction info?
sub GET_DEG_TARGET {  #($$)
    my ($prediction, $cleavage_info) = @_;
    my $cumulative_align;   # cumulative prediction 'event' under each score
    my ($equal_site_total, $equal_site_passed_pval);    # for stat only
    unless (-e $plot_dir) {
        mkdir $plot_dir or die "Can not mkdir: $plot_dir.\n";
    }
    for my $srna_id (sort keys %$prediction) {    # $srna
        my $all_mrna_info = $prediction->{$srna_id};
        for my $each_mrna_info (@$all_mrna_info) {
            my $mrna_total_score    = $each_mrna_info->[5]->{'total_score'};
            for (my $i = $mrna_total_score; $i <= 4.5; $i += 0.5) {
                $cumulative_align->{$srna_id}->{$i}++;
            }
        }
    }
    open (TARGET, ">", $target_out) or die "Can not open $target_out!\n";
    print TARGET "#sRNA_ID\tmRNA_ID\tdeg_read_ID\tcleavage_site\tdeg_hit_raw\tdeg_hit_norm\tp_value\n\n";
    my $i = 1;
    for my $srna_id (sort keys %$prediction) {    # $srna
        my $all_mrna_info = $prediction->{$srna_id};
        for my $each_mrna_info (@$all_mrna_info) {
            my $srna_match_seq      = $each_mrna_info->[0];
            my $mrna_id             = $each_mrna_info->[1];
            my $mrna_match_seq      = $each_mrna_info->[2];
            my $mrna_coord_start    = $each_mrna_info->[3];
            my $mrna_coord_end      = $each_mrna_info->[4];
            my $mrna_total_score    = $each_mrna_info->[5]->{'total_score'};
            my $align_char          = $each_mrna_info->[5]->{'align_char'};
            my $cleavage_site       = $each_mrna_info->[5]->{'cleavage_site'};
            #-----------------------------------------------------------
            if (exists $cleavage_info->{$mrna_id}) {
                for my $each_deg_map_info (@{$cleavage_info->{$mrna_id}}) {
                    my $deg_id          = $each_deg_map_info->[0];
                    my $deg_map_start   = $each_deg_map_info->[1];
                    my $deg_hit_raw     = $each_deg_map_info->[2];
                    my $deg_hit_norm    = sprintf "%.4f", $each_deg_map_info->[3];
                    my $deg_category    = $each_deg_map_info->[4];
                    # whether the position is present in both data
                    if ($cleavage_site == $deg_map_start) {
                        $equal_site_total++;
                        my $cumulative = $cumulative_align->{$srna_id}->{$mrna_total_score};
                        my $proportion = $category_proportion->{$deg_category};
                        my $p_val      = sprintf "%.4f", 1 - pbinom0($cumulative, $proportion);
                        # print STDERR "--$proportion\t$p_val--\n";
                        if ($p_val <= $p_value) {
                            $equal_site_passed_pval++;
                            print TARGET "$srna_id\t$mrna_id\t$deg_id\t$cleavage_site\t$deg_hit_raw\t$deg_hit_norm\t$p_val\n";
                            my $j = sprintf "%03d", $i++;
                            my $plot_fn = "${srna_id}_${mrna_id}_$j";
                            t_plot($plot_fn, $mrna_id, $cleavage_site, $deg_hit_raw, $mrna_total_score, $deg_category, $p_val, $cleavage_info->{$mrna_id});
                            # last;
                        }
                    }
                    # else { print STDERR "--unequal site--:\t$cleavage_site\t$deg_map_start\n";}
                }
            }
            # else { print STDERR "No match for mrna_id: $mrna_id\n"; }
        }
        $i = 1;
    }
    close TARGET;
    print STDERR "Total mRNAs that have equal sites in degradome-seq and target prediction: $equal_site_total\n";
    print STDERR "Total targets identified: $equal_site_passed_pval\n";

    # generate plots by Rscript (for linux only)
    if (($^O =~ /linux|cygwin/i) and (qx/Rscript --version 2>&1/ =~ /^R /)) {
        print STDERR "Generating t_plots by R ...";
        chdir $plot_dir or die "Can not chdir to $plot_dir.\n";
        map {system "Rscript --vanilla $_"} grep {-f} glob "*.R";
        chdir $pwd or die "Can not chdir to $pwd.\n";
        print STDERR "Done\n";
    } else {
        print STDERR "R scripts for T-plots can be found in dir: $plot_dir.\n";
    }
}

sub pbinom0 {
    my ($trial, $probability) = @_;
    return (1 - $probability) ** $trial;
}

sub t_plot {
    my ($plot_fn, $mrna_id, $cleavage_site, $deg_hit_raw, $mrna_total_score,
        $deg_category, $p_val, $all_deg_map_info) = @_;
    my $max_x_pos   = (sort {$b <=> $a} map {$_->[1]} @$all_deg_map_info)[0];
    my $max_y_pos   = (sort {$b <=> $a} map {$_->[2]} @$all_deg_map_info)[0];
    my $xlim        = int(1.1 * $max_x_pos);
    my $ylim        = int(1.1 * $max_y_pos);
    my $Rscript_filename = "Rcmd_$plot_fn.R";
    my $R_out = File::Spec->catfile($plot_dir, $Rscript_filename);
    open(ROUT, '>', $R_out) or die "Can not open file for write: '$R_out'.\n";

    #-------------------------------------------------------------------
    print ROUT
"# R script for plotting cleaved sites of mRNA: $mrna_id

pdf('plot_$plot_fn.pdf')
plot(
    $cleavage_site,
    $deg_hit_raw,
    type='h',
    xlim=c(0, $xlim),
    ylim=c(0, $ylim),
    col='red',
    lwd=2,
    xlab=\"Position of $mrna_id (nt)\",
    ylab='Read Count (raw)'
)
mtext(\"Prediction score: $mrna_total_score; Category: $deg_category; p-value: $p_val; Cleavage site: $cleavage_site\")
";
    print ROUT for
        map {"lines(c($_->[1], $_->[1]), c(0, $_->[2]))\n"}
        grep {$_->[1] != $cleavage_site} @$all_deg_map_info;
    print ROUT "\n# $Rscript_filename ends here\n";
    #-------------------------------------------------------------------

    close ROUT;
    # NOTE: ROUT will not be removed, and thus it can be modified
    # to be replotted on user's demand and then removed manually.
    # system "Rscript --vanilla $Rscript_filename";
}

sub usage {
    (my $prog = basename $0) =~ s/\.pl$//;
    print STDERR
"
                        $prog v$VERSION
Target Prediction and Degradome-based Target Identification for miRNAs.

Usage:
$prog [--srna|-s FILE] [--mrna|-m FILE] [--fasta|-a FILE] [args...]

Options/Arguments:
==================

Required:
---------
  -s, --srna FILE         miRNA sequence(s) file, fasta format. Please ensure
                          that the header line contains [A-Za-z0-9_] only.

  -m, --mrna FILE         mRNA sequence(s) file, fasta format. Please ensure
                          that the header line contains [A-Za-z0-9_] only.

  -a, --fasta FILE        FASTA alignment result in default (-m 0) format.

Required for 'deg' mode only:
-----------------------------
  -b, --bowtie FILE       Bowtie alignment result in Tab-delimited format.

  -l, --deglen INT        Length of degradome sequencing (clean) reads.

Optional:
---------
  -f, --func MODE         Either 'predict' or 'deg' (default). If 'predict'
                          is given, then it will perform FASTA3x-based
                          target prediction ONLY. By default the program will
                          do target prediction and then degradome-based
                          target identification.

  -c, --cutoff INT/FLOAT  Total score cutoff used for target prediction,
                          default is 4.5 (for plant species).

  -r, --predictout FILE   FASTA3x-based target prediction output,
                          defalut is 'out_prediction'.

  -h, --help              Print this help message.

Optional for 'deg' mode only:
-----------------------------
  -p, --pvalue FLOAT      p-value used to calculate 'random' alignment
                          probability under binomial distribution,
                          default is 0.05.

  -d, --degout FILE       Processed degradome/mRNA mapping information,
                          default is 'out_degradome_hit'.

  -t, --targetout FILE    Final identified targets based on both sRNA/mRNA
                          and degradome/mRNA mapping information,
                          default is 'out_targets'.
  -o, --plotdir DIR       Directory for writing plot data (R scripts and
                          outputs), default is 'out_plot'.

";
    exit 1;
}
