#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil floor);
use Cwd qw(abs_path getcwd);
use File::Path qw(make_path);

my $in_tgtdir;
my $in_fasta;
my $in_use_rnafold;
my $in_use_rnaplfold;
my $in_switch_plfold_length;
my $in_new_only;
my $in_verbose;

my $rnafold_params   = " --noLP ";
my $rnaplfold_params = " --noLP -c 0.0005 -L 100 -W 200 ";

my $tmp = "/var/tmp/foldfasta_$$";
my $in_vrna_path;

GetOptions(
    "tgtdir=s"         => \$in_tgtdir,
    "fasta=s"          => \$in_fasta,
    "rnafold"          => \$in_use_rnafold,
    "rnaplfold"        => \$in_use_rnaplfold,
    "switch-len=i"     => \$in_switch_plfold_length,
    "new-only"         => \$in_new_only,
    "rnaplfold-opts=s" => \$rnaplfold_params,
    "rnafold-opts=s"   => \$rnafold_params,
    "verbose"          => \$in_verbose,
    "vrna-path=s"      => \$in_vrna_path,
);

die
"Please use either --rnafold or --rnaplfold or provide a length where to switch. Exit...\n\n"
  if ( !$in_switch_plfold_length && !$in_use_rnafold && !$in_use_rnaplfold );

die "Please use either a fold method or a switch length. Exit...\n\n"
  if ( $in_switch_plfold_length && ( $in_use_rnafold || $in_use_rnaplfold ) );

die "Please use either --rnafold or --rnaplfold. Exit...\n\n"
  if ( $in_use_rnafold && $in_use_rnaplfold );

my $vrna_path = "";
$vrna_path = $in_vrna_path . "/" if ($in_vrna_path);

if ( !$vrna_path || !-e "$vrna_path" . "RNAfold" ) {
    my $loc = `which RNAfold`;
    chomp($loc);
    die "\nCannot find RNAfold binary! Exit...\n\n" if ( !$loc );
    die "\nCannot find RNAfold binary! Exit...\n\n" if ( !-e $loc );
    $loc =~ /(.*)\/[^\/]+$/;
    $vrna_path = $1 . "/";
}

my $CURRDIR = getcwd;
$in_tgtdir = $CURRDIR if ( !$in_tgtdir );

make_path($in_tgtdir);
$in_tgtdir = abs_path($in_tgtdir);

my @fa = read_fasta_file($in_fasta);

foreach my $id ( @{ $fa[1] } ) {

    my $seq         = $fa[0]->{$id};
    my $dp_filename = "$in_tgtdir/$id";
    my $len         = length($seq);

    ## check for existing structure constraints
    ## use always RNAfold if we have structure
    my $constr;
    $constr = $fa[3]->{$id}->{"#S"}  if ( exists $fa[3]->{$id}->{"#S"} );
    $constr = $fa[3]->{$id}->{"#FS"} if ( exists $fa[3]->{$id}->{"#FS"} );

    undef($constr) if ( $constr && length($constr) != length($seq) );

    #print "len:$len dp:$dp_filename seq:$seq\n";

    next if ( $in_new_only && ( -e $dp_filename ) );

    if ( ( !$in_switch_plfold_length && $in_use_rnafold ) || $constr ) {
        ## RNAfold
        fold_RNAfold( $seq, $dp_filename, $constr );

    }
    elsif ( !$in_switch_plfold_length && $in_use_rnaplfold ) {
        ## RNAplfold
        fold_RNAplfold( $seq, $dp_filename );

    }
    elsif (
        ( $in_switch_plfold_length && length($seq) < $in_switch_plfold_length )
        || $constr )
    {
        ## RNAfold
        fold_RNAfold( $seq, $dp_filename, $constr );
    }
    else {
        ## RNAplfold
        fold_RNAplfold( $seq, $dp_filename );
    }

}

################################################################################

sub fold_RNAfold {
    my $seq       = $_[0];
    my $file_name = $_[1];
    my $str       = $_[2];

    mkdir($tmp);

    my $fa_tmp = "$tmp/_tmpseq_.fa";
    open( FA, ">$fa_tmp" );
    print FA ">tmp\n$seq";
    print FA "\n$str" if ($str);
    close(FA);

    chdir($tmp);

    my $call_params = "-p --noPS $rnafold_params";
    $call_params .= " -C " if ($str);

 #print "cat $fa_tmp | $vrna_path/RNAfold $call_params 2>/dev/null 1>/dev/null";

    system(
        "cat $fa_tmp | $vrna_path/RNAfold $call_params 2>/dev/null 1>/dev/null"
    );

    die "Error! RNAfold call failed! Exit...n\n" if ( !-e "tmp_dp.ps" );

    system("mv tmp_dp.ps $file_name");

    chdir($CURRDIR);
    system("rm -R $tmp");
}

sub fold_RNAplfold {
    my $seq       = $_[0];
    my $file_name = $_[1];

    mkdir($tmp);
    my $fa_tmp = "$tmp/_tmpseq_.fa";
    open( FA, ">$fa_tmp" );
    print FA ">tmp\n$seq\n";
    close(FA);

    chdir($tmp);

#print "cat $fa_tmp | $vrna_path/RNAplfold $rnaplfold_params 2>/dev/null 1>/dev/null";
    system(
"cat $fa_tmp | $vrna_path/RNAplfold $rnaplfold_params 2>/dev/null 1>/dev/null"
    );
    die "Error! RNAplfold call failed! Exit...n\n" if ( !-e "tmp_dp.ps" );

    system("mv tmp_dp.ps $file_name");

    chdir($CURRDIR);
    system("rm -R $tmp");
}

sub read_fasta_file {
    my ( $file, $make_unique ) = @_;
    my $FUNCTION = "read_fasta_file in Sequences.pm";

    my $id         = "";
    my $seqstring  = "";
    my %fasta      = ();
    my %header     = ();
    my @order      = ();
    my $line       = "";
    my %meta       = ();
    my %seq_meta   = ();
    my $uniq_count = 0;

    open( IN_HANDLE, "<$file" )
      || die "ERROR in $FUNCTION: "
      . "Couldn't open the following file in package Tool,"
      . " sub read_fasta_file: $file\n";

    while ( $line = <IN_HANDLE> ) {
        chomp($line);

        # header (can contain one space after > symbol)
        if ( $line =~ /^\>\s?(\S+)\s*([\S*\s*]*)/ ) {
            if ($id) {
                if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {

                    #          $uniq_count++;
                    #          $id .= "_$uniq_count";
                    #          print "Warning! Make Seq-id unique! now $id\n";
                    die "ERROR in $FUNCTION: "
                      . "multiple sequence id '$id', consider using function "
                      . "read_fasta_with_nonunique_headers instead of read_fasta_file";
                }
                $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
                $fasta{$id} = $seqstring;
                $meta{$id}  = {%seq_meta};    ## anonymous hash reference

                $seqstring = "";
                undef(%seq_meta);
            }

            if ($make_unique) {
                $uniq_count++;
                $id = $1 . "_$uniq_count";
            }
            else {
                $id = $1;
            }

            my $head = $2;
            $head = "" if ( !$head );
            $header{$id} = $head;
            push( @order, $id );
        }
        elsif ( $line =~ /(.+)\s+(#\S+)\s*$/ && $id ) {

            if ( exists $seq_meta{$2} ) {
                $seq_meta{$2} .= $1;
            }
            else {
                $seq_meta{$2} = $1;
            }

        }
        else {
            $seqstring .= $line if ($id);
        }
    }

    if ($id) {
        if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {

            #$uniq_count++;
            #$id .= "_$uniq_count";
            #print "Warning! Make Seq-id unique! now $id\n";
            die "ERROR in $FUNCTION: "
              . "multiple sequence id '$id', consider using function "
              . "read_fasta_with_nonunique_headers instead of read_fasta_file";
        }
        $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
        $fasta{$id} = $seqstring;
        $meta{$id}  = +{%seq_meta};    ## anonymous hash reference
        $seqstring  = "";
        undef(%seq_meta);

    }

    return ( \%fasta, \@order, \%header, \%meta );
}
