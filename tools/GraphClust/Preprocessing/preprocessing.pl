##!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil floor);
use Cwd qw(abs_path getcwd);
use List::Util 'shuffle';
use List::Util qw/ min max /;

my $tgtdir        = 'FASTA';
my $in_fasta      = $ARGV[0];
my @fa_in         = read_fasta_file( $in_fasta, 1 );
my $in_prefix     = "data";
my $SEQPREFIX     = "SEQ";
my $max_N_stretch = 15;
my $max_length    = $ARGV[1];
my $in_winShift    = $ARGV[2];         ### winshift value from config
my $fa_name        = $in_prefix;
my $min_seq_length = $ARGV[3];
my @split_wins     = ($max_length);
my $splitWin       = $split_wins[0];
my $frags_splitN = splitMasked( \@fa_in, $max_N_stretch );
my $frags_shift =
  splitLengthShift( $frags_splitN, $splitWin, $in_winShift, $max_N_stretch );
my $frags_keep  = [ keys %{$frags_shift} ];
my $genome_locs = genomeLocations( \@fa_in );

system("mkdir -p $tgtdir");

writeFrags( \@fa_in, $frags_shift, $frags_keep, "$tgtdir/$fa_name" )
  ;    ##skip grey_href
writeFiles( \@fa_in, $frags_splitN, $genome_locs );    ## skip grey_href

system("zip -r FASTA.zip FASTA > /dev/null");

print
"\n#########################################################################\n";
print " STATS\n\n";
print "input seqs                                    : "
  . ( @{ $fa_in[1] } ) . "\n";
print "    fragments (after removing masked regions) : "
  . ( map { @{$_} } @{$frags_splitN} ) . "\n";
print "    fragments (after using windows/shift)     : "
  . ( keys %{$frags_shift} ) . "\n";
print
  "#########################################################################\n";

sub writeFrags {
    my $fa        = $_[0];    ## fastan file
    my $frags     = $_[1];
    my $fragsKeep = $_[2];
    my $prefix    = $_[3];

    ## files with kept fragments
    open( NA,   ">$prefix.names" );
    open( SEQS, ">$prefix.fasta" );
    open( LEN,  ">$prefix.lens" );
    open( MAP,  ">$prefix.map" );

    my $seq_idx     = 1;
    my @keep_sorted = sort {
        ( $a =~ /^$SEQPREFIX(\d+)#/ )[0] <=> ( $b =~ /^$SEQPREFIX(\d+)#/ )[0]
          || ( $a =~ /[^#]+#(\d+)/ )[0] <=> ( $b =~ /[^#]+#(\d+)/ )[0]
          || $a cmp $b
    } @{$fragsKeep};

    foreach my $key (@keep_sorted) {

        my @ent = split( "#", $key );
        $ent[0] =~ /^$SEQPREFIX(\d+)/;
        my $sIdx    = $1;
        my $orig_id = $fa->[1]->[ $sIdx - 1 ];

        print NA "$seq_idx $key ORIGID $orig_id ORIGHEAD "
          . $fa->[2]->{$orig_id} . "\n";
        print SEQS ">$seq_idx $key ORIGID $orig_id ORIGHEAD "
          . $fa->[2]->{$orig_id} . "\n";
        print SEQS $frags->{$key} . "\n";

        ## add meta information, but avoid for splitted seqs
        map {
            print SEQS $fa->[3]->{$orig_id}->{$_} . " $_\n"
              if (
                length( $frags->{$key} ) ==
                length( $fa->[3]->{$orig_id}->{$_} ) )
        } keys %{ $fa->[3]->{$orig_id} };

        print LEN length( $frags->{$key} ) . "\n";
        print MAP $seq_idx . " " . $key . "\n";

        # print GL $seq_idx . "\n" if ( $gl->{$sIdx} == 1 );

        $seq_idx++;
    }

    close(NA);
    close(SEQS);
    close(LEN);
    close(MAP);

    #close(GL) if ($write_gl);

    ## file with all frags, not blastclust filtered
    my @all_sorted = sort {
        ( $a =~ /^$SEQPREFIX(\d+)#/ )[0] <=> ( $b =~ /^$SEQPREFIX(\d+)#/ )[0]
          || ( $a =~ /[^#]+#(\d+)/ )[0] <=> ( $b =~ /[^#]+#(\d+)/ )[0]
          || $a cmp $b
    } keys %{$frags};
    open( IM, ">$prefix.fasta.all_frags" );

    foreach my $key (@all_sorted) {
        $key =~ /^$SEQPREFIX(\d+)#/;
        my $idx = $1;
        my $id  = $fa->[1]->[ $idx - 1 ];
        print IM ">$key ORIGID $id ORIGHEAD " . $fa->[2]->{$id} . "\n";
        print IM $frags->{$key} . "\n";
    }
    close(IM);
}

sub splitMasked {
    my $fa    = $_[0];
    my $max_N = $_[1];

    my @all_frag = ();

    foreach my $idx ( 0 .. $#{ $fa->[1] } ) {

        #print " next seq $idx\n";

        my @seq_a = split( "", $fa->[0]->{ $fa->[1]->[$idx] } );

        my $curr_frag   = "";
        my $curr_N_frag = "";
        my @seq_frags   = ();

        while ( @seq_a > 0 ) {

            my $curr_nt = shift(@seq_a);

            while ( $curr_nt eq "N" ) {
                $curr_N_frag .= $curr_nt;
                $curr_nt = shift(@seq_a);
            }

            while ( $curr_nt ne "N" && @seq_a ) {
                $curr_frag .= $curr_nt;
                $curr_nt = shift(@seq_a);
            }

            ## use always last nt, also "N"
            $curr_frag .= $curr_nt if ( !@seq_a );

            if ( length($curr_N_frag) <= $max_N ) {

                ## append at last fragment, get last idx first
                $seq_frags[ 0 < $#seq_frags ? $#seq_frags : 0 ] .=
                  ( $curr_N_frag . $curr_frag );

            }
            elsif ( length($curr_N_frag) > $max_N ) {
                push( @seq_frags, $curr_N_frag );
                push( @seq_frags, $curr_frag ) if ( length($curr_frag) > 0 );
            }

            $curr_N_frag = "";
            $curr_N_frag = $curr_nt if ( $curr_nt eq "N" );

            $curr_frag = "";

        }    ## while curr seq

        my $len = 0;
        map { $len += length($_) } @seq_frags;

        die "split error! Different length after splitting! Exit..\n\n"
          if ( $len != length( $fa->[0]->{ $fa->[1]->[$idx] } ) );

        #    print join("\n",@seq_frags)."\n";
        push( @all_frag, \@seq_frags );

    }    ## foreach seq

    return \@all_frag;
}

sub splitLengthShift {
    my $seqs     = $_[0];
    my $winSize  = $_[1];
    my $winShift = $_[2];
    my $max_N    = $_[3];

    my %usedSeqs  = ();
    my $shift_abs = 0;

    if ( $winSize > 0 ) {
        $shift_abs = floor( ( $winSize / 100 ) * $winShift );
        $shift_abs = 1        if ( $shift_abs < 1 );
        $shift_abs = $winSize if ( $shift_abs > $winSize );
    }
    foreach my $seq_idx ( 1 .. @{$seqs} ) {

        my $frags = $seqs->[ $seq_idx - 1 ];

        my $pos_abs_end   = 0;
        my $pos_abs_start = 0;

        #my $seqKey = "$SEQPREFIX$seq_idx";

        foreach my $frag_idx ( 1 .. @{$frags} ) {

            my $seq_frag = $frags->[ $frag_idx - 1 ];

            $pos_abs_end += length($seq_frag);
            $pos_abs_start = $pos_abs_end - length($seq_frag) + 1;

            next if ( length($seq_frag) < $min_seq_length );
            next if ( $seq_frag =~ /N{$max_N,}/ );

            ## N's at start end end we can ignore always
            ## we need only adjust $pos_abs_start by $start_cue later
            ## end cueing can be ignored
            #$seq_frag =~ s/^([N]+)//;
            my $start_cue = 0;
            $start_cue = length($1)
              if ( $seq_frag =~ s/^([N]+)// and length($1) > 0 );
            $seq_frag =~ s/[N]+$//;

            my $end_pos = length($seq_frag);
            $end_pos = $winSize
              if ( $winSize > 0 && length($seq_frag) > $winSize );
            my $start_idx = 0;
            my $used_len  = length($seq_frag);

            do {

                ## check if current frag is last, if then use winSize length from end
                ## if too small overhang then shrink frag a bit

                if ($winSize) {

                    ## overlap of two fragments, defined by $winShift
                    my $overlap = $winSize - $shift_abs;
                    ## tolerance for certain fragmets (within 10% length diff)
                    my $special_len = $winSize * 0.1;

                    ## idea: take always frags of $winSize, only for last frags
                    ## use special length and startidx
                    ## debug
                    my $last_len = length($seq_frag) - $start_idx;
                    if ( length($seq_frag) <= $winSize ) {
                        ## case: frag is between $min_seq_length and $winSize
                        $used_len = length($seq_frag);

                        # print "case 0 last $last_len used $used_len\n";
                    }
                    elsif ( length($seq_frag) - $start_idx <= $winSize * 0.5 ) {
                        ## case: last is very short, take at least $min_seq_length, or 60% of $winSize
                        $used_len = max( $min_seq_length, $winSize * 0.6 );
                        $start_idx = length($seq_frag) - $used_len;

                        #print "case 1 last $last_len used $used_len\n";
                    }
                    elsif (
                        length($seq_frag) - $start_idx - $overlap <=
                        $shift_abs + $special_len )
                    {
                        ## case: last frag is in tolerance range +- 10% $winSize
                        $used_len = max( $winSize,
                            length($seq_frag) - $start_idx - $overlap );
                        $start_idx = length($seq_frag) - $used_len;

                   #print "case 2 last $last_len used $used_len seq $seq_idx\n";
                    } #elsif ( length($seq_frag) - $start_idx <= $winSize + $special_len ) {
                      #$used_len = length($seq_frag) - $start_idx;
                      #print "case 3 last $last_len used $used_len\n";
                      #}
                    else {
                        ## case: if last frag is long enough then take just $winSize
                        $used_len = $winSize;

                        #print "case 5 last $last_len used $used_len\n";
                    }
                }

                my $curr_seq = substr( $seq_frag, $start_idx, $used_len );

                my $abs_pos = ( $pos_abs_start + $start_cue + $start_idx ) . "#"
                  . ( $pos_abs_start + $start_idx + length($curr_seq) - 1 );
                my $key = $SEQPREFIX . $seq_idx . "#" . $abs_pos;
                $usedSeqs{ $key . "#+" } = $curr_seq;

                $start_idx += $shift_abs;

#print "$key len frag ".length($seq_frag)." start_idx $start_idx ".($start_idx-$shift_abs)." end_pos $end_pos spec $used_len\n";
              } while (
                $start_idx - $shift_abs + $used_len < length($seq_frag) );
        }
    }

    return \%usedSeqs;
}

sub writeFiles {
    my $fa       = $_[0];
    my $f_splitN = $_[1];
    my $locas    = $_[2];

    # my $gl_href  = $_[3];

    open( FA, ">$tgtdir/orig.fasta" );
    open( SC, ">$tgtdir/data.fasta.scan" );
    open( LO, ">$tgtdir/data.locations" );
    open( FR, ">$tgtdir/fragments.splitMasked.list" );

    foreach my $idx ( 1 .. @{ $fa->[1] } ) {

        my $id = $fa->[1]->[ $idx - 1 ];

        print FA ">$id " . $fa->[2]->{$id} . "\n" . $fa->[0]->{$id} . "\n";

        ## meta information from fasta
        map { print FA $fa->[3]->{$id}->{$_} . " $_\n" }
          keys %{ $fa->[3]->{$id} };

        print SC ">$SEQPREFIX$idx ORIGID $id ORIGHEAD "
          . $fa->[2]->{$id} . "\n"
          . $fa->[0]->{$id} . "\n";
        print LO "$SEQPREFIX$idx $locas->[$idx-1]\n";

        my $pos_abs_end   = 0;
        my $pos_abs_start = 0;

        foreach my $seq_frag ( @{ $f_splitN->[ $idx - 1 ] } ) {
            $pos_abs_end += length($seq_frag);
            $pos_abs_start = $pos_abs_end - length($seq_frag) + 1;

            if ( length($seq_frag) < $min_seq_length ) {
                print FR
"$SEQPREFIX$idx#$pos_abs_start#$pos_abs_end SMALLER_THAN_$min_seq_length\_NT\n";
                next;
            }

            if ( $seq_frag =~ /N{$max_N_stretch,}/ ) {
                print FR
"$SEQPREFIX$idx#$pos_abs_start#$pos_abs_end MASKED_FRAGMENT\n";
                next;
            }

            print FR
              "$SEQPREFIX$idx#$pos_abs_start#$pos_abs_end USED_FRAGMENT\n";

        }
    }
    close(FA);
    close(SC);
    close(LO);
    close(FR);

    makeCol("$tgtdir/fragments.splitMasked.list");
}

sub genomeLocations {
    my $fa = $_[0];

    my @locs = ();

    foreach my $idx ( 1 .. @{ $fa->[1] } ) {

        if ( $fa->[2]->{ $fa->[1]->[ $idx - 1 ] } =~
            /GENOME_LOC (\S+\.chr\S+:\d+-\d+:\S+)/ )
        {
            push( @locs, $1 );
            next;
        }

        ## str = id." ".header
        my $str =
          $fa->[1]->[ $idx - 1 ] . " " . $fa->[2]->{ $fa->[1]->[ $idx - 1 ] };
        $str = lc($str);

        my $genome = "undef";
        my $chr;
        my $start;
        my $end;
        my $strand;
        my $loc = "MISS";

        if ( $str =~ /chr(\S+)[:_](\d+)[-_](\d+)/ ) {
            $chr   = $1;
            $start = $2;
            $end   = $3;

            $strand = "+" if ( $start <= $end );
            if ( $start > $end ) {
                $strand = "-";
                ## we use always start before end, even for minus strand
                my $tmp = $start;
                $start = $end;
                $end   = $tmp;
            }

            if ( $str =~ /chr\S+[:_]\d+[-_]\d+[:_](\S+)/ ) {
                $strand = $1;
            }

            if ( $str =~ /strand=([+-])/ && $start <= $end ) {
                $strand = $1;
            }

            if ( $str =~ /(\S+)[\._]chr\S+[:_]\d+[-_]\d+/ ) {
                $genome = $1;
            }

            $loc = "$genome.chr$chr:$start-$end:$strand";

        }

        push( @locs, $loc );

    }

    return \@locs;

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

sub makeCol {
    my $col_file = $_[0];

    #  my $tmp_dir = "_col_tmp_" . $$;
    #  mkdir($tmp_dir);
    #  system("column -t $col_file > $tmp_dir/t");
    #  system("mv $tmp_dir/t $col_file");
    #  system("rm -R -f $tmp_dir");

    ## read content and spit into columns
    open( my $in_file, $col_file );
    my @cont    = ();
    my $max_col = 0;
    while ( my $line = <$in_file> ) {
        chomp($line);
        my @t = split( " ", $line );
        $max_col = scalar(@t) if ( $max_col < scalar(@t) );
        push( @cont, \@t );
    }
    close($in_file);

    ## get max length of all columns
    my @len = ();
    map { push( @len, 0 ) } 1 .. $max_col;

    foreach my $ci ( 0 .. $max_col - 1 ) {
        foreach my $l (@cont) {
            if ( @{$l} >= $ci + 1 ) {
                my $clen = length( $l->[$ci] );
                $len[$ci] = $clen if ( $len[$ci] < $clen );
            }
        }
    }

    ## print out cols with fixed width
    my $col_spc = 2;

    open( my $out, ">$col_file.col_$$" );

    foreach my $l (@cont) {
        map { printf $out ( "%-*s", $len[$_] + $col_spc, $l->[$_] ) }
          0 .. @{$l} - 1;
        print $out "\n";
    }
    close($out);

    system("mv $col_file.col_$$ $col_file");
}
