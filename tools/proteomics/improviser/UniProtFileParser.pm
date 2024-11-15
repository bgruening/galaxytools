#!/usr/bin/perl
package UniProtFileParser;
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Math::Round;
use lib '/home/videmp/lib/';
use Pod::Usage;


sub new {
    my $class = shift;
    my $self = {
    	_file=>shift,
        _id	=> shift,
        _len	=> shift,
        _rec_name	=> shift,
        _sub_name	=> shift,
        _seq	=> shift, 
        _mol_wt	=> shift,        
        _features=> [],
        _identifiedPeptides => [],
    };
    bless $self, $class;
    return $self;
}

#accessor method for seq
sub file {
	my( $self, $file ) = @_;
    $self->{_file} = $file if defined($file);
    return $self->{_file};
}

#accessor method for uniprot id
sub id {
	my( $self, $id) = @_;
	$self->{_id} = $id if defined($id);
    return $self->{_id};
}

#accessor method for length of sequence
sub len {
	my( $self ) = @_;
	return length($self->{_seq});
}

#accessor method for seq
sub seq {
	my( $self, $seq ) = @_;
    $self->{_seq} = $seq if defined($seq);
    return $self->{_seq};
}

#accessor method for mol_wt
sub mol_wt {
	my( $self, $mol_wt ) = @_;
	if(defined($mol_wt)){
		$self->{_mol_wt} = $mol_wt ;
	}
    return $self->{_mol_wt};
}

#accessor method for rec_name
sub rec_name {
	my( $self, $rec_name ) = @_;
	if(defined($rec_name)){
		$self->{_rec_name} = $rec_name ;
	}
    return $self->{_rec_name};
}

#accessor method for sub_name
sub sub_name {
	my( $self, $sub_name ) = @_;
	if(defined($sub_name)){
		$self->{_sub_name} = $sub_name ;
	}
    return $self->{_sub_name};
}

#accessor method for features
sub features{
    my ( $self, @features ) = @_;
    $self->{_features} = @features if (@features);
    return @{$self->{_features}};
}

#accessor method for identifiedPeptides
sub identifiedPeptides{
    my ( $self, @identifiedPeptides ) = @_;
    $self->{_identifiedPeptides} = @identifiedPeptides if(@identifiedPeptides);
    return @{$self->{_identifiedPeptides}};
}

##############################################

sub parseContent{
	my $check;
	my ($self, $file) = @_;
#	open IN, $file or die $!,"\n";
	my @feat_array = ();
	my @lines = split(/^/, $file);
	while(@lines){
		if($lines[0]=~ m/^DT/){
			$check = $lines[0];
		}
		if($lines[0] =~ m/^DE\s+RecName:\s+Full=(.+)/){
			if($check ne ""){
				$self->{_rec_name} = $1;
			}
		}
		if($lines[0] =~ m/^DE\s+SubName:\s+Full=(.+)/){
			$self->{_sub_name} = $1;
		}
		if($lines[0] =~ m/^FT\s+(\w+)\s+(\d+)\s+(\d+)\s+(.+)/){
			my %feat = ();
			my $feat_name = $1;
			my $feat_start = $2;
			my $feat_end = $3;
			my $feat_desc = $4;
			$feat{'name'} = $feat_name;
			$feat{'start'} = $feat_start;
			$feat{'end'} = $feat_end;
			$feat{'desc'} = $feat_desc;
			push (@feat_array, \%feat);
		}
		if($lines[0]=~ m/^SQ/){
			my @mol_arr = split(';',$lines[0]);
			my $mol_wt_pre = $mol_arr[1];
			$mol_wt_pre =~ m/\s+(\d+)\s+(.+)/;
			$self->{_mol_wt} = $1;
		}
		if($lines[0]=~ m/^SQ/){
			shift(@lines);
			my $sequence = join('',@lines);
			$sequence =~ s/\s|\/|\n//g;
			$self->{_seq} = $sequence;
		}
		if($lines[0]=~ m/^DT/){
			$check = $lines[0];
		}
		else{
			$check = "";
		}
		shift(@lines);
	}
	$self->{_features} =\@feat_array;
#	close IN;
	return;
}

1;
