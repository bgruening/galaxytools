package Part;

use strict;
use warnings;

use Error;

###########################################################################
# intern functions
###########################################################################

my $_set_indices = sub {
  my $self = shift(@_);

  if (@_) {
    Error::error('Wrong argument number');
  }

  my $seq = $self->{'seq'};

  my @index_seq_to_align = ();
  my @index_align_to_seq = ();

  my $seq_counter = 0;
  for (my $i = 0; $i < length($seq); ++$i) {
    my $char = substr($seq, $i, 1);
    if ($char =~ m/^[ACGURYSWKMBDHVN]$/) {
      # seq char
      push(@index_seq_to_align, $i);
      push(@index_align_to_seq, $seq_counter);
      ++ $seq_counter;
    }
    elsif ($char =~ m/^[-]$/) {
      # indel
      push(@index_align_to_seq, -1);
    }
    elsif ($char =~ m/^[~]$/) {
      # indel
      push(@index_align_to_seq, -2);
    }
    else {
      Error::error('Unkown alignment character "'.$char.'" found');
    }
  }

  $self->{'index_seq_to_align'} = \@index_seq_to_align;
  $self->{'index_align_to_seq'} = \@index_align_to_seq;
};

my $_set_bps = sub {
  my $self = shift(@_);

  if (@_) {
    Error::error('Wrong argument number');
  }

  my $struc = $self->{'struc'};

  my @bps_begin_to_end = ();
  my @bps_end_to_begin = ();
  my @bps = ();

  my @bps_stack = ();
  for (my $i = 0; $i < length($struc); ++$i) {
    my $index = $self->{'index_align_to_seq'}->[$i];
    my $char = substr($struc, $i, 1);
    if ($char eq '(') {
      my $begin = $index;

      my $size = @bps;
      push(@bps_stack, $size);

      my %bp = ();
      $bp{'begin'} = $index;
      push(@bps, \%bp);

      push(@bps_begin_to_end, $begin);
      push(@bps_end_to_begin, -1);
    }
    elsif ($char eq ')') {
      my $end = $index;

      if (!@bps_stack) {
        Error::error('No nested structure (missing "(")');
      }
      my $bps_index = pop(@bps_stack);
      $bps[$bps_index]->{'end'} = $end;
      my $begin = $bps[$bps_index]->{'begin'};

      push(@bps_begin_to_end, -1);
      push(@bps_end_to_begin, $begin);
      $bps_begin_to_end[$begin] = $end;
    }
    else {
      if ($index >= 0) {
        push(@bps_begin_to_end, -1);
        push(@bps_end_to_begin, -1);
      }
    }
  }
  if (@bps_stack) {
    Error::error('No nested structure (missing ")")');
  }
   
  $self->{'bps_begin_to_end'} = \@bps_begin_to_end;
  $self->{'bps_end_to_begin'} = \@bps_end_to_begin;
  $self->{'bps'} = \@bps;
};

###########################################################################
# constructor
###########################################################################

sub new {
  my $class = shift(@_);
  $class = ref($class) || $class;

  if (@_ < 2 || @_ > 3) {
    Error::error('Wrong argument number');
  }

  my $id = $_[0];
  my $seq = $_[1];
  my $struc = undef;
  if (@_ == 3) {
    $struc = $_[2];
  }

  my $self = {
    'id' => $id,
    'seq' => $seq,
    'struc' => $struc,
    'index_seq_to_align' => undef,
    'index_align_to_seq' => undef,
    'bps_begin_to_end' => undef,
    'bps_end_to_begin' => undef,
    'bps' => undef
  };

  bless($self, $class);

  &{$_set_indices}($self);
  if (defined($struc)) {
    &{$_set_bps}($self);
  }

  return($self);
}


###########################################################################
# methods
###########################################################################

sub str {
  my $self = shift(@_);

  if (@_ > 1) {
    Error::error('Wrong argument number');
  }

  my $indent = 0;
  if (@_ == 1) {
    $indent = shift(@_);
  }
  my $ind = " " x $indent;
  my $string = '';

  $string .= $ind.'id: '.$self->{'id'}."\n";
  $string .= $ind.'sequence:'."\n";
  $string .= $ind.'  '.$self->{'seq'}."\n";
  if (defined($self->{'struc'})) {
    $string .= $ind.'structure:'."\n";
    $string .= $ind.'  '.$self->{'struc'}."\n";
  }

  return ($string);
}

sub id {
  my $self = shift(@_);
  
  if (@_) {
    Error::error('Wrong argument number');
  }

  return($self->{'id'});
}

sub seq_size {
  my $self = shift(@_);
  
  if (@_) {
    Error::error('Wrong argument number');
  }

  return(length($self->{'seq'}));
}

sub sequence {
  my $self = shift(@_);
  
  if (@_) {
    Error::error('Wrong argument number');
  }

  return($self->{'seq'});
}

sub struc_size {
  my $self = shift(@_);
  
  if (@_) {
    Error::error('Wrong argument number');
  }

  return(length($self->{'struc'}));
}

sub structure {
  my $self = shift(@_);
  
  if (@_) {
    Error::error('Wrong argument number');
  }

  return($self->{'struc'});
}

sub pos_align_to_seq {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  if ($index >= @{$self->{'index_align_to_seq'}}) {
    Error::debug('No alignment entry with index '.$index.' exists');
  }

  return($self->{'index_align_to_seq'}->[$index]);
}

sub pos_seq_to_align {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  if ($index >= @{$self->{'index_seq_to_align'}}) {
    Error::debug('No sequence entry with index '.$index.' exists');
  }

  return($self->{'index_seq_to_align'}->[$index]);
}

sub bps_size {
  my $self = shift(@_);
  
  if (@_) {
    Error::error('Wrong argument number');
  }

  return(@{$self->{'bps'}});
}

sub bp_begin {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  if ($index >= @{$self->{'bps'}}) {
    Error::error('No bp with index '.$index.'exists');
  }

  return($self->{'bps'}->[$index]->{'begin'});
}

sub bp_end {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  if ($index >= @{$self->{'bps'}}) {
    Error::error('No bp with index '.$index.'exists');
  }

  return($self->{'bps'}->[$index]->{'end'});
}

sub bp_begin_to_end {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  if ($index >= @{$self->{'bps_begin_to_end'}}) {
    Error::error('No bp pos with index '.$index.'exists');
  }

  return($self->{'bps_begin_to_end'}->[$index]);
}

sub bp_end_to_begin {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  if ($index >= @{$self->{'bps_end_to_begin'}}) {
    Error::error('No bp pos with index '.$index.'exists');
  }

  return($self->{'bps_end_to_begin'}->[$index]);
}

sub set_gap {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  # sequence
  my $seq = $self->{'seq'};
  $self->{'seq'} = substr($seq, 0 , $index).'-'.
      substr($seq, $index + 1);

  # structure
  my $struc = $self->{'struc'};
  if (defined($struc) && $struc ne '') {
    $self->{'struc'} = substr($struc, 0 , $index).'-'.
        substr($struc, $index + 1);
  }

  # index_align_to_seq
  $self->{'index_align_to_seq'}->[$index] = -1;
}

sub set_exclusion {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  # sequence
  my $seq = $self->{'seq'};
  $self->{'seq'} = substr($seq, 0 , $index).'~'.
      substr($seq, $index + 1);

  # structure
  my $struc = $self->{'struc'};
  if (defined($struc) && $struc ne '') {
    $self->{'struc'} = substr($struc, 0 , $index).'~'.
        substr($struc, $index + 1);
  }

  # index_align_to_seq
  $self->{'index_align_to_seq'}->[$index] = -2;
}

sub set_unkown {
  my $self = shift(@_);
  
  if (@_ != 1) {
    Error::error('Wrong argument number');
  }

  my $index = $_[0];

  # sequence
  my $seq = $self->{'seq'};
  $self->{'seq'} = substr($seq, 0 , $index).'='.
      substr($seq, $index + 1);

  # structure
  my $struc = $self->{'struc'};
  if (defined($struc) && $struc ne '') {
    $self->{'struc'} = substr($struc, 0 , $index).'='.
        substr($struc, $index + 1);
  }

  # index_align_to_seq
  $self->{'index_align_to_seq'}->[$index] = -3;
}

1;
