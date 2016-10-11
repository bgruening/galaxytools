###########################################################################
# Copyright (C) 2009 by Wolfgang Otto                      
# All Rights Reserved.                                     
#  
# Permission to use, copy, modify, and distribute this
# software and its documentation for NON-COMMERCIAL purposes
# and without fee is hereby granted provided that this     
# copyright notice appears in all copies.                     
#                                                             
# THE AUTHOR AND PUBLISHE MAKE NO REPRESENTATIONS OR          
# WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER    
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE        
# IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A        
# PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE AUTHORS        
# AND PUBLISHER SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED  
# BY LICENSEE AS A RESULT OF USING, MODIFYING OR DISTRIBUTING 
# THIS SOFTWARE OR ITS DERIVATIVES.                           
#                                                             
###########################################################################

package locarnate::Alignment;

use strict;
use warnings;

use Error;
use Structure;

###########################################################################
# intern functions
###########################################################################

my $_append = sub {
  my $self = shift(@_);
  if (@_ != 3) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_ref = $_[0];
  my $rep_tag = $_[1];
  my $string = $_[2];

  if (!defined($string)) {
    locarnate::Error::debug('Undefined string');
  }

  if (defined($rep_tag)) {
    if ($rep_tag =~ m/^\d+$/) {
      locarnate::Error::debug('Numbers can\'t be used as representation tags');
    }
    else {
      if (exists($seq_ref->[2]->{$rep_tag})) {
        # representation already exists
        my $rep_index = $seq_ref->[2]->{$rep_tag};
        $seq_ref->[1]->[$rep_index]->[1] .= $string;
        $seq_ref->[1]->[$rep_index]->[2] = undef;
        $seq_ref->[1]->[$rep_index]->[3] = undef;
      }
      else {
        # new representation with given tag
        $seq_ref->[2]->{$rep_tag} = scalar(@{$seq_ref->[1]});
        push(@{$seq_ref->[1]}, [$rep_tag, $string, undef, undef]);
      }
    }
  }
  else {
    # new representation without given tag
    push(@{$seq_ref->[1]}, [undef, $string, undef, undef]);
  }
};


my $_process_representation = sub {
  my $self = shift(@_);
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_ref = $_[0];
  my $representation = $_[1];

  if (!ref($representation)) {
    # non-tagged representation
    my $rep_tag = undef;
    my $string = $representation;
    &{$_append}($self, $seq_ref, $rep_tag, $string);
  }
  elsif (ref($representation) eq 'HASH') {
    # tagged representation
    # check for correct hash size
    if (keys(%{$representation}) != 1) {
      locarnate::Error::debug('Wrong representation syntax, '
                        .'hash must have exactly one element');
    }
    foreach my $key (keys(%{$representation})) {
      my $rep_tag = $key;
      my $string = $representation->{$key};
      &{$_append}($self, $seq_ref, $rep_tag, $string);
    }
  }
  else {
    locarnate::Error::debug('Wrong representation syntax: '
                      .'representation have to be scalar or hash');
  }
};

my $_sequence_indices = sub {
  my $self = shift(@_);
  
  if (@_ > 3) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_index = undef;
  if (@_ > 0) {
    $seq_index = $_[0];
  }
  my $rep_index = undef;
  if (@_ > 1) {
    $rep_index = $_[1];
  }
  my $str_index = undef;
  if (@_ > 2) {
    $str_index = $_[2];
  }

  # check existence of seq_index
  if (defined($seq_index) && defined(${$seq_index})) { 
    if (${$seq_index} =~ m/^\d+$/) {
      if (${$seq_index} >= @{$self->{'sequences'}->[0]}) {
        locarnate::Error::debug('Sequence index '.${$seq_index}
                          .' is higher then sequence number '
                          .@{$self->{'sequences'}->[0]});
      }
    }
    else {
      if (!exists($self->{'sequences'}->[1]->{${$seq_index}})) {
        locarnate::Error::debug('Sequence tag "'.${$seq_index}.'" does not exists');
      }
      else {
        ${$seq_index} = $self->{'sequences'}->[1]->{${$seq_index}};
      }
    }
  }

  # check existence of rep_index
  if (defined($rep_index) && defined(${$rep_index})) { 
    if (${$rep_index} =~ m/^\d+$/) {
      if (${$rep_index} >= 
          @{$self->{'sequences'}->[0]->[${$seq_index}]->[1]}) {
        locarnate::Error::debug('Representation index '.${$rep_index}
                          .' is higher then representation number '
                          .@{$self->{'sequences'}->[0]->[${$seq_index}]->[1]});
      }
    }
    else {
      if (!exists($self->{'sequences'}->[0]
                  ->[${$seq_index}]->[2]->{${$rep_index}})) {
        locarnate::Error::debug('Representation tag "'.${$rep_index}.
                          '" does not exists');
      }
      else {
        ${$rep_index} = $self->{'sequences'}->[0]
            ->[${$seq_index}]->[2]->{${$rep_index}};
      }
    }
  }

  # check existence of str_index
  if (defined($str_index) && defined(${$str_index})) { 
    if (ref($self->{'sequences'}->[0]->[${$seq_index}]
            ->[1]->[${$rep_index}]->[1])) {
      locarnate::Error::debug('Representation is no string');
    }
    elsif (${$str_index} >= length($self->{'sequences'}->[0]->
                                   [${$seq_index}]->[1]->
                                   [${$rep_index}]->[1])) {
      locarnate::Error::debug('String index '.${$str_index}
                        .' is higher then string size '
                        .length($self->{'sequences'}->[0]->[${$seq_index}]
                                ->[1]->[${$rep_index}]->[1]));
    }
  }
};


my $_consensus_indices = sub {
  my $self = shift(@_);
  
  if (@_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $rep_index = undef;
  if (@_ > 0) {
    $rep_index = $_[0];
  }
  my $str_index = undef;
  if (@_ > 1) {
    $str_index = $_[1];
  }

  # check existence of rep_index
  if (defined($rep_index) && defined(${$rep_index})) { 
    if (${$rep_index} =~ m/^\d+$/) {
      if (${$rep_index} >= @{$self->{'consensus'}->[1]}) {
        locarnate::Error::debug('Representation index '.${$rep_index}
                          .' is higher then representation number '
                          .@{$self->{'consensus'}->[1]});
      }
    }
    else {
      if (!exists($self->{'consensus'}->[2]->{${$rep_index}})) {
        locarnate::Error::debug('Representation tag "'.${$rep_index}.
                          '" does not exists');
      }
      else {
        ${$rep_index} = $self->{'consensus'}->[2]->{${$rep_index}};
      }
    }
  }

  # check existence of str_index
  if (defined($str_index) && defined(${$str_index})) { 
    if (ref($self->{'consensus'}->[1]->[${$rep_index}]->[1])) {
      locarnate::Error::debug('Representation is no string');
    }
    elsif (length($self->{'consensus'}->[1]->[${$rep_index}]->[1])) {
      locarnate::Error::debug('Str index '.${$str_index}
                        .' is higher then string size '
                        .length($self->{'consensus'}->
                                [1]->[${$rep_index}]->[1]));
    }
  }
};


###########################################################################
# constructor
###########################################################################

sub new {
  my $class = shift(@_);
  $class = ref($class) || $class;

  if (@_ != 0) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $self = {
    'columns'     => undef,
    'sequences'   => [[], {}],
    'consensus'   => [undef, [], {}],
    'score'       => ''};

  # \$self
  # {'columns'}      $columns;
  # {'sequences'}    \@data
  #                  [0]  \@sequences
  #                       []  \@sequence
  #                           [0]  $seq_tag
  #                           [1]  \@representations
  #                                []  \@representation
  #                                    [0]  $rep_tag
  #                                    [1]  $string
  #                                    [2]  \@align_to_seq
  #                                         []  $seq
  #                                    [3]  \@seq_to_align
  #                                         []  $pos
  #                           [2]  \%rep_indices
  #                  [1]  \%seq_indices
  # {'consensus'}    \@consensus
  #                  [0]  undef
  #                  [1]  \@representations
  #                       []  \@representation
  #                           [0]  $rep_tag
  #                           [1]  $string
  #                           [2]  \@align_to_seq
  #                                []  $seq
  #                           [3]  \@seq_to_align
  #                                []  $pos
  #                  [2]  \%rep_indices
  # {'score'}        $score;

  bless($self, $class);

  return($self);
}


###########################################################################
# methods
###########################################################################

sub str {
  my $self = shift(@_);

  if (@_ > 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $indent = 0;
  if (@_ == 1) {
    $indent = shift(@_);
  }
  my $ind = " " x $indent;
  my $string = '';

  # columns
  $string .= $ind.'columns: '.$self->columns()."\n";
  # score
  $string .= $ind.'score:   '.$self->score()."\n";

  $string .= $ind.'alignment:'."\n";

  # max_tag
  my $max_tag = 0;
  for (my $i = 0; $i < $self->sequence_size(); ++$i) {
    for (my $j = 0; $j < $self->sequence_size($i); ++$j) {
      my $tag = length($self->sequence_tag($i).'.'.
                       $self->sequence_tag($i, $j));
      if ($tag > $max_tag) {
        $max_tag = $tag;
      }
    }
  }
  for (my $i = 0; $i < $self->consensus_size(); ++$i) {
    my $tag = length($self->consensus_tag($i));
    if ($tag > $max_tag) {
      $max_tag = $tag;
    }
  }
  my $format_str = '%-'.($max_tag + 3).'s';

  # sequences
  for (my $i = 0; $i < $self->sequence_size(); ++$i) {
    for (my $j = 0; $j < $self->sequence_size($i); ++$j) {
      my $tag = sprintf($format_str, ($self->sequence_tag($i).'.'.
                                      $self->sequence_tag($i, $j)));
      $string .= $ind.'  '.$tag.$self->sequence($i, $j)."\n";
    }
  }
  # consensus
  $string .= "\n";
  for (my $i = 0; $i < $self->consensus_size(); ++$i) {
    my $tag = sprintf($format_str, $self->consensus_tag($i));
    $string .= $ind.'  '.$tag.$self->consensus($i)."\n";
  }

  return ($string);
}

sub append_sequence {
  my $self = shift(@_);
  
  if (@_ < 1 || @_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $sequence = $_[0];
  my $seq_tag = undef;
  if (@_ > 1) {
    $seq_tag = $_[1];
  }

  # get reference for sequence depending on sequence tag
  my $seq_ref = undef;
  if (defined($seq_tag)) {
    if ($seq_tag =~ m/^\d+$/) {
      locarnate::Error::debug('Numbers can\'t be used as sequence tags');
    }
    else {
      if (exists($self->{'sequences'}->[1]->{$seq_tag})) {
        # sequence already exists
        my $seq_index = $self->{'sequences'}->[1]->{$seq_tag};
        $seq_ref = $self->{'sequences'}->[0]->[$seq_index];
      }
      else {
        # new sequence with given tag
        $self->{'sequences'}->[1]->{$seq_tag} =
            scalar(@{$self->{'sequences'}->[0]});
        my $seq_index = @{$self->{'sequences'}->[0]};
        push(@{$self->{'sequences'}->[0]}, [$seq_tag, [], {}]);
        $seq_ref = $self->{'sequences'}->[0]->[$seq_index];
      }
    }
  }
  else {
    # new sequence without given tag
    my $seq_index = @{$self->{'sequences'}->[0]};
    push(@{$self->{'sequences'}->[0]}, [undef, [], {}]);
    $seq_ref = $self->{'sequences'}->[0]->[$seq_index];
  }

  # process representations
  if (ref($sequence) ne 'ARRAY') {
    # sequence with one representation
    my $representation = $sequence;
    &{$_process_representation}($self, $seq_ref, $representation);
  }
  else {
    # sequence with multiple representation
    for (my $i = 0; $i < @{$sequence}; ++$i) {
      my $representation = $sequence->[$i];
      &{$_process_representation}($self, $seq_ref, $representation);
    }
  }
}

sub set_consensus {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $consensus = $_[0];

  my $seq_ref = $self->{'consensus'};

  # process representations
  if (ref($consensus) ne 'ARRAY') {
    # consensus with one representation
    my $representation = $consensus;
    &{$_process_representation}($self, $seq_ref, $representation);
  }
  else {
    # consensus with multiple representation
    for (my $i = 0; $i < @{$consensus}; ++$i) {
      my $representation = $consensus->[$i];
      &{$_process_representation}($self, $seq_ref, $representation);
    }
  }
}

sub set_columns {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $columns = $_[0];
  $self->{'columns'} = $columns;
}

sub set_score {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $score = $_[0];
  $self->{'score'} = $score;
}

sub erase_sequence {
  my $self = shift(@_);
  
  if (@_ < 1 || @_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_index = $_[0];
  my $rep_index = undef;
  if (@_ > 1) {
    $rep_index = $_[1];
  }

  # check existence of indices
  &{$_sequence_indices}($self, \$seq_index, \$rep_index);

  my $seq_tag = $self->{'sequences'}->[0]->[$seq_index]->[0];
  my $rep_tag = undef;
  if (defined($rep_index)) { 
    $rep_tag = $self->{'sequences'}->[0]->[$seq_index]->
        [1]->[$rep_index]->[0];
  }

  if (defined($rep_index)) {
    splice(@{$self->{'sequences'}->[0]->[$seq_index]->[1]},
           $rep_index, 1);
    foreach my $key (keys(%{$self->{'sequences'}
                            ->[0]->[$seq_index]->[2]})) {
      if ($self->{'sequences'}->[0]->[$seq_index]->[2]->{$key} > 
          $rep_index) {
        --$self->{'sequences'}->[0]->[$seq_index]->[2]->{$key};
      }
    }
    if (defined($rep_tag)) {
      delete($self->{'sequences'}->[0]->[$seq_index]
             ->[2]->{$rep_tag});
    }
  }
  if (!defined($rep_index) || !$self->sequence_size($seq_index)) {
    splice(@{$self->{'sequences'}->[0]}, $seq_index, 1);
    foreach my $key (keys(%{$self->{'sequences'}->[1]})) {
      if ($self->{'sequences'}->[1]->{$key} > $seq_index) {
        --$self->{'sequences'}->[1]->{$key};
      }
    }
    if (defined($seq_tag)) {
      delete($self->{'sequences'}->[1]->{$seq_tag});
    }
  }
}

sub erase_consensus {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $rep_index = $_[0];

  # check existence of indices
  &{$_sequence_indices}($self, \$rep_index);

  my $rep_tag = $self->{'consensus'}->[1]->[$rep_index]->[0];

  splice(@{$self->{'consensus'}->[1]}, $rep_index, 1);
  foreach my $key (keys(%{$self->{'consensus'}->[2]})) {
    if ($self->{'consensus'}->[2]->{$key} > $rep_index) {
      --$self->{'consensus'}->[2]->{$key};
    }
  }
  if (defined($rep_tag)) {
    delete($self->{'consensus'}->[2]->{$rep_tag});
  }
}

sub columns {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  if (defined($self->{'columns'})) {
    return($self->{'columns'});
  }
  else {
    return(0);
  }
}

sub sequence_size {
  my $self = shift(@_);
  
  if (@_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_index = undef;
  if (@_ > 0) {
    $seq_index = $_[0];
  }
  my $rep_index = undef;
  if (@_ > 1) {
    $rep_index = $_[1];
  }

  # check existence of indices
  &{$_sequence_indices}($self, \$seq_index, \$rep_index);

  if (!defined($seq_index)) {
    return(scalar(@{$self->{'sequences'}->[0]}));
  }
  elsif (!defined($rep_index)) {
    return(scalar(@{$self->{'sequences'}->[0]->[$seq_index]->[1]}));
  }
  else {
    return(length($self->{'sequences'}->[0]->[$seq_index]->
                  [1]->[$rep_index]->[1]));
  }
}

sub consensus_size {
  my $self = shift(@_);

  if (@_ > 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $rep_index = undef;
  if (@_ > 0) {
    $rep_index = $_[0];
  }

  # check existence of indices
  &{$_consensus_indices}($self, \$rep_index);
  
  if (!defined($rep_index)) {
    return(scalar(@{$self->{'consensus'}->[1]}));
  }
  else {
    return(length($self->{'consensus'}->[1]->[$rep_index]->[1]));
  }
}

sub sequence_tag {
  my $self = shift(@_);
  
  if (@_ < 1 || @_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_index = $_[0];
  my $rep_index = undef;
  if (@_ > 1) {
    $rep_index = $_[1];
  }

  # check existence of indices
  &{$_sequence_indices}($self, \$seq_index, \$rep_index);

  if (defined($rep_index)) {
    my $tag = $self->{'sequences'}->[0]->[$seq_index]
        ->[1]->[$rep_index]->[0];
    if (defined($tag)) {
      return($tag);
    }
    else {
      return($rep_index);
    }
  }
  else {
    my $tag = $self->{'sequences'}->[0]->[$seq_index]->[0];
    if (defined($tag)) {
      return($tag);
    }
    else {
      return($seq_index);
    }
  }
}

sub exists_sequence_tag {
  my $self = shift(@_);
  
  if (@_ < 1 || @_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_index = $_[0];
  my $rep_index = undef;
  if (@_ > 1) {
    $rep_index = $_[1];
  }

  # check existence of seq_index
  if ($seq_index =~ m/^\d+$/) {
    if ($seq_index >= @{$self->{'sequences'}->[0]}) {
      return(0);
    }
  }
  else {
    if (!exists($self->{'sequences'}->[1]->{$seq_index})) {
      return(0);
    }
    else {
      $seq_index = $self->{'sequences'}->[1]->{$seq_index};
    }
  }

  # check existence of rep_index
  if (defined($rep_index)) { 
    if ($rep_index =~ m/^\d+$/) {
      if ($rep_index >= @{$self->{'sequences'}->[0]->[$seq_index]->[1]}) {
        return(0);
      }
    }
    else {
      if (!exists($self->{'sequences'}->[0]
                  ->[$seq_index]->[2]->{$rep_index})) {
        return(0);
      }
      else {
        return(1);
      }
    }
  }
  else {
    return(1);
  }
}

sub consensus_tag {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $rep_index = $_[0];

  # check existence of indices
  &{$_consensus_indices}($self, \$rep_index);

  my $tag = $self->{'consensus'}->[1]->[$rep_index]->[0];
  if (defined($tag)) {
    return($tag);
  }
  else {
    return($rep_index);
  }
}

sub exists_consensus_tag {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $rep_index = $_[0];

  # check existence of rep_index
  if ($rep_index =~ m/^\d+$/) {
    if ($rep_index >= @{$self->{'consensus'}->[1]}) {
      return(0);
    }
  }
  else {
    if (!exists($self->{'consensus'}->[2]->{$rep_index})) {
      return(0);
    }
    else {
      return(1);
    }
  }
}

sub sequence {
  my $self = shift(@_);
  
  if (@_ < 2 || @_ > 3) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_index = $_[0];
  my $rep_index = $_[1];
  my $str_index = undef;
  if (@_ > 2) {
    $str_index = $_[2];
  }

  # check existence of indices
  &{$_sequence_indices}($self, \$seq_index, \$rep_index, \$str_index);

  if (defined($str_index)) { 
    return(substr($self->{'sequences'}->[0]->[$seq_index]
                  ->[1]->[$rep_index]->[1], $str_index, 1));
  }
  else {
    return($self->{'sequences'}->[0]->[$seq_index]
           ->[1]->[$rep_index]->[1]);
  }
}

sub consensus {
  my $self = shift(@_);
  
  if (@_ < 1 || @_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $rep_index = $_[0];
  my $str_index = undef;
  if (@_ > 1) {
    $str_index = $_[1];
  }

  # check existence of indices
  &{$_consensus_indices}($self, \$rep_index, \$str_index);

  if (defined($str_index)) { 
    return(substr($self->{'consensus'}->[1]->[$rep_index]->[1], 
                  $str_index, 1));
  }
  else {
    return($self->{'consensus'}->[1]->[$rep_index]->[1]);
  }
}

sub score {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  return($self->{'score'});
}

sub check {
  my $self = shift(@_);

  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_tags = $_[0];
  my $cons_tags = $_[1];

  # check sequences
  for (my $i = 0; $i < $self->sequence_size(); ++$i) {
    foreach my $tag (@{$seq_tags}) {
      my $length = $self->sequence_size($i, $tag);
      if (!defined($self->{'columns'})) {
        $self->{'columns'} = $length;
      }
      elsif ($self->columns() != $length) {
        locarnate::Error::debug('Representation "'.$tag.'" of sequence "'.$i.
                          '" has wrong length '.$length.
                          ' (alignment length is '.$self->columns().')');
      }
    }
  }
  # check consensus
  foreach my $tag (@{$cons_tags}) {
    my $length = $self->consensus_size($tag);
    if (!defined($self->{'columns'})) {
      $self->{'columns'} = $length;
    }
    elsif ($self->columns() != $length) {
      locarnate::Error::debug('Consensus representation "'.$tag.
                        '" has wrong length '.$length.
                        ' (alignment length is '.$self->columns().')');
    }
  }
}

sub sequence_align_to_seq {
  my $self = shift(@_);

  if (@_ != 3) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_tag = $_[0];
  my $rep_tag = $_[1];
  my $align_pos = $_[2];

  &{$_sequence_indices}($self, \$seq_tag, \$rep_tag);

  if (!defined($self->{'sequences'}->[0]->[$seq_tag]
               ->[1]->[$rep_tag]->[2])) {
    my @align_to_seq = ();
    my @seq_to_align = ();
    my $seq_index = 0;
    for (my $i = 0; $i < $self->sequence_size($seq_tag, $rep_tag); ++$i) {
      my $char = $self->sequence($seq_tag, $rep_tag, $i);
      if ($char =~ /^[A-Za-z]$/) {
        push(@align_to_seq, $seq_index);
        push(@seq_to_align, $i);
        ++$seq_index;
      }
      else {
        push(@align_to_seq, -1);
      }
    }
    $self->{'sequences'}->[0]->[$seq_tag]
        ->[1]->[$rep_tag]->[2] = \@align_to_seq;
    $self->{'sequences'}->[0]->[$seq_tag]
        ->[1]->[$rep_tag]->[3] = \@seq_to_align;
  }

  # check pos
  if ($align_pos > @{$self->{'sequences'}->[0]->[$seq_tag]
                         ->[1]->[$rep_tag]->[2]}) {
    locarnate::Error::debug('Alignment index '.$align_pos.
                      ' is higher then alignment size '.
                      @{$self->{'sequences'}->[0]->[$seq_tag]
                            ->[1]->[$rep_tag]->[2]});
  }
  return($self->{'sequences'}->[0]->[$seq_tag]
         ->[1]->[$rep_tag]->[2]->[$align_pos]);
}

sub sequence_seq_to_align {
  my $self = shift(@_);

  if (@_ != 3) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_tag = $_[0];
  my $rep_tag = $_[1];
  my $seq_pos = $_[2];

  &{$_sequence_indices}($self, \$seq_tag, \$rep_tag);

  if (!defined($self->{'sequences'}->[0]->[$seq_tag]
               ->[1]->[$rep_tag]->[3])) {
    my @align_to_seq = ();
    my @seq_to_align = ();
    my $seq_index = 0;
    for (my $i = 0; $i < $self->sequence_size($seq_tag, $rep_tag); ++$i) {
      my $char = $self->sequence($seq_tag, $rep_tag, $i);
      if ($char =~ /^[A-Za-z]$/) {
        push(@align_to_seq, $seq_index);
        push(@seq_to_align, $i);
        ++$seq_index;
      }
      else {
        push(@align_to_seq, -1);
      }
    }
    $self->{'sequences'}->[0]->[$seq_tag]
        ->[1]->[$rep_tag]->[2] = \@align_to_seq;
    $self->{'sequences'}->[0]->[$seq_tag]
        ->[1]->[$rep_tag]->[3] = \@seq_to_align;
  }

  # check pos
  if ($seq_pos > @{$self->{'sequences'}->[0]->[$seq_tag]
                         ->[1]->[$rep_tag]->[3]}) {
    locarnate::Error::debug('Alignment index '.$seq_pos.
                      ' is higher then alignment size '.
                      @{$self->{'sequences'}->[0]->[$seq_tag]
                            ->[1]->[$rep_tag]->[3]});
  }

  return($self->{'sequences'}->[0]->[$seq_tag]
         ->[1]->[$rep_tag]->[3]->[$seq_pos]);
}

sub erase_empty_columns {
  my $self = shift(@_);
  
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $seq_tags = $_[0];
  my $cons_tags = $_[1];

  # check tags
  my @seq_tag_list = keys(%{$seq_tags});
  my @cons_tag_list = keys(%{$cons_tags});
  $self->check(\@seq_tag_list, \@cons_tag_list);
  foreach my $tag (keys(%{$seq_tags})) {
    if (length($seq_tags->{$tag}) != 1 
        && length($seq_tags->{$tag}) != 3) {
      locarnate::Error::debug('Tag description have to be single gap char for'.
                        ' sequences or gap, open and close char for'.
                        ' structures');
    }
  }
  foreach my $tag (keys(%{$cons_tags})) {
    if (length($cons_tags->{$tag}) != 1 
        && length($cons_tags->{$tag}) != 3) {
      locarnate::Error::debug('Tag description have to be single gap char for'.
                        ' sequences or gap, open and close char for'.
                        ' structures');
    }
  }

  # determine empty columns
  my @empty_columns = ();
  for (my $i = 0; $i < $self->columns(); ++$i) {
    # *STDERR->print("\r".'determine empty columns: '.$i.'/'.
    #                $self->columns());
    my $is_empty = 1;
    for (my $j = 0; $j < $self->sequence_size(); ++$j) {
      foreach my $tag (keys(%{$seq_tags})) {
        if (length($seq_tags->{$tag}) == 1) {
          my $gap_char = $seq_tags->{$tag};
          if ($self->sequence($j, $tag, $i) ne $gap_char) {
            $is_empty = 0;
            last;
          }
        } 
      }
      if (!$is_empty) {
        last;
      }
    }
    if ($is_empty) {
      push(@empty_columns, $i);
    }
  }
  # *STDERR->print("\r".(' 'x75)."\r");


  # remove ghost structures
  for (my $i = 0; $i < $self->sequence_size(); ++$i) {
    foreach my $tag (keys(%{$seq_tags})) {
      my $syntax = $seq_tags->{$tag};
      if (length($seq_tags->{$tag}) == 3) { 
        &{$_sequence_indices}($self, \$i, \$tag);
        my $structure = new locarnate::Structure($self->{'sequences'}->[0]->
                                           [$i]->[1]->[$tag]->[1],
                                           $syntax);
        for (my $j = 0; $j < @empty_columns; ++$j) {
          my $pos = $empty_columns[$j];
          if (defined($structure->begin_to_end($pos))) {
            my $end = $structure->begin_to_end($pos);
            $structure->delete_bp($pos, $end);
          }
          elsif (defined($structure->end_to_begin($pos))) {
            my $begin = $structure->end_to_begin($pos);
            $structure->delete_bp($begin, $pos);
          }
        }
        $self->{'sequences'}->[0]->[$i]->[1]->[$tag]->[1] = 
            $structure->string();
      }
    }
  }
  foreach my $tag (keys(%{$cons_tags})) {
    my $syntax = $cons_tags->{$tag};
    if (length($cons_tags->{$tag}) == 3) { 
      &{$_consensus_indices}($self, \$tag);
      my $structure = new locarnate::Structure($self->{'consensus'}->[1]->
                                         [$tag]->[1], $syntax);
      for (my $j = 0; $j < @empty_columns; ++$j) {
        my $pos = $empty_columns[$j];
        if (defined($structure->begin_to_end($pos))) {
          my $end = $structure->begin_to_end($pos);
          $structure->delete_bp($pos, $end);
        }
        elsif (defined($structure->end_to_begin($pos))) {
          my $begin = $structure->end_to_begin($pos);
          $structure->delete_bp($begin, $pos);
        }
      }
      $self->{'consensus'}->[1]->[$tag]->[1] = $structure->string();
    }
  }

  # remove empty columns
  for (my $i = @empty_columns - 1; $i >= 0; --$i) {
    my $pos = $empty_columns[$i];

    # seq
    for (my $i = 0; $i < $self->sequence_size(); ++$i) {
      foreach my $tag (keys(%{$seq_tags})) {
        &{$_sequence_indices}($self, \$i, \$tag);
        my $string = $self->{'sequences'}->[0]->
            [$i]->[1]->[$tag]->[1];
        $self->{'sequences'}->[0]->[$i]->[1]->[$tag]->[1] = 
            substr($string, 0, $pos).substr($string, $pos + 1);
      }
    }
    # cons
    foreach my $tag (keys(%{$cons_tags})) {
      &{$_consensus_indices}($self, \$tag);
      my $string = $self->{'consensus'}->[1]->[$tag]->[1];
      $self->{'consensus'}->[1]->[$tag]->[1] =
          substr($string, 0, $pos).substr($string, $pos + 1);
    }
  }
  $self->{'columns'} -= @empty_columns;
}

1;

__END__
