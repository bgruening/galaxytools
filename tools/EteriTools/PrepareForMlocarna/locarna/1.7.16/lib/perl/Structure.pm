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

package locarnate::Structure;

use strict;
use warnings;

use Error;

###########################################################################
# intern functions
###########################################################################

my $_parse_structure = sub {
  my $self = shift(@_);
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  # check for unknown characters
  my $allowed_chars = $self->{'neutral'}.$self->{'open'}.$self->{'close'};
  my $pattern = '([^'.$allowed_chars.'])';
  if ($self->{'string'} =~ m/$pattern/) {
    my $char = $1;
    locarnate::Error::debug('Unkown structure character "'.$char.'" found');
  }

  # determine structure
  my @bp_stack = ();

  for (my $i = 0; $i < length($self->{'string'}); ++$i) {
    push(@{$self->{'begin_to_end'}}, undef);
    push(@{$self->{'end_to_begin'}}, undef);
    my $char = substr($self->{'string'}, $i, 1);
    if ($char eq $self->{'open'}) {
      push(@bp_stack, $i);
    }
    elsif ($char eq $self->{'close'}) {
      if (!@bp_stack) {
        locarnate::Error::debug('Unmatched structure closing at pos '.$i);
      }
      my $start = pop(@bp_stack);
      $self->{'begin_to_end'}->[$start] = $i;
      $self->{'end_to_begin'}->[$i] = $start;
      push(@{$self->{'begin_indices'}}, $start);
    }
  }
  if (@bp_stack) {
    locarnate::Error::debug('Unmatched structure opening at pos '.$bp_stack[0]);
  }
};

###########################################################################
# constructor
###########################################################################

sub new {
  my $class = shift(@_);
  $class = ref($class) || $class;

  if (@_ < 1 || @_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $string = $_[0];
  my $syntax = '.()';
  if (@_ > 1) {
    $syntax = $_[1];
  }

  my @chars = split(//,$syntax);
  if (@chars != 3) {
  }

  my $self = {
    'string'  => $string,
    'neutral' => $chars[0],
    'open'    => $chars[1],
    'close'   => $chars[2],
    'begin_to_end' => [],
    'end_to_begin' => [],
    'begin_indices' => [],
  };

  # \$self
  # {'string'}         $string;
  # {'neutral'}        $neutral_char,
  # {'open'}           $open_char,
  # {'close'}          $close_char,
  # {'begin_to_end'}   \@begin_to_end
  #                    []  $end_pos
  # {'end_to_begin'}   \@end_to_begin
  #                    []  $begin_pos
  # {'begin_indices'}  \@begin_indices
  #                    []  $begin_index

  &{$_parse_structure}($self);
  
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

  # string
  $string .= $ind.'string: '.$self->string()."\n";

  # basepairs
  $string .= $ind.'basepairs:'."\n";
  for (my $i = 0; $i < $self->size(); ++$i) {
    $string .= $ind.'  basepair '.$i.': '
        .$self->begin_of($i).', '.$self->end_of($i)."\n";
  }

  return ($string);
}

sub string {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }
  return($self->{'string'});
}

sub size {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  return(scalar(@{$self->{'begin_indices'}}));
}

sub is_begin {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $begin = $_[0];

  # check boarders
  if ($begin >= @{$self->{'begin_to_end'}}) {
    locarnate::Error::debug('Begin is bigger then structure length');
  }

  return(defined($self->{'begin_to_end'}->[$begin]));
}

sub is_end {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $end = $_[0];

  # check boarders
  if ($end >= @{$self->{'end_to_begin'}}) {
    locarnate::Error::debug('End is bigger then structure length');
  }

  return(defined($self->{'end_to_begin'}->[$end]));
}

sub is_bp {
  my $self = shift(@_);
  
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $begin = $_[0];
  my $end = $_[1];

  # check boarders
  if ($begin >= @{$self->{'begin_to_end'}}) {
    locarnate::Error::debug('Begin is bigger then structure length');
  }
  elsif ($end >= @{$self->{'end_to_begin'}}) {
    locarnate::Error::debug('End is bigger then structure length');
  }
  elsif (!defined($self->{'begin_to_end'}->[$begin])) {
    locarnate::Error::debug('Basepair does not exist');
  }
  elsif (!defined($self->{'end_to_begin'}->[$end])) {
    locarnate::Error::debug('Basepair does not exist');
  }

  return($self->{'begin_to_end'}->[$begin] == $end);
}

sub begin_of {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $index = $_[0];

  # check boarders
  if ($index >= @{$self->{'begin_indices'}}) {
    locarnate::Error::debug('Index is higher then basepair number');
  }

  return($self->{'begin_indices'}->[$index]);
}

sub end_of {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $index = $_[0];

  # check boarders
  if ($index >= @{$self->{'begin_indices'}}) {
    locarnate::Error::debug('Index is higher then basepair number');
  }

  return($self->{'begin_to_end'}->[$self->{'begin_indices'}->[$index]]);
}

sub begin_to_end {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $begin = $_[0];

  # check boarders
  if ($begin >= @{$self->{'begin_to_end'}}) {
    locarnate::Error::debug('Pos is higher then structure length');
  }

  return($self->{'begin_to_end'}->[$begin]);
}

sub end_to_begin {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $end = $_[0];

  # check boarders
  if ($end >= @{$self->{'end_to_begin'}}) {
    locarnate::Error::debug('Pos is higher then structure length');
  }

  return($self->{'end_to_begin'}->[$end]);
}

sub add_bp {
  my $self = shift(@_);
  
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $begin = $_[0];
  my $end = $_[1];

  # check boarders
  if ($begin >= $end) {
    locarnate::Error::debug('Begin must be smaller then end');
  }
  elsif ($end >= length($self->{'string'})) {
    locarnate::Error::debug('End is bigger then string size');
  }
  elsif (defined($self->{'begin_to_end'}->[$begin])) {
    locarnate::Error::debug('Begin pos is allready involved in basepair');
  }
  elsif (defined($self->{'end_to_begin'}->[$begin])) {
    locarnate::Error::debug('Begin pos is allready involved in basepair');
  }
  elsif (defined($self->{'begin_to_end'}->[$end])) {
    locarnate::Error::debug('End pos is allready involved in basepair');
  }
  elsif (defined($self->{'end_to_begin'}->[$end])) {
    locarnate::Error::debug('End pos is allready involved in basepair');
  }

  $self->{'string'} = substr($self->{'string'}, 0, $begin).
      $self->{'open'}.
      substr($self->{'string'}, $begin + 1, $end - $begin - 1).
      $self->{'close'}.
      substr($self->{'string'}, $end + 1);
  for (my $i = $begin + 1; $i < $end; ++$i) {
    if (defined($self->{'begin_to_end'}->[$i]) 
        && $self->{'begin_to_end'}->[$i] > $end) {
      locarnate::Error::debug('Structure becomes unnested');
    }
  }
  
  $self->{'begin_to_end'}->[$begin] = $end;
  $self->{'end_to_begin'}->[$end] = $begin;

  my $pos = $self->size();
  for (my $i = 0; $i < $self->size(); ++$i) {
    if ($self->{'begin_to_end'}->[$self->{'begin_indices'}->[$i]] >
        $end) {
      $pos = $i;
      last;
    }
  }
  splice(@{$self->{'begin_indices'}}, $pos, 0, $begin);
}

sub delete_bp {
  my $self = shift(@_);
  
  if (@_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $begin = $_[0];
  my $end = $_[1];

  # check boarders
  if ($begin >= @{$self->{'begin_to_end'}}) {
    locarnate::Error::debug('Begin is bigger then structure length');
  }
  elsif ($end >= @{$self->{'end_to_begin'}}) {
    locarnate::Error::debug('End is bigger then structure length');
  }
  elsif (!defined($self->{'begin_to_end'}->[$begin])) {
    locarnate::Error::debug('Basepair does not exist');
  }
  elsif (!defined($self->{'end_to_begin'}->[$end])) {
    locarnate::Error::debug('Basepair does not exist');
  }

  $self->{'string'} = substr($self->{'string'}, 0, $begin).
      $self->{'neutral'}.
      substr($self->{'string'}, $begin + 1, $end - $begin - 1).
      $self->{'neutral'}.
      substr($self->{'string'}, $end + 1);
  
  $self->{'begin_to_end'}->[$begin] = undef;
  $self->{'end_to_begin'}->[$end] = undef;

  my $pos = undef;
  for (my $i = 0; $i < $self->size(); ++$i) {
    if ($self->{'begin_indices'}->[$i] == $begin) {
      $pos = $i;
      last;
    }
  }
  splice(@{$self->{'begin_indices'}}, $pos, 1);
}

1;

__END__
