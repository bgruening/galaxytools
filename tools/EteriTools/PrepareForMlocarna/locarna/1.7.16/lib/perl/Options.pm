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

package locarnate::Options;

use strict;
use warnings;

use FileHandle;

use Error;

###########################################################################
# definitions
###########################################################################

my $message_indent = 30;
my $message_columns = 75;
my $message_defaults = 1;

###########################################################################
# intern functions
###########################################################################

my $_find_index = sub {
  my $self = shift(@_);
  if (@_ < 3 || @_ > 4) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];
  my $name = $_[1];
  my $counter = $_[2];
  my $remote_file = '';
  if (@_ > 3) {
    $remote_file = $_[3];
  }

  my @names = ();
  my $definition_index = @{$definition};

  # search options that match name
  for (my $index = 0; $index < @{$definition}; ++$index) {
    my $long_name = $definition->[$index]->[0];
    my $short_name = $definition->[$index]->[1];

    if ($name eq $short_name || $name eq $long_name) {
      # exact match (there can be no other)
      @names = ();
      if ($name eq $short_name) {
        push(@names, $short_name);
      }
      else {
        push(@names, $long_name);
      }
      $definition_index = $index;
      last;
    }
    elsif ($name eq substr($long_name, 0, length($name))) {
      # name matchs prefix of option long name
      push(@names, $long_name);
      $definition_index = $index;
    }
    # (since short names consists of a single char they always 
    # match exact or they don't match so there is no need for check)
  }

  # test if argument exists exactly one time
  if (@names == 0) {
    # no definition for argument
    my $message = '';
    if ($remote_file ne '') {
      $message .= 'Error: Unkown option \''.$name.'\' in remote file \''
          .$remote_file.'\', line '.$counter.'!'."\n";
    }
    else {
      $message .= 'Error: Unkown named command line argument \''.
          $name.'\' at position '.$counter.'!'."\n";
    }
    $self->{'error_message'} .= $message;
  }
  elsif (@names > 1) {
    # ambiguous argument
    my $message = '';
    if ($remote_file ne '') {
      $message .= 'Error: Argument \''.$name
          .'\' in remote file \''.$remote_file.'\', line '
          .$counter.', is ambiguous (\''.$names[0];
    }
    else {
      $message .= 'Error: Argument \''.$name
          .'\' at position '.$counter.', is ambiguous (\''.$names[0];
    }
    for (my $i = 1; $i < @names - 1; ++$i) {
      $message .= '\', \''.$names[$i];
    }
    $message .= '\' and \''.$names[@names - 1].'\')!'."\n";
    $self->{'error_message'} .= $message;
  }
  return($definition_index);
};


my $_check_definition = sub {
  my $self = shift(@_);
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];

  my %used_names = ();

  my $unnamed_default = 0;
  my $unnamed_default_index = 0;
  for (my $index = 0; $index < @{$definition}; ++$index) {
    my $long_name = $definition->[$index]->[0];
    my $short_name = $definition->[$index]->[1];
    my $default_value = $definition->[$index]->[2];
    my $description = $definition->[$index]->[3];

    if ($long_name eq '' && $short_name eq '' && $description ne '') {
      # this definition is a header
      if ($default_value ne '') {
        # headers can't have a default value
        my $message = 'Syntax error in definition '.($index + 1) 
            .': headers can\'t have a default value';
        locarnate::Error::error($message);
      }
    }
    elsif ($long_name ne '' || $short_name ne '') {
      # this definition is an option
      if ($short_name eq 'NN') {
        # this definition is an unnamed option
        if ($long_name eq '') {
          if ($default_value ne '' || $description eq '') {
            # unnamed options need a long name as key
            my $message = 'Syntax error in definition '.($index + 1) 
                .': unnamed options need a long name as key';
            locarnate::Error::error($message);
          }
        }
        elsif ($description ne '') {
          # unnamed options can't have a description
          my $message = 'Syntax error in definition '.($index + 1) 
              .': unnamed options can\'t have a description';
          locarnate::Error::error($message);
        }
        elsif ($unnamed_default && $default_value eq 'REQ') {
          # no default value for unnamed option after a
          # former option has one
          my $message = 'Syntax error in definition '
              .($unnamed_default_index + 1) 
              .': unnamed options can only have default values'
              .' at the end';
          locarnate::Error::error($message);
        }
        if ($long_name ne '' && $default_value ne 'REQ') {
          # found default value for unnamed option
          $unnamed_default = 1;
          $unnamed_default_index = $index;
        }

      }
      elsif (length($short_name) > 1) {
        # short names can only consists of a single char
        my $message = 'Syntax error in definition '.($index + 1) 
            .': short names can only consists of a single char';
        locarnate::Error::error($message);
      }

      # check for multiple names
      if ($long_name ne '') {
        if (exists($used_names{$long_name})) {
          my $message = 'Syntax error in definition '.($index + 1) 
              .': options name \''.$long_name 
              .'\' is already used in definition ' 
              .($used_names{$long_name} + 1);
          locarnate::Error::error($message);
        }
        $used_names{$long_name} = $index;
      }
      if ($short_name ne '' && $short_name ne 'NN') {
        if (exists($used_names{$short_name})) {
          my $message = 'Syntax error in definition '.($index + 1) 
              .': options name \''.$short_name 
              .'\' is already used in definition ' 
              .($used_names{$short_name} + 1);
          locarnate::Error::error($message);
        }
        $used_names{$short_name} = $index;
      }
    }
    else {
      # only default value is set (because of no end definition)
      my $message = 'Syntax error in definition '.($index + 1) 
          .': only default value is set';
      locarnate::Error::error($message);
    }
  }
};


my $_create_version_meassage = sub {
  my $self = shift(@_);
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $program_name = $_[0];
  my $version = $_[1];

  # remove path information
  $program_name =~ s/^.*\///;

  my $message = $program_name.' version '.$version."\n";
  $message .= "\n";
  $message .= "Written by Wolfgang Otto.\n";

  $self->{'version_message'} = $message;
};


my $_create_help_meassage = sub {
  my $self = shift(@_);
  if (@_ != 3) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];
  my $program_description = $_[1];
  my $program_name = $_[2];

  # remove path information
  $program_name =~ s/^.*\///;

  my $message = 'Usage: '.$program_name.' [OPTIONS]';

  # unnamed options
  for (my $index = 0; $index < @{$definition}; ++$index) {
    my $long_name = $definition->[$index]->[0];
    my $short_name = $definition->[$index]->[1];
    my $default_value = $definition->[$index]->[2];
    my $description = $definition->[$index]->[3];
    if ($short_name eq 'NN') {
      if ($long_name eq '') {
        $message .= ' '.$description;
      }
      elsif ($default_value eq "REQ") {
        $message .= ' '.$long_name;
      }
      else {
        $message .= ' ['.$long_name.']';
      }
    }
  }
  $message .= "\n";

  # program_description
  my $length = 0;
  for (my $i = 0; $i < length($program_description); ++$i) {
    if ($length > $message_columns 
        && substr($program_description, $i, 1) eq ' ') {
      $message .= "\n";
      $length = 0;
    }
    else {
      $message .= substr($program_description, $i, 1);
      ++$length;
    }
  }
  $message .= "\n\n";

  # header
  $message .= 'Options:'."\n";

  # named options
  for (my $index = 0; $index < @{$definition}; ++$index) {
    my $long_name = $definition->[$index]->[0];
    my $short_name = $definition->[$index]->[1];
    my $default_value = $definition->[$index]->[2];
    my $description = $definition->[$index]->[3];

    if ($short_name eq 'NN' || $description eq '') {
      # unnamed or secret options 
      next;
    }
    elsif ($long_name eq '' && $short_name eq '') {
      # header
      my $length = 0;
      for (my $i = 0; $i < length($description); ++$i) {
        if ($length > $message_columns 
            && substr($description, $i, 1) eq ' ') {
          $message .= "\n";
          $length = 0;
        }
        else {
          $message .= substr($description, $i, 1);
          ++$length;
        }
      }
      $message .= "\n";
    }
    else {
      # indent
      my $length = 0;
      $message .= '  ';
      $length += 2;

      # short name
      if ($short_name ne '') {
        $message .= '-'.$short_name;
        $length += 1 + length($short_name);
        if ($long_name ne '') {
          $message .= ', ';
          $length += 2;
        }
      }
      # long name
      if ($long_name ne '') {
        $message .= '--'.$long_name;
        $length += 2 + length($long_name);
      }
      # indent
      $message .= '  ';
      $length += 2;
      while ($length < $message_indent) {
        $message .= ' ';
        ++$length;
      }
      # description
      for (my $i = 0; $i < length($description); ++$i) {
        if ($length > $message_columns 
            && substr($description, $i, 1) eq ' ') {
          $message .= "\n".' 'x($message_indent + 2);
          $length = $message_indent + 2;
        }
        else {
          $message .= substr($description, $i, 1);
          ++$length;
        }
      }
      # default
      if ($message_defaults && $default_value ne 'FLAG' 
          && $default_value ne 'REQ' && $default_value ne '') {
        my $default_string = '['.$default_value.']';
        if ($length + 1 + length($default_string) > $message_columns) {
          $message .= "\n".' 'x($message_indent + 2);
        }
        else {
          $message .= ' ';
        }
        $message .= $default_string;
      }

      $message .= "\n";
    }
  }

  $self->{'help_message'} = $message;
};


my $_set_defaults = sub {
  my $self = shift(@_);
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];

  for (my $index = 0; $index < @{$definition}; ++$index) {
    my $long_name = $definition->[$index]->[0];
    my $short_name = $definition->[$index]->[1];
    my $default_value = $definition->[$index]->[2];

    if ($long_name ne '' || ($short_name ne '' && $short_name ne 'NN')) {
      # determine key string (primary: long_name, secondary: short_name)
      my $key = $long_name;
      if ($long_name eq '') {
        $key = $short_name;
      }
      if ($default_value eq 'FLAG') {
        $self->{'options'}->{$key} = 0;
      }
      else {
        $self->{'options'}->{$key} = $default_value;
      }
    }
  }
};


my $_parse_file = sub {
  my $self = shift(@_);
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];
  my $remote_file = $_[1];

  # check if remote file is used
  if (!defined($remote_file) || $remote_file eq '') {
    return;
  }

  # open file
  my $file = new FileHandle($remote_file, 'r');
  unless ($file) {
    return;
    # my $message = 'Can not open file \''.$remote_file.'\'';
    # locarnate::Error::error($message);
  }

  my $line_counter = 0;
  my %used_names = ();

  # check each line in remaote file
  while (my $line = $file->getline()) {
    ++$line_counter;
    chomp($line);

    # remove comments, trim line
    for (my $i = 0; $i < length($line); ++$i) {
      if (substr($line, $i, 1) eq '#' 
          && ($i == 0 || substr($line, $i - 1, 1) ne '\\')) {
        $line = substr($line, 0, $i);
        last;
      }
    }
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;

    # skip empty lines
    if ($line eq '') {
      next;
    }

    # divide line in name and value
    my $name = $line;
    my $value = 'FLAG';
    if ($line =~ m/^([^ \t]+)[ \t](.*)$/) {
      $name = $1;
      $value = $2;
      $name =~ s/^\s+//;
      $name =~ s/\s+$//;
      $value =~ s/^\s+//;
      $value =~ s/\s+$//;
      # test for special values
      if ($value eq 'REQ' || $value eq 'FLAG') {
        $value = '\\'.$value;
      }
    }

    # test for string
    if ((substr($value, 0, 1) eq '"' &&
         substr($value, length($value) - 1, 1) eq '"') || 
        (substr($value, 0, 1) eq '\'' &&
         substr($value, length($value) - 1, 1) eq '\'')) {
      $value = substr($value, 1, length($value) - 2);
    }

    # search options that match name
    my $index = &{$_find_index}($self, $definition, $name, 
                                $line_counter, $remote_file);

    if ($index < @{$definition}) {
      my $long_name = $definition->[$index]->[0];
      my $short_name = $definition->[$index]->[1];
      my $default_value = $definition->[$index]->[2];

      # determine key string (primary: long_name, secondary: short_name)
      my $key = $long_name;
      if ($long_name eq '') {
        $key = $short_name;
      }

      # test if argument was already used
      if (exists($used_names{$key})) { 
        my $message = 'Error: Option \''.$name.'\' in remote file \''.
            $remote_file.'\'", line '.$line_counter.
            ', is already set in line '.$used_names{$key}.'!'."\n";
        $self->{'error_message'} .= $message;
        next;
      }
      $used_names{$key} = $line_counter;
      
      # test syntax of default values and set options
      if ($default_value eq 'FLAG') {
        if ($value eq 'FLAG') {
          $self->{'options'}->{$key} = 1;
        }
        else {
          my $message = 'Error: Option \''.$name.'\' in remote file \''.
              $remote_file.'\'", line '.$line_counter.
              ', does not take an argument!'."\n";
          $self->{'error_message'} .= $message;
          $self->{'options'}->{$key} = 'undef';
        }
      }
      else {
        if ($value eq 'FLAG') {
          my $message = 'Error: Option \''.$name.'\' in remote file \''.
              $remote_file.'\'", line '.$line_counter.
              ', takes an argument!'."\n";
          $self->{'error_message'} .= $message;
          $self->{'options'}->{$key} = 'undef';
        }
        else {
          $self->{'options'}->{$key} = $value;
        }
      }
    }
  }
  $file->close();
};


my $_parse_named = sub {
  my $self = shift(@_);
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];
  my $arguments = $_[1];

  my $pos_counter = 0;
  my %used_names = ();

  my $i = 0;
  while ($i < @{$arguments}) {
    my $arg = $arguments->[$i];

    # check for unnamed argument

    if (substr($arg, 0, 1) ne '-') {
      ++$i;
      next;
    }
    
    ++$pos_counter;

    # argument is named -> remove it from list
    splice(@{$arguments}, $i, 1);

    my $name = '';
    my $value = 'FLAG';

    # determine argument name
    if (length($arg) > 1) {
      if (substr($arg, 1, 1) ne '-') {
        # short name
        $name = substr($arg, 1, 1);
        if (length($arg) > 2) {
          if (substr($arg, 2, 1) eq '=') {
            # value given per '='
            $value = substr($arg, 3);
            # test for special values
            if ($value eq 'REQ' || $value eq 'FLAG') {
              $value = '\\'.$value;
            }
          }
          else {
            # concatenation of short names (FLAGs) 
            # -> push remaining short names at begin
            my $corrected_arg = '-'.substr($arg, 2);
            splice(@{$arguments}, $i, 0, $corrected_arg);
          }
        }
      }
      else {
        # long name
        if ($arg =~ /^--([^=]+)=(.*)/) {
          $name = $1;
          $value = $2;
          # test for special values
          if ($value eq 'REQ' || $value eq 'FLAG') {
            $value = '\\'.$value;
          }
        }
        else {
          $name = substr($arg, 2);
        }
      }
    }
    if ($name eq '') {
      my $message = 'Error: Named command line argument at position '
          .$pos_counter.' has no name!'."\n";
      $self->{'error_message'} .= $message;
      next;
    }

    # search options that match name
    my $index = &{$_find_index}($self, $definition, $name, $pos_counter);

    if ($index < @{$definition}) {
      my $long_name = $definition->[$index]->[0];
      my $short_name = $definition->[$index]->[1];
      my $default_value = $definition->[$index]->[2];

      # determine key string (primary: long_name, secondary: short_name)
      my $key = $long_name;
      if ($long_name eq '') {
        $key = $short_name;
      }

      # test if argument was already used
      if (exists($used_names{$key})) { 
        my $message = 'Error: Named command line argument \''.$name
            .'\' at position '.$pos_counter.
            ', is already used at position '.$used_names{$key}.'!'."\n";
        $self->{'error_message'} .= $message;
        next;
      }
      $used_names{$key} = $pos_counter;
      
      # test syntax of default values and set options
      if ($default_value eq 'FLAG') {
        if ($value eq 'FLAG') {
          $self->{'options'}->{$key} = 1;
        }
        else {
          my $message = 'Error: Named command line argument \''.$name
              .'\' at position '.$pos_counter.
              ', does not take an argument!'."\n";
          $self->{'error_message'} .= $message;
          $self->{'options'}->{$key} = 'undef';
        }
      }
      else {
        # check if default_value was given in argument by '='
        if ($value eq 'FLAG') {
          # take next parameter as default_value
          if ($i == @{$arguments}) {
            my $message = 'Error: Named command line argument \''.$name
                .'\' at position '.$pos_counter.
                ', takes an argument!'."\n";
            $self->{'error_message'} .= $message;
            $value = 'undef';
          }
          else {
            $value = $arguments->[$i];
            if (substr($value, 0, 1) eq '-') {
              my $message = 'Error: Named command line argument \''.$name
                  .'\' at position '.$pos_counter.
                  ', takes an argument!'."\n";
              $self->{'error_message'} .= $message;
              $value = 'undef';
            }
            else {
              # remove value from list
              splice(@{$arguments}, $i, 1);
              # test for special values
              if ($value eq 'REQ' || $value eq 'FLAG') {
                $value = '\\'.$value;
              }
            }
          }
        }
        $self->{'options'}->{$key} = $value;
      }
    }
  }
};


my $_parse_unnamed = sub {
  my $self = shift(@_);
  if (@_ != 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];
  my $arguments = $_[1];

  for (my $index = 0; $index < @{$definition}; ++$index) {
    my $long_name = $definition->[$index]->[0];
    my $short_name = $definition->[$index]->[1];

    if ($short_name eq 'NN' && $long_name ne '') {
      # this definition is an unnamed option
      
      if (@{$arguments} == 0) {
        last;
      }
      else {
        my $value = shift(@{$arguments});
        # test for special values
        if ($value eq 'REQ' || $value eq 'FLAG') {
          $value = '\\'.$value;
        }

        # long name is key string
        $self->{'options'}->{$long_name} = $value;
      }
    }
  }

  # test for unparsed arguments
  if (@{$arguments} != 0) {
    $self->{'unparsed_arguments'} = $arguments;

    # my $message = 'Error: Unparsed ';
    # if (@{$arguments} == 1) {
    #   $message .= 'argument (\''.$arguments->[0].'\')!';
    # }
    # else {
    #   $message .= 'arguments (\''.$arguments->[0];
    #   for (my $i = 1; $i < @{$arguments} - 1; ++$i) {
    #     $message .= '\', \''.$arguments->[$i];
    #   }
    #   $message .= '\' and \''.$arguments->[@{$arguments} - 1].'\')!';
    # }
    # $message .= "\n";
    # 
    # $self->{'error_message'} .= $message;
  }
};


my $_check_completeness = sub {
  my $self = shift(@_);
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $definition = $_[0];

  my @names = ();

  for (my $index = 0; $index < @{$definition}; ++$index) {
    my $long_name = $definition->[$index]->[0];
    my $short_name = $definition->[$index]->[1];
    my $default_value = $definition->[$index]->[2];

    if ($long_name ne '' || ($short_name ne '' && $short_name ne 'NN')) {
      # determine key string (primary: long_name, secondary: short_name)
      my $key = $long_name;
      if ($long_name eq '') {
        $key = $short_name;
      }
      my $value = $self->{'options'}->{$key};

      if ($value eq 'REQ' || $value eq 'NN') {
        # value was not read yet
        push(@names, $key);
      }
    }
  }

  if (@names != 0) {
    my $message = 'Error: Missing ';
    if (@names == 1) {
      $message .= 'argument for option (\''.$names[0].'\')!';
    }
    else {
      $message .= 'arguments for options (\''.$names[0];
      for (my $i = 1; $i < @names - 1; ++$i) {
        $message .= '\', \''.$names[$i];
      }
      $message .= '\' and \''.$names[@names - 1].'\')!';
    }
    $message .= "\n";

    $self->{'error_message'} .= $message;
  }
};

###########################################################################
# constructor
###########################################################################

sub new {
  my $class = shift(@_);
  $class = ref($class) || $class;

  if (@_ < 3 || @_ > 4) {
    locarnate::Error::debug('Wrong argument number');
  }

  my %options = ();
  my @unparsed_arguments = ();

  my $self = {
    'version_message' => '',
    'help_message' => '',
    'error_message' => '',
    'options' => \%options,
    'unparsed_arguments' => \@unparsed_arguments
  };

  bless($self, $class);

  my $version = $_[0];
  my $description = $_[1];
  my $definition = $_[2];
  my $remote_file = '';
  if (@_ > 3) {
    $remote_file = $_[3];
  }

  my $program_name = $0;

  my @arguments = @ARGV;

  # translate empty strings in definition to undef
  for (my $index = 0; $index < @{$definition}; ++$index) {
    for (my $column = 0; $column < 4; ++$column) {
      if (!defined($definition->[$index]->[$column])) { 
        $definition->[$index]->[$column] = '';
      }
    }
  }

  &{$_check_definition}($self, $definition);
  &{$_create_version_meassage}($self, $program_name, $version);
  &{$_create_help_meassage}($self, $definition, $description, 
                            $program_name);
  &{$_set_defaults}($self, $definition);
  &{$_parse_file}($self, $definition, $remote_file);
  &{$_parse_named}($self, $definition, \@arguments);
  &{$_parse_unnamed}($self, $definition, \@arguments);
  &{$_check_completeness}($self, $definition);

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

  # $string .= $ind.'version message:'."\n";
  # $string .= '>'.$self->{'version_message'}.'<'."\n";
  # $string .= $ind.'help message:'."\n";
  # $string .= '>'.$self->{'help_message'}.'<'."\n";
  # $string .= $ind.'error message:'."\n";
  # $string .= '>'.$self->{'error_message'}.'<'."\n";

  $string .= $ind.'values:'."\n";
  foreach my $key (sort(keys(%{$self->{'options'}}))) {
    my $value = value($self, $key);

    $string .= $ind.'  '.$key.': ';
    if ($self->{'options'}->{$key} eq 'REQ' 
        || $self->{'options'}->{$key} eq 'NN') {
      $string .= 'undef';
    }
    else {
      $string .= $value;
    }
    $string .= "\n";
  }
  $string .= $ind.'arguments:'."\n";
  for (my $i = 0; $i < @{$self->{'unparsed_arguments'}}; ++$i) {
    my $value = argument($self, $i);
    $string .= $ind.'  '.$i.': '.$value."\n";
  }

  return ($string);
}


sub version_message {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  return($self->{'version_message'});
}

sub help_message {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  return($self->{'help_message'});
}

sub error {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  return($self->{'error_message'} ne '');
}

sub error_message {
  my $self = shift(@_);
  
  if (@_) {
    locarnate::Error::debug('Wrong argument number');
  }

  return($self->{'error_message'});
}

sub value {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $key = $_[0];
  
  unless (exists($self->{'options'}->{$key})) {
    my $message = 'No option with key \''.$key.'\' exists';
    locarnate::Error::debug($message);
  }

  my $value = $self->{'options'}->{$key};
  # test for escape
  if (substr($value, 0, 1) eq '\\') {
    $value = substr($value, 1);
  }
  return($value);
}

sub arguments_size {
  my $self = shift(@_);
  
  if (@_ != 0) {
    locarnate::Error::debug('Wrong argument number');
  }

  return(@{$self->{'unparsed_arguments'}});
}

sub argument {
  my $self = shift(@_);
  
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $index = $_[0];
  
  unless ($index < @{$self->{'unparsed_arguments'}}) {
    my $message = 'No argument with index \''.$index.'\' exists';
    locarnate::Error::debug($message);
  }

  my $value = $self->{'unparsed_arguments'}->[$index];
  # test for escape
  if (substr($value, 0, 1) eq '\\') {
    $value = substr($value, 1);
  }
  return($value);
}

1;
