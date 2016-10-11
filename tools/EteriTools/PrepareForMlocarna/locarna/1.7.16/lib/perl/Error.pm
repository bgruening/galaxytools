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

package locarnate::Error;

use strict;
use warnings;

use FileHandle;

###########################################################################
# methods
###########################################################################

sub warning {
  if (@_ > 2) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $message = '';
  if (@_ > 0) {
    $message = $_[0];
  }
  my $break = 0;
  if (@_ > 1) {
    $break = $_[1];
  }
  
  my $what = 'Warning';

  if ($message ne '') {
    $what .= ': '.$message;
  }
  $what .= '!';

  if ($break) {
    $what .= "\n";
    $what .= 'Press "Enter" to continue or "Control-C" for breakup...';
    *STDERR->print($what);
    getc(*STDIN);
  }
  else {
    *STDERR->print($what."\n");
  }
}

sub error {
  if (@_ > 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $message = '';
  if (@_ > 0) {
    $message = $_[0];
  }

  my $what = 'Error';
  if ($message ne '') {
    $what .= ': '.$message;
  }
  $what .= '!';

  *STDERR->print($what."\n");
  exit(1);
}

sub debug {
  if (@_ > 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $message = '';
  if (@_ > 0) {
    $message = $_[0];
  }

  my @calls = ();
  my $i = 0;
  while (my @caller = caller($i)) {
    push(@calls, \@caller);
    ++$i
  }
  push(@calls, ['', '', 0, 'main']);

  my $what = 'Abort';
  if ($message ne '') {
    $what .= ': '.$message;
  }
  $what .= '!';
  *STDERR->print($what."\n");
  for (my $i = 0; $i < @calls - 1; ++$i) {
    my $file_name = $calls[$i]->[1];
    my $line = $calls[$i]->[2];
    my $function_name = $calls[$i + 1]->[3];
    my $location = $file_name.': '.$function_name.': '.$line;
    *STDERR->print('('.$i.') '.$location."\n");
  }
  exit(1);
}

1;
