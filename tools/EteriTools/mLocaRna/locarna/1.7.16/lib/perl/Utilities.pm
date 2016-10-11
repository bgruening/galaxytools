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

package locarnate::Utilities;

use strict;
use warnings;

use FileHandle;
use FindBin;
use Cwd;

use Error;


###########################################################################
# methods
###########################################################################

### string ################################################################

sub trim {
  # check correct arguments number and parse them
  if (@_ > 1) {
    locarnate::Error::debug('Wrong argument number');
  }
  my $string = $_[0];

  if (ref($string) eq 'SCALAR') {
    ${$string} =~ s/^\s+//;
    ${$string} =~ s/\s+$//;
    return(1);
  }
  else {
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return($string);
  }
}

sub uncomment {
  # check correct arguments number and parse them
  if (@_ > 1) {
    locarnate::Error::debug('Wrong argument number');
  }
  my $string = $_[0];

  if (ref($string) eq 'SCALAR') {
    if (${$string} =~ m/^\#/) {
      ${$string} = '';
    }
    elsif (${$string} =~ m/^(.*[^\\])\#/) {
      ${$string} = $1;
    }
    ${$string} =~ s/\\\#/\#/g;
    return(1);
  }
  else {
    if ($string =~ m/^\#/) {
      $string = '';
    }
    elsif ($string =~ m/^(.*[^\\])\#/) {
      $string = $1;
    }
    $string =~ s/\\\#/\#/g;
    return($string);
  }
}


### path ##################################################################

sub bin_dir {
  # check correct arguments number and parse them
  if (@_ > 0) {
    locarnate::Error::debug('Wrong argument number');
  }

  return($FindBin::RealBin);
}

sub cw_dir {
  # check correct arguments number and parse them
  if (@_ > 0) {
    locarnate::Error::debug('Wrong argument number');
  }

  return(cwd());
}

sub is_abs_dir {
  # check correct arguments number and parse them
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $path = $_[0];

  if (substr($path, 0, 1) eq '/') {
    return(1);
  }
  else {
    return(0);
  }
}

sub abs_dir {
  # check correct arguments number and parse them
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $path = $_[0];
  my $abs_dir = $path;
  if (!is_abs_dir($path)) {
    $abs_dir = cw_dir().'/'.$path;
  }
  $abs_dir =~ s/\/\.\//\//g;
  $abs_dir =~ s/\/+/\//g;

  return($abs_dir);
}

### file ##################################################################

sub dir_exists {
  # check correct arguments number and parse them
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $path = $_[0];

  if (-d $path) {
    return(1);
  }
  else {
    return(0);
  }
}

sub mk_dir {
  # check correct arguments number and parse them
  if (@_ != 1) {
    locarnate::Error::debug('Wrong argument number');
  }

  my $path = $_[0];

  my $subpath = '';
  if (substr($path, 0, 1) eq '/') {
    $subpath .= '/';
  }

  while ($path =~ m/([^\/]+)/g) {
    $subpath .= $1.'/';
    if (!dir_exists($subpath)) {
      my $success = mkdir($subpath);
      if (!$success) {
        locarnate::Error::error('Can\'t create directory \''.$subpath.'\'');
      }
    }
  }
}


1;
