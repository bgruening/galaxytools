package Cwd;
use 5.006;

=head1 NAME

Cwd - get pathname of current working directory

=head1 SYNOPSIS

    use Cwd;
    my $dir = getcwd;

    use Cwd 'abs_path';
    my $abs_path = abs_path($file);

=head1 DESCRIPTION

This module provides functions for determining the pathname of the
current working directory.  It is recommended that getcwd (or another
*cwd() function) be used in I<all> code to ensure portability.

By default, it exports the functions cwd(), getcwd(), fastcwd(), and
fastgetcwd() into the caller's namespace.  


=head2 getcwd and friends

Each of these functions are called without arguments and return the
absolute path of the current working directory.

=over 4

=item getcwd

    my $cwd = getcwd();

Returns the current working directory.

Re-implements the getcwd(3) (or getwd(3)) functions in Perl.

Taint-safe.

=item cwd

    my $cwd = cwd();

The cwd() is the most natural form for the current architecture. For
most systems it is identical to `pwd` (but without the trailing line
terminator).

Taint-safe.

=item fastcwd

    my $cwd = fastcwd();

A more dangerous version of getcwd(), but potentially faster.

It might conceivably chdir() you out of a directory that it can't
chdir() you back into.  If fastcwd encounters a problem it will return
undef but will probably leave you in a different directory.  For a
measure of extra security, if everything appears to have worked, the
fastcwd() function will check that it leaves you in the same directory
that it started in. If it has changed it will C<die> with the message
"Unstable directory path, current directory changed
unexpectedly". That should never happen.

=item fastgetcwd

  my $cwd = fastgetcwd();

The fastgetcwd() function is provided as a synonym for cwd().

=back


=head2 abs_path and friends

These functions are exported only on request.  They each take a single
argument and return the absolute pathname for it.

=over 4

=item abs_path

  my $abs_path = abs_path($file);

Uses the same algorithm as getcwd().  Symbolic links and relative-path
components ("." and "..") are resolved to return the canonical
pathname, just like realpath(3).

Taint-safe.

=item realpath

  my $abs_path = realpath($file);

A synonym for abs_path().

Taint-safe.

=item fast_abs_path

  my $abs_path = fast_abs_path($file);

A more dangerous, but potentially faster version of abs_path.

This function is B<Not> taint-safe : you can't use it in programs
that work under taint mode.

=back

=head2 $ENV{PWD}

If you ask to override your chdir() built-in function, 

  use Cwd qw(chdir);

then your PWD environment variable will be kept up to date.  Note that
it will only be kept up to date if all packages which use chdir import
it from Cwd.


=head1 NOTES

=over 4

=item *

Since the path seperators are different on some operating systems ('/'
on Unix, ':' on MacPerl, etc...) we recommend you use the File::Spec
modules wherever portability is a concern.

=item *

Actually, on Mac OS, the C<getcwd()>, C<fastgetcwd()> and C<fastcwd()>
functions  are all aliases for the C<cwd()> function, which, on Mac OS,
calls `pwd`. Likewise, the C<abs_path()> function is an alias for
C<fast_abs_path()>.

=back

=head1 SEE ALSO

L<File::chdir>

=cut

use strict;

use Carp;

our $VERSION = '2.06';

use base qw/ Exporter /;
our @EXPORT = qw(cwd getcwd fastcwd fastgetcwd);
our @EXPORT_OK = qw(chdir abs_path fast_abs_path realpath fast_realpath);

# sys_cwd may keep the builtin command

# All the functionality of this module may provided by builtins,
# there is no sense to process the rest of the file.
# The best choice may be to have this in BEGIN, but how to return from BEGIN?

if ($^O eq 'os2' && defined &sys_cwd && defined &sys_abspath) {
    local $^W = 0;
    *cwd		= \&sys_cwd;
    *getcwd		= \&cwd;
    *fastgetcwd		= \&cwd;
    *fastcwd		= \&cwd;
    *abs_path		= \&sys_abspath;
    *fast_abs_path	= \&abs_path;
    *realpath		= \&abs_path;
    *fast_realpath	= \&abs_path;
    return 1;
}

eval {
    require XSLoader;
    undef *Cwd::fastcwd; # avoid redefinition warning
    XSLoader::load('Cwd');
};


# Find the pwd command in the expected locations.  We assume these
# are safe.  This prevents _backtick_pwd() consulting $ENV{PATH}
# so everything works under taint mode.
my $pwd_cmd;
foreach my $try (qw(/bin/pwd /usr/bin/pwd)) {
    if( -x $try ) {
        $pwd_cmd = $try;
        last;
    }
}
$pwd_cmd ||= 'pwd';

# The 'natural and safe form' for UNIX (pwd may be setuid root)
sub _backtick_pwd {
    local @ENV{qw(PATH IFS CDPATH ENV BASH_ENV)};
    my $cwd = `$pwd_cmd`;
    # Belt-and-suspenders in case someone said "undef $/".
    local $/ = "\n";
    # `pwd` may fail e.g. if the disk is full
    chomp($cwd) if defined $cwd;
    $cwd;
}

# Since some ports may predefine cwd internally (e.g., NT)
# we take care not to override an existing definition for cwd().

unless(defined &cwd) {
    # The pwd command is not available in some chroot(2)'ed environments
    if( $^O eq 'MacOS' || (defined $ENV{PATH} && 
                           grep { -x "$_/pwd" } split(':', $ENV{PATH})) ) 
    {
	*cwd = \&_backtick_pwd;
    }
    else {
	*cwd = \&getcwd;
    }
}

# set a reasonable (and very safe) default for fastgetcwd, in case it
# isn't redefined later (20001212 rspier)
*fastgetcwd = \&cwd;

# By Brandon S. Allbery
#
# Usage: $cwd = getcwd();

sub getcwd
{
    abs_path('.');
}


# By John Bazik
#
# Usage: $cwd = &fastcwd;
#
# This is a faster version of getcwd.  It's also more dangerous because
# you might chdir out of a directory that you can't chdir back into.
    
sub fastcwd {
    my($odev, $oino, $cdev, $cino, $tdev, $tino);
    my(@path, $path);
    local(*DIR);

    my($orig_cdev, $orig_cino) = stat('.');
    ($cdev, $cino) = ($orig_cdev, $orig_cino);
    for (;;) {
	my $direntry;
	($odev, $oino) = ($cdev, $cino);
	CORE::chdir('..') || return undef;
	($cdev, $cino) = stat('.');
	last if $odev == $cdev && $oino == $cino;
	opendir(DIR, '.') || return undef;
	for (;;) {
	    $direntry = readdir(DIR);
	    last unless defined $direntry;
	    next if $direntry eq '.';
	    next if $direntry eq '..';

	    ($tdev, $tino) = lstat($direntry);
	    last unless $tdev != $odev || $tino != $oino;
	}
	closedir(DIR);
	return undef unless defined $direntry; # should never happen
	unshift(@path, $direntry);
    }
    $path = '/' . join('/', @path);
    if ($^O eq 'apollo') { $path = "/".$path; }
    # At this point $path may be tainted (if tainting) and chdir would fail.
    # Untaint it then check that we landed where we started.
    $path =~ /^(.*)\z/s		# untaint
	&& CORE::chdir($1) or return undef;
    ($cdev, $cino) = stat('.');
    die "Unstable directory path, current directory changed unexpectedly"
	if $cdev != $orig_cdev || $cino != $orig_cino;
    $path;
}


# Keeps track of current working directory in PWD environment var
# Usage:
#	use Cwd 'chdir';
#	chdir $newdir;

my $chdir_init = 0;

sub chdir_init {
    if ($ENV{'PWD'} and $^O ne 'os2' and $^O ne 'dos' and $^O ne 'MSWin32') {
	my($dd,$di) = stat('.');
	my($pd,$pi) = stat($ENV{'PWD'});
	if (!defined $dd or !defined $pd or $di != $pi or $dd != $pd) {
	    $ENV{'PWD'} = cwd();
	}
    }
    else {
	my $wd = cwd();
	$wd = Win32::GetFullPathName($wd) if $^O eq 'MSWin32';
	$ENV{'PWD'} = $wd;
    }
    # Strip an automounter prefix (where /tmp_mnt/foo/bar == /foo/bar)
    if ($^O ne 'MSWin32' and $ENV{'PWD'} =~ m|(/[^/]+(/[^/]+/[^/]+))(.*)|s) {
	my($pd,$pi) = stat($2);
	my($dd,$di) = stat($1);
	if (defined $pd and defined $dd and $di == $pi and $dd == $pd) {
	    $ENV{'PWD'}="$2$3";
	}
    }
    $chdir_init = 1;
}

sub chdir {
    my $newdir = @_ ? shift : '';	# allow for no arg (chdir to HOME dir)
    $newdir =~ s|///*|/|g unless $^O eq 'MSWin32';
    chdir_init() unless $chdir_init;
    my $newpwd;
    if ($^O eq 'MSWin32') {
	# get the full path name *before* the chdir()
	$newpwd = Win32::GetFullPathName($newdir);
    }

    return 0 unless CORE::chdir $newdir;

    if ($^O eq 'VMS') {
	return $ENV{'PWD'} = $ENV{'DEFAULT'}
    }
    elsif ($^O eq 'MacOS') {
	return $ENV{'PWD'} = cwd();
    }
    elsif ($^O eq 'MSWin32') {
	$ENV{'PWD'} = $newpwd;
	return 1;
    }

    if ($newdir =~ m#^/#s) {
	$ENV{'PWD'} = $newdir;
    } else {
	my @curdir = split(m#/#,$ENV{'PWD'});
	@curdir = ('') unless @curdir;
	my $component;
	foreach $component (split(m#/#, $newdir)) {
	    next if $component eq '.';
	    pop(@curdir),next if $component eq '..';
	    push(@curdir,$component);
	}
	$ENV{'PWD'} = join('/',@curdir) || '/';
    }
    1;
}


# In case the XS version doesn't load.
*abs_path = \&_perl_abs_path unless defined &abs_path;
sub _perl_abs_path
{
    my $start = @_ ? shift : '.';
    my($dotdots, $cwd, @pst, @cst, $dir, @tst);

    unless (@cst = stat( $start ))
    {
	carp "stat($start): $!";
	return '';
    }
    $cwd = '';
    $dotdots = $start;
    do
    {
	$dotdots .= '/..';
	@pst = @cst;
	unless (opendir(PARENT, $dotdots))
	{
	    carp "opendir($dotdots): $!";
	    return '';
	}
	unless (@cst = stat($dotdots))
	{
	    carp "stat($dotdots): $!";
	    closedir(PARENT);
	    return '';
	}
	if ($pst[0] == $cst[0] && $pst[1] == $cst[1])
	{
	    $dir = undef;
	}
	else
	{
	    do
	    {
		unless (defined ($dir = readdir(PARENT)))
	        {
		    carp "readdir($dotdots): $!";
		    closedir(PARENT);
		    return '';
		}
		$tst[0] = $pst[0]+1 unless (@tst = lstat("$dotdots/$dir"))
	    }
	    while ($dir eq '.' || $dir eq '..' || $tst[0] != $pst[0] ||
		   $tst[1] != $pst[1]);
	}
	$cwd = (defined $dir ? "$dir" : "" ) . "/$cwd" ;
	closedir(PARENT);
    } while (defined $dir);
    chop($cwd) unless $cwd eq '/'; # drop the trailing /
    $cwd;
}


# added function alias for those of us more
# used to the libc function.  --tchrist 27-Jan-00
*realpath = \&abs_path;

sub fast_abs_path {
    my $cwd = getcwd();
    require File::Spec;
    my $path = @_ ? shift : File::Spec->curdir;
    CORE::chdir($path) || croak "Cannot chdir to $path: $!";
    my $realpath = getcwd();
    -d $cwd && CORE::chdir($cwd) ||
	croak "Cannot chdir back to $cwd: $!";
    $realpath;
}

# added function alias to follow principle of least surprise
# based on previous aliasing.  --tchrist 27-Jan-00
*fast_realpath = \&fast_abs_path;


# --- PORTING SECTION ---

# VMS: $ENV{'DEFAULT'} points to default directory at all times
# 06-Mar-1996  Charles Bailey  bailey@newman.upenn.edu
# Note: Use of Cwd::chdir() causes the logical name PWD to be defined
#   in the process logical name table as the default device and directory
#   seen by Perl. This may not be the same as the default device
#   and directory seen by DCL after Perl exits, since the effects
#   the CRTL chdir() function persist only until Perl exits.

sub _vms_cwd {
    return $ENV{'DEFAULT'};
}

sub _vms_abs_path {
    return $ENV{'DEFAULT'} unless @_;
    my $path = VMS::Filespec::pathify($_[0]);
    croak("Invalid path name $_[0]") unless defined $path;
    return VMS::Filespec::rmsexpand($path);
}

sub _os2_cwd {
    $ENV{'PWD'} = `cmd /c cd`;
    chop $ENV{'PWD'};
    $ENV{'PWD'} =~ s:\\:/:g ;
    return $ENV{'PWD'};
}

sub _win32_cwd {
    $ENV{'PWD'} = Win32::GetCwd();
    $ENV{'PWD'} =~ s:\\:/:g ;
    return $ENV{'PWD'};
}

*_NT_cwd = \&_win32_cwd if (!defined &_NT_cwd && 
                            defined &Win32::GetCwd);

*_NT_cwd = \&_os2_cwd unless defined &_NT_cwd;

sub _dos_cwd {
    if (!defined &Dos::GetCwd) {
        $ENV{'PWD'} = `command /c cd`;
        chop $ENV{'PWD'};
        $ENV{'PWD'} =~ s:\\:/:g ;
    } else {
        $ENV{'PWD'} = Dos::GetCwd();
    }
    return $ENV{'PWD'};
}

sub _qnx_cwd {
	local $ENV{PATH} = '';
	local $ENV{CDPATH} = '';
	local $ENV{ENV} = '';
    $ENV{'PWD'} = `/usr/bin/fullpath -t`;
    chop $ENV{'PWD'};
    return $ENV{'PWD'};
}

sub _qnx_abs_path {
	local $ENV{PATH} = '';
	local $ENV{CDPATH} = '';
	local $ENV{ENV} = '';
    my $path = @_ ? shift : '.';
    my $realpath=`/usr/bin/fullpath -t $path`;
    chop $realpath;
    return $realpath;
}

sub _epoc_cwd {
    $ENV{'PWD'} = EPOC::getcwd();
    return $ENV{'PWD'};
}

{
    no warnings;	# assignments trigger 'subroutine redefined' warning

    if ($^O eq 'VMS') {
        *cwd		= \&_vms_cwd;
        *getcwd		= \&_vms_cwd;
        *fastcwd	= \&_vms_cwd;
        *fastgetcwd	= \&_vms_cwd;
        *abs_path	= \&_vms_abs_path;
        *fast_abs_path	= \&_vms_abs_path;
    }
    elsif ($^O eq 'NT' or $^O eq 'MSWin32') {
        # We assume that &_NT_cwd is defined as an XSUB or in the core.
        *cwd		= \&_NT_cwd;
        *getcwd		= \&_NT_cwd;
        *fastcwd	= \&_NT_cwd;
        *fastgetcwd	= \&_NT_cwd;
        *abs_path	= \&fast_abs_path;
        *realpath   = \&fast_abs_path;
    }
    elsif ($^O eq 'os2') {
        # sys_cwd may keep the builtin command
        *cwd		= defined &sys_cwd ? \&sys_cwd : \&_os2_cwd;
        *getcwd		= \&cwd;
        *fastgetcwd	= \&cwd;
        *fastcwd	= \&cwd;
        *abs_path	= \&fast_abs_path;
    }
    elsif ($^O eq 'dos') {
        *cwd		= \&_dos_cwd;
        *getcwd		= \&_dos_cwd;
        *fastgetcwd	= \&_dos_cwd;
        *fastcwd	= \&_dos_cwd;
        *abs_path	= \&fast_abs_path;
    }
    elsif ($^O =~ m/^(?:qnx|nto)$/ ) {
        *cwd		= \&_qnx_cwd;
        *getcwd		= \&_qnx_cwd;
        *fastgetcwd	= \&_qnx_cwd;
        *fastcwd	= \&_qnx_cwd;
        *abs_path	= \&_qnx_abs_path;
        *fast_abs_path	= \&_qnx_abs_path;
    }
    elsif ($^O eq 'cygwin') {
        *getcwd	= \&cwd;
        *fastgetcwd	= \&cwd;
        *fastcwd	= \&cwd;
        *abs_path	= \&fast_abs_path;
    }
    elsif ($^O eq 'epoc') {
        *cwd            = \&_epoc_cwd;
        *getcwd	        = \&_epoc_cwd;
        *fastgetcwd	= \&_epoc_cwd;
        *fastcwd	= \&_epoc_cwd;
        *abs_path	= \&fast_abs_path;
    }
    elsif ($^O eq 'MacOS') {
    	*getcwd     = \&cwd;
    	*fastgetcwd = \&cwd;
    	*fastcwd    = \&cwd;
    	*abs_path   = \&fast_abs_path;
    }
}


1;
