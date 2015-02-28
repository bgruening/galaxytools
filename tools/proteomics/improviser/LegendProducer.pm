#!/usr/bin/perl
use GD::Simple;
use strict;
use warnings;

sub drawLegend {
	my $zn = $_[0];
	my $zp = $_[1];
	my $rmin = $_[2];
	my $rmax = $_[3];
	my $outdir = $_[4];
	 	
	my $x1=10;
	my $y1=10;
	my $x2=40;
	my $y2=20;
	my $box_width=20;

	my $image_width= $box_width*25+50+20;
	my $img = GD::Simple->new($image_width,$y2+30);
	my $black = $img->colorAllocate(0,0,0);	
	$img->penSize(0,0);
	for (my $r=128, my $g = 0, my $b=0; $r<256;){
				$img->bgcolor($r,$g,$b);
				$img->fgcolor(undef);
				$img->rectangle($x1,$y1,$x2,$y2);    			
				$r=$r+10;
				$g=$g+20;
				$b=$b+20;
				$x1=$x1+$box_width;
				$x2=$x2+$box_width;
	}

	$img->bgcolor('white');
	$img->fgcolor('white');
	$img->rectangle($x1,$y1,$x2,$y2);    			
	
	# zn
	$zn = sprintf "%.2f", $zn;	 
	$img->moveTo($x2-60,$y2+25);
	$img->fgcolor('black');
	$img->font('Times');    
	$img->fontsize(14);    
	#$img->stringFT($black,'Times', 12, 0, 20, 20, $zn);
	$img->string($zn);	
	
	#jump to the last postion
	$x1=$image_width-30;
	$x2=$image_width-10;
	
	for (my $r=0, my $g = 128, my $b=0; $r<256;){
				$img->bgcolor($r,$g,$b);
				$img->fgcolor(undef);
				$img->rectangle($x1,$y1,$x2,$y2);    			
				$r=$r+20;
				$g=$g+10;
				$b=$b+20;
				$x1=$x1-$box_width;
				$x2=$x2-$box_width;
	}
	
	# zp
	$zp = sprintf "%.2f", $zp;
	$img->moveTo($x2,$y2+25);
	$img->fgcolor('black');
	$img->font('Times');    
	$img->fontsize(14);  
	#$img->stringFT($black,'Times', 12, 0, 20, 20, $zp);	  
	$img->string($zp);
	
	# rmin 
	$rmin = sprintf "%.2f", $rmin;
	$img->moveTo(10,$y2+25);
	$img->fgcolor('black');
	$img->font('Times');    
	$img->fontsize(14);   
	$img->string("\x{2264} $rmin");
	
	# rmax
	$rmax = sprintf "%.2f", $rmax;
	$img->moveTo($image_width-60,$y2+25);
	$img->fgcolor('black');
	$img->font('Times');    
	$img->fontsize(14); 
	$img->string("\x{2265} $rmax");
	  
	open OUT, ">$outdir/css/ratio_color_legend.png" or die $!,"\n";
	# convert into png data
	binmode(OUT);
	print OUT $img->png;
	close OUT;
}
1;
