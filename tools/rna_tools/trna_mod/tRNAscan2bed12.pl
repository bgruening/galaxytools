#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my $in = shift; 	#csv
my $out = shift; 	#bed12

open IN, "< $in" or die "can t open $in\n";
open OUT, "> $out" or die "can t open $out\n";


while(<IN>){
  chomp;
  my @F=split/\s+/,$_;

  my $chr        = $F[0];
  my $left       = '';
  my $right      = '';
  my $id         = "${chr}.tRNA${F[1]}-${F[4]}${F[5]}";
  my $score      = int($F[8]*10);
  $score      = ($score >1000)?(1000):($score);
  my $strand     = '';
  my $blockCount = '';
  my $blockSize  = '';
  my $blockStart = '';

  if($F[3]>$F[2]){
    $left = $F[2]-1;
    $right=$F[3];
    $strand = "+"
  }
  else{
    $left = $F[3]-1;
    $right=$F[2];
    $strand = "-";
  }

  if ($F[6] && $F[7]){
    $blockCount = 2;
    
    if ($strand eq "+"){
      $blockStart = join(",", 0, $F[7] - $F[2] + 1);
      $blockSize  = join(",", $F[6] - $F[2], $F[3] - $F[7]);
    }
    elsif($strand eq "-"){
      $blockStart = join(",", 0, $F[6] - $F[3] + 1);
      $blockSize  = join(",", $F[7] - $F[3], $F[2] - $F[6]);
    }
  }
  else{
    $blockCount = 1;
    $blockSize  = abs($F[3]-$F[2])+1;
    $blockStart = 0;
  }
  print OUT join("\t", $chr, $left, $right, $id, $score, $strand, $left, $right, 0, $blockCount, $blockSize, $blockStart)."\n";
}
