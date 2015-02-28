#!/usr/bin/perl
use warnings;
use strict;
use FindBin '$Bin';
#use lib "$FindBin::Bin/";
#use lib "/home/wsimp/perl5/lib/perl5";
use Getopt::Long;
use Cwd 'abs_path';
use Pod::Usage;
use Bio::Graphics;
use XML::Twig;
use Bio::SeqFeature::Generic;
use Math::Round;
use UniProtFileParser;
use List::Util qw[min max sum];
use LegendProducer;
use File::Copy::Recursive;
use LWP::Simple;

=head1 NAME

plotPeptideASAPRatios.pl

=head1 SYNOPSIS

perl plotPeptideASAPRatios.pl -in [input XML file] -out [output directory for plots] 

=head1 OPTIONS

	-help	    brief help message
	-man	    full documentation
	-in	XML     input file (interact.prot.xml)
	-out	    output directory for plots
	-threshold  ProteinProphet probability cutoff i.e. if the ProteinProphet probability for a protein is lesser than this cutoff, it is discarded. Default: 0.90. 
	-rmax	    maximum peptide ratio. i.e, color of peptide with ratio greater than or equal to rmax is greenish. Default: 3.0.
	-rmin	    minimum peptide ratio. i.e, color of peptide with ratio less than or equal to rmin is redish. Default: -3.0. 
	-zp		    positive no change zone. petide with ratio from 0 to zp is white. Default: 0.25  	
	-zn		    negetive no change zone. petide with ratio from 0 to -zp is white. Default: -0.25 

=head1 DESCRIPTION



=head1 AUTHOR

Pavankumar Videm - <videmp@informatik.uni-freiburg.de>
Deepika Gunasekaran - <mailtodeepika@gmail.com>

=cut


my ($help, $man, $in, $outdir, $invert_ratios, $validate_asap_xpress, $elaborate_peptides);

my $threshold = 0.90;
my $rmin = -3.0;
my $rmax = 3.0;
my $zp = 0.25;
my $zn = -0.25;

my $options = GetOptions ("help"	    => \$help,
                         "man"  	    => \$man,
                         "in=s"		    => \$in,
                         "out=s"	    => \$outdir,
                         "threshold=s"	=> \$threshold,
                         "rmin=s"	    => \$rmin,
                         "rmax=s"	    => \$rmax,
                         "zp=s"		    => \$zp,
                         "zn=s"		    => \$zn,
                         "invert"	    => \$invert_ratios,
                         "validate"	    => \$validate_asap_xpress,
                         "elaborate"    => \$elaborate_peptides
                         ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(-exitstatus => 1, -verbose => 1) if (!$in && !$outdir);

#$in = abs_path($in);
mkdir $outdir if(!-d $outdir);
#$outdir = abs_path($outdir);

my $considered_count = 0;
my $discarded_count = 0;
my $mol_wgt_sum;
my $molecular_weight = 0;
my $mol_wgt_avg = 0;

my @prot_asap_ratios = ();
my @prot_xpress_ratios = ();

## create required directories ##
mkdir "$outdir/images" if (! -d "$outdir/images");
mkdir "$outdir/small_images" if (! -d "$outdir/small_images");
mkdir "$outdir/index_files" if (! -d "$outdir/index_files");
mkdir "$outdir/js";
mkdir "$outdir/css";
mkdir "$outdir/temp";

## Copy the java scripts and style sheets to output directory ## 
File::Copy::Recursive::dircopy "$Bin/js/", "$outdir/js" or die "Copy failed: $!";
File::Copy::Recursive::dircopy "$Bin/css/", "$outdir/css" or die "Copy failed: $!";
my %colorHash = ();
populateColorHash();
my %green_bins_hash = ();
my %red_bins_hash = ();
computeBins($rmin, $rmax);
drawLegend($zn, $zp, $rmin, $rmax, $outdir);
my @prot_ratio_means_all = ();
my @prot_adj_ratio_means_all = ();
my @pep_ratio_means_all = ();
my @prot_ratio_means_filtered = ();
my @prot_adj_ratio_means_filtered = ();
my @pep_ratio_means_filtered = ();
## Open log file ##  
open LOG, ">$outdir/run_stats.out" or die "Could not create LOG file: $!\n";
#print "Caluculating average ratios...<br>\n";
print LOG "Following proteins are considered for average calculation:\n";
print LOG "uniprot_id\tasap ratio\txpress ratio\n";  

my $parser = XML::Twig->new();
$parser->parsefile($in);
my $root= $parser->root;
my @asap_summary_tags = $root->getElementsByTagName("ASAP_pvalue_analysis_summary");
my $correction_factor = $asap_summary_tags[0]->{'att'}->{'background_ratio_mean'};
$correction_factor = 10**$correction_factor; 

close LOG;

## Open main output HTML file. Always named as index.html in output directory ##
## Add required metadata and create main table ## 
open HTML, ">$outdir/index.html" or die $!;
print HTML "<html>\n";
print HTML "<head>\n";
print HTML "\t<link rel=\"stylesheet\" type=\"text/css\" href=\"css/style.css\"/>\n";
## print HTML "\t<script type=\"text/javascript\" src=\"js/jquery.tablesorter/jquery.tablesorter.min.js\"></script>\n"; ##
print HTML "\t<script type=\"text/javascript\" src=\"js/jquery.tablesorter/jquery-latest.js\"></script>\n";
print HTML "\t<script type=\"text/javascript\" src=\"js/jquery.tablesorter/jquery.tablesorter.js\"></script>\n";
print HTML "\t<script type=\"text/javascript\" src=\"js/jquery.tablesorter/jquery.metadata.js\"></script>\n";
print HTML "\t<script type=\"text/javascript\" src=\"js/picnet.table.filter.min.js\"></script>\n";
print HTML "\t<script type=\"text/javascript\" src=\"js/jquery.thfloat-0.3.min.js\"></script>\n";		
print HTML "\t<script type=\"text/javascript\">\n";
print HTML "\t\t\$(document).ready(function() {\n";
print HTML "\t\t\t\$(\"#asapTable\").tableFilter();\n";
print HTML "\t\t\t\$(\"table\").tablesorter({debug: true});\n";
print HTML "\$('table tr:eq(1) td:last').css({'background-image': 'url(css/ratio_color_legend.png)', 'background-repeat': 'no-repeat', 'background-position': 'center', 'height': '50px', 'background-size': '77%' });\n";
print HTML "\t\t\$(window).scroll(function(){\n";
print HTML "\t\t\t\$('thead').css('left','-'+\$(window).scrollLeft());\n";
print HTML "\t\t});\n";
print HTML "\t});\n";
print HTML "\t</script>\n";		
print HTML "</head>\n";

## Create table with required fields in index.html ##
print HTML "<body>\n";
print HTML "<table border=\"1\" cellpadding=\"10\" id=\"asapTable\" class=\"tablesorter\" style=\"border-collapse:collapse;\">\n";
print HTML "<thead>\n";
print HTML "\t<tr>\n";
print HTML "\t\t<th id=\"upid\"><div title=\"Unique UniProt id consistent with the UniProt database (sortable).\">Uniprot id</div></th>\n";
print HTML "\t\t<th id=\"recname\"><div title=\"The recommended name obtained from the UniProt database (alphabetically sortable).\">Recommended name</div></th>\n";
print HTML "\t\t<th id=\"molwt\"><div title=\"The molecular weight of the given proteins obtained from the UniProt database (numerically sortable). Additionally the sum of the molecular weights of all the proteins present in the table are provided in the average_molecular weight.txt output file.\">MW</div></th>\n";
print HTML "\t\t<th id=\"protratio\"><div title=\"The ASAP Ratios obtained from the ProtXML file (numerically sortable).\">Global ratio of protein</div></th>\n";
print HTML "\t\t<th id=\"protratiostddev\"><div title=\"Log(base2) transformed standard deviations from the ASAP Ratios obtained from the ProtXML file (numerically sortable).\">Std dev</div></th>\n";
print HTML "\t\t<th id=\"signalpep\" filter-type='ddl'><div title=\"Proteins are classified based on the presence of an N-terminal signal peptide. Information was obtained from the UniProt database (sortable).\">Signal</div></th>\n";
print HTML "\t\t<th id=\"tm\" filter-type='ddl'><div title=\"Proteins are classified based on whether they are transmembrane proteins or not. Information was obtained from the UniProt database (sortable).\">Transmem</div></th>\n";
print HTML "\t\t<th id=\"icr\"><div title=\"Intracellular ratio is the ratio of non-membrane region of a membrane-spanning (transmembrane) protein found in the cytoplasm. Information was obtained from the UniProt database (numerically sortable).\">Intra-cellular ratio</div></th>\n";
print HTML "\t\t<th id=\"ecr\"><div title=\"Extracellular ratio is the ratio of non-membrane region of a membrane-spanning (transmembrane) protein found in the extracellular region. Information was obtained from the UniProt database (numerically sortable).\">Extra-cellular ratio</div></th>\n";
print HTML "\t\t<th id=\"plot\" filter='false'><div title=\"Visualization of the peptides on the protein. Clicking on the image for each protein will give a detailed description of the peptides.\">Image</div></th>\n";
print HTML "\t</tr>\n";
print HTML "</thead>\n";

print HTML "<tbody>\n";
open LOG, ">>$outdir/run_stats.out" or die $!;
my $parser2 = XML::Twig->new(
    twig_handlers => { 
    				'protein' => \&addProteinToTable
    				});
    				
$parser2->parsefile($in);    	
#$parser2->flush;
print HTML "</tbody>\n";
print HTML "</table>\n";
print HTML "</body>\n";
print HTML "</html>\n";
close HTML;
print LOG "---------------------------------------------\n";
print LOG "Number of proteins considered: $considered_count\n";
print LOG "Number of proteins discarded: $discarded_count\n";
print LOG "---------------------------------------------\n\n";
print LOG "\n";
close LOG;
print "<br>---------------------------------------------<br>\n";
print "Number of proteins considered: $considered_count<br>\n";
print "Number of proteins discarded: $discarded_count<br>\n";
print "---------------------------------------------<br><br>\n\n"; 

system "rm -r $outdir/temp";

#plotFCValues();
#my $twig_obj=XML::Twig->new();
#$twig_obj->parsefile( $in); 
#plotNormalDistribution( $twig_obj);  

#######################################################################################
## Print output text file with compatible Browsers and operating systems
#######################################################################################

open FILE, ">$outdir/suppoted_browsers_and_os.txt" or die $!; 
print FILE "!!! WARNING !!!\n";
print FILE "Formatted output is NOT supported in Windows 2003 due to reasons unknown\n\n";
print FILE "!!! WARNING !!!\n";
print FILE "Formatted output is NOT supported in Internet Explorer due to reasons unknown\n\n";
print FILE "COMPATIBLE OPERATING SYSTEMS\n";
print FILE "Formatted output is supported in the following operating systems:\n"; 
print FILE "Mac OS\n";
print FILE "Ubuntu\n";
print FILE "WIndows XP\n\n";
print FILE "COMPATIBLE BROWSERS\n";
print FILE "Formatted output is supported in the following browsers:\n"; 
print FILE "Opera: version 10.5 and above\n";
print FILE "Google Chrome: version 5 and above\n";
print FILE "Mozilla Firefox: version 4 and above\n";
print FILE "Safari: version 5 and above\n\n";  
close FILE;   

#######################################################################################
## Print Average molecular weight to another file
#######################################################################################
$mol_wgt_avg = $mol_wgt_sum / $considered_count;
$mol_wgt_avg = sprintf "%.2f", $mol_wgt_avg;
open FILE, ">$outdir/average_molecular weight.txt" or die $!; 
print FILE "Average moleculat weight = $mol_wgt_avg\n";
close FILE;        

#######################################################################################
## Draw protein and peptide graphics using position information from uniprot database,
## label using ASAP ratios and populate in HTML table. Not all proteins are added to 
## the table, a filtering step has been taken.  
## Filtering constraints:
##		- Xpress and ASAP ratios agreement: Xpress/ASAP should be in between 0.5 and 2
## 		- probability should be greater than ProteinProphet cutoff (default = 0.90)
##		- xpress and asap ratios in between 0.33 and 3.0
#######################################################################################

sub addProteinToTable {
	my( $t, $protein_tag)= @_;
	my $content;
	my $uniprot_id;
	my $protein_name_id = $protein_tag->{'att'}->{'protein_name'};
	if($protein_name_id =~ m/^sp/ || $protein_name_id =~ m/^tr/){
		@_ = split('\|', $protein_name_id);
		$uniprot_id = $_[1];	
	}
	elsif($protein_name_id =~ m/^IPI/) {
		my @protein_annotation_tag =  $protein_tag->getElementsByTagName("annotation");
		my $uniprot_ipi = $protein_annotation_tag[0]->{'att'}->{'protein_description'};
		if($uniprot_ipi =~ /^|SWISS/){
				my @uni_ipi_name = split('\|', $uniprot_ipi);
				my $swissprot_id = $uni_ipi_name[1];
				my @uni_swiss_name = split(':', $swissprot_id);
				my $uniprot_id_pre = $uni_swiss_name[1];
				$uniprot_id = substr($uniprot_id_pre, 0, 6);
			}
		else {
			print "Error: Unidentified Swissprot/Uniprot identifier<br><br>\n\n"; 
		}
	}
	
	#print $uniprot_id,"<br>\n";
	my $probability = $protein_tag->{'att'}->{'probability'};
	my @ASAPRatio =  $protein_tag->getElementsByTagName("ASAPRatio");

	my @XPressRatio_tags =  $protein_tag->getElementsByTagName("XPressRatio");
	my $prot_xpress_ratio = $XPressRatio_tags[0]->{'att'}->{'ratio_mean'};
	my $prot_xpress_sd = $XPressRatio_tags[0]->{'att'}->{'ratio_standard_dev'};

	my @ASAPRatio_pvalue_tags =  $protein_tag->getElementsByTagName("ASAPRatio_pvalue");
	my $prot_adj_asap_ratio = $ASAPRatio_pvalue_tags[0]->{'att'}->{'adj_ratio_mean'};
	my $prot_adj_asap_std_dev = $ASAPRatio_pvalue_tags[0]->{'att'}->{'adj_ratio_standard_dev'};

	## return when there is no Xpress or ASAP tag ##
	if(!$prot_adj_asap_ratio || !$prot_xpress_ratio){
		print LOG "No ASAP or Xpress tags found for '$uniprot_id'\n";
		return;
	}
	
	## give a minimum real positive value to the ratio when it is zero ##
	$prot_adj_asap_ratio = 0.01 if($prot_adj_asap_ratio == 0);
	$prot_xpress_ratio = 0.01 if($prot_xpress_ratio == 0);
	$prot_xpress_ratio = $prot_xpress_ratio/$correction_factor;
	
	## invert ratios (H<==>L) if the invert option is active ## 			
	if($invert_ratios){
		$prot_adj_asap_ratio= 1/$prot_adj_asap_ratio;
		$prot_xpress_ratio = 1/$prot_xpress_ratio;
	}

	if(!$prot_adj_asap_ratio|| !$prot_xpress_ratio){
		$discarded_count++; 
		return;
	}
	
	## quotient for Xpress and ASAP agreement calculation ##
	my $quotient = $prot_xpress_ratio/$prot_adj_asap_ratio;
	my $prot_adj_asap_ratio_log;

	## Filter the protein when validating with Xpress Ratio##
	if($validate_asap_xpress){
	if (($probability > $threshold) && ((($quotient < 2) || ($quotient > 0.5)) ||
		(($prot_adj_asap_ratio >= 3) && ($prot_xpress_ratio >= 3)) || 
		(($prot_adj_asap_ratio <= 0.33) && ($prot_xpress_ratio <= 0.33)) 
		)){      
		my $prot_adj_asap_std_dev_ratio = $prot_adj_asap_std_dev/$prot_adj_asap_ratio;
		eval{
			$prot_adj_asap_ratio_log = log2($prot_adj_asap_ratio);
			$prot_adj_asap_std_dev = abs($prot_adj_asap_std_dev_ratio*$prot_adj_asap_ratio_log);
			push(@prot_ratio_means_filtered, $prot_adj_asap_ratio_log);			
		};
		if($@){
			print LOG "Warning: protein $uniprot_id is discarded. ", $@,"\n";
			$discarded_count++;
			return;
		}
		
		$prot_adj_asap_ratio_log = sprintf "%.2f", $prot_adj_asap_ratio_log;
		$prot_adj_asap_std_dev = sprintf "%.2f", $prot_adj_asap_std_dev;		
		
		## create the UniProtFileParser object ##
		my $uniprot_obj = new UniProtFileParser();
		$uniprot_obj->id($uniprot_id);
		## fetch and store the current protein information from uniprot database to a scalar ##    
		my $content = "";		
		eval {
			my $up_url = "http://www.uniprot.org/uniprot/$uniprot_id.txt";
			$content = get("$up_url");
		};
		if($@){
			print LOG "Warning: protein $uniprot_id is discarded.  Could not fetch information from uniprot database.\n";
			$discarded_count++;
			return;  
		}
		if($content eq ""){
			print LOG "Warning: protein $uniprot_id is discarded. It is not found in Uniprot database.\n";
			$discarded_count++;   				
   			return;
		}
		
		## parse the scalar containing uniprot protein information ## 
		$uniprot_obj->parseContent($content);
					
		my $protein_name = $uniprot_obj->rec_name;

		if(!$protein_name){
			$protein_name = $uniprot_obj->sub_name;			
			print LOG "No recommended name found for $uniprot_id\n";
		}
	    my $protein_sequence = $uniprot_obj->seq;
		my $protein_length = length($protein_sequence);
		
		## create two panels for images: small panel for the preview in the table, big one for the detailed view ##
		my $small_panel = createPanel($protein_name, $protein_length, 700);
		my $big_panel = createPanel($protein_name, $protein_length, 1200);
			    
	    my $signal = "NO"; my $transmembrane = "NO";
	    my @uniprot_features = $uniprot_obj->features;
	    foreach my $feature(@uniprot_features){
	    	if($feature->{'name'} eq "SIGNAL"){
	    		$signal = "YES";
	    	}
	    	if($feature->{'name'} eq "TRANSMEM"){
	    		$transmembrane = "YES";
	    	}
	    }
	    
	    ## get the molecular weight of the proteins
	    $molecular_weight = $uniprot_obj->mol_wt;
	    $mol_wgt_sum += $molecular_weight;
	    
	    ## add protein and peptides information to the panel ##  
		addTrackToPanel($big_panel, $uniprot_id, $protein_tag, $uniprot_obj);
		addTrackToPanel($small_panel, $uniprot_id, $protein_tag, $uniprot_obj);
		
		## add uniprot database features such as transmembrane, extra and intra cellular regions ##
	    if($transmembrane eq 'YES'){
				addUPFeaturesTrack($big_panel, $uniprot_obj);
				addUPFeaturesTrack($small_panel, $uniprot_obj);
	    }

		## create image map for the whole protein, link to the uniprot database entry ##
		my ($url,$map,$mapname) = $big_panel->image_and_map(
		                            -root => "$outdir",
	                                -url  => "temp",
	                                -link => sub {
	                                      my $feature = shift;
	                                      my $name = $feature->display_name;
	                                      return "http://www.uniprot.org/uniprot/$uniprot_id";
	                                 });

		open BIGIMG,">$outdir/images/$uniprot_id.png";
		binmode(BIGIMG);
		print BIGIMG $big_panel->png;
		close BIGIMG;	    
		$big_panel->finished;
		
		open SMALLIMG,">$outdir/small_images/$uniprot_id.png";
		binmode(SMALLIMG);
		print SMALLIMG $small_panel->png;
		close SMALLIMG;	    
		$small_panel->finished;

		open MAP, ">$outdir/index_files/$uniprot_id.html";
		print MAP "<html>\n";
		print MAP "<head>\n";
		print MAP "<style type=\"text/css\">\n";
		print MAP ".col1{width: 35%;}\n.col2{width:65%;}\n.col2 div{width:100%; overflow:auto;}";
		print MAP "</style>\n";
		print MAP "</head>\n";
		print MAP "<body style=\"font-family: 'Times New Roman' !important; font-size: 12px !important;\">\n";
		print MAP "<img id=\"image\" src=\"../images/$uniprot_id.png\" alt=\"../images/$uniprot_id.png\" USEMAP='#$mapname' BORDER='0'/>\n";
		print MAP "$map";
		print MAP "<img style=\"margin: 0px 30px;\" src=\"../css/ratio_color_legend.png\" height=\"50\" width=\"500\" />\n";
		if($transmembrane eq "YES"){
			print MAP "<img style=\"position:relative; top: 0px;margin: 0px 0px;\" src=\"../css/UP_feat_legend.png\" height=\"18\" width=\"650\" />\n";
		}
		print MAP "<br/>\n";
		print MAP "<a href=\"../index.html#$uniprot_id\">back to main page</a>\n";
		print MAP "<a href=\"#image\" style=\"position:fixed; bottom: 0px; right:0px;width: 100px; \">to the top</a>\n";
		print MAP "<br/>\n";
		print MAP "<br/>\n";
		print MAP "<div width=1300><tr>\n";
		if($invert_ratios){
			addInvertPeptideTables($protein_tag, $uniprot_obj, $uniprot_id);
		}
		else{
			addPeptideTables($protein_tag, $uniprot_obj, $uniprot_id);
		}
		print MAP "</div>\n";
		print MAP "</body>\n";
		print MAP "</html>\n";
		close MAP;
		
		my $icr = "NA";
		my $ecr = "NA";		
		if($transmembrane eq "YES"){
			$icr = calculateDomainAvgRatio($uniprot_obj, "Cytoplasmic");
			$ecr = calculateDomainAvgRatio($uniprot_obj, "Extracellular");
		}
		## add protein information to the main table ##
		print HTML "\t<tr>\n";
		print HTML "\t\t<td class=\"upid\" id=\"$uniprot_id\"><div>$uniprot_id</div></td>\n";
		print HTML "\t\t<td class=\"recname\"><div>$protein_name</div></td>\n";
		print HTML "\t\t<td class=\"molwt\"><div>$molecular_weight</div></td>\n";		
		print HTML "\t\t<td class=\"protratio\"><div>$prot_adj_asap_ratio_log</div></td>\n";
		print HTML "\t\t<td class=\"protratiostddev\"><div>$prot_adj_asap_std_dev</div></td>\n";
		print HTML "\t\t<td class=\"signalpep\"><div>$signal</div></td>\n";
		print HTML "\t\t<td class=\"tm\"><div>$transmembrane</div></td>\n";
		print HTML "\t\t<td class=\"icr\"><div>$icr</div></td>\n";
		print HTML "\t\t<td class=\"ecr\"><div>$ecr</div></td>\n";
		print HTML "\t\t<td class=\"plot\"><div><a href=\"index_files/$uniprot_id.html\" title=\"Click to view details\"><img src=\"small_images/$uniprot_id.png\" alt=\"small_images/$uniprot_id.png\"/></a></div></td>\n";
		print HTML "\t</tr>\n";
		$considered_count ++;
	}
		else{
		print LOG "Warning: protein $uniprot_id is discarded. Failed in filtering.\n";
		$discarded_count++;
	}
}
		## Without the filtering step when validation with XPress ratio is not required ##

	elsif(!$validate_asap_xpress){
		if ($probability > $threshold){
		my $prot_adj_asap_std_dev_ratio = $prot_adj_asap_std_dev/$prot_adj_asap_ratio;
			eval{
				$prot_adj_asap_ratio_log = log2($prot_adj_asap_ratio);
				$prot_adj_asap_std_dev = abs($prot_adj_asap_std_dev_ratio*$prot_adj_asap_ratio_log);
				push(@prot_ratio_means_filtered, $prot_adj_asap_ratio_log);			
			};
			if($@){
				print LOG "Warning: protein $uniprot_id is discarded. ", $@,"\n";
				$discarded_count++;
				return;
			}
		
			$prot_adj_asap_ratio_log = sprintf "%.2f", $prot_adj_asap_ratio_log;
			$prot_adj_asap_std_dev = sprintf "%.2f", $prot_adj_asap_std_dev;		
			
			## create the UniProtFileParser object ##
			my $uniprot_obj = new UniProtFileParser();
			$uniprot_obj->id($uniprot_id);
			## fetch and store the current protein information from uniprot database to a scalar ##    
			my $content = "";
			eval {
				my $up_url = "http://www.uniprot.org/uniprot/$uniprot_id.txt";
				$content = get("$up_url");
			};
			if($@){
				print LOG "Warning: protein $uniprot_id is discarded.  Could not fetch information from uniprot database.\n";
				$discarded_count++;
				return;
			}
			if($content eq ""){
				print LOG "Warning: protein $uniprot_id is discarded. It is not found in Uniprot database.\n";
				$discarded_count++;
				return;
			}
		
			## parse the scalar containing uniprot protein information ## 
				$uniprot_obj->parseContent($content);	
			my $protein_name = $uniprot_obj->rec_name;

			if(!$protein_name){
				$protein_name = $uniprot_obj->sub_name;			
				print LOG "No recommended name found for $uniprot_id\n";
			}
			my $protein_sequence = $uniprot_obj->seq;
			my $protein_length = length($protein_sequence);
		
			## create two panels for images: small panel for the preview in the table, big one for the detailed view ##
			my $small_panel = createPanel($protein_name, $protein_length, 700);
			my $big_panel = createPanel($protein_name, $protein_length, 1200);
			    
			my $signal = "NO"; my $transmembrane = "NO";
			my @uniprot_features = $uniprot_obj->features;
			foreach my $feature(@uniprot_features){
				if($feature->{'name'} eq "SIGNAL"){
					$signal = "YES";
				}
				if($feature->{'name'} eq "TRANSMEM"){
					$transmembrane = "YES";
				}
			}
	    
	    ## get the molecular weight of the proteins
	    $molecular_weight = $uniprot_obj->mol_wt;
	    $mol_wgt_sum += $molecular_weight;
	    
	    ## add protein and peptides information to the panel ##  
		addTrackToPanel($big_panel, $uniprot_id, $protein_tag, $uniprot_obj);
		addTrackToPanel($small_panel, $uniprot_id, $protein_tag, $uniprot_obj);
		
		## add uniprot database features such as transmembrane, extra and intra cellular regions ##
	    if($transmembrane eq 'YES'){
				addUPFeaturesTrack($big_panel, $uniprot_obj);
				addUPFeaturesTrack($small_panel, $uniprot_obj);
	    }

		## create image map for the whole protein, link to the uniprot database entry ##
		my ($url,$map,$mapname) = $big_panel->image_and_map(
		                            -root => "$outdir",
	                                -url  => "temp",
	                                -link => sub {
	                                      my $feature = shift;
	                                      my $name = $feature->display_name;
	                                      return "http://www.uniprot.org/uniprot/$uniprot_id";
	                                 });

		open BIGIMG,">$outdir/images/$uniprot_id.png";
		binmode(BIGIMG);
		print BIGIMG $big_panel->png;
		close BIGIMG;	    
		$big_panel->finished;
		
		open SMALLIMG,">$outdir/small_images/$uniprot_id.png";
		binmode(SMALLIMG);
		print SMALLIMG $small_panel->png;
		close SMALLIMG;	    
		$small_panel->finished;

		open MAP, ">$outdir/index_files/$uniprot_id.html";
		print MAP "<html>\n";
		print MAP "<head>\n";
		print MAP "<style type=\"text/css\">\n";
		print MAP ".col1{width: 35%;}\n.col2{width:65%;}\n.col2 div{width:100%; overflow:auto;}";
		print MAP "</style>\n";
		print MAP "</head>\n";
		print MAP "<body style=\"font-family: 'Times New Roman' !important; font-size: 12px !important;\">\n";
		print MAP "<img id=\"image\" src=\"../images/$uniprot_id.png\" alt=\"../images/$uniprot_id.png\" USEMAP='#$mapname' BORDER='0'/>\n";
		print MAP "$map";
		print MAP "<img style=\"margin: 0px 30px;\" src=\"../css/ratio_color_legend.png\" height=\"50\" width=\"500\" />\n";
		if($transmembrane eq "YES"){
			print MAP "<img style=\"position:relative; top: 0px;margin: 0px 0px;\" src=\"../css/UP_feat_legend.png\" height=\"18\" width=\"650\" />\n";
		}
		print MAP "<br/>\n";
		print MAP "<a href=\"../index.html#$uniprot_id\">back to main page</a>\n";
		print MAP "<a href=\"#image\" style=\"position:fixed; bottom: 0px; right:0px;width: 100px; \">to the top</a>\n";
		print MAP "<br/>\n";
		print MAP "<br/>\n";
		print MAP "<div width=1300><tr>\n";
		if($invert_ratios){
			addInvertPeptideTables($protein_tag, $uniprot_obj, $uniprot_id);
		}
		else{
			addPeptideTables($protein_tag, $uniprot_obj, $uniprot_id);
		}
		print MAP "</div>\n";
		print MAP "</body>\n";
		print MAP "</html>\n";
		close MAP;
		
		my $icr = "NA";
		my $ecr = "NA";		
		if($transmembrane eq "YES"){
			$icr = calculateDomainAvgRatio($uniprot_obj, "Cytoplasmic");
			$ecr = calculateDomainAvgRatio($uniprot_obj, "Extracellular");
		}
		## add protein information to the main table ##
		print HTML "\t<tr>\n";
		print HTML "\t\t<td class=\"upid\" id=\"$uniprot_id\"><div>$uniprot_id</div></td>\n";
		print HTML "\t\t<td class=\"recname\"><div>$protein_name</div></td>\n";
		print HTML "\t\t<td class=\"molwt\"><div>$molecular_weight</div></td>\n";		
		print HTML "\t\t<td class=\"protratio\"><div>$prot_adj_asap_ratio</div></td>\n";
		print HTML "\t\t<td class=\"protratiostddev\"><div>$prot_adj_asap_std_dev</div></td>\n";
		print HTML "\t\t<td class=\"signalpep\"><div>$signal</div></td>\n";
		print HTML "\t\t<td class=\"tm\"><div>$transmembrane</div></td>\n";
		print HTML "\t\t<td class=\"icr\"><div>$icr</div></td>\n";
		print HTML "\t\t<td class=\"ecr\"><div>$ecr</div></td>\n";
		print HTML "\t\t<td class=\"plot\"><div><a href=\"index_files/$uniprot_id.html\" title=\"Click to view details\"><img src=\"small_images/$uniprot_id.png\" alt=\"small_images/$uniprot_id.png\"/></a></div></td>\n";
		print HTML "\t</tr>\n";
		$considered_count ++;
	}
		else{
		print LOG "Warning: protein $uniprot_id is discarded. Failed in filtering.\n";
		$discarded_count++;
	}
}

	$t->purge;
	return;
}


################################################################################################
## This subroutine creates the tables with eaech peptide information in the detailed view page
################################################################################################
sub addPeptideTables{
	my $protein_tag = $_[0];
	my $uniprot_obj = $_[1];
	my $uniprot_id = $_[2];
	
	my @asap_seq_tags = $protein_tag->getElementsByTagName("ASAP_Seq");
	my @peptide_tags = $protein_tag->getElementsByTagName("peptide");

	my %pp_prob_scores = ();
	foreach my $peptide_tag (@peptide_tags){
		my $peptide_sequence = $peptide_tag->{'att'}->{'peptide_sequence'};
		my $pp_prob_score = $peptide_tag->{'att'}->{'initial_probability'};
		push(@{$pp_prob_scores{$peptide_sequence}},$pp_prob_score);
	}	
	my $peptide_nr = 1;
	my $discarded_pep = 0;
	my @peptides;
	foreach my $asap_seq_tag (@asap_seq_tags){
		my $include_parameter = $asap_seq_tag->{'att'}->{'include'};
		if($include_parameter == 1 && !$elaborate_peptides){ 
			my $peptide_hash;
			my $peptide_sequence = $asap_seq_tag->{'att'}->{'light_sequence'};
			my $pep_ratio_mean = $asap_seq_tag->{'att'}->{'ratio_mean'};
			my $pep_std_dev = $asap_seq_tag->{'att'}->{'ratio_standard_dev'};
			my @pp_prob_scores = @{$pp_prob_scores{$peptide_sequence}};
			my $nr_of_times_identified = scalar @pp_prob_scores;
			my $pp_prob_score = max(@pp_prob_scores);
			my ($pep_ratio_mean_log, $pep_std_dev_log);
			$pep_ratio_mean = 0.01 if($pep_ratio_mean == 0);
			$pep_ratio_mean = $pep_ratio_mean/$correction_factor;
			if($invert_ratios){
				$pep_ratio_mean= 1/$pep_ratio_mean;
			}
			eval {
				my $pep_sd_mean_ratio = $pep_std_dev/$pep_ratio_mean; 
				$pep_ratio_mean_log = log2($pep_ratio_mean);
				push(@pep_ratio_means_filtered, $pep_ratio_mean_log);
				$pep_std_dev_log = abs($pep_sd_mean_ratio*$pep_ratio_mean_log);
			};
			if($@){
				print LOG "Warning at '",$uniprot_id,"' peptide $peptide_nr: ",$@;
				$discarded_pep++;
				$peptide_nr++;
				return;
			}
			else{
				my $peptide_length = length($peptide_sequence);
				my $peptide_start = index($uniprot_obj->seq, $peptide_sequence)+1;
				my $peptide_end = $peptide_start+$peptide_length-1;
				$pep_ratio_mean_log = sprintf "%.2f", $pep_ratio_mean_log;
				$pep_std_dev_log = sprintf "%.2f", $pep_std_dev_log;
			
				$peptide_hash->{"id"} = $peptide_nr;
				$peptide_hash->{"start"} = $peptide_start;
				$peptide_hash->{"end"} = $peptide_end;
				$peptide_hash->{"ratio"} = $pep_ratio_mean;
				$peptide_hash->{"sd"} = $pep_std_dev;
				push @peptides, $peptide_hash;
				
				## create and write peptide information to table ##
				print MAP "<table id=\"peptide_$peptide_nr\" table border=\"1\"  cellpadding=\"5\" width=\"500\" style=\"border-collapse:collapse; margin:10px 5px; float:left;\">\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Peptide nr</td>\n\t\t<td class=\"col2\">$peptide_nr</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Start</td>\n\t\t<td class=\"col2\">$peptide_start</td></tr>\n"; 
				print MAP "\t<tr>\t\t<td class=\"col1\">End</td>\n\t\t<td class=\"col2\">$peptide_end</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Length</td>\n\t\t<td class=\"col2\">$peptide_length</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Ratio (log<sub>2</sub> transformed)</td>\n\t\t<td class=\"col2\">$pep_ratio_mean_log &plusmn $pep_std_dev_log</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">nr. of times identified</td>\n\t\t<td class=\"col2\">$nr_of_times_identified</td></tr>\n";
				if($nr_of_times_identified > 1){
					print MAP "\t<tr>\t\t<td class=\"col1\">Highest PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
				}
				else{
					print MAP "\t<tr>\t\t<td class=\"col1\">PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
				}
				print MAP "\t<tr>\t\t<td class=\"col1\">Sequence</td>\n\t\t<td class=\"col2\"><div>$peptide_sequence</div></td></tr>\n";
				print MAP "</table>\n";
				}
			$peptide_nr++;
		}
		if($include_parameter == 1 && $elaborate_peptides){ 
			my $peptide_hash;
			my $peptide_sequence = $asap_seq_tag->{'att'}->{'light_sequence'};
			my $pep_ratio_mean = $asap_seq_tag->{'att'}->{'ratio_mean'};
			my $pep_std_dev = $asap_seq_tag->{'att'}->{'ratio_standard_dev'};
			my @pp_prob_scores = @{$pp_prob_scores{$peptide_sequence}};
			my $nr_of_times_identified = scalar @pp_prob_scores;
			foreach my $pp_prob_score (@pp_prob_scores) {
				my ($pep_ratio_mean_log, $pep_std_dev_log);
				$pep_ratio_mean = 0.01 if($pep_ratio_mean == 0);
				$pep_ratio_mean = $pep_ratio_mean/$correction_factor;
				if($invert_ratios){
					$pep_ratio_mean= 1/$pep_ratio_mean;
				}
				eval {
					my $pep_sd_mean_ratio = $pep_std_dev/$pep_ratio_mean; 
					$pep_ratio_mean_log = log2($pep_ratio_mean);
					push(@pep_ratio_means_filtered, $pep_ratio_mean_log);
					$pep_std_dev_log = abs($pep_sd_mean_ratio*$pep_ratio_mean_log);
				};
				if($@){
					print LOG "Warning at '",$uniprot_id,"' peptide $peptide_nr: ",$@;
					$discarded_pep++;
					$peptide_nr++;
					return;
				}
				else{
					my $peptide_length = length($peptide_sequence);
					my $peptide_start = index($uniprot_obj->seq, $peptide_sequence)+1;
					my $peptide_end = $peptide_start+$peptide_length-1;
					$pep_ratio_mean_log = sprintf "%.2f", $pep_ratio_mean_log;
					$pep_std_dev_log = sprintf "%.2f", $pep_std_dev_log;
			
					$peptide_hash->{"id"} = $peptide_nr;
					$peptide_hash->{"start"} = $peptide_start;
					$peptide_hash->{"end"} = $peptide_end;
					$peptide_hash->{"ratio"} = $pep_ratio_mean;
					$peptide_hash->{"sd"} = $pep_std_dev;
					push @peptides, $peptide_hash;
				
					## create and write peptide information to table ##
					print MAP "<table id=\"peptide_$peptide_nr\" table border=\"1\"  cellpadding=\"5\" width=\"500\" style=\"border-collapse:collapse; margin:10px 5px; float:left;\">\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Peptide nr</td>\n\t\t<td class=\"col2\">$peptide_nr</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Start</td>\n\t\t<td class=\"col2\">$peptide_start</td></tr>\n"; 
					print MAP "\t<tr>\t\t<td class=\"col1\">End</td>\n\t\t<td class=\"col2\">$peptide_end</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Length</td>\n\t\t<td class=\"col2\">$peptide_length</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Ratio (log<sub>2</sub> transformed)</td>\n\t\t<td class=\"col2\">$pep_ratio_mean_log &plusmn $pep_std_dev_log</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">nr. of times identified</td>\n\t\t<td class=\"col2\">$nr_of_times_identified</td></tr>\n";
					if($nr_of_times_identified > 1){
						print MAP "\t<tr>\t\t<td class=\"col1\">Highest PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
					}
					else{
						print MAP "\t<tr>\t\t<td class=\"col1\">PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
					}
					print MAP "\t<tr>\t\t<td class=\"col1\">Sequence</td>\n\t\t<td class=\"col2\"><div>$peptide_sequence</div></td></tr>\n";
					print MAP "</table>\n";
					}
				$peptide_nr++;
			}
		}
	}
#	print LOG $uniprot_id,"\t",$peptide_nr,"\t",$discarded_pep,"\n";
	
	$uniprot_obj->{_identifiedPeptides}= \@peptides;
}

sub addInvertPeptideTables{
	my $protein_tag = $_[0];
	my $uniprot_obj = $_[1];
	my $uniprot_id = $_[2];
	
	my @asap_seq_tags = $protein_tag->getElementsByTagName("ASAP_Seq");
	my @peptide_tags = $protein_tag->getElementsByTagName("peptide");

	my %pp_prob_scores = ();
	foreach my $peptide_tag (@peptide_tags){
		my $peptide_sequence = $peptide_tag->{'att'}->{'peptide_sequence'};
		my $pp_prob_score = $peptide_tag->{'att'}->{'initial_probability'};
		push(@{$pp_prob_scores{$peptide_sequence}},$pp_prob_score);
	}	
	my $peptide_nr = 1;
	my $discarded_pep = 0;
	my @peptides;
	foreach my $asap_seq_tag (@asap_seq_tags){
		my $include_parameter = $asap_seq_tag->{'att'}->{'include'};
		if($include_parameter == 1 && !$elaborate_peptides){ 
			my $peptide_hash;
			my $peptide_sequence = $asap_seq_tag->{'att'}->{'light_sequence'};
			my $pep_ratio_mean = $asap_seq_tag->{'att'}->{'ratio_mean'};
			my $pep_std_dev = $asap_seq_tag->{'att'}->{'ratio_standard_dev'};
			my @pp_prob_scores = @{$pp_prob_scores{$peptide_sequence}};
			my $nr_of_times_identified = scalar @pp_prob_scores;
			my $pp_prob_score = max(@pp_prob_scores);
			my ($pep_ratio_mean_log, $pep_std_dev_log);
			$pep_ratio_mean = 0.01 if($pep_ratio_mean == 0);
			$pep_ratio_mean = $pep_ratio_mean/$correction_factor;
			my $pep_ratio_mean_inv= 1/$pep_ratio_mean;
			eval {
				my $pep_sd_mean_ratio = $pep_std_dev/$pep_ratio_mean_inv; 
				$pep_ratio_mean_log = log2($pep_ratio_mean_inv);
				push(@pep_ratio_means_filtered, $pep_ratio_mean_log);
				$pep_std_dev_log = abs($pep_sd_mean_ratio*$pep_ratio_mean_log);
			};
			if($@){
				print LOG "Warning at '",$uniprot_id,"' peptide $peptide_nr: ",$@;
				$discarded_pep++;
				$peptide_nr++;
				return;
			}
			else{
				my $peptide_length = length($peptide_sequence);
				my $peptide_start = index($uniprot_obj->seq, $peptide_sequence)+1;
				my $peptide_end = $peptide_start+$peptide_length-1;
				$pep_ratio_mean_log = sprintf "%.2f", $pep_ratio_mean_log;
				$pep_std_dev_log = sprintf "%.2f", $pep_std_dev_log;
			
				$peptide_hash->{"id"} = $peptide_nr;
				$peptide_hash->{"start"} = $peptide_start;
				$peptide_hash->{"end"} = $peptide_end;
				$peptide_hash->{"ratio"} = $pep_ratio_mean;
				$peptide_hash->{"sd"} = $pep_std_dev;
				push @peptides, $peptide_hash;
				
				## create and write peptide information to table ##
				print MAP "<table id=\"peptide_$peptide_nr\" table border=\"1\"  cellpadding=\"5\" width=\"500\" style=\"border-collapse:collapse; margin:10px 5px; float:left;\">\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Peptide nr</td>\n\t\t<td class=\"col2\">$peptide_nr</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Start</td>\n\t\t<td class=\"col2\">$peptide_start</td></tr>\n"; 
				print MAP "\t<tr>\t\t<td class=\"col1\">End</td>\n\t\t<td class=\"col2\">$peptide_end</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Length</td>\n\t\t<td class=\"col2\">$peptide_length</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">Ratio (log<sub>2</sub> transformed)</td>\n\t\t<td class=\"col2\">$pep_ratio_mean_log &plusmn $pep_std_dev_log</td></tr>\n";
				print MAP "\t<tr>\t\t<td class=\"col1\">nr. of times identified</td>\n\t\t<td class=\"col2\">$nr_of_times_identified</td></tr>\n";
				if($nr_of_times_identified > 1){
					print MAP "\t<tr>\t\t<td class=\"col1\">Highest PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
				}
				else{
					print MAP "\t<tr>\t\t<td class=\"col1\">PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
				}
				print MAP "\t<tr>\t\t<td class=\"col1\">Sequence</td>\n\t\t<td class=\"col2\"><div>$peptide_sequence</div></td></tr>\n";
				print MAP "</table>\n";
				}
			$peptide_nr++;
		}
		if($include_parameter == 1 && $elaborate_peptides){ 
			my $peptide_hash;
			my $peptide_sequence = $asap_seq_tag->{'att'}->{'light_sequence'};
			my $pep_ratio_mean = $asap_seq_tag->{'att'}->{'ratio_mean'};
			my $pep_std_dev = $asap_seq_tag->{'att'}->{'ratio_standard_dev'};
			my @pp_prob_scores = @{$pp_prob_scores{$peptide_sequence}};
			my $nr_of_times_identified = scalar @pp_prob_scores;
			foreach my $pp_prob_score (@pp_prob_scores) {
				my ($pep_ratio_mean_log, $pep_std_dev_log);
				$pep_ratio_mean = 0.01 if($pep_ratio_mean == 0);
				$pep_ratio_mean = $pep_ratio_mean/$correction_factor;
				my $pep_ratio_mean_inv= 1/$pep_ratio_mean;
				eval {
					my $pep_sd_mean_ratio = $pep_std_dev/$pep_ratio_mean_inv; 
					$pep_ratio_mean_log = log2($pep_ratio_mean_inv);
					push(@pep_ratio_means_filtered, $pep_ratio_mean_log);
					$pep_std_dev_log = abs($pep_sd_mean_ratio*$pep_ratio_mean_log);
				};
				if($@){
					print LOG "Warning at '",$uniprot_id,"' peptide $peptide_nr: ",$@;
					$discarded_pep++;
					$peptide_nr++;
					return;
				}
				else{
					my $peptide_length = length($peptide_sequence);
					my $peptide_start = index($uniprot_obj->seq, $peptide_sequence)+1;
					my $peptide_end = $peptide_start+$peptide_length-1;
					$pep_ratio_mean_log = sprintf "%.2f", $pep_ratio_mean_log;
					$pep_std_dev_log = sprintf "%.2f", $pep_std_dev_log;
			
					$peptide_hash->{"id"} = $peptide_nr;
					$peptide_hash->{"start"} = $peptide_start;
					$peptide_hash->{"end"} = $peptide_end;
					$peptide_hash->{"ratio"} = $pep_ratio_mean;
					$peptide_hash->{"sd"} = $pep_std_dev;
					push @peptides, $peptide_hash;
				
					## create and write peptide information to table ##
					print MAP "<table id=\"peptide_$peptide_nr\" table border=\"1\"  cellpadding=\"5\" width=\"500\" style=\"border-collapse:collapse; margin:10px 5px; float:left;\">\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Peptide nr</td>\n\t\t<td class=\"col2\">$peptide_nr</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Start</td>\n\t\t<td class=\"col2\">$peptide_start</td></tr>\n"; 
					print MAP "\t<tr>\t\t<td class=\"col1\">End</td>\n\t\t<td class=\"col2\">$peptide_end</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Length</td>\n\t\t<td class=\"col2\">$peptide_length</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">Ratio (log<sub>2</sub> transformed)</td>\n\t\t<td class=\"col2\">$pep_ratio_mean_log &plusmn $pep_std_dev_log</td></tr>\n";
					print MAP "\t<tr>\t\t<td class=\"col1\">nr. of times identified</td>\n\t\t<td class=\"col2\">$nr_of_times_identified</td></tr>\n";
					if($nr_of_times_identified > 1){
						print MAP "\t<tr>\t\t<td class=\"col1\">Highest PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
					}
					else{
						print MAP "\t<tr>\t\t<td class=\"col1\">PP prob. score	</td>\n\t\t<td class=\"col2\">$pp_prob_score</td></tr>\n";				
					}
					print MAP "\t<tr>\t\t<td class=\"col1\">Sequence</td>\n\t\t<td class=\"col2\"><div>$peptide_sequence</div></td></tr>\n";
					print MAP "</table>\n";
					}
				$peptide_nr++;
			}
		}
	}
#	print LOG $uniprot_id,"\t",$peptide_nr,"\t",$discarded_pep,"\n";
	
	$uniprot_obj->{_identifiedPeptides}= \@peptides;
}

#############################################################################################
## Creates a track with protein information 
## Adds each peptide as a Bio SeqFeature object with appropriate color (by ratio as score)  
############################################################################################# 
sub addTrackToPanel{
	my $panel = $_[0];
	my $uniprot_id = $_[1];
	my $protein_tag = $_[2];
	my $uniprot_obj = $_[3];
	
	my @asap_seq_tags = $protein_tag->getElementsByTagName("ASAP_Seq");
	my @peptide_tags = $protein_tag->getElementsByTagName("peptide");

	my %pp_prob_scores = ();
	foreach my $peptide_tag (@peptide_tags){
		my $peptide_sequence = $peptide_tag->{'att'}->{'peptide_sequence'};
		my $pp_prob_score = $peptide_tag->{'att'}->{'initial_probability'};
		push(@{$pp_prob_scores{$peptide_sequence}},$pp_prob_score);
	}	
	my $peptide_nr = 1;
	my $discarded_pep = 0;

	## create protein track with protein name and anchored to uniprot database entry 
	my $track = $panel->add_track(
							  -glyph    => 'generic',
                              -label    => 1,
                              -link		=> "#peptide_".'$name',
                              -font	=> 'gdLargeFont',
                              -font2color  => 'black',
                              -bgcolor	=> sub {
                              					getColor(shift);
                              				},
                              -description => sub {
                              					my $feature = shift;
                                				my $score   = $feature->score;
                                				return "$score";
                               				}           
                             );
              
	foreach my $asap_seq_tag (@asap_seq_tags){
		my $include_parameter = $asap_seq_tag->{'att'}->{'include'};
		if($include_parameter == 1 && !$elaborate_peptides){ 
			my $peptide_sequence = $asap_seq_tag->{'att'}->{'light_sequence'};
			my $ratio_mean = $asap_seq_tag->{'att'}->{'ratio_mean'};
			my $standard_deviation = $asap_seq_tag->{'att'}->{'ratio_standard_dev'};
			my ($ratio_mean_log, $standard_deviation_log);
			$ratio_mean = 0.01 if($ratio_mean == 0);
			$ratio_mean = $ratio_mean/$correction_factor;
			if($invert_ratios){
				$ratio_mean= 1/$ratio_mean;
			}
			eval {
				my $sd_mean_ratio = $standard_deviation/$ratio_mean;
				$ratio_mean_log = log2($ratio_mean);
				$standard_deviation_log = abs($ratio_mean_log*$sd_mean_ratio);
			};
			if($@){
				$peptide_nr++;
				return;
			}
			else{
				my $peptide_length = length($peptide_sequence);
				my $peptide_start = index($uniprot_obj->seq, $peptide_sequence)+1;
				my $peptide_end = $peptide_start+$peptide_length-1;
		   		
				$ratio_mean_log = sprintf "%.2f", $ratio_mean_log;
				$standard_deviation_log = sprintf "%.2f", $standard_deviation_log;
				my $feat;
				## create a Bio::SeqFeature object using peptide coordinates and ratios
				if($panel->width() > 1000){
					$feat = new Bio::SeqFeature::Generic(-display_name => $peptide_nr,	
														-start	=> $peptide_start,
														-end	=> $peptide_end-1,
														-score	=> $ratio_mean_log,
														);
				}
				## display only start, end and ratio in the small panel			
				else{
					$feat = new Bio::SeqFeature::Generic(
														-start	=> $peptide_start,
														-end	=> $peptide_end-1,
														-score	=> $ratio_mean_log,
														);
				}
				$track->add_feature($feat);
			}
			$peptide_nr++;	
		}
		elsif($include_parameter == 1 && $elaborate_peptides){ 
			my $peptide_sequence = $asap_seq_tag->{'att'}->{'light_sequence'};
			my $ratio_mean = $asap_seq_tag->{'att'}->{'ratio_mean'};
			my $standard_deviation = $asap_seq_tag->{'att'}->{'ratio_standard_dev'};
			my ($ratio_mean_log, $standard_deviation_log);
			my @pp_prob_scores = @{$pp_prob_scores{$peptide_sequence}};
			my $nr_of_times_identified = scalar @pp_prob_scores;
			foreach my $pp_prob_score (@pp_prob_scores) {
				$ratio_mean = 0.01 if($ratio_mean == 0);
				$ratio_mean = $ratio_mean/$correction_factor;
				if($invert_ratios){
					$ratio_mean= 1/$ratio_mean;
				}
				eval {
					my $sd_mean_ratio = $standard_deviation/$ratio_mean;
					$ratio_mean_log = log2($ratio_mean);
					$standard_deviation_log = abs($ratio_mean_log*$sd_mean_ratio);
				};
				if($@){
					$peptide_nr++;
					return;
				}
				else{
					my $peptide_length = length($peptide_sequence);
					my $peptide_start = index($uniprot_obj->seq, $peptide_sequence)+1;
					my $peptide_end = $peptide_start+$peptide_length-1;
		   		
					$ratio_mean_log = sprintf "%.2f", $ratio_mean_log;
					$standard_deviation_log = sprintf "%.2f", $standard_deviation_log;
					my $feat;
					## create a Bio::SeqFeature object using peptide coordinates and ratios
					if($panel->width() > 1000){
						$feat = new Bio::SeqFeature::Generic(-display_name => $peptide_nr,	
															-start	=> $peptide_start,
															-end	=> $peptide_end-1,
															-score	=> $ratio_mean_log,
															);
					}	
					## display only start, end and ratio in the small panel			
					else{
						$feat = new Bio::SeqFeature::Generic(
															-start	=> $peptide_start,
															-end	=> $peptide_end-1,
															-score	=> $ratio_mean_log,
															);
					}
					$track->add_feature($feat);
				}
				$peptide_nr++;
			}	
		}	
	}
}

#######################################
## Adds uniprot features to the track
#######################################
 
sub addUPFeaturesTrack {
	my $panel = $_[0];
	my $uniprot_obj = $_[1];
	
	my $uniprot_id = $uniprot_obj->id;
	my @uniprot_features = $uniprot_obj->features;
	my $track = $panel->add_track(
                              -font	=> 'gdLargeFont',
                              -glyph     => 'graded_segments',
                              -font2color  => 'red',
                              -link		=> "http://www.uniprot.org/uniprot/$uniprot_id#section_features",
                              -bgcolor	=> sub {
                              				my $feature = shift;
                              				if($feature->display_name eq "TRANSMEM"){
                              					return '#00FFFF';
                              				}
                              				else{
                              					my $color = "#FF00FF";
												if ($feature->has_tag('description')) {
										        	for my $desc ($feature->get_tag_values('description')){
									         			if($desc =~ m/Cytoplasmic/){
										         			$color = "#000000";
										         		}
										         		elsif($desc =~ m/Lumenal/){
										         			$color = "#777777";
										         		}
										         		elsif($desc =~ m/Extracellular/){
										         			$color = "#CCCCCC";
										         		}
										         		elsif($desc =~ m/Mitochondrial/){
										         			$color = "#FEA8FF";
										         		}
										         		elsif($desc =~m/Nuclear|Perinuclear space/){
										         			$color = "#A8BFFF";
										         		}
												      	else{
												      		print LOG "New feature detected for '$uniprot_id'\n";
												      	}                              					
										         	}
										      	}
                              					return $color;
                              				}
                                  		}
                             		);	    	
		                             
   	foreach my $feature(@uniprot_features){
   		if($feature->{'name'} =~ m/TRANSMEM|TOPO_DOM/){
	        my $feat = Bio::SeqFeature::Generic->new(
                                             -display_name => $feature->{'name'},
                                             -start => $feature->{'start'}+1,
                                             -end   => $feature->{'end'}-1,
                                             -tag          => {
                                                         description => $feature->{'desc'}
                                                     }
                                            );
            $track-> add_feature($feat);                               		    		    		    			
   		}
   	}
}

################################################################
## Automatic generation of color shades by altering RGB values
################################################################
sub populateColorHash{
	my $red_shade_count = 0;
	for (my $r=128, my $g = 0, my $b=0; $r<256;){
		$r=$r+10;
		$g=$g+20;
		$b=$b+20;	
		$red_shade_count++;			
	}
	for (my $r=128, my $g = 0, my $b=0; $r<256;){
		my $red_shade = rgbToHex($r,$g,$b);
		$colorHash{"red$red_shade_count"} = $red_shade;
		$r=$r+10;
		$g=$g+20;
		$b=$b+20;
		$red_shade_count--;			
	}
	my $green_shade_count = 0;
	for (my $r=0, my $g = 128, my $b=0; $g<256;){
		$r=$r+20;
		$g=$g+10;
		$b=$b+20;
		$green_shade_count++;			
	}
	for (my $r=0, my $g = 128, my $b=0; $g<256;){
		my $green_shade = rgbToHex($r,$g,$b);
		$colorHash{"green$green_shade_count"} = $green_shade;
		$r=$r+20;
		$g=$g+10;
		$b=$b+20;
		$green_shade_count--;			
	}
	$colorHash{"white"} = "#FFFFFF";
	return;
}

#########################################################################
## subroutine converts from RGB values to hexadecimal(HTML) color codes 
#########################################################################
sub rgbToHex {
    my $red=$_[0];
    my $green=$_[1];
    my $blue=$_[2];
    my $hexcolor=sprintf("#%2.2X%2.2X%2.2X",$red,$green,$blue);
    return ($hexcolor);
}

#############################################################
## Computes the in which the peptide falls using the ratio 
#############################################################
sub computeBins {
	my $global_min = $_[0];
	my $global_max = $_[1];
	my $red_bin_size = $zn-$global_min;
	$red_bin_size = $red_bin_size/12; 
	$red_bin_size = $red_bin_size + 0.00001;
	my $green_bin_size = $global_max-$zp;
	$green_bin_size = $green_bin_size/12;
	$green_bin_size = $green_bin_size+0.00001;
	
	my $bin_count = 1;
	my $green_bin_end;
	for (my $green_bin_start = $zp; $green_bin_start <= $global_max-$green_bin_size;){
		$green_bin_end = $green_bin_start+$green_bin_size;
		$green_bins_hash{"$green_bin_start\t$green_bin_end"} = "green$bin_count";
		$green_bin_start = $green_bin_start+$green_bin_size;
		$bin_count++;
	}
	my $p_inf = 9**9**9;
	$green_bins_hash{"$green_bin_end\t$p_inf"} = "green$bin_count";
	
	$bin_count = 1;
	my $red_bin_end;
	for (my $red_bin_start = $zn; $red_bin_start >= $global_min+$red_bin_size;){
		$red_bin_end = $red_bin_start-$red_bin_size;
		$red_bins_hash{"$red_bin_end\t$red_bin_start"} = "red$bin_count";
		$red_bin_start = $red_bin_start-$red_bin_size;
		$bin_count++;
	}
	my $n_inf = -9**9**9;
	$red_bins_hash{"$n_inf\t$red_bin_end"} = "red$bin_count";		
}

##########################################################################################
## Returns the color for the peptide by looking up the in which bin its score falls into  
##########################################################################################
sub getColor {
	my $feature = shift;
    my $score   = $feature->score;

	if($score >= $zn && $score <= $zp){
		return $colorHash{'white'};
	}	
	if($score > $zp){
		foreach my $bin (keys %green_bins_hash){
			my ($start, $end) = split('\t', $bin);
			if($score >= $start && $score <= $end){
				return $colorHash{$green_bins_hash{$bin}};				
			}
		}
	}
 	if($score < $zn){
		foreach my $bin (keys %red_bins_hash){
			my ($start, $end) = split('\t', $bin);
			if($score >= $start && $score <= $end){
				return $colorHash{$red_bins_hash{$bin}};				
			}
		}
	}
}

#######################################################################
## Domain specific average ratio calculator: intra and extra cellular 
#######################################################################
sub calculateDomainAvgRatio {
	my $uniprot_obj = $_[0];
	my $domain = $_[1];
	
	my @peptides = $uniprot_obj->identifiedPeptides();
	my @uniprot_features = $uniprot_obj->features();
   	
  	
   	my @domain_peptide_ratios = (); 
   	my @domain_peptide_sds = ();
   	foreach my $feature (@uniprot_features){
   		if($feature->{'desc'} =~ m/$domain/){
   			my $feat_start = $feature->{'start'};
   			my $feat_end = $feature->{'end'};
			foreach my $peptide_hash (@peptides){
				my $peptide_start = $peptide_hash->{'start'};
				my $peptide_end = $peptide_hash->{'end'};
				if($peptide_start >= $feat_start && $peptide_end <= $feat_end){
					push(@domain_peptide_ratios, $peptide_hash->{'ratio'});
					push(@domain_peptide_sds, $peptide_hash->{'sd'});
				}
		   	}		 
   		}
	}
	if($#domain_peptide_ratios >= 0){
		my $domain_ratios_avg = sum(@domain_peptide_ratios)/($#domain_peptide_ratios+1);
		my $domain_sds_avg = sum(@domain_peptide_sds)/($#domain_peptide_sds+1);
		my ($domain_ratios_avg_log, $domain_sds_avg_log);
		$domain_ratios_avg = 0.01 if($domain_ratios_avg == 0);
		if($invert_ratios){
			$domain_ratios_avg = 1/$domain_ratios_avg;
			## invert sd also?
			##$domain_sds_avg = 1/$domain_sds_avg;
		}
		my $domain_sd_mean_ratio = $domain_sds_avg/$domain_ratios_avg;
		eval {		
			$domain_ratios_avg_log = log2($domain_ratios_avg);
			$domain_sds_avg_log = abs($domain_sd_mean_ratio*$domain_ratios_avg_log);
		};
		if($@){
			print LOG "Cannot calculate $domain ratio of '",$uniprot_obj->id,"': ",$@,"<br>\n";
			return "NA";
		}
		$domain_ratios_avg_log = sprintf "%.2f", $domain_ratios_avg_log;
		$domain_sds_avg_log = sprintf "%.2f", $domain_sds_avg_log;
		
		## return "$domain_ratios_avg_log \x{B1} $domain_sds_avg_log";##
		return "$domain_ratios_avg_log";
	}
	else{
		return "NA";
	}
}

#############################################################

sub createPanel {
	my $protein_name = $_[0];
	my $protein_length = $_[1];
	my $panel_width = $_[2];
	
	my $padding = 0.025*$panel_width;
	
	my $panel = Bio::Graphics::Panel->new(
										-length    => $protein_length,
										-pad_left  => $padding,
										-pad_right => $padding,
										-pad_top	=> $padding,
										-pad_bottom  => $padding,		                                      
										-grid		 => 0,
										-width     => $panel_width,
										-key_style => 'between'
										);
	my $pos_ref_bar = Bio::SeqFeature::Generic->new(
										-display_name => $protein_name,
										-start => 1,
										-end   => $protein_length,
										-description => "$protein_name"
										);
	$panel->add_track($pos_ref_bar,
		                  -glyph   => 'arrow',
		                  -tick    => 2,
		                  -fgcolor => 'black',
		                  -double  => 1,
		                  -label   => 1,
		                  -font	=> 'gdLargeFont',
		                 );
	return $panel;
}

#############################################################

sub log2 {
	my $n = shift;
    return log($n)/log(2);
}



#############################################################

sub plotNormalDistribution{
	my $twig_obj = $_[0];
	my @asap_seq_tags = $twig_obj->getElementsByTagName("ASAP_Seq");
	
}

#############################################################

sub plotFCValues{
	mkdir "$outdir/fc_plots" if (! -d "$outdir/fc_plots");
	## open R handler in perl, load required R libraries, create pdf output and setup page styles 
	open(R, "| R --vanilla --slave");
	print R "library(\"grid\")\n";
	print R "library(\"ggplot2\")\n";
	print R "pdf(\"$outdir/fc_plots/fc_plots.pdf\", width=15, height=10, title=\"fc-values plots\")\n";
	print R "vplayout <- function(x, y) { viewport(layout.pos.row = x, layout.pos.col = y) }\n";
	print R "pushViewport(viewport(layout = grid.layout(2, 3)))\n";
	
	## Calculate the minimum and maximum fc-values, these values are the coordinate limits in all plots
	my $min_fc = min(@prot_ratio_means_all);
	my $max_fc = max(@prot_ratio_means_all);
	open FC, ">$outdir/fc_plots/protein_ratio_means_all.out" or die $!;
	print FC "fc_value\n";
	foreach (@prot_ratio_means_all){
		print FC $_,"\n";
	}
	close FC;
	print R "fc<-read.table(\"$outdir/fc_plots/protein_ratio_means_all.out\", header=T)\n";
	print R "p1<-qplot(fc_value, data=fc, geom=\"bar\", xlab=\"protein ratios all\")+coord_cartesian(xlim = c($min_fc, $max_fc))\n";
	
	## Calculate the minimum and maximum fc-values, these values are the coordinate limits in all plots
	$min_fc = min(@prot_adj_ratio_means_all);
	$max_fc = max(@prot_adj_ratio_means_all);
	open FC, ">$outdir/fc_plots/protein_adj_ratio_means_all.out" or die $!;
	print FC "fc_value\n";
	foreach (@prot_adj_ratio_means_all){
		print FC $_,"\n";
	}
	close FC;
	print R "fc<-read.table(\"$outdir/fc_plots/protein_adj_ratio_means_all.out\", header=T)\n";
	print R "p2<-qplot(fc_value, data=fc, geom=\"bar\", xlab=\"protein adj ratios all\")+coord_cartesian(xlim = c($min_fc, $max_fc))\n";
	
	## Calculate the minimum and maximum fc-values, these values are the coordinate limits in all plots
	$min_fc = min(@pep_ratio_means_all);
	$max_fc = max(@pep_ratio_means_all);
	open FC, ">$outdir/fc_plots/peptide_ratio_means_all.out" or die $!;
	print FC "fc_value\n";
	foreach (@pep_ratio_means_all){
		print FC $_,"\n";
	}
	close FC;
	print R "fc<-read.table(\"$outdir/fc_plots/peptide_ratio_means_all.out\", header=T)\n";
	print R "p3<-qplot(fc_value, data=fc, geom=\"bar\", xlab=\"peptide ratios all\")+coord_cartesian(xlim = c($min_fc, $max_fc))\n";



	## Calculate the minimum and maximum fc-values, these values are the coordinate limits in all plots
	$min_fc = min(@prot_ratio_means_filtered);
	$max_fc = max(@prot_ratio_means_filtered);
	open FC, ">$outdir/fc_plots/protein_ratio_means_filtered.out" or die $!;
	print FC "fc_value\n";
	foreach (@prot_ratio_means_filtered){
		print FC $_,"\n";
	}
	close FC;
	print R "fc<-read.table(\"$outdir/fc_plots/protein_ratio_means_filtered.out\", header=T)\n";
	print R "p4<-qplot(fc_value, data=fc, geom=\"bar\", xlab=\"protein ratios filtered\")+coord_cartesian(xlim = c($min_fc, $max_fc))\n";
	
	## Calculate the minimum and maximum fc-values, these values are the coordinate limits in all plots
	$min_fc = min(@prot_adj_ratio_means_filtered);
	$max_fc = max(@prot_adj_ratio_means_filtered);
	open FC, ">$outdir/fc_plots/protein_adj_ratio_means_filtered.out" or die $!;
	print FC "fc_value\n";
	foreach (@prot_adj_ratio_means_filtered){
		print FC $_,"\n";
	}
	close FC;
	print R "fc<-read.table(\"$outdir/fc_plots/protein_adj_ratio_means_filtered.out\", header=T)\n";
	print R "p5<-qplot(fc_value, data=fc, geom=\"bar\", xlab=\"protein adj ratios filtered\")+coord_cartesian(xlim = c($min_fc, $max_fc))\n";
	
	## Calculate the minimum and maximum fc-values, these values are the coordinate limits in all plots
	$min_fc = min(@pep_ratio_means_filtered);
	$max_fc = max(@pep_ratio_means_filtered);
	open FC, ">$outdir/fc_plots/peptide_ratio_means_filtered.out" or die $!;
	print FC "fc_value\n";
	foreach (@pep_ratio_means_filtered){
		print FC $_,"\n";
	}
	close FC;
	print R "fc<-read.table(\"$outdir/fc_plots/peptide_ratio_means_filtered.out\", header=T)\n";
	print R "p6<-qplot(fc_value, data=fc, geom=\"bar\", xlab=\"peptide ratios filtered\")+coord_cartesian(xlim = c($min_fc, $max_fc))\n";
	
	
	print R "print(p1, vp=vplayout(1,1))\n";
	print R "print(p2, vp=vplayout(1,2))\n";
	print R "print(p3, vp=vplayout(1,3))\n";
	print R "print(p4, vp=vplayout(2,1))\n";
	print R "print(p5, vp=vplayout(2,2))\n";
	print R "print(p6, vp=vplayout(2,3))\n";			
				
	print R "dev.off()\n";
	return;
}
