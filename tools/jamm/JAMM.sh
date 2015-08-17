########################################################################
# JAMMv1.0.7rev1 is a peak finder for joint analysis of NGS replicates.
# Copyright (C) 2014-2015  Mahmoud Ibrahim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: mahmoud.ibrahim@mdc-berlin.de
########################################################################


##Finding out the path
sPath="`dirname \"$0\"`"
sPath="`( cd \"$sPath\" && pwd )`"



usage()
{
cat << EOF
Welcome to JAMM v1.0.7rev1 (GNU GPLv3). Copyright (C) 2014-2015  Mahmoud Ibrahim.

This program comes with ABSOLUTELY NO WARRANTY; for details visit http://www.gnu.org/licenses/gpl.html. This is free software, and you are welcome to redistribute it under certain conditions; visit http://www.gnu.org/licenses/gpl.html for details.

OPTIONS:
   -s      directory containing Sample files (required)
   -g      Genome size file (required)
   -o      Output directory (required)
   -c      directory containing input or Control files
   -f      Fragment length(s) (default: estimated)
   -r      Resolution, peak or region or window (default: peak)
   -m      Mode, normal or narrow (default: normal)
   -i      clustering Initialization window selection, deterministic or stochastic (default: deterministic)
   -b	   Bin Size (default: estimated)
   -w      minimum Window size (default: 2 --- Note: this means minimum_window_size = bin_size x the_value_of_-w)
   -e	   window Enrichment cutoff, auto or any numeric value (default: 1 --- Set this to "auto" to estimate the window enrichment cutoff)
   -d	   keep PCR Dupicates in single-end mode, y or n (default: n --- if -t is "paired", this option has no effect) 
   -t	   Type, single or paired (default: single, requires BED files. paired requires BEDPE files)
   -p	   Number of processors used by R scripts (default: 1)

EOF
}


# ========================= 
# Process Input parameters
# =========================


#Defaults -- Change those if you want
mode="normal"
resol="peak"
cores="1"
window="2"
type="single"
windowe="1"
initModel="deterministic"
uniq="n"

#Defaults -- Do not change
sdir=""
gsize=""
out=""
binsize="ns"
fraglen="ns"
ran=$RANDOM
wdir=$(mktemp -d)
export LANG=C #locale defaults
export LC_ALL=C #locale defaults

while getopts "s:g:o:c:m:r:f:p:w:b:t:e:i:d:" OPTION
do
	case $OPTION in
	s) sdir=$OPTARG
	;;
	g) gsize=$OPTARG
	;;
	o) out=$OPTARG
	;;	
	c) bdir=$OPTARG
	;;
	m) mode=$OPTARG
	;;
	r) resol=$OPTARG
	;;
	f) fraglen=$OPTARG
	;;
	p) cores=$OPTARG
	;;
	w) window=$OPTARG
	;;
	b) binsize=$OPTARG
	;;
	t) type=$OPTARG
	;;
	e) windowe=$OPTARG
	;;
	i) initModel=$OPTARG
	;;
	d) uniq=$OPTARG
	;;
	?)
	usage
	exit
	;;
	esac
done
if [ "$mode" == "normal" ]; then
	clustno="2"
fi
if [ "$mode" == "narrow" ]; then
	clustno="3"
fi

if [[ -z $sdir ]] || [[ -z $gsize ]] || [[ -z $out ]]
then
     usage
     exit 1
fi
if [[ -d "$out/peaks" ]]; then
	printf "\n\nOutput directory $out/peaks already exists. I can't override existing results!\n\n"
	exit 0
fi
if [ $fraglen == "ns" ]; then
	if [[ -d "$out/xcorr" ]]; then
		printf "\n\nOutput directory $out/xcorr already exists. I can't override existing results!\n\n"
		exit 0
	fi
fi
#=======================> DONE!




# ============================= 
# Step One: Initial Processing
# =============================
printf "\n\n============================================\nStarted JAMM Pipeline v1.0.7rev1...Hang on!\n============================================\n\n"

if [ ! -d "$wdir" ]; then
	mkdir $wdir #make working directory
fi
if [ ! -d "$out" ]; then
	mkdir $out #make output directory
fi
mkdir $wdir/bkgd.$ran/ #directory to store background files
mkdir $wdir/sizes.$ran/ #chromosomes and sizes
mkdir $wdir/samples.$ran/ #store sample files

dupnum=$(ls -1 $sdir | wc -l) #count how many sample files


#separate chromosome sizes
printf "Loading genome size file..."
ext="$wdir/sizes.$ran/"
awk -v ext="$ext" '{ print >> ext"/size." $1 ".bed" }' $gsize
printf "Done!\n"


printf "Processing sample files..."
#load each chromosome from each sample file
for i in $sdir/*.bed; do
samplefile=$(basename $i)	
	for f in $wdir/sizes.$ran/*; do
		sizefile=$(basename $f)
		chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
		awk -v chr="$chr" -v ext="$wdir/samples.$ran/" -v samplefile="$samplefile" -F"\t" '$1 == chr { print $2"\t"$6 >> ext"sample."chr"."samplefile }' "$i" 
	done
done
printf "Done!\n"


if [ ! -z $bdir ]; then
#concatenate all background files into one file
printf "Processing control files..."
cat $bdir/*.bed > $wdir/bkgd.$ran/ctrl.bed

for f in $wdir/sizes.$ran/*; do
	sizefile=$(basename $f)
	chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
	awk -v chr="$chr" -v ext="$wdir/bkgd.$ran/" -F"\t" '$1 == chr { print $2"\t"$6 >> ext"bkgd."chr".ctrl.bed" }' "$wdir/bkgd.$ran/ctrl.bed"
done

printf "Done!\n"
fi

#determine average read lengths
printf "Getting average read lengths...\n"
readL=""
if [ ! -z $bdir ]; then
	readC=$(awk '{a=$3-$2;print a;}' "$wdir/bkgd.$ran/ctrl.bed" | perl -lane '$a+=$_;END{print $a/$.}' | awk '{a=$1+0.5;print a;}' | cut -d"." -f1)
	printf "Control: $readC\n"
fi
readL=""
for s in $sdir/*.bed; do #and for each sample file
	file=$(basename $s)
	samplefile=$(echo $file | awk -F"." '{print $1}');	
	read=$(awk '{a=$3-$2;print a;}' "$s" | perl -lane '$a+=$_;END{print $a/$.}' | awk '{a=$1+0.5;print a;}' | cut -d"." -f1)
	printf "$samplefile: $read\n"
	readL="$readL,$read"
done
readL=${readL#","}
#=======================> DONE!



# ============================= 
# Step Two: Fragment Length
# =============================
#single-end
if [ $type == "single" ]; then

if [ $fraglen == "ns" ]; then

	##Counting Where Reads Start and Calculating Cross Correlation
	mkdir $wdir/stats.$ran/ #store count files
	mkdir $out/xcorr #final xcorr results


	printf "Calculating Fragment Length(s)...\n"
	for f in $wdir/sizes.$ran/*; do #for each chromosome
		samplelist=""
		readlist=""
		sizefile=$(basename $f)
		chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
		
		#list of sample bed files and read lengths
		for s in $wdir/samples.$ran/*.bed; do #and for each sample file
			samplefile=$(basename $s)
			chr2=$(echo $samplefile | awk -F"." '{print $2}');
			if [ $chr == $chr2 ] #belonging to this chromosome
			then
				samplelist="$samplelist,$wdir/samples.$ran/$samplefile"
			fi
		done
		readlist="$readL"
		
		#list of control bed files and read lengths
		if [ ! -z $bdir ]; then
			for s in $wdir/bkgd.$ran/*.bed; do #and for each sample file
				samplefile=$(basename $s)
				chr2=$(echo $samplefile | awk -F"." '{print $2}');
				if [ $chr == $chr2 ] #belonging to this chromosome
				then
					samplelist="$samplelist,$wdir/bkgd.$ran/$samplefile"
					readlist="$readL,$readC"
				fi
			done		
		fi
		
		#remove leading comma
		samplelist=${samplelist#","}
		
		#call R script for xcorr calculation
		Rscript "$sPath/xcorr.r" -ibed="$samplelist" -s="$wdir/sizes.$ran/size.$chr.bed"  -rl="$readlist" -d="$wdir/stats.$ran" -p="$cores"	
	done

	#report xcorr results (samples)
	for f in $sdir/*.bed; do
		file=$(basename $f)
		samplefile=$(echo $file | awk -F"." '{print $1}');	
		mkdir "$out/xcorr/$samplefile" #final xcorr results
		if [ -f "$wdir/stats.$ran/xc.$samplefile.tab" ]; then
			cp $wdir/stats.$ran/xc.$samplefile.tab $out/xcorr/$samplefile/shifts.txt	
		fi
		Rscript "$sPath/xcorrhelper.r" -infile="$out/xcorr/$samplefile/shifts.txt" -out="$out/xcorr/$samplefile"
	done
	#report xcorr results (control)
	if [ ! -z $bdir ]; then
		mkdir "$out/xcorr/ctrl" #final xcorr results
		if [ -f "$wdir/stats.$ran/xc.ctrl.tab" ]; then
			cp $wdir/stats.$ran/xc.ctrl.tab $out/xcorr/ctrl/shifts.txt
		fi
		Rscript "$sPath/xcorrhelper.r" -infile="$out/xcorr/ctrl/shifts.txt" -out="$out/xcorr/ctrl"
	fi

fi
fi

#paired-end
if [ $type == "paired" ]; then
	printf "Getting Average Fragment Length(s)...\n"
	mkdir $out/xcorr #final xcorr results
	
	for f in $sdir/*.bed; do
		file=$(basename $f)
		samplefile=$(echo $file | awk -F"." '{print $1}');	
		mkdir "$out/xcorr/$samplefile"
		frag=$(awk '{a=$6-$2;print a;}' $f | perl -lane '$a+=$_;END{print $a/$.}' | awk '{a=$1+0.5;print a;}' | cut -d"." -f1)
		echo "Average_from_paired	$frag" > $out/xcorr/$samplefile/shifts.txt
		Rscript "$sPath/xcorrhelper.r" -infile="$out/xcorr/$samplefile/shifts.txt" -out="$out/xcorr/$samplefile"
	done
fi
#=======================> DONE!



# ================================= 
# Step Three: Calculating Bin Size
# =================================
if [ $binsize == "ns" ]; then
	printf "Getting Bin Size: "

	chr=$(sort -nr -k2 $gsize | head -n 1 | awk -F"\t" '{print $1}');
	samplelist=""
	frag=""
	if [ $fraglen != "ns" ]; then
		frag=$fraglen
		k=1
	fi
	
	#list of sample bed files and read lengths
	for s in $wdir/samples.$ran/*.bed; do
		samplefile=$(basename $s)
		chr2=$(echo $samplefile | awk -F"." '{print $2}');
		if [ $chr == $chr2 ]
		then
			samplelist="$samplelist,$wdir/samples.$ran/$samplefile"
			samplename=$(echo $samplefile | awk -F"." '{ print $3 }')
			samplefilename=$(echo $samplefile | cut -d'.' -f 3-)
			if [ $fraglen == "ns" ]; then
				shift=$(awk -F":" '$1 == "Fragment Length" { print $2 }' "$out/xcorr/$samplename/xcorrsummary.txt")
				frag="$frag,$shift"
			fi
		fi
	done
	#remove leading comma
	samplelist=${samplelist#","}
	frag=${frag#","}
	Rscript "$sPath/bincalculator.r" -ibed="$samplelist" -s="$gsize" -rl="$readL" -d="$wdir" -p="$cores"	-f="$frag" -type="$type"
fi
if [ $binsize != "ns" ]; then
	printf "You set a Bin Size: $binsize \n"
fi
#=======================> DONE!



# =========================== 
# Step Four: Calling Peaks
# ===========================
mkdir $wdir/peaks.$ran/ #store count files
mkdir $out/peaks #store peak files

printf "Calling Peaks...(mode: $mode, resolution: $resol)\n"


#single-end reads
if [ $type == "single" ]; then

if [ $binsize == "ns" ]; then
	binsize=$(cat "$wdir/binsize.txt")
fi

counting=1;			
for f in $wdir/sizes.$ran/*; do #for each chromosome
		samplelist=""
		frag=""
		k=1
		if [ $fraglen != "ns" ]; then
			frag=$fraglen
		fi

		
		sizefile=$(basename $f)
		chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');

		printf "Chromosome $chr: "
		
		#list of sample bed files and fragment lengths
		for s in $wdir/samples.$ran/*.bed; do #and for each sample file
			samplefile=$(basename $s)
			chr2=$(echo $samplefile | awk -F"." '{print $2}');
			if [ $chr == $chr2 ] #belonging to this chromosome
			then
				samplelist="$samplelist,$wdir/samples.$ran/ext.$samplefile"
				samplename=$(echo $samplefile | awk -F"." '{ print $3 }')
				samplefilename=$(echo $samplefile | cut -d'.' -f 3-)
				if [ $fraglen == "ns" ]; then
					shift=$(awk -F":" '$1 == "Fragment Length" { print $2 }' "$out/xcorr/$samplename/xcorrsummary.txt")
					frag="$frag,$shift"
					read=$(echo $readL | cut -f "$k" -d ",")
					k=$(($k+1))
				fi
				if [ $fraglen != "ns" ]; then
					shift=$(echo $frag | cut -f "$k" -d ",")	
					read=$(echo $readL | cut -f "$k" -d ",")
					k=$(($k+1))
				fi
				if [ $uniq == "y" ]; then
					perl "$sPath/readshifter.pl" "$wdir/samples.$ran/$samplefile" $shift $read > "$wdir/samples.$ran/ext.$samplefile"
				fi
				if [ $uniq == "n" ]; then
					perl "$sPath/readshifter.pl" "$wdir/samples.$ran/$samplefile" $shift $read | uniq > "$wdir/samples.$ran/ext.$samplefile"
				fi				
			fi
		done
		
		#control file
		bkgdfile="None"
		if [ ! -z $bdir ]; then
			if [ $fraglen == "ns" ]; then
				bshift=$(awk -F":" '$1 == "Fragment Length" { print $2 }' "$out/xcorr/ctrl/xcorrsummary.txt")
				frag="$frag,$bshift"
			fi
			if [ $fraglen != "ns" ]; then
				l=$(($dupnum+1))
				bshift=$(echo $frag | cut -f "$l" -d ",")
			fi
			
			if [ $uniq == "y" ]; then
				perl "$sPath/readshifter.pl" "$wdir/bkgd.$ran/bkgd.$chr.ctrl.bed" $bshift $readC > "$wdir/bkgd.$ran/ext.bkgd.$chr.ctrl.bed"
			fi
			if [ $uniq == "n" ]; then
				perl "$sPath/readshifter.pl" "$wdir/bkgd.$ran/bkgd.$chr.ctrl.bed" $bshift $readC | uniq > "$wdir/bkgd.$ran/ext.bkgd.$chr.ctrl.bed"
			fi
			bkgdfile="$wdir/bkgd.$ran/ext.bkgd.$chr.ctrl.bed"
		fi
		
		#remove leading comma
		samplelist=${samplelist#","}
		frag=${frag#","}
		
		#call the peak calling R script
		Rscript "$sPath/peakfinder.r" -sfile="$f" -bednames="$samplelist" -frag="$frag" -bkgd=$bkgdfile -out="$wdir/peaks.$ran/" -clustnummer="$clustno" -resolution="$resol" -window="$window" -p="$cores" -bin="$binsize" -type="$type" -chrcount="$counting" -initModel="$initModel" -windowe="$windowe"
		counting=$(($counting+1));
		cp "$wdir/peaks.$ran/$chr.peaks.bed" "$out/peaks/$chr.peaks.bed"
		rm "$wdir/peaks.$ran/$chr.peaks.bed"
done
counting=1;
fi


#paired-end reads
if [ $type == "paired" ]; then

if [ $binsize == "ns" ]; then
	binsize=$(cat "$wdir/binsize.txt")
fi
			

counting=1;
for f in $wdir/sizes.$ran/*; do #for each chromosome
	samplelist=""
	
	sizefile=$(basename $f)
	chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');

	printf "Chromosome $chr: "
		
	#list of sample bed files and fragment lengths
	for s in $wdir/samples.$ran/*.bed; do #and for each sample file
		samplefile=$(basename $s)
		chr2=$(echo $samplefile | awk -F"." '{print $2}');
		if [ $chr == $chr2 ] #belonging to this chromosome
		then
			samplelist="$samplelist,$wdir/samples.$ran/$samplefile"
			samplename=$(echo $samplefile | awk -F"." '{ print $3 }')
			samplefilename=$(echo $samplefile | cut -d'.' -f 3-)
			x="$sdir/$samplefilename"
		fi
	done
		
	#control file
	bkgdfile="None"
	if [ ! -z $bdir ]; then
		bkgdfile="$wdir/bkgd.$ran/bkgd.$chr.ctrl.bed"
	fi
		
	#remove leading comma
	samplelist=${samplelist#","}
	frag=${frag#","}
		
	#call the peak calling R script
	Rscript "$sPath/peakfinder.r" -sfile=$f -bednames=$samplelist -frag="NA" -bkgd=$bkgdfile -out="$wdir/peaks.$ran/" -clustnummer="$clustno" -resolution="$resol" -window="$window" -p="$cores" -bin="$binsize" -type="$type" -chrcount="$counting" -initModel="$initModel" -windowe="$windowe"
	counting=$(($counting+1));
	cp "$wdir/peaks.$ran/$chr.peaks.bed" "$out/peaks/$chr.peaks.bed"
	rm "$wdir/peaks.$ran/$chr.peaks.bed"
done
counting=1;
fi


cp $wdir/peaks.$ran/min.peaksize $out/peaks/min.peaksize


#concatenate, sort and filter
cat $out/peaks/*.bed > $out/peaks/all.narrowPeak 
if [[ -s  $out/peaks/all.narrowPeak ]]; then
	Rscript "$sPath/peakhelper.r" -filelist="$out/peaks/all.narrowPeak"
	perl "$sPath/peakfilter.pl" $out/peaks/all.narrowPeak | sort -nr -k7 > $out/peaks/filtered.peaks.narrowPeak
	cut -f1-10 $out/peaks/all.narrowPeak | awk -F"\t" -v j=0 '$7 > j' | sort -nr -k7 > $out/peaks/all.peaks.narrowPeak
fi
rm $out/peaks/all.narrowPeak
rm $out/peaks/*.bed
rm $out/peaks/min.peaksize
#=======================> DONE!



rm -rf $wdir



printf "\n\n========================================\nWe're done...Congratulations!\n========================================\n\n"
