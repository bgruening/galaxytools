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
Welcome to JAMM v1.0.7rev1 Signal Generator Script (GNU GPLv3). Copyright (C) 2014-2015  Mahmoud Ibrahim.

This program comes with ABSOLUTELY NO WARRANTY; for details visit http://www.gnu.org/licenses/gpl.html. This is free software, and you are welcome to redistribute it under certain conditions; visit http://www.gnu.org/licenses/gpl.html for details.

OPTIONS:
   -s      Directory containing sample files (required)
   -g      Genome size file (required)
   -o      Output Directory (required)
   -c      directory containing input or Control files
   -r 	   file with Regions to get signal for (required)
   -b      Bin size for signal generation (default: 10)
   -f      Fragment lengths (required if -t is "single")
   -p	   Number of processors used by R scripts (default: 1)
   -t	   Alignment type, paired or single (default: single)
 
EOF
}


# ========================= 
# Process Input parameters
# =========================


sdir=""
gsize=""
out=""
signal="10"
regs=""
fraglen=""
wdir=$(mktemp -d)
ran=$RANDOM
cores="1"
type="single"

while getopts "s:o:c:r:b:f:g:p:t:" OPTION
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
	r) regs=$OPTARG
	;;
	b) signal=$OPTARG
	;;
	f) fraglen=$OPTARG
	;;
	p) cores=$OPTARG
	;;
	t) type=$OPTARG
	;;
	?)
	usage
	exit
	;;
	esac
done


if [[ -z $sdir ]] || [[ -z $regs ]] || [[ -z $out ]]
then
     usage
     exit 1
fi

if [[ -z $fraglen ]]
then
	if [ $type == "single" ]
	then
     usage
     exit 1
	fi
fi
if [[ -d "$out/signal" ]]; then
	printf "\n\nOutput directory $out/signal already exists. I can't override existing results!\n\n"
	exit 0
fi
#=======================> DONE!



# ============================= 
# Step One: Initial Processing
# =============================
printf "\n\n=======================================================\nStarted JAMM1.0.7rev1 Signal Generator Pipeline...Hang on!\n=======================================================\n\n"

if [ ! -d "$wdir" ]; then
	mkdir $wdir #make working directory
fi
if [ ! -d "$out" ]; then
	mkdir $out #make working directory
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
if [ $type == "single" ]; then
printf "Getting average read lengths..."
	if [ ! -z $bdir ]; then
		readC=$(awk '{a=$3-$2;print a;}' "$wdir/bkgd.$ran/ctrl.bed" | perl -lane '$a+=$_;END{print $a/$.}' | awk '{print int($1)}')
	fi
	readL=""
	for s in $sdir/*.bed; do #and for each sample file
		read=$(awk '{a=$3-$2;print a;}' "$s" | perl -lane '$a+=$_;END{print $a/$.}' | awk '{print int($1)}')
		readL="$readL,$read"
	done
	readL=${readL#","}
printf "Done!\n"
fi
#=======================> DONE!






# =========================== 
# Step Three: Getting Signal
# ===========================
mkdir $wdir/signal.$ran/
mkdir $out/signal #store signal


printf "Generating Signal...(bin size: $signal)\n"

if [ $type == "single" ]; then
counting=1;			
for f in $wdir/sizes.$ran/*; do #for each chromosome
		samplelist=""
		frag=""
		frag=$fraglen
		k=1
		
		sizefile=$(basename $f)
		chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
		chrSize=$(cat $f | cut -f2);
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
				shift=$(echo "$frag" | cut -f "$k" -d ",")
				read=$(echo "$readL" | cut -f "$k" -d ",")	
				k=$(($k+1))
				perl "$sPath/readshifter.pl" "$wdir/samples.$ran/$samplefile" $shift $read > "$wdir/samples.$ran/ext.$samplefile"
			fi
		done
		
		#control file
		bkgdfile="None"
		if [ ! -z $bdir ]; then
			l=$(($dupnum+1))
			bshift=$(echo $frag | cut -f "$l" -d ",")
			perl "$sPath/readshifter.pl" "$wdir/bkgd.$ran/bkgd.$chr.ctrl.bed" $bshift $readC > "$wdir/bkgd.$ran/ext.bkgd.$chr.ctrl.bed"
			bkgdfile="$wdir/bkgd.$ran/ext.bkgd.$chr.ctrl.bed"
		fi
		
		#remove leading comma
		samplelist=${samplelist#","}
		frag=${frag#","}
		
		#call the peak calling R script
		Rscript "$sPath/signalmaker.r" -chromoS="$chrSize" -chromo="$chr" -bednames=$samplelist -frag=$frag -bkgd=$bkgdfile -out="$wdir/signal.$ran/" -sig="$signal" -regions="$regs" -p="$cores" -chrcount="$counting" -t="$type"
		cat $wdir/signal.$ran/*.bedSignal | awk -F"\t" -v j=0 '$4 > j' > "$out/signal/$chr.bedGraph"
		rm $wdir/signal.$ran/*.bedSignal
		counting=$(($counting+1));
done
counting=1;
fi



if [ $type == "paired" ]; then
counting=1;
for f in $wdir/sizes.$ran/*; do #for each chromosome
	samplelist=""
	
	sizefile=$(basename $f)
	chr=$(echo $sizefile | awk -F"." '{print $2}' | awk -F"." '{print $1}');
	chrSize=$(cat $f | cut -f2);
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
	Rscript "$sPath/signalmaker.r" -chromoS="$chrSize" -chromo="$chr" -bednames=$samplelist -frag=$frag -bkgd=$bkgdfile -out="$wdir/signal.$ran/" -sig="$signal" -regions="$regs" -p="$cores" -chrcount="$counting" -t="$type"
	cat $wdir/signal.$ran/*.bedSignal | awk -F"\t" -v j=0 '$4 > j' > "$out/signal/$chr.bedGraph"
	rm $wdir/signal.$ran/*.bedSignal
	counting=$(($counting+1));
done
counting=1;
fi
#=======================> DONE!

rm -rf $wdir

printf "\n\n========================================\nWe're done...Congratulations!\n========================================\n\n"
