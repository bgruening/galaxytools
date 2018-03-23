#!/usr/bin/env bash
# RNAcode sometimes fails because of bugs. Since the manual suggests
# to call RNAcode on splitted alignments it is feasible to run 
# RNAcode separately on the parts. This is implemented here. Command 
# line parameters just passed on to RNAcode. 
#
# - the script ensures that the region ids are continuous (otherwise
# the results for each block would start with 0)
# - also eps file names are corrected accordingly
# if RNAcode fails for one part it just outputs the part (for bug reporting)
# and continues

# for splitting the alignment you can use breakMAF.pl from the RNAcode 
# github (it seems to be absent from the 0.3 release) and feed the output 
# with the desired RNAcode parameters into this shell script, e.g.: 
# 
# breakMAF.pl < chrM.maf | ./processMAF.sh --tabular --eps --eps-dir eps2/ 

# parse the command line options 
# - removes --outfile, --help, --version, and the input file from the arguments
# - outfile and infile are stored in variables 
declare -a args=()
while [[ $# -gt 0 ]]
do
key="$1"
	case $key in
		-d|--eps-dir)
		epsdir=$2
		args+=($1 $2)
		shift # past argument
		shift # past value
		;;
		-e|--eps)
		eps=1
		args+=($1)
		shift # past argument
		;;
		-t|--tabular)
		tabular=1
		args+=($1)
		shift # past argument
		;;
		-o|--outfile)
		outfile=$2
		shift # past argument
		shift # past value
		;;
		-g|--gtf|-b|--best-only|-r|--best-region|-s|--stop-early)
		args+=($1)
		shift # past argument
		;;
		-n|--num-samples|-p|--cutoff|-c|--pars|-i|--eps-cutoff)
		args+=($1 $2)
		shift # past argument
		;;
		-h|--help|-v|--version)
		shift # past argument
		;;
		*)    # unknown option
		file=$1 
		shift # past argument
		;;
	esac
done

# fix output (renumber blocks)
# and move eps files (if present) to tmpdir
function fix_output {
 	if [[ -z "$last" ]]; then
 		last=0
 	fi
	while read line
 	do
 		i=`echo "$line" | awk '{print $1}'`
 		j=`echo "$i+$last" | bc`
		echo $line | awk -v n=$j '{printf("%d\t", n); for(i=2; i<=NF; i++){printf("%s", $(i)); if(i==NF){printf("\n")}else{printf("\t")}}}'
		if [[ ! -z "$eps" && -f ${epsdir:-eps}/hss-$i.eps ]]; then
			mv ${epsdir:-eps}/hss-$i.eps $tmpd/hss-$j.eps
		fi
	done
 	if [[ ! -z "$j" ]]; then
 		last=`echo "$j+1" | bc` 
 		unset j
 	fi
}
 
# run RNAcode for $tempfile if >= 3 sequences
function run_rnacode {
        >&2 echo -e "processing" `cat ${tmpif} | grep $ref | cut -d" " -f1-6`
	nl=`cat ${tmpif} | grep "^s" | wc -l`
	if [[ "$nl" -ge "3" ]]; then
		streams=$(mktemp -p '.')
		RNAcode $@ &> $streams
		if [[ "$?" != "0" ]]; then
			ef=$(mktemp -u -p '.')
			cat ${tmpif} > ${ef}.maf
			>&2 echo "RNAcode failed for the alignmentblock \""`cat ${tmpif} | grep $ref | cut -d" " -f 1-6`"\" (${ef}.maf)"
		fi
		# - filter stdout for lines containing the ref and redirect everything to stderr
		#   https://github.com/wash/rnacode/issues/9
		# - we can not pipe stdout | ... | fix_output since then $last can not be used as global variable
		cat $streams | grep -v $ref | egrep -v "^ HSS #|^====|^$" 1>&2 
		sed -i -n "/$ref/p" $streams
		fix_output < $streams 
		rm $streams
	else
		>&2 echo "less than 3 sequences in the alignment block \""`cat ${tmpif} | grep $ref | cut -d" " -f 1-6`"\""
	fi
}

ref=""
last=0

if [[ ! -z "$tabular" ]]; then
	if [[ ! -z "$outfile"  ]]; then
		echo -e "HSS #\tFrame\tLength\tFrom\tTo\tName\tStart\tEnd\tScore\tP" > "$outfile"
	else
		echo -e "HSS #\tFrame\tLength\tFrom\tTo\tName\tStart\tEnd\tScore\tP"
	fi
fi

tmpif=$(mktemp -p '.')
tmpof=$(mktemp -p '.')
tmpd=$(mktemp -d -p '.')

# process lines of the alignment 
# - save lines to tmpif 
# - empty lines: process tmpif (ie. last alignment block) with RNAcode, clear tmpif
# - on the go the name of the reference species is determined from the 1st line 
#   of the alignment this is used then for filtering the RNAcode output
while read line
do
	if [[ "$line" =~ ^# ]]; then
		echo -n > ${tmpif}
	elif [[ "$line" =~ ^$ ]]; then
		run_rnacode ${args[@]} ${tmpif}
		echo -n > ${tmpif}
	else
                if [[ -z $ref && "$line" =~ ^s ]]; then
			ref=`echo $line | sed 's/\./ /g' | cut -d" " -f 2`
		fi
		echo $line >> ${tmpif}
	fi
done < ${file:-/dev/stdin}
# if there is something left -> process it
if [[ "`cat ${tmpif} | wc -l`" -gt "0" ]]; then
	run_rnacode ${args[@]} ${tmpif}
fi

if [[ ! -z "$outfile" ]]; then
	cat $tmpof > "$outfile"
fi

if [[ ! -z "$eps" ]]; then
	mv ${tmpd}/*eps ${epsdir:-eps}/
fi

rm ${tmpif}
rm ${tmpof}
rmdir ${tmpd}
