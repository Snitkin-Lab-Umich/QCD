#!/usr/bin/env bash
# A simple shell Script to summarize MLST reports
# The script takes following inputs:
#
#    -Path to ARIBA MLST directory
#
####################################################################
usage="$(basename "$0 output_directory") [-h] 
-- A simple shell Script to summarize MLST reports

where:
    -h  show this help text"

while getopts ':hs:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done

######################################################################

#printf "\nSearching MLST reports in directory: $1\n\n"

for i in `find $1 -name "mlst_report.tsv" -print`; do filename=`echo $i | sed 's/\/mlst_report.tsv//g'`; header=`head -n1 $i`; mlst=`head -n2 $i | tail -n1`; printf "\t$header\n$filename\t$mlst\n" >> MLST_temp.txt; done
cat MLST_temp.txt | sort | uniq
rm MLST_temp.txt
