#!/bin/bash

# source argparse to get argparse python functionality in bash
source /standard/dremel_lab/scripts/argparse.bash || exit 1

# set required conda and conda env
source /standard/dremel_lab/conda/miniconda/etc/profile.d/conda.sh || exit 1

# define inputs and outputs
DEBUG=""

ARGPARSE_DESCRIPTION="Ramdomly subset fastq file"
argparse "$@" <<EOF || exit 1
parser.add_argument('-i', '--infastq', required=True,
                    help='input fastq file')
parser.add_argument('-o', '--outfastq', required=False, 
                    help='output fastq file .. substring "ss" added to input filename if not provided')
parser.add_argument('-n', '--nreads', default=10000, type=int,
                    help='number of reads in output [default %(default)s]')
parser.add_argument('--verbose', action='store_true',
                    default=False, help='verbose mode [default %(default)s]')
parser.add_argument('--debug', action='store_true',
                    default=False, help='enter debug mode [default %(default)s]')
EOF

if [[ ! -z $VERBOSE ]];then
# set exo
set -exo pipefail
fi


if [[ -z $OUTFASTQ ]];then
OUTFASTQ=$(echo $INFASTQ | sed "s/.fastq/.ss.fastq/g")
fi

echo $INFASTQ
echo $OUTFASTQ

if [[ ! -z "$DEBUG" ]];then
exit 0
fi

# run seqtk
conda activate biocondatools || exit 1
seqtk sample $INFASTQ $NREADS | gzip -c - > ${OUTFASTQ}
conda deactivate || exit 1