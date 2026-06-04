#!/bin/bash - 
#===============================================================================
#
#          FILE: quaak.sh
# 
#         USAGE: bash ./quaak.sh  [OPTIONS] value
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Ying Zhou
#  ORGANIZATION: 
#       CREATED: 04/27/2026 15:32
#      REVISION: v0.01
#===============================================================================

set -o nounset                              # Treat unset variables as an error
set -o pipefail                             # Pipeline fails if any stage fails
shopt -s inherit_errexit 2>/dev/null || true   # bash 4.4+; propagate failures into $(...)

SCRIPTPATH=$(dirname $0)

#input files
KMER=unset
QUE=unset
REF=unset
OUT=output
CYTOBANDS=unset

#optional parameters
BLOCKCUT=500000,10000 # kmer distance cut, indel cut
TAUCUT=0.9 # select contig with confident strand direction
SVCUT=500000,100

#for devlopers
alltrue=true

run_path=false
run_block=false
run_sel_ctgs=false
run_annot_sv=false
run_brief_summary=false
run_plot_sv=true

if $alltrue; then
  run_path=true
  run_block=true
  run_sel_ctgs=true
  run_annot_sv=true
  run_brief_summary=true
  run_plot_sv=true
fi

usage()
{
  echo "
  Usage: bash ${SCRIPTPATH}/quaak.sh  [OPTIONS] value
      [ -k | --kmer <string> references kmer set (.fa, .fa.gz)               ]
      [ -r | --reference <string> reference asm (.fa, .fa.gz, path, path.gz) ]
      [ -q | --query <string> query asm (.fa, .fa.gz, path, path.gz)         ]
      [ -o | --out <string> outpref (optional, default ./output)            ]
      [--cytobands <string> bedfile (optional)                              ]
      [--blockcut <int,int> (optional, default 10000000,100000)
         first value is distance cut between kmers, second value is distance 
         difference of kmer pair between ref and que                        ]
      [--taucut <float>  (optional, default 0.9)
         abs(tau) cutoff to select contif with confident strand direction   ]
      [--svcut <int,int> (optional, default 100000,100)
         length and kmer count cutoff to report large SV and asm errors     ]
  Example :
     bash quaak.sh -r ref.path.gz -q query.path.gz -o test.out
     bash quaak.sh -k ref.kmer.fa.gz -r ref.fa.gz -q query.fa.gz -o test.out
  Note: '-k' is required if any of -r and -q use fa/fa.gz as input
  "
  exit 2
}


PARSED_ARGUMENTS=$(getopt -a -n quaak \
  -o k:r:q:o: \
  --long kmer:,reference:,query:,out:,blockcut:,taucut:,svcut:,cytobands:\
  -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
  usage
  exit 1
fi

#echo "PARSED_ARGUMENTS is $PARSED_ARGUMENTS"
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -k | --kmer)  KMER="$2"    ; shift 2 ;;
    -r | --reference)  REF="$2"    ; shift 2 ;;
    -q | --query)  QUE="$2"    ; shift 2 ;;
    -o | --out)  OUT="$2"    ; shift 2 ;;
    --blockcut)  BLOCKCUT="$2"    ; shift 2 ;;
    --taucut)  TAUCUT="$2"    ; shift 2 ;;
    --svcut)  SVCUT="$2"    ; shift 2 ;;
    --cytobands)  CYTOBANDS="$2"    ; shift 2 ;;
    --) shift; break ;;
    *) echo "Unexpected option: $1."
      usage
      ;;
  esac
done


out_dir=$(dirname -- "${OUT}")
[[ -d "${out_dir}" && -w "${out_dir}" ]] || {
  echo "Error: output directory '${out_dir}' missing or not writable"; 
  exit 1;
}


if [[ $QUE == unset || $REF == unset ]]; then
  echo "Error: both '-r/--reference, -q/--query' are required"
  usage
  exit 1
fi

run_que_path=false
if [[ -f ${QUE} && ( ${QUE} == *.fa || ${QUE} == *.fa.gz || ${QUE} == *.fasta || ${QUE} == *.fasta.gz) ]]; then
  if [ $KMER == unset ]; then
    echo "Error: '-k/--kmer' is required for generating path file of ${QUE}"
    exit 1
  fi
  run_que_path=true
fi

run_ref_path=false
if [[ -f ${REF} && ( ${REF} == *.fa || ${REF} == *.fa.gz || ${REF} == *.fasta || ${REF} == *.fasta.gz) ]]; then
  if [ $KMER == unset ]; then
    echo "Error: '-k/--kmer' is required for generating path file of ${REF}"
    exit 1
  fi
  run_ref_path=true
fi

if $alltrue; then

  start_second=`date +%s`
  start=`date +%D-%H:%M:%S`

  echo "########################################"
  echo "##Welcome###############################"
  echo "########################################"
  echo ""
  echo "Starting time: ${start}"
  echo "#####parameters:########################"
  if [[ ${run_ref_path} || ${run_que_path} ]]; then
    echo "KMER(-k/--kmer)                  : $KMER"
  fi
  echo "REFERENCE(-r/--reference)        : $REF"
  echo "QUERY(-q/--query)                : $QUE"
  echo "OUT(-o/--out)                    : $OUT"
  echo "CYTOBANDS(--cytobands)           : $CYTOBANDS"
  echo "BLOCKCUT(--blockcut)             : $BLOCKCUT"
  echo "TAUCUT(--taucut)                 : $TAUCUT"
  echo "SVCUT(--svcut)                   : $SVCUT"

  echo ""
fi
out_pref=${OUT}

kmap=${SCRIPTPATH}/kmer-C-ult/kmer-map

ref_path=${REF}
que_path=${QUE}

if ${run_path} ; then

  if ${run_ref_path} ; then
    echo "## Calculate reference path"
    ref_path=${out_pref}.ref.path.gz
    ${kmap} ${KMER} ${REF} | gzip -c > ${ref_path}
    [[ -s "${ref_path}" ]] || { echo "ERROR: empty ref path"; exit 1; }
  fi

  if ${run_que_path} ; then
    echo "## Calculate query path"
    que_path=${out_pref}.que.path.gz
    ${kmap} ${KMER} ${QUE} | gzip -c > ${que_path}
    [[ -s "${que_path}" ]] || { echo "ERROR: empty que path"; exit 1; }
  fi
fi

kb=${SCRIPTPATH}/kmer-C-ult/kmer-block

if ${run_block} ; then
  if [[ -f ${ref_path} && -f ${que_path} ]] ; then
    echo "## Build blocks"
    ${kb} --cut ${BLOCKCUT} ${ref_path} ${que_path} ${out_pref}
  else :
    echo 'Incomplete input to process blocking, check input of -r/--reference, -t/--query'
  fi
fi

if ${run_sel_ctgs} ; then
  echo "## Separate ctgs with confident strand direction"
  sc="python3 ${SCRIPTPATH}/kmer-py-ult/sel_ctg.py"
  sb="python3 ${SCRIPTPATH}/kmer-py-ult/sample_blocks.py"
  pp="python3 ${SCRIPTPATH}/block-panel-plot/panel-plot.py"
  ## input
  strand_file=${out_pref}.ctg-strand.tsv
  block_file=${out_pref}.block.tsv.gz
  ## output
  susp_ctgs=${out_pref}.suspect-ctgs.txt
  susp_blks=${out_pref}.suspect-block.tsv.gz
  susp_blks_pdf=${out_pref}.suspect-block.pdf
  conf_ctgs=${out_pref}.conf-ctgs.txt
  rescue_ctg=${out_pref}.rescue-ctgs.txt
  sel_ctgs=${out_pref}.sel-ctgs.txt

  cat ${strand_file} | ${sc} ${TAUCUT} ${conf_ctgs} ${susp_ctgs}
  # extract uncertain ctgs and use plot to rescue some
  if [ ! -f "$rescue_ctg" ]; then
    touch ${rescue_ctg}
  fi
  if [ -s ${susp_ctgs} ]; then
    zcat ${block_file} | ${sb} ${susp_ctgs} | gzip -c > ${susp_blks}
    ${pp} ${susp_blks} ${susp_blks_pdf} --highlight ${rescue_ctg} --compress 1
  fi
  cat ${conf_ctgs} ${rescue_ctg} > ${sel_ctgs}

fi

if ${run_annot_sv} ; then
  echo "## Annotate SV"
  sel_ctgs=${out_pref}.sel-ctgs.txt
  block_file=${out_pref}.block.tsv.gz
  #missing_kmer=${out_pref}.missing.tsv.gz
  ab="python3 ${SCRIPTPATH}/kmer-py-ult/annot_block.py"
  annot_file=${out_pref}.annot.tsv.gz
  ${ab} ${sel_ctgs} ${block_file} ${ref_path}\
    | gzip -c > ${annot_file}
fi


if ${run_brief_summary}; then
  echo "## Generate Summary"
  strand_file=${out_pref}.ctg-strand.tsv
  cnv_kmer=${out_pref}.cnv.tsv.gz
  kmer_count=${out_pref}.kmer-summary.tsv
  annotation_file=${out_pref}.annot.tsv.gz
  conf_ctgs=${out_pref}.conf-ctgs.txt
  susp_ctgs=${out_pref}.suspect-ctgs.txt
  resc_ctgs=${out_pref}.rescue-ctgs.txt
  # what to summary
  # 1) large missing segments
  # 2) large SV
  # 3) inter-chromosome insertion, large
  # 4) total kmer summary
  out_summary=${out_pref}.summary.txt

  br="python3 ${SCRIPTPATH}/kmer-py-ult/brief.py"

  n_conf_ctgs=$(wc -l < ${out_pref}.conf-ctgs.txt)
  n_susp_ctgs=$(wc -l < ${out_pref}.suspect-ctgs.txt)
  n_resuce_ctgs=$(wc -l < ${out_pref}.rescue-ctgs.txt)

  ${br} ${kmer_count} ${TAUCUT} ${conf_ctgs} ${susp_ctgs} ${resc_ctgs} ${SVCUT} ${strand_file} ${cnv_kmer} ${annotation_file} > ${out_summary}

fi


if ${run_plot_sv}; then
  block_file=${out_pref}.block.tsv.gz
  annotation_file=${out_pref}.annot.tsv.gz
  strand=${out_pref}.sel-ctgs.txt

  sls="python3 ${SCRIPTPATH}/kmer-py-ult/sel_large_sv.py"
  mpp="python3 ${SCRIPTPATH}/block-panel-plot/m-panel-plot.py"
  agp="python3 ${SCRIPTPATH}/block-panel-plot/anchor-gap-plot.py"
  esb="python3 ${SCRIPTPATH}/kmer-py-ult/ext_SV_block.py"

  overall_ref_pdf=${out_pref}.allref.pdf
  overall_ctg_pdf=${out_pref}.allctg.pdf
  sel_gid=${out_pref}.SV.gid.txt
  sel_sv_block=${out_pref}.tmp-block.tsv
  sel_sv_pdf=${out_pref}.sel-sv.pdf
  sel_ctgs=${out_pref}.sel-ctgs.txt

  ## 1, overall plot
  echo "## Panel plot for all blocks"
  if [[ ${CYTOBANDS} == unset ]]; then
    ${mpp} ${block_file} ${overall_ref_pdf} --ref-centric --one-page --strand ${sel_ctgs}
    ${mpp} ${block_file} ${overall_ctg_pdf} --ref-centric --one-plot --strand ${sel_ctgs}
  else
    ${mpp} ${block_file} ${overall_ref_pdf} --ref-centric --one-page --strand ${sel_ctgs} --cytobands ${CYTOBANDS}
    ${mpp} ${block_file} ${overall_ctg_pdf} --ref-centric --one-plot --strand ${sel_ctgs} --cytobands ${CYTOBANDS}
  fi

  ## 2, large SV plot
  echo "## Panel plot for large SV only"
  ${sls} ${SVCUT} ${annotation_file} | uniq \
    | while read -r line;do
  gid=$(echo ${line} | cut -f1 -d' ')
  ctg=$(echo ${line} | cut -f2 -d' ')
  zgrep -w "^${ctg}" ${annotation_file} | ${esb} ${gid} ${sel_ctgs} 
done | uniq > ${sel_sv_block}

  if [ -s "$sel_sv_block" ]; then
    ${agp} ${sel_sv_block} ${sel_sv_pdf} --compress 0.05
  fi

fi


if $alltrue ; then
  end_second=`date +%s`
end=`date +%D-%H:%M:%S`
runtime=$((end_second - start_second))
echo -e "\n\nEnding time: ${end}, Wallclock time :${runtime} seconds"
fi
echo "####done####"
