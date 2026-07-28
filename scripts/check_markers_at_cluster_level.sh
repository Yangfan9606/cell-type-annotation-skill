#!/usr/bin/env bash
# check_markers_at_cluster_level.sh — one-shot cluster-level marker check.
#
# Runs sc_006.find_marker.DEG.Heter.R --do_find_marker, then aggregates the
# per-cluster top50 markers with a strict->sensitive fallback (see
# references/scripts_reference.md §4a for why the fallback is needed) plus a
# specificity flag: any positive marker whose pct.2 (background) exceeds
# HIGHBG is tagged "high_bg_check" for manual review instead of being
# silently trusted.
#
# Usage:
#   bash check_markers_at_cluster_level.sh <rds> [cluster_col] [thresh] [high_bg_pct2]
#   nohup bash check_markers_at_cluster_level.sh <rds> > check_markers.log 2>&1 &   # to background it
#
# Output: results/sc006_top50_per_cluster.csv — hand this (plus the compute
# log, if anything looks off) to the skill for a per-cluster judgment.

set -euo pipefail

RDS="${1:?usage: check_markers_at_cluster_level.sh <rds> [cluster_col] [thresh] [high_bg_pct2]}"
CLU="${2:-seurat_clusters}"
THRESH="${3:-10}"
HIGHBG="${4:-0.5}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[1/2] sc_006.find_marker.DEG.Heter.R --do_find_marker (rds=$RDS, ident=$CLU) ..."
Rscript "${SCRIPT_DIR}/sc_006.find_marker.DEG.Heter.R" -i "$RDS" \
  -c "$CLU" --do_find_marker > sc_006.find_marker.DEG.Heter.R.log 2>&1

cd results

echo "[2/2] aggregating per-cluster top50 (strict->sensitive fallback, high-bg flag, thresh=$THRESH, high_bg_pct2=$HIGHBG) ..."
awk -F',' -v th="$THRESH" 'NR>1 && $8!="NA"{cnt[$6]++} END{for (c in cnt) if (cnt[c]>=th) print c}' \
  sc006_find_marker_strict_final.csv > .strict_ok_clusters.tmp
{
  head -n1 sc006_find_marker_strict_final.csv | sed 's/$/,"tier","pct_diff","specificity_flag"/'
  awk -F',' -v hb="$HIGHBG" 'NR==FNR{ok[$1]=1; next} FNR==1{next} $8!="NA" && ($6 in ok){
    diff=$3-$4
    flag=($2+0>0 && $4+0>hb)?"\"high_bg_check\"":"\"ok\""
    print $0",\"strict\","diff","flag
  }' .strict_ok_clusters.tmp sc006_find_marker_strict_final.csv \
    | sort -t',' -k6,6 -k8,8gr | awk -F',' '{c=$6; if(c!=prev){n=0;prev=c}; n++; if(n<=50) print}'
  awk -F',' -v hb="$HIGHBG" 'NR==FNR{ok[$1]=1; next} FNR==1{next} $8!="NA" && !($6 in ok){
    diff=$3-$4
    flag=($2+0>0 && $4+0>hb)?"\"high_bg_check\"":"\"ok\""
    print $0",\"sensitive\","diff","flag
  }' .strict_ok_clusters.tmp sc006_find_marker_sensitive_final.csv \
    | sort -t',' -k6,6 -k8,8gr | awk -F',' '{c=$6; if(c!=prev){n=0;prev=c}; n++; if(n<=50) print}'
} > sc006_top50_per_cluster.csv
rm -f .strict_ok_clusters.tmp

echo "done -> results/sc006_top50_per_cluster.csv"
echo "hand this file (and sc_006.find_marker.DEG.Heter.R.log if anything looks off) to the skill for judgment."
