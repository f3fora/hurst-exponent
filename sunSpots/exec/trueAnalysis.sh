#!/usr/bin/env bash

# I got sunSpots from
# http://sidc.be/silso/datafiles

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $DIR

detrend=1 # hurst exponent <=> detrend=1
fluctuation=2 # hurst exponent <=> flutctuation=2
minGroup=4 # should be >=4
nGroup=300
numberOfWindows=10

# Paths to find data and to save plots
sunSpotOriginalData="../data/raw/SN_m_tot_V2.0.csv"
sunSpotEditData="../data/processed/kdata.dat"
imgPath="../img/true/"
sunSpotPlot="${imgPath}dataPlot.tex"
sunProfilePlot="${imgPath}profilePlot.tex"
sunDFAPlot="${imgPath}DFAPlot.tex"

k1Path="../data/processed/k1/"
k1profileName="${k1Path}profile.dat"
k1sunSpotEditData="${k1Path}editedSN.dat"
k1DFAData="${k1Path}DFA.dat"
k1params="${k1Path}params.dat"
k1detrend=1

k3Path="../data/processed/k3/"
k3profileName="${k3Path}profile.dat"
k3sunSpotEditData="${k3Path}editedSN.dat"
k3DFAData="${k3Path}DFA.dat"
k3params="${k3Path}params.dat"
k3detrend=2
k8Path="../data/processed/k8/"
k8profileName="${k8Path}profile.dat"
k8sunSpotEditData="${k8Path}editedSN.dat"
k8DFAData="${k8Path}DFA.dat"
k8params="${k8Path}params.dat"
k8detrend=4

k20Path="../data/processed/k20/"
k20profileName="${k20Path}profile.dat"
k20sunSpotEditData="${k20Path}editedSN.dat"
k20DFAData="${k20Path}DFA.dat"
k20params="${k20Path}params.dat"
k20detrend=7

# Store in a new file: index (eq. to days pass since 01.01.1818), value, error
awk -F';' '$4 != -1 {print NR-1 " " $4 " "$5 " " $1}' $sunSpotOriginalData > $sunSpotEditData
numberOfColumns=4

# Choose some years to show as xticlabel 
numberOfPoints=$(wc -l < $sunSpotEditData)

k1max=$((numberOfPoints/(k1detrend+3)))
k3max=$((numberOfPoints/(k3detrend+3)))
k8max=$((numberOfPoints/(k8detrend+3)))
k20max=$((numberOfPoints/(k20detrend+3)))

./../../exec/DFA $numberOfPoints $numberOfColumns $k1detrend $fluctuation $minGroup $k1max $nGroup $numberOfWindows $sunSpotEditData $k1Path
./../../exec/DFA $numberOfPoints $numberOfColumns $k3detrend $fluctuation $minGroup $k3max $nGroup $numberOfWindows $sunSpotEditData $k3Path
./../../exec/DFA $numberOfPoints $numberOfColumns $k8detrend $fluctuation $minGroup $k8max $nGroup $numberOfWindows $sunSpotEditData $k8Path
./../../exec/DFA $numberOfPoints $numberOfColumns $k20detrend $fluctuation $minGroup $k20max $nGroup $numberOfWindows $sunSpotEditData $k20Path

QsunDFAPlot="'${sunDFAPlot}'"
Qk1="'$k1DFAData'"
Qk3="'$k3DFAData'"
Qk8="'$k8DFAData'"
Qk20="'$k20DFAData'"

get_data_for_plot ()
{
    local min=$1
    local max=$2
    local params=$4
    local DFAData=$3
    ./../../exec/HE $nGroup $min $max $DFAData $params
    local c0=$(awk 'FNR == 1 {print $1}' $params)
    local c1=$(awk 'FNR == 2 {print $1}' $params)
    local chi2=$(awk 'FNR == 3 {print $1}' $params)
    local mmax=$(awk -F' ' -v j=$((min + 1)) 'FNR == j {print $1}' $DFAData)
    local mmin=$(awk -F' ' -v j=$((max + min)) 'FNR == j {print $1}' $DFAData)
    echo "$c0 $c1 $chi2 $mmax $mmin"
}

#k1=($(get_data_for_plot 218 $((nGroup - 218)) $k1DFAData $k1params))
#k3=($(get_data_for_plot 185 $((nGroup - 185)) $k3DFAData $k3params))
#k8=($(get_data_for_plot 165 $((nGroup - 165)) $k8DFAData $k8params))
#k20=($(get_data_for_plot 125 $((nGroup - 125)) $k20DFAData $k20params))

k1=($(get_data_for_plot 0 $nGroup $k1DFAData $k1params))
k3=($(get_data_for_plot 0 $nGroup $k3DFAData $k3params))
k8=($(get_data_for_plot 0 $nGroup $k8DFAData $k8params))
k20=($(get_data_for_plot 0 $nGroup $k20DFAData $k20params))

gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm
set output $QsunDFAPlot

set logscale xy
set key right bottom 
set title 'DFAk for k=$k1detrend, $k3detrend, $k8detrend, $k20detrend analysis for q=$fluctuation'
set ylabel 'F(s)'
set xlabel 'time range'
set format x '$ 10^{%L}$'
set format y '$ 10^{%L}$'

set xtics ("4 months" 4, "1 year" 12, "4 years" 48, "16 years" 192, "48 years" 576)

k1(x) = (${k1[4]} <= x && x <= ${k1[3]}) ? ${k1[0]} * ( x ** ${k1[1]}) : 1/0
k3(x) = (${k3[4]} <= x && x <= ${k3[3]}) ? ${k3[0]} * ( x ** ${k3[1]}) : 1/0
k8(x) = (${k8[4]} <= x && x <= ${k8[3]}) ? ${k8[0]} * ( x ** ${k8[1]}) : 1/0
k20(x) = (${k20[4]} <= x && x <= ${k20[3]}) ? ${k20[0]} * ( x ** ${k20[1]}) : 1/0

set arrow from 4,graph(0,0) to 4,graph(1,1) nohead
set arrow from 17,graph(0,0) to 17,graph(1,1) nohead
set arrow from 132,graph(0,0) to 132,graph(1,1) nohead
set arrow from 450,graph(0,0) to 450,graph(1,1) nohead

plot $Qk1 using 1:2 t 'DFA$k1detrend' pt 1,\
    $Qk3 using 1:2 t 'DFA$k3detrend' pt 1,\
    $Qk8 using 1:2 t 'DFA$k8detrend' pt 1,\
    $Qk20 using 1:2 t 'DFA$k20detrend' pt 1, \
	k1(x) t '$\alpha$ = $(printf "%.2f" ${k1[1]})' lw 4, \
	k3(x) t '$\alpha$ = $(printf "%.2f" ${k3[1]})' lw 4, \
	k8(x) t '$\alpha$ = $(printf "%.2f" ${k8[1]})' lw 4, \
	k20(x) t '$\alpha$ = $(printf "%.2f" ${k20[1]})' lw 4

EOF

gnuplot <<EOF

EOF

