#!/usr/bin/env bash

# I got sunSpots from
# http://sidc.be/silso/datafiles

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $DIR

detrend=1 # hurst exponent <=> detrend=1
fluctuation=2 # hurst exponent <=> flutctuation=2
minGroup=4 # should be >=4
maxGroup=1000 # shouldn't be too high 
nGroup=300
numberOfWindows=10

# Paths to find data and to save plots
sunSpotOriginalData="../data/raw/SN_d_tot_V2.0.csv"
imgPath="../img/daily/"
sunSpotPlot="${imgPath}dataPlot.tex"
sunProfilePlot="${imgPath}profilePlot.tex"
sunDFAPlot="${imgPath}DFAPlot.tex"
outputPath="../data/processed/daily/"
profileName="${outputPath}profile.dat"
sunSpotEditData="${outputPath}editedSN.dat"
DFAData="${outputPath}DFA.dat"
params="${outputPath}params.dat"
s0params="${outputPath}s0params.dat"
s1params="${outputPath}s1params.dat"
s2params="${outputPath}s2params.dat"


# Store in a new file: index (eq. to days pass since 01.01.1818), value, error
awk -F';' '$5 != -1 {print NR-1 " " $5 " "$6 " " $1}' $sunSpotOriginalData > $sunSpotEditData
numberOfColumns=4

# Choose some years to show as xticlabel 
numberOfPoints=$(wc -l < $sunSpotEditData)
xTicLabel=$(awk -v range=$(((numberOfPoints)/9)) 'BEGIN {print "("} !((NR-1)%range) {print "\""$4"\" " $1 ", "} END {print  ")"}' $sunSpotEditData | tr -d '\n' ) 

# Plot the data
QsunSpotPlot="'${sunSpotPlot}'"
QsunSpotEditData="'${sunSpotEditData}'"
gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm
set output $QsunSpotPlot

set title 'Daily observed sun spots number'
set xlabel 'Time [year]'
set ylabel 'Count' 
set bars small
stats $QsunSpotEditData using 1 nooutput name 'X_'
stats $QsunSpotEditData using 2 nooutput name 'Y_'

set xrange [X_min:X_max]
set yrange [Y_min:Y_max*1.1]
set xtics $xTicLabel 

plot $QsunSpotEditData using 1:2:3 notitle with yerrorbar pt 0 
EOF

aus=$((numberOfPoints/(detrend+3)))
maxGroup=$(( maxGroup < aus ?  maxGroup : aus))

./../../exec/DFA $numberOfPoints $numberOfColumns $detrend $fluctuation $minGroup $maxGroup $nGroup $numberOfWindows $sunSpotEditData $outputPath

QsunProfilePlot="'${sunProfilePlot}'"
QprofileData="'${profileName}'"
gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm
set output $QsunProfilePlot

set title 'profile function for sun spots number'
set xlabel 'Time [year]'
stats $QprofileData using 1 nooutput name 'X_'
stats $QprofileData using 2 nooutput name 'Y_'
set xtics $xTicLabel

set format y '$ %2.0t \times10^{%L}$'

set xrange [X_min:X_max]
set yrange [Y_min:Y_max*1.1]

plot $QprofileData using 1:2 notitle pt 0
EOF

QsunDFAPlot="'${sunDFAPlot}'"
QDFAData="'$DFAData'"

get_data_for_plot ()
{
    local min=$1
    local max=$2
    local params=$3
    ./../../exec/HE $nGroup $min $max $DFAData $params
    local c0=$(awk 'FNR == 1 {print $1}' $params)
    local c1=$(awk 'FNR == 2 {print $1}' $params)
    local chi2=$(awk 'FNR == 3 {print $1}' $params)
    local mmax=$(awk -F' ' -v j=$((min + 1)) 'FNR == j {print $1}' $DFAData)
    local mmin=$(awk -F' ' -v j=$((max + min)) 'FNR == j {print $1}' $DFAData)
    echo "$c0 $c1 $chi2 $mmax $mmin"
}

#total=($(get_data_for_plot 0 $nGroup $params))
total=($(get_data_for_plot 0 17 $params))
s0=($(get_data_for_plot 18 $((83 - 18)) $s0params))
s1=($(get_data_for_plot 84 $((194  - 84)) $s1params))
s2=($(get_data_for_plot 195 $((nGroup - 195)) $s2params))

gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm
set output $QsunDFAPlot

set logscale xy
set key right bottom 
set title 'DFA$detrend analysis for q=$fluctuation'
set ylabel 'F(s)'
set xlabel 'time range'
stats $QDFAData using 1 nooutput name 'X_'
stats $QDFAData using 2 nooutput name 'Y_'
set format x '$ 10^{%L}$'
set format y '$ 10^{%L}$'

set xrange [X_min:X_max]
set yrange [Y_min:Y_max*1.1]
set xtics ("4 months" 120, "1 years" 365, "4 years" 1460, "16 years" 5840, "48 years" 17520)

tf(x) = (${total[4]} <= x && x <= ${total[3]}) ? ${total[0]} * ( x ** ${total[1]}) : 1/0
s0f(x) = (${s0[4]} <= x && x<= ${s0[3]}) ? ${s0[0]} * ( x ** ${s0[1]}) : 1/0
s1f(x) = (${s1[4]}  <= x && x <= ${s1[3]}) ? ${s1[0]} * ( x ** ${s1[1]}) : 1/0
s2f(x) = (${s2[4]} <= x && x <= ${s2[3]}) ? ${s2[0]} * ( x ** ${s2[1]}) : 1/0

set arrow from 120,graph(0,0) to 120,graph(1,1) nohead
set arrow from 510,graph(0,0) to 510,graph(1,1) nohead
set arrow from 3960,graph(0,0) to 3960,graph(1,1) nohead
set arrow from 13500,graph(0,0) to 13500,graph(1,1) nohead

plot $QDFAData using 1:2 t 'DFA' pt 0, \
	tf(x) t '$\alpha$ = $(printf "%.2f" ${total[1]})' lw 4, \
 	s0f(x) t '$\alpha$ = $(printf "%.2f" ${s0[1]})'lw 4, \
	s1f(x) t '$\alpha$ = $(printf "%.2f" ${s1[1]})' lw 4, \
	s2f(x) t '$\alpha$ = $(printf "%.2f" ${s2[1]})' lw 4

EOF

