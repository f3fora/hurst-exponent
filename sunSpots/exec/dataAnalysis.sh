#!/usr/bin/env bash

# I got sunSpots from
# http://sidc.be/silso/datafiles

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $DIR

detrend=3
fluctuation=2 # hurst exponent <=> flutctuation=2
minGroup=4 # should be >=4
maxGroup=1000 # shouldn't be too high 
nGroup=300
numberOfWindows=10

# Paths to find data and to save plots
sunSpotOriginalData="../data/raw/SN_d_tot_V2.0.csv"
sunSpotPlot="../img/dataPlot.tex"
sunProfilePlot="../img/profilePlot.tex"
sunDFAPlot="../img/DFAPlot.tex"
outputPath="../data/processed/"
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

plot $QprofileData using 1:2 notitle pt 1
EOF

QsunDFAPlot="'${sunDFAPlot}'"
QDFAData="'$DFAData'"

min=0
max=$nGroup
./../../exec/HE $nGroup $min $max $DFAData $params
tc0=$(awk 'FNR == 1 {print $1}' $params)
tc1=$(awk 'FNR == 2 {print $1}' $params)
tH=$(awk 'FNR == 3 {print $1}' $params)
tMax=$(awk -F' ' -v j=$((min + 1)) 'FNR == j {print $1}' $DFAData)
tMin=$(awk -F' ' -v j=$((max + min)) 'FNR == j {print $1}' $DFAData)

min=0
max=80
./../../exec/HE $nGroup $min $max $DFAData $s0params
s0c0=$(awk 'FNR == 1 {print $1}' $s0params)
s0c1=$(awk 'FNR == 2 {print $1}' $s0params)
s0H=$(awk 'FNR == 3 {print $1}' $s0params)
s0Max=$(awk -F' ' -v j=$((min + 1)) 'FNR == j {print $1}' $DFAData)
s0Min=$(awk -F' ' -v j=$((max + min + 1)) 'FNR == j {print $1}' $DFAData)

min=80
max=70
./../../exec/HE $nGroup $min $max $DFAData $s1params
s1c0=$(awk 'FNR == 1 {print $1}' $s1params)
s1c1=$(awk 'FNR == 2 {print $1}' $s1params)
s1H=$(awk 'FNR == 3 {print $1}' $s1params)
s1Max=$(awk -F' ' -v j=$((min + 1)) 'FNR == j {print $1}' $DFAData)
s1Min=$(awk -F' ' -v j=$((max + min + 1)) 'FNR == j {print $1}' $DFAData)

min=150
max=$((nGroup - min))
./../../exec/HE $nGroup $min $max $DFAData $s2params
s2c0=$(awk 'FNR == 1 {print $1}' $s2params)
s2c1=$(awk 'FNR == 2 {print $1}' $s2params)
s2H=$(awk 'FNR == 3 {print $1}' $s2params)
s2Max=$(awk -F' ' -v j=$((min + 1)) 'FNR == j {print $1}' $DFAData)
s2Min=$(awk -F' ' -v j=$((max + min)) 'FNR == j {print $1}' $DFAData)

gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm
set output $QsunDFAPlot

set logscale xy
set key right bottom 
set title 'DFA$detrend analysis for q=$fluctuation'
set ylabel 'F(s)'
set xlabel 's'
stats $QDFAData using 1 nooutput name 'X_'
stats $QDFAData using 2 nooutput name 'Y_'
set format x '$ 10^{%L}$'
set format y '$ 10^{%L}$'

set xrange [X_min:X_max]

tf(x) = ($tMin <= x && x <= $tMax) ? $tc0 * ( x ** $tc1) : 1/0
s0f(x) = ($s0Min <= x && x<= $s0Max) ? $s0c0 * ( x ** $s0c1) : 1/0
s1f(x) = ($s1Min  <= x && x <= $s1Max) ? $s1c0 * ( x ** $s1c1) : 1/0
s2f(x) = ($s2Min <= x && x <= $s2Max) ? $s2c0 * ( x ** $s2c1) : 1/0

plot $QDFAData using 1:2 t 'DFA' pt 6, \
	tf(x) t 'H = $tH' , \
 	s0f(x) t 'H = $s0H', \
	s1f(x) t 'H = $s1H', \
	s2f(x) t 'H = $s2H'


EOF

