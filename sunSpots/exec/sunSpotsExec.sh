#!/usr/bin/env bash

# I got sunSpots from
# http://sidc.be/silso/datafiles

detrend=2
fluctuation=2 # hurst exponent <=> flutctuation=2
minGroup=4 # should be >=4
maxGroup=1500 # shouldn't be too high 
numberOfWindows=30

# Paths to find data and to save plots
sunSpotOriginalData="../data/SN_d_tot_V2.0.csv"
sunSpotEditData="../data/editedSN.dat"
sunSpotPlot="../img/dataPlot.tex"
sunProfilePlot="../img/profilePlot.tex"
sunDFAPlot="../img/DFAPlot.tex"


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

./../../src/hurst-exponent $numberOfPoints $numberOfColumns $detrend $fluctuation $minGroup $maxGroup $numberOfWindows $sunSpotEditData

QsunProfilePlot="'${sunProfilePlot}'"
QprofileData="'profile.dat'"
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
QDFAData="'DFA.dat'"

c0=$(awk "FNR == 1 {print $1}" params.dat)
c1=$(awk "FNR == 2 {print $1}" params.dat)
H=$(awk "FNR == 3 {print $1}" params.dat)

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

f(x) = $c0 * ( x ** $c1)

plot $QDFAData using 1:2 t 'DFA' pt 6, \
	f(x) t 'H = $H'

EOF

