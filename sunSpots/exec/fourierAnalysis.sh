#!/usr/bin/env bash

# I got sunSpots from
# http://sidc.be/silso/datafiles

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

cd $DIR

detrend=1 # hurst exponent <=> detrend=1
fluctuation=2 # hurst exponent <=> flutctuation=2
minGroup=4 # should be >=4
maxGroup=900 # shouldn't be too high 
nGroup=300
numberOfWindows=10

# Paths to find data and to save plots
sunSpotOriginalData="../data/raw/SN_m_tot_V2.0.csv"
imgPath="../img/fourier/"
sunSpotPlot="${imgPath}dataPlot.tex"
sunDFAPlot="${imgPath}DFAPlot.tex"
outputPath="../data/processed/fou1/"
profileName="${outputPath}profile.dat"
sunSpotEditData="${outputPath}editedSN.dat"

f1Path="../data/processed/fou1/"
f1ProfilePlot="${imgPath}f1profilePlot.tex"
f1Name="${f1Path}profile.dat"
f1DFAData="${f1Path}DFA.dat"
f1s0params="${f1Path}s0params.dat"
f1s1params="${f1Path}s1params.dat"
f1s2params="${f1Path}s2params.dat"
f1Detrend=50

f2Path="../data/processed/fou2/"
f2ProfilePlot="${imgPath}f2profilePlot.tex"
f2Name="${f2Path}profile.dat"
f2DFAData="${f2Path}DFA.dat"
f2s0params="${f2Path}s0params.dat"
f2s1params="${f2Path}s1params.dat"
f2Detrend=100

# Store in a new file: index (eq. to days pass since 01.01.1818), value, error
awk -F';' '$4 != -1 {print NR-1 " " $4 " "$5 " " $1}' $sunSpotOriginalData > $sunSpotEditData
numberOfColumns=4

# Choose some years to show as xticlabel 
numberOfPoints=$(wc -l < $sunSpotEditData)
xTicLabel=$(awk -v range=$(((numberOfPoints)/9)) 'BEGIN {print "("} !((NR-1)%range) {print "\""$4"\" " $1 ", "} END {print  ")"}' $sunSpotEditData | tr -d '\n' ) 

aus=$((numberOfPoints/(detrend+3)))
maxGroup=$(( maxGroup < aus ?  maxGroup : aus))

./../../exec/PR $numberOfPoints $numberOfColumns $sunSpotEditData $outputPath
# fourier remove trend
./../../exec/FO $numberOfPoints 2 $f1Detrend $profileName $f1Name
./../../exec/DFA $numberOfPoints 2 $detrend $fluctuation $minGroup $maxGroup $nGroup $numberOfWindows $f1Name $f1Path

./../../exec/FO $numberOfPoints 2 $f2Detrend $profileName $f2Name
./../../exec/DFA $numberOfPoints 2 $detrend $fluctuation $minGroup $maxGroup $nGroup $numberOfWindows $f2Name $f2Path

QsunProfilePlot="'${f1ProfilePlot}'"
QprofileData="'${f1Name}'"
gnuplot << EOF
set terminal cairolatex pdf transparent size 8cm, 4.5cm font ',8'
set output $QsunProfilePlot
load 'pal1.pal'

stats $QprofileData using 1 nooutput name 'X_'
stats $QprofileData using 2 nooutput name 'Y_'
set xtics $xTicLabel rotate by 45 center offset 0,-1

set format y ' '

set xrange [X_min:X_max]
set yrange [Y_min:Y_max*1.1]

plot $QprofileData using 1:2 notitle pt 0
EOF

QsunProfilePlot="'${f2ProfilePlot}'"
QprofileData="'${f2Name}'"
gnuplot << EOF
set terminal cairolatex pdf transparent size 8cm, 4.5cm font ',8'
set output $QsunProfilePlot
load 'pal1.pal'

stats $QprofileData using 1 nooutput name 'X_'
stats $QprofileData using 2 nooutput name 'Y_'
set xtics $xTicLabel rotate by 45 center offset 0,-1

set format y ' '

set xrange [X_min:X_max]
set yrange [Y_min:Y_max*1.1]

plot $QprofileData using 1:2 notitle pt 0
EOF

QsunDFAPlot="'${sunDFAPlot}'"
Q1DFAData="'$f1DFAData'"
Q2DFAData="'$f2DFAData'"

get_data_for_plot ()
{
    local min=$1
    local max=$2
    local params=$3
    local DFAData=$4
    ./../../exec/HE $nGroup $min $max $DFAData $params
    local c0=$(awk 'FNR == 1 {print $1}' $params)
    local c1=$(awk 'FNR == 2 {print $1}' $params)
    local chi2=$(awk 'FNR == 3 {print $1}' $params)
    local mmax=$(awk -F' ' -v j=$((min + 1)) 'FNR == j {print $1}' $DFAData)
    local mmin=$(awk -F' ' -v j=$((max + min)) 'FNR == j {print $1}' $DFAData)
    echo "$c0 $c1 $chi2 $mmax $mmin"
}

s10=($(get_data_for_plot 0 103 $f1s0params $f1DFAData))
s11=($(get_data_for_plot 103 $((218 - 103)) $f1s0params $f1DFAData))
s12=($(get_data_for_plot 218 $((nGroup - 218))  $f1s2params $f1DFAData))

s20=($(get_data_for_plot 0 150 $f2s0params $f2DFAData))
s21=($(get_data_for_plot 150 $((nGroup- 150)) $f2s1params $f2DFAData))


gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm 
set output $QsunDFAPlot
load 'pal1.pal'

set logscale xy
set key inside bottom right horizontal maxcols 2 
stats $Q1DFAData using 1 nooutput name 'X_'
stats $Q1DFAData using 2 nooutput name 'Y_'
set ylabel '$ F^{k}_$fluctuation (\tau) $'
set xlabel '$\tau$ [months]'
set format x '$ 10^{%L}$'
set format y '$ 10^{%L}$'

set xrange [X_min:X_max]
set yrange [Y_min:Y_max*1.1]

s10f(x) = (${s10[4]} <= x && x<= ${s10[3]}) ? ${s10[0]} * ( x ** ${s10[1]}) : 1/0
s11f(x) = (${s11[4]}  <= x && x <= ${s11[3]}) ? ${s11[0]} * ( x ** ${s11[1]}) : 1/0
s12f(x) = (${s12[4]}  <= x && x <= ${s12[3]}) ? ${s12[0]} * ( x ** ${s12[1]}) : 1/0

s20f(x) = (${s20[4]} <= x && x<= ${s20[3]}) ? ${s20[0]} * ( x ** ${s20[1]}) : 1/0
s21f(x) = (${s21[4]}  <= x && x <= ${s21[3]}) ? ${s21[0]} * ( x ** ${s21[1]}) : 1/0


set arrow from 4,graph(0,0) to 4,graph(1,1) nohead
set arrow from 17,graph(0,0) to 17,graph(1,1) nohead
set arrow from 132,graph(0,0) to 132,graph(1,1) nohead
set arrow from 450,graph(0,0) to 450,graph(1,1) nohead

plot $Q1DFAData using 1:2 t'F$f1Detrend-DFA1' pt 6 lc 2, \
 	s10f(x) t '$(printf "%.2f" ${s10[1]})' lw 4 lc 10, \
	s11f(x) t '$(printf "%.2f" ${s11[1]})' lw 4 lc 9, \
	s12f(x) t '$(printf "%.2f" ${s12[1]})' lw 4 lc 8, \
	$Q2DFAData using 1:2 t 'F$f2Detrend-DFA1' pt 6 lc 1, \
 	s20f(x) t '$(printf "%.2f" ${s20[1]})' lw 4 lc 7, \
	s21f(x) t '$(printf "%.2f" ${s21[1]})' lw 4 lc 6 ,\
	0 t " " lc rgb '#ffffffff'
EOF

