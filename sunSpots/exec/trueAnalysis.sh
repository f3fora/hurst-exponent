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
profileName="${k1Path}profile.dat"
k1sunSpotEditData="${k1Path}editedSN.dat"
k1DFAData="${k1Path}DFA.dat"
k1params="${k1Path}params.dat"
dk1params="${k1Path}dparams.dat"
k1detrend=1

k3Path="../data/processed/k3/"
k3sunSpotEditData="${k3Path}editedSN.dat"
k3DFAData="${k3Path}DFA.dat"
k3params="${k3Path}params.dat"
dk3params="${k3Path}dparams.dat"
k3detrend=2

k8Path="../data/processed/k8/"
k8sunSpotEditData="${k8Path}editedSN.dat"
k8DFAData="${k8Path}DFA.dat"
k8params="${k8Path}params.dat"
dk8params="${k8Path}pdarams.dat"
k8detrend=4

k20Path="../data/processed/k20/"
k20sunSpotEditData="${k20Path}editedSN.dat"
k20DFAData="${k20Path}DFA.dat"
k20params="${k20Path}params.dat"
dk20params="${k20Path}dparams.dat"
uk20params="${k20Path}uparams.dat"
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


./../../exec/PR $numberOfPoints $numberOfColumns $sunSpotEditData $k1Path
./../../exec/DFA $numberOfPoints 2 $k1detrend $fluctuation $minGroup $k1max $nGroup $numberOfWindows $profileName $k1Path

./../../exec/DFA $numberOfPoints 2 $k3detrend $fluctuation $minGroup $k3max $nGroup $numberOfWindows $profileName $k3Path

./../../exec/DFA $numberOfPoints 2 $k8detrend $fluctuation $minGroup $k8max $nGroup $numberOfWindows $profileName $k8Path

./../../exec/DFA $numberOfPoints 2 $k20detrend $fluctuation $minGroup $k20max $nGroup $numberOfWindows $profileName $k20Path

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
range0=0
range1=$((range0 + numberOfPoints/4))
range2=$((range1 + 576))
range3=$((range2 + 450))
range4=$((range3 + 400))
range5=$((range4 + 311))
range6=$((range5 + 280))
range7=$((range6 + 220))
range8=$((range7 + 132))
range9=$((range8 + 60))
range10=$((range9 + 17))

dk1=($(get_data_for_plot 218 $((nGroup - 218)) $k1DFAData $dk1params))
dk3=($(get_data_for_plot 185 $((nGroup - 185)) $k3DFAData $dk3params))
dk8=($(get_data_for_plot 165 $((nGroup - 165)) $k8DFAData $dk8params))
dk20=($(get_data_for_plot 125 $((nGroup - 125)) $k20DFAData $dk20params))
uk20=($(get_data_for_plot 75 50 $k20DFAData $uk20params))

k1=($(get_data_for_plot 0 $nGroup $k1DFAData $k1params))
k3=($(get_data_for_plot 0 $nGroup $k3DFAData $k3params))
k8=($(get_data_for_plot 0 $nGroup $k8DFAData $k8params))
k20=($(get_data_for_plot 0 $nGroup $k20DFAData $k20params))

gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm
set output $QsunDFAPlot

set logscale xy
set key inside bottom right horizontal maxcols 4 
set ylabel '$ F^{k}_$fluctuation (\tau) $'
set xlabel '$\tau$ [months]'
set format x '$ 10^{%L}$'
set format y '$ 10^{%L}$'

k1(x) = (${k1[4]} <= x && x <= ${k1[3]}) ? ${k1[0]} * ( x ** ${k1[1]}) : 1/0
k3(x) = (${k3[4]} <= x && x <= ${k3[3]}) ? ${k3[0]} * ( x ** ${k3[1]}) : 1/0
k8(x) = (${k8[4]} <= x && x <= ${k8[3]}) ? ${k8[0]} * ( x ** ${k8[1]}) : 1/0
k20(x) = (${k20[4]} <= x && x <= ${k20[3]}) ? ${k20[0]} * ( x ** ${k20[1]}) : 1/0

dk1(x) = (${dk1[4]} <= x && x <= ${dk1[3]}) ? ${dk1[0]} * ( x ** ${dk1[1]}) : 1/0
dk3(x) = (${dk3[4]} <= x && x <= ${dk3[3]}) ? ${dk3[0]} * ( x ** ${dk3[1]}) : 1/0
dk8(x) = (${dk8[4]} <= x && x <= ${dk8[3]}) ? ${dk8[0]} * ( x ** ${dk8[1]}) : 1/0
dk20(x) = (${dk20[4]} <= x && x <= ${dk20[3]}) ? ${dk20[0]} * ( x ** ${dk20[1]}) : 1/0
uk20(x) = (${uk20[4]} <= x && x <= ${uk20[3]}) ? ${uk20[0]} * ( x ** ${uk20[1]}) : 1/0

load 'pal1.pal'
set arrow from $((range1 - range0)), graph(0,0) to $((range1 - range0)), graph(1,1) nohead lc 1 dt 9 lw 4
set arrow from $((range2 - range1)), graph(0,0) to $((range2 - range1)), graph(1,1) nohead lc 2 dt 9 lw 4
set arrow from $((range3 - range2)), graph(0,0) to $((range3 - range2)), graph(1,1) nohead lc 3 dt 9 lw 4
set arrow from $((range4 - range3)), graph(0,0) to $((range4 - range3)), graph(1,1) nohead lc 4 dt 9 lw 4
set arrow from $((range5 - range4)), graph(0,0) to $((range5 - range4)), graph(1,1) nohead lc 5 dt 9 lw 4
set arrow from $((range6 - range5)), graph(0,0) to $((range6 - range5)), graph(1,1) nohead lc 6 dt 9 lw 4
set arrow from $((range7 - range6)), graph(0,0) to $((range7 - range6)), graph(1,1) nohead lc 7 dt 9 lw 4
set arrow from $((range8 - range7)), graph(0,0) to $((range8 - range7)), graph(1,1) nohead lc 8 dt 9 lw 4
set arrow from $((range9 - range8)), graph(0,0) to $((range9 - range8)), graph(1,1) nohead lc 9 dt 9 lw 4
set arrow from $((range10 - range9)), graph(0,0) to $((range10 - range9)), graph(1,1) nohead lc 10 dt 9 lw 4

load 'pal2.pal'

plot $Qk1 using 1:2 t 'DFA$k1detrend' pt 6 lc 1,\
	k1(x) t '$(printf "%.2f" ${k1[1]})' lw 4 lc 2,\
	dk1(x) t '$(printf "%.2f" ${dk1[1]})' lw 4 lc 2 dt 7,\
	0 t " " lc rgb '#ffffffff' ,\
	$Qk3 using 1:2 t 'DFA$k3detrend' pt 6 lc 3,\
	k3(x) t '$(printf "%.2f" ${k3[1]})' lw 4 lc 4,\
	dk3(x) t '$(printf "%.2f" ${dk3[1]})' lw 4 lc 4 dt 7,\
	0 t " " lc rgb '#ffffffff' ,\
	$Qk8 using 1:2 t 'DFA$k8detrend' pt 6 lc 7,\
	k8(x) t '$(printf "%.2f" ${k8[1]})' lw 4 lc 8,\
	dk8(x) t '$(printf "%.2f" ${dk8[1]})' lw 4 lc 8 dt 7,\
	0 t " " lc rgb '#ffffffff' ,\
	$Qk20 using 1:2 t 'DFA$k20detrend' pt 6 lc 9,\
	k20(x) t '$(printf "%.2f" ${k20[1]})' lw 4 lc 10,\
	dk20(x) t '$(printf "%.2f" ${dk20[1]})' lw 4 lc 10 dt 7,\
	uk20(x) t '$(printf "%.2f" ${uk20[1]})' lw 4 lc 10 dt 9

EOF

QsunProfilePlot="'${sunProfilePlot}'"
QprofileData="'${profileName}'"

xTicLabel=$(awk -v range=$(((numberOfPoints)/9)) 'BEGIN {print "("} !((NR-1)%range) {print "\""$4"\" " $1 ", "} END {print  ")"}' $sunSpotEditData | tr -d '\n' ) 

echo $range10 $numberOfPoints
gnuplot << EOF
set terminal cairolatex pdf transparent size 16cm, 9cm
set output $QsunProfilePlot
set fit quiet
set fit logfile '/dev/null'
load 'pal1.pal'

set key inside bottom center horizontal

set ylabel '$ Y_n$' offset -1,0
set xlabel 't [year]'
stats $QprofileData using 1 nooutput name 'X_'
stats $QprofileData using 2 nooutput name 'Y_'
set xtics $xTicLabel
set ytics rotate by 45 center offset -2,0

set format y '$ %2.0t \times10^{%L}$'

set xrange [X_min:X_max]
set yrange [1.1*Y_min:Y_max*1.1]
N=11
array aa[N]
array bb[N]
array cc[N]
array dd[N]
array ee[N]
array ff[N]
array gg[N]
array hh[N]
array r[N]
r[1] = $range0
r[2] = $range1
r[3] = $range2
r[4] = $range3
r[5] = $range4
r[6] = $range5
r[7] = $range6
r[8] = $range7
r[9] = $range8
r[10] = $range9
r[11] = $range10

z(x,a,b,c,d,e,f,g,h) = a + b*x + c*x**2 +  d*x**3 + e*x**4 + f*x**5 + g*x**6 + h*x**7
do for [i=1:N-1]{
	set arrow from r[i+1],graph(0,0) to r[i+1],graph(1,1) nohead
	fit z(x,tmpa,tmpb,tmpc,tmpd,tmpe,tmpf,tmpg,tmph) $QprofileData using 1:2 every ::r[i]::r[i+1] via tmpa,tmpb,tmpc,tmpd,tmpe,tmpf,tmpg,tmph
	aa[i]=tmpa
	bb[i]=tmpb
	cc[i]=tmpc
	dd[i]=tmpd
	ee[i]=tmpe
	ff[i]=tmpf
	gg[i]=tmpg
	hh[i]=tmph
}


plot $QprofileData using 1:2 notitle pt 6 lc 11,\
	for [i=1:N-1] [r[i]:r[i+1]] z(x, aa[i], bb[i], cc[i], dd[i], ee[i], ff[i], gg[i], hh[i]) title sprintf("%d", r[i+1]-r[i] ) lw 4 lc i 

EOF

