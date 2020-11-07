# An investigation of HURST EXPONENT

Clone this repository

```
git clone https://github.com/f3fora/hurst-exponent.git
```

## Requirements

- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [Gnuplot](http://www.gnuplot.info/)

## Instructions

Compile the C code
```
make
```
Eventually, remove the executables
```
make clean
```

Process data and generate plots
```
make monthlySunSpots
make dailySunSpots
make trueSunSpots
make fourierSunSpots
```

Get the essay as PDF, using latex with `sunSpots/tex/sunSpots.tex`
