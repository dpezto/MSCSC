set ter cairolatex size 10cm, 7cm
set datafile separator ','

dplt = 'Plots/' # relative path to plots folder
ddat = 'Data/' # relative path to data folder
gnupal = '~/Documents/Gnuplot/Palettes/' # path to gnuplot palettes

# out files *.py
file(n) = word('ttilde.csv utilde.csv', n)

# Number of different N,M
stats ddat.file(1) nooutput
N = STATS_blocks - 1 # -1 *.py adds extra \n\n at the end

# ttilde, utilde
vs(n) = word('tt uu', n)

graph(n) = word('fc Dn2 Ek', n)
labelx(n) = word('t/U U/t', n)
labely(n) = word('$f_c$ $\Delta(\hat{n})^2$ $E_\text{kin}$', n)

load gnupal.'mathematica.pal'
do for [i=2:4]{ # columns fir fc, Δn2 y E_kin (2, 3 y 4)
  do for [j=1:2]{ # files t/U y U/t
    set k tm center horizontal samplen 2 width 2
    set dataf sep ','
    # With Gnuplot 5.5
    if (i == 2 && j == 1){set logscale x}
    else if (j == 2){set xran [0:20]}
    else {set xran[0:0.5]}
    # With Gnuplot 5.4
    # if (i == 2 && j == 1){set logscale x
    # }else{set xran[0:0.5]}
    # if (j == 2){set xran [0:20]}
    set xlab labelx(j); set ylab labely(i-1)
    set o dplt.'plot_'.vs(j).'_vs_'.graph(i-1).'.tex'
    p for[k=1:N] ddat.file(j) i (k-1) u 1:i w l lw 3 t 'M='.(k+2)
    reset
  }
}
load 'paired.pal'
do for [i=1:2]{
  set dataf sep ','
  set k tm center horizontal samplen 2 width 2
  # graph gap
  if (i == 1){set xran[0:0.5]}else{set xran[0:20]}
  set xlab labelx(i); set ylab '$\Delta\mu$'
  set o dplt.'plot_'.vs(i).'_vs_mgap.tex'
  p for[k=1:N] ddat.file(i) i (k-1) u 1:($5-$6) w l lt 2*k lw 3 t 'M='.(k+2)
  # graph chem pot
  set ylab '$\mu^\pm$'
  set o dplt.'plot_'.vs(i).'_vs_mpm.tex'
  p for[k=1:N] ddat.file(i) i (k-1) u 1:5 w l lt (2*k) lw 3 t 'M='.(k+2), \
    for[k=1:N] ddat.file(i) i (k-1) u 1:6 w l lt (2*k-1) lw 3 notit
  reset
}

# hamiltonian graph
M(n) = word('M6 M10', n)
label(n) = word('\num{1e2} \num{1e4}', n)
r = 5e-2
cir = 5e-3

do for [i=1:2]{
  set dataf sep ','
  unset k
  set style fill transparent solid 0.8 noborder
  set size ratio 1
  H = 'H_'.M(i).'.csv' # file

  stats ddat.H u 1:3 name M(i) nooutput
  eval 'min = '.M(i).'_min_x-'.M(i).'_max_x*r'
  eval 'max = '.M(i).'_max_x+'.M(i).'_max_x*r'
  eval 'cirmax = sqrt(abs('.M(i).'_max_y))'
  mul = cir/(cirmax/(max-min))

  set format xy '%.0t'
  set xran [min:max]; set yran [max:min]
  set mxtics 2; set mytics 2
  set xlab label(i); set ylab label(i)
  load gnupal.'spectral.pal'
  set palette negative maxcolors 200
  set o dplt.'plot_H_'.M(i).'.tex'
  p ddat.H u 1:2:(mul*sqrt(abs($3))):3 w circles lc palette

  reset
}
