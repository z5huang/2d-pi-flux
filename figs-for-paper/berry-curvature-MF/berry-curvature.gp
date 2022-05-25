unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='berry-curvature'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 1in,1in  #fontscale 0.6
  set output output.'.tex'
}

fn='berry_curvature'
dat=fn.'.dat'
par=fn.'.par'
load par

###
# do stuff

set xlabel '$\kappa_x/\pi$'
set ylabel '$\kappa_y/\pi$'

set xrange [0:2]
set yrange [0:2]
set xtics 0,1
set mxtics 2
set ytics 0,1
set mytics 2

set grid front

set cbrange [1:3]
set cbtics 1,1
set mcbtics 2
unset colorbox
p dat matrix u (($1+0.5)*2/ntx):(($2+0.5)*2/ntx):(2*$3) w image pixels notit

#
###

if (tex == 1){
  unset output
  set term wxt
  build = buildtex(output)
  print '"eval build" to build and preview'
} else {
  #print "press enter to continue"
  #pause -1
}
