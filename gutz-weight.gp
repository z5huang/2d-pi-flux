unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='gutz-weight'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 2in,2in  #fontscale 0.6
  set output output.'.tex'
}

fn='gutz_weight'
dat=fn.'.dat'
par=fn.'.par'
load par

###
# do stuff

set xlabel '$\vartheta_x/\pi$'
set ylabel '$\vartheta_y/\pi$'

dx = 2.0/ntx
dy = 2.0/nty
set xrange [0-dx:2]
set yrange [0-dy:2]
set xtics 0,1
set mxtics 2
set ytics 0,1
set mytics 2

set cbrange [0.00015:0.00019]
set cbtics 0.00001

set grid front

p dat matrix u (($1)*2/ntx):(($2)*2/ntx):3 w image pixels notit

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
