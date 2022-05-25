unset multiplot
reset
load 'moreland.gp'
eval setdef('tex','0')

output='sv-flow-log'
if (tex == 1) {
  set lmargin at screen 0
  set rmargin at screen 1
  set bmargin at screen 0
  set tmargin at screen 1
  set term lua tikz standalone createstyle size 2in,2in  #fontscale 0.6
  set output output.'.tex'
}

#fn='svflow_boson'
fn='svflow_nonint_half'
#fn='svflow_fermion'
dat=fn.'.dat'
par=fn.'.par'
load par

###
# do stuff
set xrange [0:2]
set xlabel '$\vartheta_x / \pi$'
set xtics 0, 1
set mxtics 2

set yrange [-40:0]
set ylabel '$\log(w)$'
set ytics -80,20
set mytics 4

set grid xtics ytics  mxtics mytics

p dat u ($2 * 2 / ntx):(log($3**2 / nty)) w linesp pt 1 ps 0.4  pointinterval 1 notit
#p dat u ($2 * 2 / ntx):(log($3**2 / nty)) w l notit
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
