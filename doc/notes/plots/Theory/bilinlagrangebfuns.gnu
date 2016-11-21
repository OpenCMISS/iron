#set title "Bilinear Lagrange basis functions"
set nokey
set samples 21
set isosample 11
set xlabel "$\xi_1$"
set ylabel "$\xi_2$"
set label 1 "$\lbfn{1}{\xi_1,\xi_2}$" at 0.00, 0.00, 1.20 centre
set label 2 "$\lbfn{2}{\xi_1,\xi_2}$" at 1.00, 0.00, 1.20 centre
set label 3 "$\lbfn{3}{\xi_1,\xi_2}$" at 0.00, 1.00, 1.20 centre
set label 4 "$\lbfn{4}{\xi_1,\xi_2}$" at 1.00, 1.00, 1.20 centre
#set xtics  0.00,0.25,1
#set ytics -0.25,0.25,1
psi1(x,y)=(1.0-x)*(1.0-y)
psi2(x,y)=x*(1.0-y)
psi3(x,y)=(1.0-x)*y
psi4(x,y)=x*y
set xrange [0:1]
set yrange [0:1]
set zrange [0:1]
splot psi1(x,y),psi2(x,y),psi3(x,y),psi4(x,y)
 
