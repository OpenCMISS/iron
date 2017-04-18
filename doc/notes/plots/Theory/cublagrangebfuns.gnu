#set title "Cubic Lagrange basis functions"
set nokey
set xlabel "$\xi$"
set label 1 "$\lbfn{1}{\xi}$" at 0.19, 0.480,0 centre
set label 2 "$\lbfn{2}{\xi}$" at 0.26, 0.990,0 centre
set label 3 "$\lbfn{3}{\xi}$" at 0.74, 0.990,0 centre
set label 4 "$\lbfn{4}{\xi}$" at 0.83, 0.480,0 centre
#set xtics  0.00,0.25,1
#set ytics -0.25,0.25,1
psi1(x)=0.5*(3.0*x-1.0)*(3.0*x-2.0)*(1.0-x)
psi2(x)=4.5*x*(3.0*x-2.0)*(x-1.0)
psi3(x)=4.5*x*(3.0*x-1.0)*(1.0-x)
psi4(x)=0.5*x*(3.0*x-1.0)*(3.0*x-2.0)
plot[0:1] psi1(x),psi2(x),psi3(x),psi4(x)
 
