#set title "Quadratic Lagrange basis functions"
set nokey
set xlabel "$\xi$"
set label 1 "$\lbfn{1}{\xi}$" at 0.29, 0.45, 0 centre
set label 2 "$\lbfn{2}{\xi}$" at 0.50, 0.96, 0 centre
set label 3 "$\lbfn{3}{\xi}$" at 0.73, 0.45, 0 centre
#set xtics  0.00,0.25,1
#set ytics -0.25,0.25,1
psi1(x)=2.0*(x-0.5)*(x-1.0)
psi2(x)=4.0*x*(1.0-x)
psi3(x)=2.0*x*(x-0.5)
plot[0:1] psi1(x),psi2(x),psi3(x)
 
