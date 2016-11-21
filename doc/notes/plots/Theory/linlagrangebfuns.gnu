#set title "Linear Lagrange basis functions"
set nokey
set xlabel "$\xi$"
set label 1 "$\lbfn{1}{\xi}$" at 0.30, 0.800,0 centre
set label 2 "$\lbfn{2}{\xi}$" at 0.30, 0.200,0 centre
#set xtics  0.00,0.25,1
#set ytics -0.25,0.25,1
psi1(x)=1-x
psi2(x)=x
plot[0:1] psi1(x),psi2(x)
 
