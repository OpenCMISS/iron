A=2.0
L=2.0*pi

a=1.5
b=3.5
#a=L/6.0
#b=5.0*L/6.0
#a=L/3.0
#b=2.0*L/3.0

S1=a
S2=b-a
S3=L-b

U1=A*cos(0)
U2=A*cos(a)
U3=A*cos(b)
U4=A*cos(L)
dU1ds=-A*sin(0)
dU2ds=-A*sin(a)
dU3ds=-A*sin(b)
dU4ds=-A*sin(L)

max(x,y) = x>=y ? x : y
min(x,y) = x<=y ? x : y

psi10(x)=1.0-3.0*x*x+2.0*x*x*x
psi11(x)=x*x*x-2.0*x*x+x
psi20(x)=3.0*x*x-2.0*x*x*x
psi21(x)=x*x*x-x*x
dpsi10(x)=-6.0*x+6.0*x*x
dpsi11(x)=3.0*x*x-4.0*x+1.0
dpsi20(x)=6.0*x-6.0*x*x
dpsi21(x)=3.0*x*x-2.0*x
d2psi10(x)=-6.0+12.0*x
d2psi11(x)=6.0*x-4.0
d2psi20(x)=6.0-12.0*x
d2psi21(x)=6.0*x-2.0

uanal(x) = A*cos(x)
duanal(x) = -A*sin(x)
d2uanal(x) = -A*cos(x)

u1interp(x,s1,s2)=psi10(x)*U1+psi20(x)*U2+psi11(x)*dU1ds*s1+psi21(x)*dU2ds*s2
u2interp(x,s1,s2)=psi10(x)*U2+psi20(x)*U3+psi11(x)*dU2ds*s1+psi21(x)*dU3ds*s2
u3interp(x,s1,s2)=psi10(x)*U3+psi20(x)*U4+psi11(x)*dU3ds*s1+psi21(x)*dU4ds*s2
du1interp(x,s1,s2)=dpsi10(x)*U1+dpsi20(x)*U2+dpsi11(x)*dU1ds*s1+dpsi21(x)*dU2ds*s2
du2interp(x,s1,s2)=dpsi10(x)*U2+dpsi20(x)*U3+dpsi11(x)*dU2ds*s1+dpsi21(x)*dU3ds*s2
du3interp(x,s1,s2)=dpsi10(x)*U3+dpsi20(x)*U4+dpsi11(x)*dU3ds*s1+dpsi21(x)*dU4ds*s2
d2u1interp(x,s1,s2)=d2psi10(x)*U1+d2psi20(x)*U2+d2psi11(x)*dU1ds*s1+d2psi21(x)*dU2ds*s2
d2u2interp(x,s1,s2)=d2psi10(x)*U2+d2psi20(x)*U3+d2psi11(x)*dU2ds*s1+d2psi21(x)*dU3ds*s2
d2u3interp(x,s1,s2)=d2psi10(x)*U3+d2psi20(x)*U4+d2psi11(x)*dU3ds*s1+d2psi21(x)*dU4ds*s2

ual(x)= x<=a ? u1interp(x/S1,S1,S1) : x<=b ? u2interp((x-a)/S2,S2,S2) : u3interp((x-b)/S3,S3,S3)
uam(x)= x<=a ? u1interp(x/S1,S1,(S1+S2)/2.0) : x<=b ? u2interp((x-a)/S2,(S1+S2)/2.0,(S2+S3)/2.0) : u3interp((x-b)/S3,(S2+S3)/2.0,S3)
ugm(x)= x<=a ? u1interp(x/S1,S1,sqrt(S1*S2)) : x<=b ? u2interp((x-a)/S2,sqrt(S1*S2),sqrt(S2*S3)) : u3interp((x-b)/S3,sqrt(S2*S3),S3)
uhm(x)= x<=a ? u1interp(x/S1,S1,2.0*S1*S2/(S1+S2)) : x<=b ? u2interp((x-a)/S2,2*S1*S2/(S1+S2),2.0*S2*S3/(S2+S3)) : u3interp((x-b)/S3,2.0*S2*S3/(S2+S3),S3)
umax(x)= x<=a ? u1interp(x/S1,S1,max(S1,S2)) : x<=b ? u2interp((x-a)/S2,max(S1,S2),max(S2,S3)) : u3interp((x-b)/S3,max(S2,S3),S3)
umin(x)= x<=a ? u1interp(x/S1,S1,min(S1,S2)) : x<=b ? u2interp((x-a)/S2,min(S1,S2),min(S2,S3)) : u3interp((x-b)/S3,min(S2,S3),S3)
urms(x)= x<=a ? u1interp(x/S1,S1,sqrt((S1*S1+S2*S2)/2.0)) : x<=b ? u2interp((x-a)/S2,sqrt((S1*S1+S2*S2)/2.0),sqrt((S2*S2+S3*S3)/2.0)) : u3interp((x-b)/S3,sqrt((S2*S2+S3*S3)/2.0),S3)
dual(x)= x<=a ? du1interp(x/S1,S1,S1) : x<=b ? du2interp((x-a)/S2,S2,S2) : du3interp((x-b)/S3,S3,S3)
duam(x)= x<=a ? du1interp(x/S1,S1,(S1+S2)/2.0) : x<=b ? du2interp((x-a)/S2,(S1+S2)/2.0,(S2+S3)/2.0) : du3interp((x-b)/S3,(S2+S3)/2.0,S3)
dugm(x)= x<=a ? du1interp(x/S1,S1,sqrt(S1*S2)) : x<=b ? du2interp((x-a)/S2,sqrt(S1*S2),sqrt(S2*S3)) : du3interp((x-b)/S3,sqrt(S2*S3),S3)
duhm(x)= x<=a ? du1interp(x/S1,S1,2.0*S1*S2/(S1+S2)) : x<=b ? du2interp((x-a)/S2,2*S1*S2/(S1+S2),2.0*S2*S3/(S2+S3)) : du3interp((x-b)/S3,2.0*S2*S3/(S2+S3),S3)
dumax(x)= x<=a ? du1interp(x/S1,S1,max(S1,S2)) : x<=b ? du2interp((x-a)/S2,max(S1,S2),max(S2,S3)) : du3interp((x-b)/S3,max(S2,S3),S3)
dumin(x)= x<=a ? du1interp(x/S1,S1,min(S1,S2)) : x<=b ? du2interp((x-a)/S2,min(S1,S2),min(S2,S3)) : du3interp((x-b)/S3,min(S2,S3),S3)
durms(x)= x<=a ? du1interp(x/S1,S1,sqrt((S1*S1+S2*S2)/2.0)) : x<=b ? du2interp((x-a)/S2,sqrt((S1*S1+S2*S2)/2.0),sqrt((S2*S2+S3*S3)/2.0)) : du3interp((x-b)/S3,sqrt((S2*S2+S3*S3)/2.0),S3)
d2ual(x)= x<=a ? d2u1interp(x/S1,S1,S1) : x<=b ? d2u2interp((x-a)/S2,S2,S2) : d2u3interp((x-b)/S3,S3,S3)
d2uam(x)= x<=a ? d2u1interp(x/S1,S1,(S1+S2)/2.0) : x<=b ? d2u2interp((x-a)/S2,(S1+S2)/2.0,(S2+S3)/2.0) : d2u3interp((x-b)/S3,(S2+S3)/2.0,S3)
d2ugm(x)= x<=a ? d2u1interp(x/S1,S1,sqrt(S1*S2)) : x<=b ? d2u2interp((x-a)/S2,sqrt(S1*S2),sqrt(S2*S3)) : d2u3interp((x-b)/S3,sqrt(S2*S3),S3)
d2uhm(x)= x<=a ? d2u1interp(x/S1,S1,2.0*S1*S2/(S1+S2)) : x<=b ? d2u2interp((x-a)/S2,2*S1*S2/(S1+S2),2.0*S2*S3/(S2+S3)) : d2u3interp((x-b)/S3,2.0*S2*S3/(S2+S3),S3)
d2umax(x)= x<=a ? d2u1interp(x/S1,S1,max(S1,S2)) : x<=b ? d2u2interp((x-a)/S2,max(S1,S2),max(S2,S3)) : d2u3interp((x-b)/S3,max(S2,S3),S3)
d2umin(x)= x<=a ? d2u1interp(x/S1,S1,min(S1,S2)) : x<=b ? d2u2interp((x-a)/S2,min(S1,S2),min(S2,S3)) : d2u3interp((x-b)/S3,min(S2,S3),S3)
d2urms(x)= x<=a ? d2u1interp(x/S1,S1,sqrt((S1*S1+S2*S2)/2.0)) : x<=b ? d2u2interp((x-a)/S2,sqrt((S1*S1+S2*S2)/2.0),sqrt((S2*S2+S3*S3)/2.0)) : d2u3interp((x-b)/S3,sqrt((S2*S2+S3*S3)/2.0),S3)

set samples 1000
set linetype 10 dt 2
set xlabel "$\xi$"
set label 1 "$a$" at a,0.1,0 centre
set label 2 "$b$" at b,0.1,0 centre
set label 3 "$2\pi$" at L,0.1,0 centre

set title "Second derivative interpolation"
set yrange [-25:20]
set key at 6.0,17.5,0
set arrow 1 from a,-25 to a,20 nohead linetype 10
set arrow 2 from b,-25 to b,20 nohead linetype 10
plot[0:L] d2uanal(x) title "Analytic", d2ual(x) title "Arc length", d2umin(x) title "Minimum", d2uhm(x) title "Harmonic mean", d2ugm(x) title "Geometric mean", d2uam(x) title "Arithmetic mean", d2urms(x) title "Root Mean Square", d2umax(x) title "Maximum"
