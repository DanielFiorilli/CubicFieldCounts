\\ This is the pari/gp code used for the computations in the paper https://arxiv.org/abs/2102.08077 . It produces vectors containing 5000 values of the function (N^+_p(X,T) - A_p^+(T) X - B_p^+(X,T) X^{5/6}  ) / X^{1/2}. Here, N^+_p(X,T) denotes the number of cubic fields of Galois group S_3 of (positive) discriminant at most X in which p (set to 5 in this code) has splitting type T, where T is one of the 5 possible splitting types in cubic fields. For a more precise description of the splitting types as well as the constants A_p^+(T) and B_p^+(X,T), see the paper https://arxiv.org/abs/2102.08077 .
\\ This code works with development version 2.14 of pari/gp, which can be found here : https://pari.math.u-bordeaux.fr/Events/PARI2022/talks/sources.pdf
\\ The paramters which can be changed are as follows: 
\\ lim is the maximal value of X.
\\ step is the distance between every two points of the function which will be computed. 
\\ p can be set to any prime.
\\ Finally, if one wants to compute values of (N^-_p(X,T) - A_p^-(T) X - B_p^-(X,T) X^{5/6}  ) / X^{1/2} (e.g. negative discriminants), then A=nflist("S3",[1,lim],0) should be replaced with A=nflist("S3",[1,lim],1), and d=nfdisc(A[i]) should be replaced with d=abs(nfdisc(A[i])) .  



lim=10^8; 
step=lim/5000;
c1=1/12/zeta(3); c2=4*zeta(1/3)/5/gamma(2/3)^3/zeta(5/3); p = 5; S1=parvector(lim/step, i,[0,0,0,0,0]);  
{my(A=nflist("S3",[1,lim],0)); my(p=p);
c11=c1*(1+1/p+p^(-2))^(-1) /6; c21=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p) *(1+p^(-1/3))^3/6; 
c13=c1*(1+1/p+p^(-2))^(-1) /3; c23=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p) /3; 
c15=c1*(1+1/p+p^(-2))^(-1) /p^2; c25=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p^(1/3)) /p^2; 
c12=c1*(1+1/p+p^(-2))^(-1)/2; c22=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p^(1/3)) *(1+1/p^(2/3)) /2; 
c14=c1*(1+1/p+p^(-2))^(-1)/p; c24=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p^(1/3))^2 /p;
for(i=1,#A, 
my(c=idealprimedec(nfinit(A[i]),p));my(d=nfdisc(A[i]));my(d=(d-d%step)/step+1);
if(#c==3,S1[d][1]++, 
if(#c==1 && c[1][4]==3 , S1[d][3]++,  
if(#c==1 && c[1][4]==1,S1[d][5]++,
if(#c==2 && c[1].e==1 && c[2].e==1,S1[d][2]++;, 
if( #c==2 && c[1].f==1 && c[2].f==1,  S1[d][4]++;,print("ERROR!!")))));))}
S2= parvector(lim/step,i,[0,0,0,0,0]); S2[1]=S1[1];
for(i=2,lim/step,S2[i]=S2[i-1]+S1[i] )
ET1=vector(lim/step,i,0); ET2=vector(lim/step,i,0); ET3=vector(lim/step,i,0); ET4=vector(lim/step,i,0); ET5=vector(lim/step,i,0);
for(i=1,lim/step,ET1[i]=(S2[i][1]-c11* i*step - c21* (i*step)^(5/6) )/(i*step)^(.5) );
for(i=1,lim/step,ET2[i]=(S2[i][2]-c12* i*step - c22* (i*step)^(5/6) )/(i*step)^(.5) );
for(i=1,lim/step,ET3[i]=(S2[i][3]-c13* i*step - c23* (i*step)^(5/6) )/(i*step)^(.5) );
for(i=1,lim/step,ET4[i]=(S2[i][4]-c14* i*step - c24* (i*step)^(5/6) )/(i*step)^(.5) );
for(i=1,lim/step,ET5[i]=(S2[i][5]-c15* i*step - c25* (i*step)^(5/6) )/(i*step)^(.5) );


\\Finally, one can add lines of the form write("**your path**/ET1.txt",ET1); write("**your path**/ET2.txt",ET2);... to store the computations. 

