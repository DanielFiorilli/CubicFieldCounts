\\ This is the pari/gp code used for the computations in the paper https://arxiv.org/abs/2102.08077 . It produces vectors containing values of the function  f_p(X,T):= max_{x\leq X} |N^+_p(x,T) - A_p^+(T) x - B_p^+(x,T) x^{5/6}  | / x^{1/2}. Here, N^+_p(x,T) denotes the number of cubic fields of Galois group S_3 of (positive) discriminant at most x in which p has splitting type T, where T is one of the 5 possible splitting types in cubic fields. For a more precise description of the splitting types as well as the constants A_p^+(T) and B_p^+(X,T), see the paper https://arxiv.org/abs/2102.08077 .
\\ This code works with development version 2.14 of pari/gp, which can be found here : https://pari.math.u-bordeaux.fr/Events/PARI2022/talks/sources.pdf
\\ The paramters which can be changed are as follows: 
\\ lim is the value of X.
\\ limp is the size of the maximial prime p to be tested.
\\ Finally, if one wants to compute values of (N^-_p(X,T) - A_p^-(T) X - B_p^-(X,T) X^{5/6}  ) / X^{1/2} (e.g. negative discriminants), then A=nflist("S3",[1,lim],0) should be replaced with A=nflist("S3",[1,lim],1), and nfdisc(A[i]) should be replaced with abs(nfdisc(A[i])) .  
\\ This code is made for relatively modest values of lim (and substantial values of limp). If one wants to compute very large values of lim with limp of small (say fixed) size, then the memory usage in the for loop should be optimized (as in the companion file cubicnfcounts.gp). 




lim=10^5;
limp=10^4;
st=3;
A=nflist("S3",[1,lim],0);
c1=1/12/zeta(3);
c2=4*zeta(1/3)/5/gamma(2/3)^3/zeta(5/3);
P=primes(primepi(limp));
limN=#P;
step=1;
N=vector(floor(limN/step),i,0);
{
for(i=1,limN/step, 
	my(p=P[i*step]);    
	my(C=vector(#A,k,idealprimedec(nfinit(A[k]),p))); my(T1=vector(#A,k,if( (st==1 && length(C[k])==3) || (st==3 && length(C[k])==1  && C[k][1][4]==3) || (st==5 && length(C[k])==1  && C[k][1][4]==1 ) ||  (st==2 && length(C[k])==2  && C[k][1].e==1 && C[k][2].e==1 ) || (st==4 && length(C[k])==2  && C[k][1].f==1 && C[k][2].f==1 ) ,1,0)));
	my(TT1= vector(lim,k,0)); for(k=1,#A,if(T1[k],TT1[nfdisc(A[k])]++;,));
	my(ET1= vector(lim,k,0)); ET1[1]=TT1[1]; for(k=2,lim,ET1[k]=ET1[k-1]+TT1[k];);
	c11=c1*(1+1/p+p^(-2))^(-1) /6; c21=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p) *(1+p^(-1/3))^3/6; 
	c13=c1*(1+1/p+p^(-2))^(-1) /3; c23=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p) /3; 
	c15=c1*(1+1/p+p^(-2))^(-1) /p^2; c25=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p^(1/3)) /p^2; 
	c12=c1*(1+1/p+p^(-2))^(-1)/2; c22=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p^(1/3)) *(1+1/p^(2/3)) /2; 
	c14=c1*(1+1/p+p^(-2))^(-1)/p; c24=c2*(1-p^(-1/3))/(1-p^(-5/3))/(1+1/p)*(1+1/p^(1/3))^2 /p;
	for(k=1,lim,ET1[k]=abs(ET1[k] -if(st==1,c11* k + c21* (k)^(5/6) ,0)- if(st==2,c12* k + c22* (k)^(5/6) ,0) -if(st==3,c13* k + c23* (k)^(5/6) ,0) -if(st==4,c14* k + c24* (k)^(5/6) ,0)  -if(st==5,c15* k + c25* (k)^(5/6) ,0)    )/(k)^(.5) ); 
	N[i]=vecmax(ET1);  );
}



\\Finally, one can add a line of the form write("**your path**/uniformity.txt",N); to store the computations. 

