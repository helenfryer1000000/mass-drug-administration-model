function [dx] = MDAode2(t,x,inputStruc)

eta=inputStruc.eta;
mu=inputStruc.mu;
  alpha=inputStruc.alpha;
  sigma=inputStruc.sigma;
      p=inputStruc.p;
      w=inputStruc.w;
     vL=inputStruc.vL;
     vB=inputStruc.vB;
     vI=inputStruc.vI;    
     vBt=inputStruc.vBt;
     vIt=inputStruc.vIt;   
     vLb=inputStruc.vLb;
     vBb=inputStruc.vBb;
     vIb=inputStruc.vIb;
     f=inputStruc.f;
     z=inputStruc.z;
     g=inputStruc.g;
theta=inputStruc.theta;

%%%%  x(1)   x(2)   x(3)   x(4)   x(5)   x(6)   x(7)   x(8)   x(9)   x(10)
%%%%  x(11)  x(12)  x(13)  x(14)  x(15)  x(16)  x(17)  x(18)  x(19)  x(20)

%     S      L      B      Bt     I      I~     R      Lb     Bb     Ib    

% Note that the twiddles are before the dashes
I=x(5)+x(6)+x(14)+ eta*(x(10)+x(18)); % infectious 
N=sum(x); %total population

%%%%%%%%%%%%%%%%%%% Transmission parameter

q=mod(t-0.05,1);
%beta=f*(g*cos(2*pi*q)+1); %original
beta=f;

dx=[(1-theta)*mu*N-beta*x(1)*I/N+(1-z)*(vL(1)*x(2)+vB(1)*x(3)+vBt(1)*x(4)+vI(1)*x(5)+vIt(1)*x(6))+w*x(7)-mu*x(1);...% clinical treatment + MDA (vs for both)
beta*x(1)*I/N-(alpha+vL(1)+mu)*x(2);...
p(1)*alpha*x(2)-(sigma+vB(1)+mu)*x(3);...
(1-p(1))*alpha*x(2)-(sigma+vBt(1)+mu)*x(4);...
sigma*x(3)-(vI(1)+mu)*x(5);...
sigma*x(4)-(vIt(1)+mu)*x(6);...
-beta*x(7)*I/N+z*(vL(1)*x(2)+vB(1)*x(3)+vBt(1)*x(4)+vI(1)*x(5)+vIt(1)*x(6))+vLb(1)*x(8)+vBb(1)*x(9)+vIb(1)*x(10)-(mu+w)*x(7)
beta*x(7)*I/N-(alpha+vLb(1)+mu)*x(8);...
alpha*x(8)-(sigma+vBb(1)+mu)*x(9);...
sigma*x(9)-(vIb(1)+mu)*x(10);...
theta*mu*N-beta*x(11)*I/N+(1-z)*(vL(2)*x(12)+vB(2)*x(13)+vI(2)*x(14))+w*x(15)-mu*x(11);...% clinical treatment + MDA (vs for both)
beta*x(11)*I/N-(alpha+vL(2)+mu)*x(12);...
alpha*x(12)-(sigma+vB(2)+mu)*x(13);...
sigma*x(13)-(vI(2)+mu)*x(14);...
-beta*x(15)*I/N+z*(vL(2)*x(12)+vB(2)*x(13)+vI(2)*x(14))+vLb(2)*x(16)+vBb(2)*x(17)+vIb(2)*x(18)-(mu+w)*x(15)
beta*x(15)*I/N-(alpha+vLb(2)+mu)*x(16);...
alpha*x(16)-(sigma+vBb(2)+mu)*x(17);...
sigma*x(17)-(vIb(2)+mu)*x(18)];
end
