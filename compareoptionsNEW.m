clear all
clc

D=1;

mumu=15/1000;
alphaalpha=365/5;
sigmasigma=365/15;
gg=0.67;

inputStruc.eta=1;
% inputStruc.eta=0.2;
inputStruc.vL= 365*[0    D];
inputStruc.vB= 365*[D    D];
inputStruc.vI= 365*[D    D];   
inputStruc.vBt=365*[0    0];
inputStruc.vIt=365*[1/60 0];  
inputStruc.vLb=365*[0,   D];
inputStruc.vBb=365*[0,   D];
inputStruc.vIb=365*[1/60 D];
inputStruc.mu=mumu;
inputStruc.alpha=alphaalpha;
inputStruc.sigma=sigmasigma;
inputStruc.g=gg;

inputStruc2.eta=1;
% inputStruc2.eta=0.2;
inputStruc2.vL= 365*[0    0];
inputStruc2.vB= 365*[D    1/3];
inputStruc2.vI= 365*[D    1/21];   
inputStruc2.vBt=365*[0    0];
inputStruc2.vIt=365*[1/60 0];  
inputStruc2.vLb=365*[0,   0];
inputStruc2.vBb=365*[0,   1/3];
inputStruc2.vIb=365*[1/60 1/21];
inputStruc2.mu=mumu;
inputStruc2.alpha=alphaalpha;
inputStruc2.sigma=sigmasigma;
inputStruc2.g=gg;

option=1; %  CT stable (no mda)
%option=2;  %3 rounds of mda CT stable
% option=3;  %4 rounds of MDA CT stable

justrun=0;   % this is for doing the clines picture
%justrun=1; % just run a clinical treatment figure
% % justrun=2; % just run a mda (3 rounds) figure
% justrun=3; % just run a mda (4 rounds) figure

%use for clines (quite slow)
fVEC=[6:0.025:20];
pVEC=[0.01:0.01:0.04,0.05:0.05:0.95,0.99]; %CT coverage
% mdaVEC=0.6; (%or 0)
mdaVEC=0;

%testing fasterclines
fVEC=[6:0.05:70];
fVEC=[6:0.05:20];
pVEC=[0.01:0.01:0.04,0.05:0.05:0.95,0.99]; %CT coverage
% mdaVEC=0.6; (%or 0)
mdaVEC=0;

%use for ct comparison
% fVEC=9.17; % 30% prev for oscillating
% pVEC=0.2;
% pVEC=0.4;  %CT 
% pVEC=0.6;
% mdaVEC=0;

%for mda comparisons
% fVEC=8.14;
% fVEC=8.72;% 
% pVEC=0.4;
% pVEC=0.3;
% pVEC=0.2;

% mdaVEC=0;
%  mdaVEC=0.5;
% % mdaVEC=0.9;

%-------------
% wVEC=1.07; % thing says 1
% zVEC=0.5;

wVEC=0.8; % thing says 1
zVEC=0.8;

if option==1  %just clinical treatment
    
time0=-50;               % in years
time(1)=0;               % start first mda
time(2)=35;              % just treat

elseif (option==2)  %mda 3 rounds changed this from just option 2
    
addon=0.0694;
time0=-50;               % in years
time(1)=0;               % start first mda
time(2)=addon;           % just treat
time(3)=1/12;            % start second mda
time(4)=time(3)+addon;   % just treat
time(5)=2/12;            % start third mda
time(6)=time(5)+addon;   % just treat
time(7)=35;

elseif (option==3) % mda=4 rounds % have changed this 

addon=0.0694;
time0=-50;               % in years
time(1)=0;               % start first mda
time(2)=addon;           % just treat
time(3)=1/12;            % start second mda
time(4)=time(3)+addon;   % just treat
time(5)=2/12;            % start third mda
time(6)=time(5)+addon;   % just treat
time(7)=3/12;            % start fourth mda
time(8)=time(7)+addon;   % just treat
time(9)=8;
end

%%%%  x(1)   x(2)   x(3)   x(4)   x(5)   x(6)   x(7)   x(8)   x(9)   x(10)
%%%%  x(11)  x(12)  x(13)  x(14)  x(15)  x(16)  x(17)  x(18)  x(19)  x(20)

%     S      L      B     Bt       I     It     R      Lb     Bb     Ib    
x0=40*[10,    5,     5,    5,    5,    5,    35,    10,    10,   10];

limitthing=0*pVEC;

howmany=(length(fVEC))*(length(wVEC))*(length(zVEC))*(length(pVEC))*(length(mdaVEC));
prevend=zeros(length(fVEC),length(wVEC),length(zVEC),length(pVEC),length(mdaVEC));
k=0;
for i1=1:length(fVEC)
    inputStruc.f=fVEC(i1);
    inputStruc2.f=fVEC(i1);
    limit=0;
    for i2=1:length(wVEC)
        inputStruc.w=wVEC(i2);
        inputStruc2.w=wVEC(i2);
        for i3=1:length(zVEC)
            
            inputStruc.z=zVEC(i3);
            inputStruc2.z=zVEC(i3);
            
            inputStruc.p=0;
            inputStruc2.p=0;
            theta0=0; % no treat
            
            inputStruc.theta=theta0;
            inputStruc2.theta=theta0;            
            x0start=newstart([x0,x0],theta0,inputStruc.p);
            [TT0,ZZ0]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time0,time(1)], x0start);
            
            %%%%  x(1)   x(2)   x(3)   x(4)   x(5)   x(6)   x(7)   x(8)   x(9)   x(10)
           %%%%  x(11)  x(12)  x(13)  x(14)  x(15)  x(16)  x(17)  x(18)  x(19)  x(20)
             %     S      L      B     Bt       I     It     R      Lb     Bb     Ib    
            
            Uinf0=ZZ0(:,1)+ZZ0(:,7)+ZZ0(:,11)+ZZ0(:,15);
            Ainf0=ZZ0(:,2)+ZZ0(:,8)+ZZ0(:,9)+ZZ0(:,10)+ZZ0(:,12)+ZZ0(:,16)+ZZ0(:,17)+ZZ0(:,18);
            Sinf0=ZZ0(:,3)+ZZ0(:,4)+ZZ0(:,5)+ZZ0(:,6)+ZZ0(:,13)+ZZ0(:,14);
            
            fracASYstart0=Ainf0./(Sinf0+Ainf0);
            fracASYstart=fracASYstart0(length(TT0));
            
            total0=sum(ZZ0,2);
            prev0=(Sinf0+Ainf0)./total0;
            prevstart=prev0(length(TT0));
            
            for i4=1:length(mdaVEC)
                inputStruc.mda=mdaVEC(i4);
                inputStruc2.mda=mdaVEC(i4);
                for i5=1:length(pVEC)
                    inputStruc.p=pVEC(i5); %p=proportion who get clinical treatment
                    inputStruc2.p=pVEC(i5); %p=proportion who get clinical treatment                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if option==1
                        
                    theta1=0; % mda
                    
                    inputStruc.theta=theta1;
                    inputStruc2.theta=theta1;
                    
                    x1start=newstart(ZZ0(end,:),theta1,inputStruc.p);
                    [TT1,ZZ1]= ode15s(@(t,x) MDAode2(t,x,inputStruc), [time(1),time(2)], x1start);
                    
                    TT=[TT0;TT1];
                    ZZ=[ZZ0;ZZ1];
                    
                    elseif (option==2)||(option==4)||(option==6) %three rounds of MDA

                    theta1=inputStruc.mda; % MDA
                    theta2=0;              % no MDA
                    theta3=inputStruc.mda; % MDA+clinical
                    theta4=0;              % no MDA
                    theta5=inputStruc.mda; % MDA+clinical
                    theta6=0;              % no MDA
                                  
                    %%%%%%%
                    inputStruc.theta=theta1;
                    x1start=newstart(ZZ0(end,:),theta1,inputStruc.p);
                    [TT0b,ZZ0b]= ode15s(@(t,x) MDAode2(t,x,inputStruc), time(1)+[0,3/365], x1start);
                    x1startb=ZZ0b(end,:);
                    [TT1,ZZ1]= ode15s(@(t,x) MDAode2(t,x,inputStruc2), [time(1)+3/365,time(2)], x1startb);
                                       
                    inputStruc.theta=theta2;
                    x2start=newstart(ZZ1(end,:),theta2,inputStruc.p);
                    [TT2,ZZ2]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(2),time(3)], x2start);
                    
                    inputStruc.theta=theta3;
                    x3start=newstart(ZZ2(end,:),theta3,inputStruc.p);
                    [TT2b,ZZ2b]= ode15s(@(t, x) MDAode2(t, x, inputStruc), time(3)+[0,3/365], x3start);
                    x3startb=ZZ2b(end,:);
                    [TT3,ZZ3]= ode15s(@(t, x) MDAode2(t, x, inputStruc2), [time(3)+3/365,time(4)], x3startb);
                    
                    inputStruc.theta=theta4;
                    x4start=newstart(ZZ3(end,:),theta4,inputStruc.p);
                    [TT4,ZZ4]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(4),time(5)], x4start);
                       
                    inputStruc.theta=theta5;
                    x5start=newstart(ZZ4(end,:),theta5,inputStruc.p);
                    [TT4b,ZZ4b]= ode15s(@(t, x) MDAode2(t, x, inputStruc), time(5)+[0,3/365], x5start);
                    x5startb=ZZ4b(end,:);
                    [TT5,ZZ5]= ode15s(@(t, x) MDAode2(t, x, inputStruc2), [time(5)+3/365,time(6)], x5startb);
                    
                    inputStruc.theta=theta6;
                    x6start=newstart(ZZ5(end,:),theta6,inputStruc.p);
                    [TT6,ZZ6]= ode15s(@(t, x) MDAode2(t, x, inputStruc),[time(6),time(7)], x6start);
                                        
                    TT=[TT0;TT0b;TT1;TT2;TT2b;TT3;TT4;TT4b;TT5;TT6];
                    ZZ=[ZZ0;ZZ0b;ZZ1;ZZ2;ZZ2b;ZZ3;ZZ4;ZZ4b;ZZ5;ZZ6];
                    
                    elseif (option==3) % 4 rounds of MDA
                        
                    theta1=inputStruc.mda; % MDA round 1
                    theta2=0  ; % clinical
                    theta3=inputStruc.mda; % MDA round 2
                    theta4=0  ; % clinical
                    theta5=inputStruc.mda; % MDA round 3
                    theta6=0  ; % clinical
                    theta7=inputStruc.mda; % MDA round 4
                    theta8=0  ; % clinical                    
                    
                    inputStruc.theta=theta1;
                    x1start=newstart(ZZ0(end,:),theta1,inputStruc.p);
                    [TT0b,ZZ0b]= ode15s(@(t,x) MDAode2(t,x,inputStruc), time(1)+[0,7/365], x1start);
                    x1startb=ZZ0b(end,:);
                    [TT1,ZZ1]= ode15s(@(t,x) MDAode2(t,x,inputStruc2), [time(1)+7/365,time(2)], x1startb);
                                       
                    inputStruc.theta=theta2;
                    x2start=newstart(ZZ1(end,:),theta2,inputStruc.p);
                    [TT2,ZZ2]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(2),time(3)], x2start);

                    inputStruc.theta=theta3;
                    x3start=newstart(ZZ2(end,:),theta3,inputStruc.p);
                    [TT2b,ZZ2b]= ode15s(@(t, x) MDAode2(t, x, inputStruc), time(3)+[0,7/365], x3start);
                    x3startb=ZZ2b(end,:);
                    [TT3,ZZ3]= ode15s(@(t, x) MDAode2(t, x, inputStruc2), [time(3)+7/365,time(4)], x3startb);
                    
                    inputStruc.theta=theta4;
                    x4start=newstart(ZZ3(end,:),theta4,inputStruc.p);
                    [TT4,ZZ4]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(4),time(5)], x4start);
                    
                    inputStruc.theta=theta5;
                    x5start=newstart(ZZ4(end,:),theta5,inputStruc.p);
                    [TT4b,ZZ4b]= ode15s(@(t, x) MDAode2(t, x, inputStruc), time(5)+[0,7/365], x5start);
                    x5startb=ZZ4b(end,:);
                    [TT5,ZZ5]= ode15s(@(t, x) MDAode2(t, x, inputStruc2), [time(5)+7/365,time(6)], x5startb);
                    
                    inputStruc.theta=theta6;
                    x6start=newstart(ZZ5(end,:),theta6,inputStruc.p);
                    [TT6,ZZ6]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(6),time(7)], x6start);
                    
                    inputStruc.theta=theta7;
                    x7start=newstart(ZZ6(end,:),theta7,inputStruc.p);
                    [TT6b,ZZ6b]= ode15s(@(t, x) MDAode2(t, x, inputStruc), time(7)+[0,7/365], x7start);
                    x7startb=ZZ6b(end,:);
                    [TT7,ZZ7]= ode15s(@(t, x) MDAode2(t, x, inputStruc2), [time(7)+7/365,time(8)], x7startb);
                    
                    inputStruc.theta=theta8;
                    x8start=newstart(ZZ7(end,:),theta8,inputStruc.p);
                    [TT8,ZZ8]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(8),time(9)], x8start);

                    
                    TT=[TT0;TT0b;TT1;TT2;TT2b;TT3;TT4;TT4b;TT5;TT6;TT6b;TT7;TT8];
                    ZZ=[ZZ0;ZZ0b;ZZ1;ZZ2;ZZ2b;ZZ3;ZZ4;ZZ4b;ZZ5;ZZ6;ZZ6b;ZZ7;ZZ8];
                    end
   
                                        
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
                    total=sum(ZZ,2);
                   
                    Uinf=ZZ(:,1)+ZZ(:,7)+ZZ(:,11)+ZZ(:,15);
                    Ainf=ZZ(:,2)+ZZ(:,8)+ZZ(:,9) +ZZ(:,10)+ZZ(:,12)+ZZ(:,16)+ZZ(:,17)+ZZ(:,18);
                    Sinf=ZZ(:,3)+ZZ(:,4)+ZZ(:,5) +ZZ(:,6) +ZZ(:,13)+ZZ(:,14);
            

                    prev=(Sinf+Ainf)./total; 
                    
                    if justrun==1; %just run a single clinical treatment figure
%                         figure(20)
%                         hold on; plot(TT,100*Ainf./total,'k-')
%                         hold on; plot(TT,100*Sinf./total,'r-')
%                         hold on; area([time0,  time(1),time(1),time0  ],100*[0,0,1,1],'facecolor','g','LineStyle','none')
%                         hold on; area([time(1),time(2),time(2),time(1)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
%                         hold on; plot(TT,100*Ainf./total,'k-')
%                         hold on; plot(TT,100*Sinf./total,'r-')
%                         legend('asymptomatic','symptomatic')
%                         axis([-4,10,0,100])
%                         xlabel('Time (years)')
%                         ylabel('Time (years)')
                        
                        figure(1)
%                         hold on; area([time0,  time(1),time(1),time0  ],100*[0,0,1,1],'facecolor','g','LineStyle','none')
%                         hold on; area([time(1),time(2),time(2),time(1)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; plot(TT,100*prev,'k-')
                        axis([-5,35,0,100])
%                         axis([11,12,0.05,0.15])
%                          hold on; plot([-4,16],[0.01,0.01])
                        xlabel('Time (years)')
                        ylabel('Malaria prevalence (%)')
                        
                    elseif justrun==2; %just run a 3 round mda figure
                        figure(1)   
% 
%                         hold on; area([time0,  time(1),time(1),time0  ],100*[0,0,1,1],'facecolor','g','LineStyle','none')
%                         hold on; area([time(1),time(2),time(2),time(1)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
%                         hold on; area([time(2),time(3),time(3),time(2)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
%                         hold on; area([time(3),time(4),time(4),time(3)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
%                         hold on; area([time(4),time(5),time(5),time(4)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
%                         hold on; area([time(5),time(6),time(6),time(5)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
%                         hold on; area([time(6),time(7),time(7),time(6)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                             hold on; plot(TT,100*prev,'r-')
%                         axis([-0.2,1,0,50])
                         axis([-5,35,0,100])
%                          axis([-0.1,0.3,0,0.3])
                         axis([7,8,0.05,0.15])
                        
                        xlabel('Time (years)')
                        ylabel('Malaria prevalence (%)')
                        
                    elseif justrun==3
                        figure(7)
                        hold on; area([time0,  time(1),time(1),time0  ],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(1),time(2),time(2),time(1)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; area([time(2),time(3),time(3),time(2)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(3),time(4),time(4),time(3)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; area([time(4),time(5),time(5),time(4)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(5),time(6),time(6),time(5)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; area([time(6),time(7),time(7),time(6)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(7),time(8),time(8),time(7)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; area([time(8),time(9),time(9),time(8)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; plot(TT,100*prev,'k-')
                        axis([-4,30,0,100])
                        xlabel('Time (years)')
                        ylabel('Malaria prevalence (%)')
                     end
%                     
                    Lveryend=length(TT);
                    if (option==1)||(option==2)
                       if(prev(Lveryend)<0.001)    
                        limit=prevstart;
                        limitthing(i5)=limit;
                       end
                    end
                    k=k+1;
                    farthrough=100*k/howmany
                 end
            end
        end
    end
end

if justrun==0
    figure(1)
    xlabel('Malaria prevalence (%)')
    ylabel('Clinical treatment coverage (%)')
    axis([0,100,0,100])
    hold on; plot(limitthing*100,pVEC*100,'k-')
end
