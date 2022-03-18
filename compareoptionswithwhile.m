clear all
clc

vvL=365/5;
vvB=365/5;
% vvI=365/1;
vvI=365/5*2;
I00=365/60; %(check this as it might be the sum of two)

inputStruc.vL= [0,vvL];
inputStruc.vB= [vvB,vvB];
inputStruc.vI= [vvI,vvI]+I00;    
inputStruc.vBt=[0,vvB];
inputStruc.vIt=[0,vvI]+I00;   
inputStruc.vLb=[0,vvL];
inputStruc.vBb=[0,vvB];
inputStruc.vIb=[0,vvI]+I00;
inputStruc.mu=15/1000;
inputStruc.alpha=365/5;
inputStruc.sigma=365/15;
inputStruc.g=0.67;

option=1; %CT
option=2; %3 mda
% option=3 %4 mda

 justrun=0;
% justrun=1;  % clinical
justrun=2;  % mda 3
% justrun=3;  %mda 4


% fVEC=[6:0.1:7];
% fVEC=[5:0.1:15];

fVEC=[5:0.05:10]
% fVEC=8;
wVEC=[0.2];
zVEC=[0.95];
% mdaeffVEC=0.02:0.02:0.98 %original
mdaVEC=0.7; %0.7;
% pVEC=[0.04:0.02:0.94]; %CT coverage
pVEC=[0.025,0.05:0.05:0.95,0.99];
%pVEC=0.9

if option==1  %just clinical treatment
    
time0=-50;               % in years
time(1)=0;               % start first mda
time(2)=15;              % just treat

elseif option==2   %mda 3 rounds
    
addon=0.0694;
time0=-50;               % in years
time(1)=0;               % start first mda
time(2)=addon;           % just treat
time(3)=1/12;            % start second mda
time(4)=time(3)+addon;   % just treat
time(5)=2/12;            % start third mda
time(6)=time(5)+addon;   % just treat
time(7)=15;

else % mda=4 rounds

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
time(9)=15;
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
    limit=0;
    for i2=1:length(wVEC)
        inputStruc.w=wVEC(i2);
        for i3=1:length(zVEC)
            
            inputStruc.z=zVEC(i3);
            
            inputStruc.p=[0,0];
            theta0=0; % no treat
            
            inputStruc.theta=theta0;
            x0start=newstart([x0,x0],theta0,inputStruc.p);
            [TT0,ZZ0]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time0,time(1)], x0start);
            
            %%%%  x(1)   x(2)   x(3)   x(4)   x(5)   x(6)   x(7)   x(8)   x(9)   x(10)
           %%%%  x(11)  x(12)  x(13)  x(14)  x(15)  x(16)  x(17)  x(18)  x(19)  x(20)
             %     S      L      B     Bt       I     It     R      Lb     Bb     Ib    
            
            Uinf0=ZZ0(:,1)+ZZ0(:,7)+ZZ0(:,11)+ZZ0(:,17);
            Ainf0=ZZ0(:,2)+ZZ0(:,8)+ZZ0(:,9)+ZZ0(:,10)+ZZ0(:,12)+ZZ0(:,18)+ZZ0(:,19)+ZZ0(:,20);
            Sinf0=ZZ0(:,3)+ZZ0(:,4)+ZZ0(:,5)+ZZ0(:,6)+ZZ0(:,13)+ZZ0(:,14)+ZZ0(:,15)+ZZ0(:,16);
            
            fracASYstart0=Ainf0./(Sinf0+Ainf0);
            fracASYstart=fracASYstart0(length(TT0));
            
            total0=sum(ZZ0,2);
            prev0=(Sinf0+Ainf0)./total0;
            prevstart=prev0(length(TT0));
            
            for i4=1:length(mdaVEC)
                inputStruc.mda=mdaVEC(i4);
                
                for i5=1:length(pVEC)
                    inputStruc.p=[pVEC(i5),1]; %p=proportion who get clinical treatment
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if option==1
                        
                    theta1=0; % mda
                    
                    inputStruc.theta=theta1;
                    x1start=newstart(ZZ0(end,:),theta1,inputStruc.p);
                    [TT1,ZZ1]= ode15s(@(t,x) MDAode2(t,x,inputStruc), [time(1),time(2)], x1start);
                    
                    TT=[TT0;TT1];
                    ZZ=[ZZ0;ZZ1];
                    
                    elseif option==2

                    theta1=inputStruc.mda; % MDA
                    theta2=0;              % no MDA
                    theta3=inputStruc.mda; % MDA+clinical
                    theta4=0;              % no MDA
                    theta5=inputStruc.mda; % MDA+clinical
                    theta6=0;              % no MDA
                
                    inputStruc.theta=theta1;
                    x1start=newstart(ZZ0(end,:),theta1,inputStruc.p);
                    [TT1,ZZ1]= ode15s(@(t,x) MDAode2(t,x,inputStruc), [time(1),time(2)], x1start);
                    
                    inputStruc.theta=theta2;
                    x2start=newstart(ZZ1(end,:),theta2,inputStruc.p);
                    [TT2,ZZ2]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(2),time(3)], x2start);
                    
                    inputStruc.theta=theta3;
                    x3start=newstart(ZZ2(end,:),theta3,inputStruc.p);
                    [TT3,ZZ3]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(3),time(4)], x3start);
                    
                    inputStruc.theta=theta4;
                    x4start=newstart(ZZ3(end,:),theta4,inputStruc.p);
                    [TT4,ZZ4]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(4),time(5)], x4start);
                    
                    inputStruc.theta=theta5;
                    x5start=newstart(ZZ4(end,:),theta5,inputStruc.p);
                    [TT5,ZZ5]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(5),time(6)], x5start);
                    
                    inputStruc.theta=theta6;
                    x6start=newstart(ZZ5(end,:),theta6,inputStruc.p);
                    [TT6,ZZ6]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(6),time(7)], x6start);
                                        
                    TT=[TT0;TT1;TT2;TT3;TT4;TT5;TT6];
                    ZZ=[ZZ0;ZZ1;ZZ2;ZZ3;ZZ4;ZZ5;ZZ6];
                    
                    elseif option==3
                        
                    theta1=inputStruc.mda; % MDA round 1
                    theta2=0  ; % clinical
                    theta3=inputStruc.mda; % MDA round 2
                    theta4=0  ; % clinical
                    theta5=inputStruc.mda; % MDA round 3
                    theta6=0  ; % clinical
                    theta7=inputStruc.mda; % MDA round 3
                    theta8=0  ; % clinical
                
                    inputStruc.theta=theta1;
                    x1start=newstart(ZZ0(end,:),theta1,inputStruc.p);
                    [TT1,ZZ1]= ode15s(@(t,x) MDAode2(t,x,inputStruc), [time(1),time(2)], x1start);
                    
                    inputStruc.theta=theta2;
                    x2start=newstart(ZZ1(end,:),theta2,inputStruc.p);
                    [TT2,ZZ2]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(2),time(3)], x2start);
                    
                    inputStruc.theta=theta3;
                    x3start=newstart(ZZ2(end,:),theta3,inputStruc.p);
                    [TT3,ZZ3]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(3),time(4)], x3start);
                    
                    inputStruc.theta=theta4;
                    x4start=newstart(ZZ3(end,:),theta4,inputStruc.p);
                    [TT4,ZZ4]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(4),time(5)], x4start);
                    
                    inputStruc.theta=theta5;
                    x5start=newstart(ZZ4(end,:),theta5,inputStruc.p);
                    [TT5,ZZ5]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(5),time(6)], x5start);
                    
                    inputStruc.theta=theta6;
                    x6start=newstart(ZZ5(end,:),theta6,inputStruc.p);
                    [TT6,ZZ6]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(6),time(7)], x6start);
                    
                    inputStruc.theta=theta7;
                    x7start=newstart(ZZ6(end,:),theta7,inputStruc.p);
                    [TT7,ZZ7]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(7),time(8)], x5start);
                    
                    inputStruc.theta=theta8;
                    x8start=newstart(ZZ7(end,:),theta8,inputStruc.p);
                    [TT8,ZZ8]= ode15s(@(t, x) MDAode2(t, x, inputStruc), [time(8),time(9)], x6start);
                    
                    TT=[TT0;TT1;TT2;TT3;TT4;TT5;TT6;TT7;TT8];
                    ZZ=[ZZ0;ZZ1;ZZ2;ZZ3;ZZ4;ZZ5;ZZ6;ZZ7;ZZ8];
                    
                    end
                                        
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
                    total=sum(ZZ,2);

                    Uinf=ZZ(:,1)+ZZ(:,7)+ZZ(:,11)+ZZ(:,17);
                    Ainf=ZZ(:,2)+ZZ(:,8)+ZZ(:,9) +ZZ(:,10)+ZZ(:,12)+ZZ(:,18)+ZZ(:,19)+ZZ(:,20);
                    Sinf=ZZ(:,3)+ZZ(:,4)+ZZ(:,5) +ZZ(:,6) +ZZ(:,13)+ZZ(:,14)+ZZ(:,15)+ZZ(:,16);
                    prev=(Sinf+Ainf)./total; 
                    
                    if justrun==1;
                        figure(20)
                        hold on; plot(TT,100*Ainf./total,'k-')
                        hold on; plot(TT,100*Sinf./total,'r-')
                        hold on; area([time0,  time(1),time(1),time0  ],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(1),time(2),time(2),time(1)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; plot(TT,100*Ainf./total,'k-')
                        hold on; plot(TT,100*Sinf./total,'r-')
                        legend('asymptomatic','symptomatic')
                        axis([-4,10,0,100])
                        xlabel('Time (years)')
                        ylabel('Time (years)')
                        
                        figure(21)
                        hold on; area([time0,  time(1),time(1),time0  ],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(1),time(2),time(2),time(1)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; plot(TT,100*prev,'k-')
                        axis([-4,10,0,100])
                        xlabel('Time (years)')
                        ylabel('Malaria prevalence (%)')
                        
                    elseif justrun==2;
                        figure(6)   

                        hold on; area([time0,  time(1),time(1),time0  ],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(1),time(2),time(2),time(1)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; area([time(2),time(3),time(3),time(2)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(3),time(4),time(4),time(3)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; area([time(4),time(5),time(5),time(4)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; area([time(5),time(6),time(6),time(5)],100*[0,0,1,1],'facecolor','y','LineStyle','none')
                        hold on; area([time(6),time(7),time(7),time(6)],100*[0,0,1,1],'facecolor','g','LineStyle','none')
                        hold on; plot(TT,100*prev,'k-')
                        axis([-4,10,0,100])
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
                        axis([-4,10,0,100])
                        xlabel('Time (years)')
                        ylabel('Malaria prevalence (%)')
                    end
                    
                        
                    Lveryend=length(TT);
                    
                    if(prev(Lveryend)<0.001)
                        
                         limit=prevstart;
                         limitthing(i5)=limit;
                         %   hold on; plot(100*fracASYstart,100*ct,'kx')
                         %                         figure(2); hold on; plot(100*fracASYstart,100*mda,'kx')
                         %                         figure(4); hold on; plot(100*prevstart,100*mda,'kx')
                    else
                        
                        %   hold on; plot(100*fracASYstart,100*ct,'rx')
                        %                         figure(2); hold on; plot(100*fracASYstart,100*mda,'rx')
                        %                         figure(4); hold on; plot(100*prevstart,100*mda,'rx')
                    end
                    k=k+1;
                    farthrough=100*k/howmany
                end
            end
        end
    end
end

% figure(2)
% xlabel('% aysomptomatic')
% ylabel('% MDA coverage')
% axis([0,100,0,100])
% % hold on; plot(TT,()./total,'r')

if justrun==0
    figure(4)
    xlabel('Malaria prevalence (%)')
    ylabel('Clinical treatment coverage (%)')
    axis([0,100,0,100])
    % hold on; plot(TT,()./total,'r')
    hold on; plot(limitthing*100,pVEC*100,'r-')
end
