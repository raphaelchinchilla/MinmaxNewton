close all
clear all
addpath('Algorithm')
delete('temp*');rc=rmdir('@temp*','s');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This algorithm implements the classic homicidal chauffeur problem.
% The pursuer is a driving Dubins vehicle and the evader is an
% pedestrian, represented by a simple accumulator.

% System dynamics:
% pursuer (minimizer)
% x1(k+1)    = x1(k) + v cos theta(k) 
% x2(k+1)    = x2(k) + v sin theta(k) 
% theta(k+1) = theta(k) + u(k)         u(k) in [-uMax,uMax]

% evader (maximizer)
% y1(k+1)    = y1(k) + d1(k)           d1(k), d2(k) in [-dMax,dMax]
% y2(k+1)    = y2(k) + d2(k)  

% controlled output
% z1 = y1 - x1
% z2 = y2 - x2

T=50;% forward horizon  (T=50 gives nicer results)

% Building symbolic problem

nXmin=3;
nXmax=2;
nU=1;
nD=2;

Tvariable v [];  % pursuer's velocity

dXfun_min=@(x,u,v)[x(1,:)+v*cos(x(3,:));
               x(2,:)+v*sin(x(3,:));
               x(3,:)+u];
           
dXfun_max=@(x,d)x+d;         

Tvariable x0min    [nXmin,1];     % xmin(t)
Tvariable x0max    [nXmax,1];     % xmax(t)

Tvariable x1min [nXmin,T] % xmin(t+1), ... , xmin(t+T)
Tvariable x1max [nXmax,T] % xmax(t+1), ... , xmax(t+T)

Tvariable d     [nD,T];   % d(t), ... , d(t+T-1)
Tvariable u  [nU,T];     % u(t),   ... , u(t+T-1);

xmin=[x0min,x1min]; % xmin(t), ... , xmin(t+T)
xmax=[x0max,x1max]; % xmax(t), ... , xmax(t+T)


dynamicsMin = x1min==dXfun_min(xmin(:,1:end-1),u,v);
dynamicsMax = x1max==dXfun_max(xmax(:,1:end-1),d);

% parameters for state and input constraints
Tvariable uMax [];
Tvariable dMax [];

% Criterion

JJ=[norm2(x1min(1:2,:)-x1max)/(T^2); % distance from pursuer to evader
    norm2(u)/(T^2);              % pursuer control
    -norm2(d)/(T^2)];               % disturbance
Tvariable cc size(JJ);

J=cc*JJ;

% Warm start
tol=.9;
uWarm=[u(:,2:end),zeros(nU,1)];
uWarm=min(uWarm,tol*uMax);
uWarm=max(uWarm,-tol*uMax);
dWarm=[d(:,2:end),zeros(nD,1)];
dWarm=min(dWarm,tol*dMax);
dWarm=max(dWarm,-tol*dMax);

% Generating solver

classname=['temp_pursuit_evasion',num2str(T)];

objective=J;
minimizationVariables={u,x1min};
maximizationVariables={d,x1max};

minimizationConstraints={u<=uMax,u>=-uMax,dynamicsMin};
maximizationConstraints={sum(d.*d,1)<=dMax*dMax,dynamicsMax};
parameters={uMax,dMax,cc,v,x0min,x0max};
outputExpressions={J,JJ,u(:,1),d(:,1),xmin,xmax,full(uWarm),full(dWarm)};
code_type='c';


generate_tens_functions(classname,objective,minimizationVariables,maximizationVariables,minimizationConstraints,maximizationConstraints,outputExpressions,parameters,code_type)
obj=feval(classname);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v=.1;
setP_v(obj,v);
uMax=.5;
setP_uMax(obj,uMax);
dMax=.06;
setP_dMax(obj,dMax);



SigmaMin=.00;
SigmaMax=.00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Weights associated with each term of the cost
uWeight=80;
dWeight=1000;
cc=[1;uWeight;dWeight];
setP_cc(obj,cc);


params_optim.delta_L=1e-7;
params_optim.delta_G=1e-8;
params_optim.delta_F=1e-8;
params_optim.delta_gap=1e-7;
params_optim.gamma_agressive=0.7;
params_optim.gamma_conservative=0.95;
mu0=1;
params_optim.max_iter=500;
params_optim.decouple_alpha=true;
params_optim.gamma_py=1e-4;



s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% set initial condition

t0=0;
x0min=[0;0;0];
x0max=[.5;.5];

uWarm=zeros(nU,T);
dWarm=zeros(nD,T);
lastu=uWarm;
lastuT=0;

 setV_u(obj,uWarm);
 setV_d(obj,dWarm);
 % fix xWarm
 xWarmMin(:,1)=dXfun_min(x0min,uWarm(:,1),v);
 xWarmMax(:,1)=dXfun_max(x0max,dWarm(:,1));     

 for i=2:T
     xWarmMin(:,i)=dXfun_min(xWarmMin(:,i-1),uWarm(:,i),v);
     xWarmMax(:,i)=dXfun_max(xWarmMax(:,i-1),dWarm(:,i));
 end
 setP_x0min(obj,x0min);
 setV_x1min(obj,xWarmMin);

 setP_x0max(obj,x0max);
 setV_x1max(obj,xWarmMax);     




xWarmMin=zeros(nXmin,T);
xWarmMax=zeros(nXmax,T);


closedloop.t=t0;
closedloop.xmin=x0min;
closedloop.xmax=x0max;
closedloop.status=[];
closedloop.iter=[];
closedloop.stime=[];
closedloop.J=[];
closedloop.JJ=zeros(length(JJ),0);
closedloop.u=[];
closedloop.d=[];

figure(1);clf;
set(1,'Name','Trajectory');
figure(2);clf;
set(2,'Name','MPC solutions');
figure(3);clf;
set(3,'Name','Solver');

for k=1:200
    

    % MPC optimization    

    [closedloop.status(end+1,1),...
     closedloop.iter(end+1,1),...
     closedloop.stime(end+1,1)]=ip_newton_minmax(obj,mu0,params_optim);

    [closedloop.J(end+1,1),...
     closedloop.JJ(:,end+1),...
     u0,d0,hatxmin,hatxmax,uWarm,dWarm]=getOutputs(obj);

    fprintf('t=%g, status=%g, J=%g computed in %g iterations & %g ms\n',...
            closedloop.t(end),closedloop.status(end),closedloop.J(end),closedloop.iter(end),1e3*closedloop.stime(end));

    if closedloop.status(end)>0
        fprintf('solver failed at time %g \n',...
                closedloop.t(end));
    end
    
    
    % Apply control
     closedloop.u(:,end+1)=u0;     
     % Apply disturbance 
     if closedloop.t(end)<=55
         closedloop.d(:,end+1)=[.05;0];
     else
         closedloop.d(:,end+1)=d0;
     end
    % Current state with perturbations
     closedloop.xmin(:,end+1)=dXfun_min(closedloop.xmin(:,end),closedloop.u(:,end),v)+[SigmaMin*randn(2,1);0];
     closedloop.xmax(:,end+1)=dXfun_max(closedloop.xmax(:,end),closedloop.d(:,end))+SigmaMax*randn(2,1);
     closedloop.t(end+1)=closedloop.t(end)+1;
     
     % apply warm start for next iteration     
     setV_u(obj,uWarm);
     setV_d(obj,dWarm);
     % fix xWarm
     xWarmMin(:,1)=dXfun_min(closedloop.xmin(:,end),uWarm(:,1),v);
     xWarmMax(:,1)=dXfun_max(closedloop.xmax(:,end),dWarm(:,1));     
     
     for i=2:T
         xWarmMin(:,i)=dXfun_min(xWarmMin(:,i-1),uWarm(:,i),v);
         xWarmMax(:,i)=dXfun_max(xWarmMax(:,i-1),dWarm(:,i));
     end
     setP_x0min(obj,closedloop.xmin(:,end));
     setV_x1min(obj,xWarmMin);
     
     setP_x0max(obj,closedloop.xmax(:,end));
     setV_x1max(obj,xWarmMax);     

     
     if  true && mod(k,5)==0
         figure(1);
         %subplot(4,4,mod(k,16)+1);
         plot(closedloop.xmin(1,:),closedloop.xmin(2,:),'g.-',...
              closedloop.xmax(1,:),closedloop.xmax(2,:),'b.-',...
              hatxmin(1,:),hatxmin(2,:),'g:',...
              hatxmax(1,:),hatxmax(2,:),'b:');grid on
         legend('pursuer','evader','pursuer prediction','evader prediction');
         axis equal
         figure(2);
         subplot(3,1,1);
         plot(closedloop.t(1:end-1),closedloop.J,'-',closedloop.t(1:end-1),closedloop.JJ,'.');grid on
         legend('J','J_x','J_u','J_d');
         subplot(3,1,2);
         plot(closedloop.t(1:end-1),closedloop.u,'.-');grid on;
         legend('pursuer \omega','location','best');
         subplot(3,1,3);
         plot(closedloop.t(1:end-1),closedloop.d,'.-');grid on
         legend('evader v_x','evader v_y','location','best');
         figure(3);
         subplot(3,1,1)
         plot(closedloop.t(1:end-1),closedloop.iter,'.');grid on
         ylabel('# iter');
         subplot(3,1,2)
         plot(closedloop.t(1:end-1),1000*closedloop.stime,'.');grid on
         ylabel('[ms]');
         subplot(3,1,3)
         plot(closedloop.t(1:end-1),closedloop.status,'*');grid on
         ylabel('status');

         drawnow;
     end
end

