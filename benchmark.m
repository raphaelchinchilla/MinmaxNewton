clear all
addpath('Algorithm')
delete('temp*');rc=rmdir('@temp*','s');

% This code implements the benchmark examples described in our paper

%% Defining minmax problem

example=2;

Tvariable x [];
Tvariable y [];

if example==1 % from Adolphs
    cost = 2*x.^2 -y.^2 + 4*x*y + 4/3*y.^3 - 1/4*y.^4;
    delta_l=1;
    box=1;
    gamma=sqrt(eps());
elseif example==2 % from Wang
    cost = (4*x.^2-(y-3*x+0.05*x.^3).^2-0.1*y.^4)*exp(-0.01*(x.^2+y.^2));
    delta_l=1;
    box=20;
    gamma=sqrt(eps());
elseif example==3 % from Mertikopoulo
    cost = (x-0.5)*(y-0.5)+exp(-(x-0.25).^2-(y-0.75).^2);
    delta_l=1;
    box=20;
    gamma=sqrt(eps());
elseif example==4 % from Mertikopoulo 10x^2
    cost = (x-0.5)*(y-0.5)+exp(-(x-0.25).^2-(y-0.75).^2)+20*x.^2;
    delta_l=0.5;
    box=2;
    gamma=0;
end


classname='temp_';
objective=cost;
minimizationVariables={x};
maximizationVariables={y};
minimizationConstraints={};
maximizationConstraints={};
parameters={};
outputExpressions={x,y};
code_type='c';

generate_tens_functions(classname,objective,minimizationVariables,maximizationVariables,minimizationConstraints,maximizationConstraints,outputExpressions,parameters,code_type)

obj=feval(classname);                               

%% Solving examples

params_optim.max_iter=500;

pure_newton.converged=0;
pure_newton.converged_minmax=0;
pure_newton.avg_iter=0;
pure_newton.sol_x=[];
pure_newton.sol_y=[];

delta_zero.converged=0;
delta_zero.converged_minmax=0;
delta_zero.avg_iter=0;
delta_zero.sol_x=[];
delta_zero.sol_y=[];

delta_inf.converged=0;
delta_inf.converged_minmax=0;
delta_inf.avg_iter=0;
delta_inf.sol_x=[];
delta_inf.sol_y=[];

mixed.converged=0;
mixed.converged_minmax=0;
mixed.avg_iter=0;
mixed.sol_x=[];
mixed.sol_y=[];

% Finding the stationary points of the system, independent of being local
% minmax or not
stationary_x=[];
stationary_y=[];


for count=1:1e3
    xinit=box*2*(rand()-0.5); 
    yinit=box*2*(rand()-0.5);
    count
   
    % Pure Newton method
    setV_x(obj,xinit);
    setV_y(obj,yinit);    
    params_optim.adjust_eps=false;
    params_optim.gamma_px=0;
    params_optim.gamma_py=0;
    [status,iter]=ip_newton_minmax(obj,1,params_optim);
    [xsol,ysol]= getOutputs(obj);
    
    if status<1
        pure_newton.converged=pure_newton.converged+1;
        stationary_x(pure_newton.converged)=xsol;
        stationary_y(pure_newton.converged)=ysol;
        if status==0
            pure_newton.converged_minmax=pure_newton.converged_minmax+1;
            pure_newton.avg_iter=pure_newton.avg_iter+(iter-pure_newton.avg_iter)/pure_newton.converged_minmax;
            pure_newton.sol_x(pure_newton.converged_minmax)=xsol;
            pure_newton.sol_y(pure_newton.converged_minmax)=ysol;
        end
    end
    
    params_optim=rmfield(params_optim,'adjust_eps');
    
    params_optim.gamma_px=gamma;
    params_optim.gamma_py=gamma;
 
   
    
    % Delta=0
    setV_x(obj,xinit);
    setV_y(obj,yinit);    
    params_optim.adjust_eps=true;
    params_optim.delta_l=0;
    [status,iter]=ip_newton_minmax(obj,1,params_optim);
    [xsol,ysol]= getOutputs(obj);
    
    if status<1
        delta_zero.converged=delta_zero.converged+1;
        if status==0
            delta_zero.converged_minmax=delta_zero.converged_minmax+1;
            delta_zero.avg_iter=delta_zero.avg_iter+(iter-delta_zero.avg_iter)/delta_zero.converged_minmax;
            delta_zero.sol_x(delta_zero.converged_minmax)=xsol;
            delta_zero.sol_y(delta_zero.converged_minmax)=ysol;
        end
    end
   
    
    % Delta=Inf
    setV_x(obj,xinit);
    setV_y(obj,yinit);    
    params_optim.delta_l=Inf;
    [status,iter]=ip_newton_minmax(obj,1,params_optim);
    [xsol,ysol]= getOutputs(obj);
    
    if status<1
        delta_inf.converged=delta_inf.converged+1;
        if status==0
            delta_inf.converged_minmax=delta_inf.converged_minmax+1;
            delta_inf.avg_iter=delta_inf.avg_iter+(iter-delta_inf.avg_iter)/delta_inf.converged_minmax;
            delta_inf.sol_x(delta_inf.converged_minmax)=xsol;
            delta_inf.sol_y(delta_inf.converged_minmax)=ysol;
        end
    end
    
    % Mixed (different Delta depending on the problem)
    setV_x(obj,xinit);
    setV_y(obj,yinit);    
    params_optim.delta_l=delta_l;
    [status,iter]=ip_newton_minmax(obj,1,params_optim);
    [xsol,ysol]= getOutputs(obj);
    
    if status<1
        mixed.converged=mixed.converged+1;
        if status==0
            mixed.converged_minmax=mixed.converged_minmax+1;
            mixed.avg_iter=mixed.avg_iter+(iter-mixed.avg_iter)/mixed.converged_minmax;
            mixed.sol_x(mixed.converged_minmax)=xsol;
            mixed.sol_y(mixed.converged_minmax)=ysol;
        end
    end
    
end

fprintf('Pure Newton: converged %d, converged to minmax %d, avg num iter to converge to minmax %4f\n',pure_newton.converged,pure_newton.converged_minmax,pure_newton.avg_iter)
fprintf('Delta=0 : converged %d, converged to minmax %d, avg num iter to converge to minmax %4f\n',delta_zero.converged,delta_zero.converged_minmax,delta_zero.avg_iter)
fprintf('Delta=Inf : converged %d, converged to minmax %d, avg num iter to converge to minmax %4f\n',delta_inf.converged,delta_inf.converged_minmax,delta_inf.avg_iter)
fprintf('Mixed: converged %d, converged to minmax %d, avg num iter to converge to minmax %4f\n',mixed.converged,mixed.converged_minmax,mixed.avg_iter)