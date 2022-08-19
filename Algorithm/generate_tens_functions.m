function generate_tens_functions(classname,objective,minimizationVariables,maximizationVariables,minimizationConstraints,maximizationConstraints,outputExpressions,parameters,code_type)
% This functions uses TensCalc symbolic backend to construct the gradient
% and Hessian that will be used in the solver.

    if ~exist('code_type','var')
        code_type='matlab';
    end

    
    if strcmp(code_type,'matlab')
        tprod2matlab=true;
    elseif strcmp(code_type,'c')
        tprod2matlab=false;
    else
        error('Code type needs to be either matlab or c')
    end
    
    debug=false;
    scratchbookType='double';    
    code=csparse(scratchbookType,debug,tprod2matlab); % using instructionsTable.c
    classhelp={'Create object';
                   sprintf('obj=%s();',classname)};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Packing variables and creating necessary dual and auxiliary variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Check input parameters

    parameters=checkParameters(parameters);


    if isstruct(minimizationVariables)
        minimizationVariables=struct2cell(minimizationVariables);
    end

    if ~iscell(minimizationVariables)
        minimizationVariables
        error('optimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(minimizationVariables)
        if ~isequal(class(minimizationVariables{i}),'Tcalculus')
            minimizationVariables{i}
            error('all optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(minimizationVariables{i}));
        end
        if ~isequal(type(minimizationVariables{i}),'variable')
            minimizationVariables{i}
            error('all optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(minimizationVariables{i}));
        end
        for j=1:length(parameters)
            if isequal(name(minimizationVariables{i}),name(parameters{j}))
                minimizationVariables{i}
                error('optimization variable ''%s'' cannot also be a parameter\n',name(minimizationVariables{i}));
            end
        end
    end


    if isstruct(maximizationVariables)
        maximizationVariables=struct2cell(maximizationVariables);
    end

    if ~iscell(maximizationVariables)
        maximizationVariables
        error('optimizationVariables must be a cell array of Tcalculus variables');
    end

    for i=1:length(maximizationVariables)
        if ~isequal(class(maximizationVariables{i}),'Tcalculus')
            maximizationVariables{i}
            error('all optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,class(maximizationVariables{i}));
        end
        if ~isequal(type(maximizationVariables{i}),'variable')
            maximizationVariables{i}
            error('all optimizationVariables must be of the type ''variable'' (%dth is of type ''%s'')\n',...
                  i,type(maximizationVariables{i}));
        end
        for j=1:length(parameters)
            if isequal(name(maximizationVariables{i}),name(parameters{j}))
                maximizationVariables{i}
                error('optimization variable ''%s'' cannot also be a parameter\n',name(maximizationVariables{i}));
            end
        end
    end

    for i=1:length(maximizationVariables)
        for j=1:length(minimizationVariables)
            if isequal(name(maximizationVariables{i}),name(minimizationVariables{j}))
                maximizationVariables{i}
                error('maximization variable ''%s'' cannot also be a minimization variable\n',name(maximizationVariables{i}));
            end
        end
    end

    if ~isempty(size(objective))
        error('Minimization criterion must be scalar (not [%s])',index2str(size(objective)));
    end

    if ~isempty(minimizationConstraints) && ~iscell(minimizationConstraints)
        error('Minimization constraints parameter must be a cell array\n');
    end

    if ~isempty(maximizationConstraints) && ~iscell(maximizationConstraints)
        error('Maximization constraints parameter must be a cell array\n');
    end

    [outputExpressions,outputNames]=checkOutputExpressions(outputExpressions);

    %% Declare 'sets' for initializing parameters
    if length(parameters)>0
        classhelp{end+1}='% Set parameters';
    end
    for i=1:length(parameters)
        declareSet(code,parameters{i},sprintf('setP_%s',name(parameters{i})));
        msize=size(parameters{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setP_%s(obj,{[%s] matrix});',...
                                name(parameters{i}),index2str(msize));
    end

    %% Declare 'sets' for initializing primal minimization variables
    classhelp{end+1}='Initialize primal minimization variables';
    for i=1:length(minimizationVariables)
        declareSet(code,minimizationVariables{i},...
                   sprintf('setV_%s',name(minimizationVariables{i})));
        msize=size(minimizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                name(minimizationVariables{i}),index2str(msize));
    end

    %% Declare 'sets' for initializing primal maximization variables
    classhelp{end+1}='Initialize primal maximization variables';
    for i=1:length(maximizationVariables)
        declareSet(code,maximizationVariables{i},...
                   sprintf('setV_%s',name(maximizationVariables{i})));
        msize=size(maximizationVariables{i});while(length(msize)<2),msize(end+1)=1;end
        classhelp{end+1}=sprintf('setV_%s(obj,{[%s] matrix});',...
                                name(maximizationVariables{i}),index2str(msize));
    end

    %% Parse the constraints
    
    Tvariable t_ [] % barrier parameter. Declaring here for initialization lambda_x,lambda_y

    
    % minimization constraints


    [G_x,F_x,~,~,outputExpressions]=parseConstraints(code,classname,minimizationConstraints,outputExpressions,'min');
    
    m_x=length(F_x);
    l_x=length(G_x);    
    
    src={};
    dst={};

    if size(F_x,1)>0
        F_x=-F_x; % Deductions in the article were made with constraints such that F<0. parseConstraints returns for F>0.
        lambda_x_=Tvariable('lambda_x_',size(F_x));
        s_x_=Tvariable('s_x_',size(F_x));   
        
        s_x_new=max(-F_x,0.001);
        lambda_x_new=t_*ones(m_x,1)./s_x_new;
        src={s_x_new,lambda_x_new};
        dst={s_x_,lambda_x_};
    else
        lambda_x_=Tzeros(0);
        s_x_=Tzeros(0);
    end

    if l_x>0
        nu_x_=Tvariable('nu_x_',size(G_x));
        
        src={src{:},full(Tzeros(l_x))};
        dst={dst{:},nu_x_};
        
    else
        nu_x_=Tzeros(0);
    end
   

    % maximization constraints

    [G_y,F_y,~,~,outputExpressions]=parseConstraints(code,classname,maximizationConstraints,outputExpressions,'max');

    m_y=length(F_y);
    l_y=length(G_y);    
    
    
    if m_y>0
        F_y=-F_y; % Deductions in the article were made with constraints such that F<0. parseConstraints returns for F>0.
        lambda_y_=Tvariable('lambda_y_',size(F_y));
        s_y_=Tvariable('s_y_',size(F_y));
        
        s_y_new=max(-F_y,0.001);
        lambda_y_new=t_*ones(m_y,1)./s_y_new;
        src={src{:},s_y_new,lambda_y_new};
        dst={dst{:},s_y_,lambda_y_};
    else
        lambda_y_=Tzeros(0);
        s_y_=Tzeros(0);
    end

    if l_y>0
        nu_y_=Tvariable('nu_y',size(G_y));
        
        src={src{:},full(Tzeros(l_y))};
        dst={dst{:},nu_y_};
        
    else
        nu_y_=Tzeros(0);
    end
    
    if l_x+m_x+l_y+m_y>0
        declareCopy(code,dst,src,'initDual__');
    end
    

    %% Pack primal variables
    % minimization
    [x_,whereVariables_x,~,~,objective,outputExpressions,F_x,G_x,F_y,G_y]...
        =packVariables(minimizationVariables,'x_',objective,outputExpressions,F_x,G_x,F_y,G_y);
    x0=packExpressions(minimizationVariables);

    % maximization
    [y_,whereVariables_y,~,~,objective,outputExpressions,F_y,G_y]...
        =packVariables(maximizationVariables,'y_',objective,outputExpressions,F_y,G_y);
    y0=packExpressions(maximizationVariables);

    src={x0,y0};
    dst={x_,y_};
    declareCopy(code,dst,src,'initPrimal__');
    
    n_x=length(x_);
    n_y=length(y_);



    %% Declare 'gets' for output expressions
    classhelp{end+1}='% Get outputs';
    classhelp{end+1}='';
    for i=1:length(outputExpressions)
        classhelp{end}=[classhelp{end},outputNames{i},','];
    end
    classhelp{end}=sprintf('[%s]=getOutputs(obj);',classhelp{end}(1:end-1));
    classhelp{end+1}=sprintf('[y (struct)]=getOutputs(obj);',classhelp{end}(1:end-1));
    
   declareGet(code,cell2struct(outputExpressions,outputNames),'getOutputs');




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Constructing the derivatives
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Defining epsilons to be adjusted
    Tvariable eps_x_ []
    Tvariable eps_y_ []
    % Defining gamma, used to stabilize the LDL decomposition p for primal,
    % d for dual
    Tvariable gamma_px_ []    
    Tvariable gamma_dx_ []
    Tvariable gamma_py_ []
    Tvariable gamma_dy_ []
    
    
    declareSet(code,t_,'setA_t_')
    declareSet(code,eps_x_,'setA_eps_x_')
    declareSet(code,eps_y_,'setA_eps_y_')
    declareSet(code,gamma_px_,'setA_gamma_px_')
    declareSet(code,gamma_dx_,'setA_gamma_dx_')
    declareSet(code,gamma_py_,'setA_gamma_py_')
    declareSet(code,gamma_dy_,'setA_gamma_dy_')
    
    % Constructing lagrangian

    L=objective+nu_x_*G_x+lambda_x_*(F_x+s_x_)+nu_y_*G_y-lambda_y_*(F_y+s_y_);

    % Constructing gradient and extended hessian
    g__=[gradient(L,x_);lambda_x_.*s_x_-t_;gradient(L,y_);-lambda_y_.*s_y_+t_;G_y;-F_y-s_y_;G_x;F_x+s_x_];
    z_={x_,s_x_,y_,s_y_,nu_y_,lambda_y_,nu_x_,lambda_x_}; 
    H__=gradientVector(g__,z_);
    
    % Normalizing with respect to the slack variables
    Sinv=[Tones(n_x);s_x_.\Tones(m_x);Tones(n_y);s_y_.\Tones(m_y);Tones(l_y+m_y+l_x+m_x)];
    
    g_=Sinv.*g__;
%     H=diag(Sinv)*H_;
    H_=tprod(Sinv,1,H__,[1,2]);
    
    % defining block used in verifying conditions
    Hyy_=H_(n_x+m_x+(1:n_y+2*m_y+l_y),n_x+m_x+(1:n_y+2*m_y+l_y));
    
    Hxy_=H_(1:(n_x+m_x),n_x+m_x+(1:n_y+2*m_y+l_y));
    
    Jx_=H_(1:(n_x+m_x),n_x+m_x+n_y+2*m_y+l_y+1:end);

    % perturbation matrices
    Ex=diag([eps_x_*Tones(n_x);Tzeros(m_x)]);   
    Ey=diag([eps_y_*Tones(n_y);Tzeros(l_y+2*m_y)]);   
    E=diag([eps_x_*Tones(n_x);Tzeros(m_x);-eps_y_*Tones(n_y+m_y);Tzeros(l_y+m_y+l_x+m_x)]);
    
    
    Gammay=diag([-gamma_py_*Tones(n_y+m_y);gamma_dy_*Tones(l_y+m_y)]);
    Gamma=diag([gamma_px_*Tones(n_x+m_x);-gamma_py_*Tones(n_y+m_y);gamma_dy_*Tones(l_y+m_y);-gamma_dx_*Tones(l_x+m_x)]);
    
    % ldl decompositions
    ldlHyy=ldl(Hyy_-Ey+Gammay);
    ldlH=ldl(H_+E+Gamma);
    
    Eyinv=diag([(1/eps_y_)*Tones(n_y);(1/gamma_dy_)*Tones(l_y+2*m_y)]); 
    Haux_=[Ex,Hxy_,Jx_;
        Hxy_',Hyy_*Eyinv*Hyy_,Tzeros(n_y+2*m_y+l_y,l_x+m_x);
        Jx_',Tzeros(l_x+m_x,n_y+2*m_y+l_y),Tzeros(l_x+m_x,l_x+m_x)];
    
    ldlHaux=ldl(Haux_+Gamma);
    
    
    % newton direction
    dz_=-(ldlH\g_);
    
    % parsing newton direction for each subvariable
    track_=0;
    wherez_=cell(length(z_),1);
    for i=1:length(z_)
        len=numel(z_{i});
        wherez_{i}=track_+1:track_+len;
        track_=track_+len;
    end
    
    
    dx=dz_(wherez_{1});
    ds_x=dz_(wherez_{2});
    dy=dz_(wherez_{3});
    ds_y=dz_(wherez_{4});
    dnu_y=dz_(wherez_{5});
    dlambda_y=dz_(wherez_{6});
    dnu_x=dz_(wherez_{7});
    dlambda_x=dz_(wherez_{8});        
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Computing steps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    Tvariable tau []
    declareSet(code,tau,'setA_tau')
    alpha_s_x=clp(s_x_,ds_x/tau);
    alpha_s_y=clp(s_y_,ds_y/tau);
    alpha_lambda_x=clp(lambda_x_,dlambda_x/tau); 
    alpha_lambda_y=clp(lambda_y_,dlambda_y/tau);    
    
    declareGet(code,{alpha_s_x,alpha_s_y,alpha_lambda_x,alpha_lambda_y},'get_alphas')
    
    % primal and dual step sizes that will be applied
    Tvariable alpha_p_x_ [] 
    Tvariable alpha_p_y_ []
    Tvariable alpha_d_x_ []
    Tvariable alpha_d_y_ []
    
    declareSet(code,alpha_p_x_,'setA_alpha_p_x_')
    declareSet(code,alpha_p_y_,'setA_alpha_p_y_')
    declareSet(code,alpha_d_x_,'setA_alpha_d_x_')
    declareSet(code,alpha_d_y_,'setA_alpha_d_y_')
    
    
    % Computing new steps
    new_x_=x_+alpha_p_x_*dx;
    new_s_x_=s_x_+alpha_p_x_*ds_x;
    new_lambda_x_=lambda_x_+alpha_d_x_*dlambda_x;
    new_nu_x_=nu_x_+alpha_d_x_*dnu_x;
    
    new_y_=y_+alpha_p_y_*dy;
    new_s_y_=s_y_+alpha_p_y_*ds_y;
    new_lambda_y_=lambda_y_+alpha_d_y_*dlambda_y;
    new_nu_y_=nu_y_+alpha_d_y_*dnu_y;    
    

    dst={x_,s_x_,lambda_x_,nu_x_,y_,s_y_,lambda_y_,nu_y_};
    src={new_x_,new_s_x_,new_lambda_x_,new_nu_x_,new_y_,new_s_y_,new_lambda_y_,new_nu_y_};
    declareCopy(code,dst,src,'applyStepsUpdateVariables__')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% declaring gets used in optimization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
    
    declareGet(code,{Tconstant(n_x),Tconstant(l_x),Tconstant(m_x),Tconstant(n_y),Tconstant(l_y),Tconstant(m_y)},'get_sizes')
    
    declareGet(code,{x_,y_},'get_xy')
    declareGet(code,{G_x,F_x,G_y,F_y},'get_constr')
    
    declareGet(code,{norminf(g__)},'get_norm_g')
    
    % norm of individual sections of g_
    norm_gx=norminf(g__(wherez_{1}));
    norm_gy=norminf(g__(wherez_{3}));
    
    if m_x>0
        norm_gs_x=norminf(g__(wherez_{2}));
        norm_glambda_x=norminf(g__(wherez_{8}));
    else
        norm_gs_x=full(Tzeros(1));
        norm_glambda_x=full(Tzeros(1));
    end

    if l_x>0
        norm_gnu_x=norminf(g__(wherez_{7}));
    else
        norm_gnu_x=full(Tzeros(1));
    end

    if m_y>0
        norm_gs_y=norminf(g__(wherez_{4}));
        norm_glambda_y=norminf(g__(wherez_{6}));
    else
        norm_gs_y=full(Tzeros(1));
        norm_glambda_y=full(Tzeros(1));
    end

    if l_y>0
        norm_gnu_y=norminf(g__(wherez_{5}));        
    else
        norm_gnu_y=full(Tzeros(1));
    end   
    
    declareGet(code,{norm_gx,norm_gy,norm_gnu_x,norm_gnu_y,norm_glambda_x,norm_glambda_y,norm_gs_x,norm_gs_y},'get_norms');
    
     
    
    declareGet(code,ldl_d(ldlH),'get_ldlH')
    declareGet(code,ldl_d(ldlHyy),'get_ldlHyy')
    declareGet(code,ldl_d(ldlHaux),'get_ldlHaux')
    
    % Helpers get for debugging
    declareGet(code,Hyy_,'get_Hyy')
    declareGet(code,H_,'get_H')
    declareGet(code,g_,'get_g')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Temporary, not needed for ip_newton
    declareSet(code,x_,'set_x')
    declareSet(code,y_,'set_y')
    
    Hxx_=H_(1:n_x,1:n_x);
    gx_=g_(1:n_x);
    gy_=g_(n_x+(1:n_y));
    
    declareGet(code,{Hxx_,Hyy_,Hxy_,gx_,gy_},'get_TR_param')    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    if strcmp(code_type,'matlab')
        class2compute('csparseObject',code,'classname',classname);
    elseif strcmp(code_type,'c')
        cmex2compute('csparseObject',code,'classname',classname)
    else
        error('Code type needs to be either matlab or c')
    end

end

