function [dxdt,Diagnostics_vec,Diagnostics_vec_protists,reproduction,detrit_flux,light] =function_model_NPZF_diags(t,x,param,theta_P,theta_cop,theta_cop_P,theta_cop_cop,theta_cop_D,Ls,zmld2,...
    dzdt2,Diff2,temp,seasonal_switch,diagnostics_switch)
N=x(1);
P=x(2:param.nbr_P+1);
C=x(param.nbr_P+2:param.nbr_P+param.nbr_cops*param.nbr_stages+1);
F=x(param.nbr_P+param.nbr_cops*param.nbr_stages+2:end);

diagnostics_switch=1;

if seasonal_switch==0 %constant environment

        Is=Ls(1);      % Light at surface
        T=temp(1);     % Temperature
        Diff=param.D;
        L=Is;
        MLD_P=0;
        MLD_C=0;
        MLD_D= - (param.sink_fp'./40).*F;% - Diff.*F;
        seeding_P=0;
        seeding_C=0;

        if mod(ceil(t),365*5) == 0
             disp(['years',num2str(ceil(t)/365)]);
        end

else % Time-varying environmental forcings
    if t==0
        Is=Ls(1);      % Light at surface
        dzdt=dzdt2(1); % Rate of change of ML
        zmld=zmld2(1); % MLD
        T=temp(1);     % Temperature
        Diff=Diff2(1); % Input of nitrogen
    else
        Is=Ls(ceil(t));
        dzdt=dzdt2(ceil(t));
        zmld=zmld2(ceil(t));
        T=temp(ceil(t));
%         Diff=Diff2(ceil(t));
        Diff=((param.Diff_min + max(0,dzdt))./(zmld));

    end
   
    
    ktot = param.kw + param.kchl .* sum(P);
    L = (Is./(ktot*zmld)).*(1-exp(-ktot.*zmld)); % Mean PAR in surface layer
%     L =(Is.*(exp(-ktot.*(zmld/2)))); % Mean PAR in surface layer
    MLD_P= - ((param.Diff_min + max(0,dzdt))./zmld).*P; %effects of ML dynamics on protists
    MLD_C= - (dzdt/zmld) .* C; %effects of ML dynamics on copepods
    MLD_D= - ((param.Diff_min + max(0,dzdt) + param.sink_fp' )./zmld).*F; %+ param.sink_fp'
    seeding_P =param.flow.*(param.inputP' - P);
    seeding_C =param.flowC.*(param.inputC' - C(param.ind_b));
    
         if mod(ceil(t),365) == 0
                disp(['years',num2str(ceil(t)/365)]);
         end
         
%             if Is>180
%                 Diff=10^(-1.5);
%             end
%          
         
end

%Q10 temperatures
I=param.I.*param.Q_I.^((T-param.Tref)/10);
alpha_N=param.alpha_N.*param.Q_N.^((T-param.Tref)/10);
R=param.R.*param.Q_k.^((T-param.Tref)/10);
k=param.k.*param.Q_k.^((T-param.Tref)/10);
mu_max=param.mu_max.*param.Q_k.^((T-param.Tref)/10);
remin=param.remin.*param.Q_k.^((T-param.Tref)/10);
alpha_c=param.alpha_c;%.*2.^((T-param.Tref)/10);
alpha_F=param.alpha_F;%.*2.^((T-param.Tref)/10);

%     %%%%%%%%%%%%%%%%%
%     Diff=0;
%     seeding_P(:) =0;
%     seeding_C(:) =0;
%     param.sink_fp=0;
%     param.d_c(:)=0;
%     k(:)=0; %este te da problemas
%     param.rep_eff(:)=1;
% %     R(:)=0;
%     theta_cop_D(:)=0;
%     theta_cop(:)=0;
% %     theta_P(:)=0;
%     theta_cop_cop(:)=0;
%     theta_cop_P(:)=0;
%     MLD_D(:)=0;
%     param.m=0;%!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%

% %food encountered
prey_vec_biom_gen=P';
prey_vec_biom=[P',C',F'];

% if L<5
%     theta_cop(param.ind_act,:)=0; %rows are predators, columns are prey   
% end

% %total food for each predator
E_P_tot=(theta_P*prey_vec_biom_gen');
% E_cop_tot=(theta_cop*prey_vec_biom'); 

    E_cop_totP=theta_cop_P*P;
    E_cop_totC=theta_cop_cop*C;
    E_cop_totF=theta_cop_D*F;
E_cop_tot=E_cop_totP+E_cop_totC+E_cop_totF;
%uptake rates by protists
J_N=(mu_max.*param.Qcn.*alpha_N.*N)./(param.Qcn.*alpha_N.*N+mu_max);
J_L=(mu_max.*param.alpha_L.*L)./(param.alpha_L.*L+mu_max);
J_F=(mu_max.*alpha_F.*E_P_tot')./(alpha_F.*E_P_tot'+mu_max);

% J_L(:)=0;%!!!!!!!!!!!!!!!!!!!!!

J_Ntot=sum(J_N'.*P); %total uptake from the N pool

J_R=R; %respiration

%total flow
Jtot=min(J_L + param.eff(1).*J_F - J_R , J_N+param.eff(1)*J_F);
mu=Jtot;
fg=(alpha_F.*E_P_tot')./(alpha_F.*E_P_tot'+mu_max);
leaks=max(0,J_N-J_L + J_R);


%mortality from protists on generalists
% pred_P=sum(theta_P(:,1:param.nbr_P).*param.alpha_F'.*(1-fg').*P,1);
pred_P=theta_P(:,1:param.nbr_P)'*(alpha_F'.*(1-fg').*P);

%Feeding level of copepods
F_lvl=(alpha_c.*E_cop_tot)./(alpha_c.*E_cop_tot+I);
% F_lvl_P=(alpha_c.*E_cop_totP)./(alpha_c.*E_cop_totP+I);
% F_lvl_C=(alpha_c.*E_cop_totC)./(alpha_c.*E_cop_totC+I);
% F_lvl_F=(alpha_c.*E_cop_totF)./(alpha_c.*E_cop_totF+I);
% F_lvl=F_lvl_P+F_lvl_C+F_lvl_F;

%     P_frac=E_cop_totP./E_cop_tot;
%     C_frac=E_cop_totC./E_cop_tot;
%     D_frac=E_cop_totF./E_cop_tot;
% 
% F_lvl_P=P_frac.*F_lvl;
% F_lvl_C=C_frac.*F_lvl;
% F_lvl_F=D_frac.*F_lvl;


%mortality from copepods on protists
% pred_C_on_gen=sum(param.migration.*theta_cop_P.*param.alpha_c.*(1-F_lvl).*C,1);
pred_C_on_gen=theta_cop_P'*(param.migration.*alpha_c.*(1-F_lvl).*C);
%  pred_C_on_gen=sum(theta_cop_P./E_cop_tot.*I.*F_lvl.*C);
% pred_C_on_gen2=theta_cop_P'*(param.migration.*alpha_c.*(1-F_lvl_P).*C);
% pred_C_on_gen3=theta_cop_P'*(I.*F_lvl_P.*C./E_cop_totP);

mortP=param.m'.*P.^2 + 0.0001.*param.V'.^(-0.25).*P;
% Protists ODE
dGdt=mu'.*P - mortP - pred_P.*P - pred_C_on_gen.*P + seeding_P + MLD_P;% - param.ml'.*P; %!!!!!!!!!!!!!!!!!!
% dGdt=mu'.*P - param.m'.*P - pred_P'.*P - pred_C_on_gen'.*P + param.flow.*(param.inputP' - P);

if any(pred_P<-1e-2)
    disp('negative predations P')
end

if any(pred_C_on_gen<-1e-2)
    disp('negative predations C')
end

if any(P<-1e-2)
    disp('negative P!')
end

% if t>100
% 
% end


%Copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dCdt=zeros(size(C)); %vector to fill

% Cddm=zeros(1,param.nbr_ddm);
for i=1:param.nbr_ddm
    range_ddm=param.vec_idx_ddm(param.idx_ddm_range1(i):param.idx_ddm_range2(i));
    dCdt(range_ddm)=(sum(C(range_ddm)));
end

dd_mort_C=param.mmax.*dCdt;%.*Diff;%.*dCdt; %uncomment if we want it density dependent
dCdt= - dd_mort_C.* C;
% dCdt(:)=0;%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Net energy gain for all stages
fc=k./(param.eff.*I); 
nu=param.migration.*param.eff.*I.*(F_lvl-fc); 
if seasonal_switch==1 %dormancy
    dorm_idx=find(nu<0);
    nu(dorm_idx)=param.eff(dorm_idx).*I(dorm_idx).*(F_lvl(dorm_idx)-fc(dorm_idx)./3); %switch to dormancy
%     if L<5
%          nu(param.ind_act)=param.eff(param.ind_act).*I(param.ind_act).*(0-param.fc(param.ind_act)./4);
%     end
end
nu_pos=max(0,nu); %positive growth
nu_neg=-min(0,nu); %starvation mortalityt, we put the minus to make it positive
nu_ca_neg=min(0,nu); %starvation of the adults

%Predation from copepods on copepods
% pred_C_on_C=sum(param.migration.*theta_cop_cop.*param.alpha_c.*(1-F_lvl).*C,1)';
pred_C_on_C=theta_cop_cop'*(param.migration.*alpha_c.*(1-F_lvl).*C);

%Predation from copepods on detritus
pred_C_on_D=theta_cop_D'*(param.migration.*alpha_c.*(1-F_lvl).*C);

%maturation rates
d_c=param.d_c+pred_C_on_C+nu_neg+dd_mort_C;
gamma=zeros(size(C)); %for convinience in later calculations we make it in this form
gamma(param.ind_j)=(nu_pos(param.ind_j)-d_c(param.ind_j))./(1-param.z(param.ind_j).^(1-(d_c(param.ind_j)./nu_pos(param.ind_j))));
% gamma(:)=0;%!!!!!!!!!!!!!!!
%ODE copepod
% Size at birth
dCdt(param.ind_b)=dCdt(param.ind_b) + param.rep_eff(param.ind_a).* nu_pos(param.ind_a).*C(param.ind_a) + nu(param.ind_b).*C(param.ind_b) ...
                 - gamma(param.ind_b).*C(param.ind_b) - pred_C_on_C(param.ind_b).*C(param.ind_b) ...
                 - param.d_c(param.ind_b).*C(param.ind_b) + seeding_C;
              
% In-between size-classes
dCdt(param.ind_rest)=dCdt(param.ind_rest) +  gamma(param.ind_rest-1).*C(param.ind_rest-1) + nu(param.ind_rest).*C(param.ind_rest) -...
                      gamma(param.ind_rest).*C(param.ind_rest) - pred_C_on_C(param.ind_rest).*C(param.ind_rest)...
                      -param.d_c(param.ind_rest).*C(param.ind_rest);

% Adults                  
dCdt(param.ind_a) = dCdt(param.ind_a) + gamma(param.ind_a-1).*C(param.ind_a-1) + nu_ca_neg(param.ind_a).*C(param.ind_a) ...
                    - param.d_c(param.ind_a).*C(param.ind_a) - pred_C_on_C(param.ind_a).*C(param.ind_a);
                
dCdt= dCdt + MLD_C; %effects from ML dynamics   


% Fecal pellet production                
fpp=(1-param.eff).*I.*F_lvl;

% ODE detritus
FPP=fpp.*C;

%Copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dFdt=zeros(size(F)); %vector to fill

% Cddm=zeros(1,param.nbr_ddm);
for i=1:param.nbr_fp
    range_fp=param.vec_idx_fp(param.idx_fp_range1(i):param.idx_fp_range2(i));
    dFdt(i)=(sum(FPP(range_fp)));
end

dFdt = dFdt - remin*F  - pred_C_on_D.*F + MLD_D;

stuff_dying=sum(param.d_c.*C) + sum(mortP);  

%ode Nutrients %0.1*!!!!!!!!!!!!!!!!!!!!
dNdt=Diff.*(param.No-N) + (-J_Ntot + sum(leaks'.*P) + sum(remin*F) + 0.1*(stuff_dying)...
        + sum((1-param.rep_eff(param.ind_a)).* nu_pos(param.ind_a).*C(param.ind_a))+ sum((1-param.eff(1)).*J_F'.*P))./param.Qcn; %

if seasonal_switch==0    
    dNdt=dNdt + (sum(k.*C))./param.Qcn; 
elseif seasonal_switch==1
    dNdt=dNdt+ (sum(k(dorm_idx)/3.*C(dorm_idx)))./param.Qcn; 
end    
    
    
    % if L<5
%     dCdt(param.ind_act)=0;   
% end

%%%%%%%
%checks
% sum(pred_C_on_gen.*P)-sum(E_cop_totP./E_cop_tot.*I.*F_lvl.*C)
% sum(pred_C_on_C.*C)-sum(E_cop_totC./E_cop_tot.*I.*F_lvl.*C)
% sum(pred_C_on_D.*F)-sum(E_cop_totF./E_cop_tot.*I.*F_lvl.*C)

if diagnostics_switch==0
    
    dxdt=[dNdt; dGdt; dCdt; dFdt];
    Diagnostics_vec=0;
    Diagnostics_vec_protists=0;
    reproduction=0;
    detrit_flux=0;
    
elseif diagnostics_switch==1
    
    % Diagnostics--------------------------------------------------
    dxdt=0;
    P_frac=theta_cop_P*P;%./E_cop_tot;
    C_frac=theta_cop_cop*C;%./E_cop_tot;
    D_frac=theta_cop_D*F;%./E_cop_tot;
    
    P_frac=((alpha_c.*P_frac)./(alpha_c.*P_frac+I))./F_lvl;
    C_frac=((alpha_c.*C_frac)./(alpha_c.*C_frac+I))./F_lvl;
    D_frac=((alpha_c.*D_frac)./(alpha_c.*D_frac+I))./F_lvl;
    
    m_prot=mortP./P; %specific background mortality
    reproduction=param.rep_eff(param.ind_a).* nu_pos(param.ind_a);%.*C(param.ind_a); 
    NPP=min(J_L - J_R , J_N)'.*P;
    Diagnostics_vec=[F_lvl, d_c, dCdt, d_c-nu_neg,pred_C_on_C, nu, dd_mort_C, nu_neg, P_frac, C_frac, D_frac];
    Diagnostics_vec_protists=[fg', pred_P,mu',pred_C_on_gen, dGdt,NPP,m_prot,J_N',J_L',J_F'];
    
    if seasonal_switch==1
       detrit_flux=(param.Diff_min + max(0,dzdt)+ param.sink_fp').*F; %total flux of carbon out of the ML
    elseif seasonal_switch==0
       detrit_flux=param.sink_fp'.*F;
    end
    
%     light=L;
% mass balance check
dNtotdt=dNdt + (sum(dGdt) + sum(dCdt) + sum(dFdt))./param.Qcn;
dNindt=Diff.*param.No;
dNoutdt=Diff.*N + (1/param.Qcn)*((1-0.1)*(stuff_dying) + sum(dd_mort_C.* C) ...
        + sum((param.sink_fp'./40).*F));%+ sum(k.*C) + sum((1-param.rep_eff(param.ind_a)).* nu_pos(param.ind_a).*C(param.ind_a)));
 light=dNtotdt-dNindt+dNoutdt;

   %sum(knormal.*C)+ sum((1-param.rep_eff(param.ind_a)).* nu_pos(param.ind_a).*C(param.ind_a))
end

end