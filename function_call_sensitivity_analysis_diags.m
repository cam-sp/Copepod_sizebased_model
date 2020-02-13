function [Ncell,Pcell,Ccell,Dcell,NPPcell,fluxcell]=function_call_sensitivity_analysis_diags(seasonal_switch,linear_mortality_switch,dd_mortality_switch,ontogeny_switch,...
    sensit_switch,sensit_param,diagnostics_switch,parallerun_mode)

months = ['J';'F';'M';'A';'M';'J';'J';'A';'S';'O';'N';'D'];

param=parameters_copepod_model(linear_mortality_switch,dd_mortality_switch,ontogeny_switch,sensit_switch,sensit_param);

yearsrun=150;
param.No=140;
param.lat=55;
param.D=10^(-1.5);

tic

[temp,Diff,zmld,Ls,dzdt]=func_seasonal_forcing(param, yearsrun, seasonal_switch);

%this is a fix
    param.WD=param.WF;
    param.nbr_D=param.nbr_fp;
    param.D_dw=param.F_dw;
    param.D_up=param.F_up;


[theta_P_P,theta_P_D, theta_P, theta_cop_P, theta_cop_cop, theta_cop_D, theta_cop]=func_feeding_kernels(param);

Ncell=cell(1,length(sensit_param));
Dcell=cell(1,length(sensit_param));
Pcell=cell(1,length(sensit_param));
Ccell=cell(1,length(sensit_param));
NPPcell=cell(1,length(sensit_param));
fluxcell=cell(1,length(sensit_param));

%Initial conditions
No=5;
Po=zeros(1,param.nbr_P);
Po(:)=1;
Co=zeros(1,length(param.Wvec));
Co(:)=1;
Fo=zeros(1,param.nbr_fp);
Fo(:)=0;

    options1 = odeset('Refine',1);
    options2 = odeset(options1,'NonNegative',1:1+param.nbr_P+param.nbr_cops*param.nbr_stages+param.nbr_fp);
    
for i=1:length(sensit_param) %here you can switch to parallel running (just change the "for" for "parfor")

Diff=sensit_param(i);


    [t,x]= ode23(@(t,x) function_model_NPZF(t,x,param,theta_P_P,theta_cop,theta_cop_P,theta_cop_cop,...
        theta_cop_D,Ls,zmld,dzdt,Diff,temp,seasonal_switch,diagnostics_switch,sensit_switch),0:1:365*yearsrun, [No, Po, Co, Fo],options2); 
    
    N=x(:,1);
    P=x(:,2:param.nbr_P+1);
    C=x(:,param.nbr_P+2:param.nbr_P+param.nbr_cops*param.nbr_stages+1);
    D=x(:,param.nbr_P+param.nbr_cops*param.nbr_stages+2:end);
    
% end
Ncell{i}=N(ceil(t(end)*3/4):end,:);
Pcell{i}=P(ceil(t(end)*3/4):end,:);
Ccell{i}=C(ceil(t(end)*3/4):end,:);
Dcell{i}=D(ceil(t(end)*3/4):end,:);


NPP_diag=zeros(length(t),length(param.V));
detrit_flux_diag=zeros(length(t),param.nbr_fp);

        for jij=ceil(t(end)*3/4):length(t)

            [dxdt, output_diagnostic, Diagnostics_vec_protists, reproduction,detrit_flux,light]= function_model_NPZF_diags(t(jij),x(jij,:)',param,theta_P_P,...
                        theta_cop,theta_cop_P,theta_cop_cop, theta_cop_D,Ls,zmld,dzdt,Diff,temp,seasonal_switch,diagnostics_switch); 


            NPP_diag(jij,:)=Diagnostics_vec_protists(:,6);
            detrit_flux_diag(jij,:)=detrit_flux';

        end

NPPcell{i}=NPP_diag;
fluxcell{i}=detrit_flux_diag;

end
toc


end