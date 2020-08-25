function [theta_P_P, theta_cop_P, theta_cop_cop, theta_cop_F, theta_cop]=...
    z_function_feeding_kernels(param)
% create the matrixes for the preference function
% in all: rows are predators, columns are prey
% OUTPUTS:
    % theta_P_P     --> preference of protists for protists
    % theta_cop_P   --> preference of copepods for protists
    % theta_cop_cop --> preference of copepods for copepods
    % theta_cop_D   --> preference of copepods for detritus
    % theta_cop     --> total preference funtion of copepods (for protists, copepods, and detritus together)

%Note that due to the adult stages of copepods being discrete, we solve the
%prefernce matrix of copepods for copepods in two steps: one for all the
%juveniles where we use the error function and one for the adults where we
%use the normal function.
%------------------------------------------------------------------------------------------------------

%%%%%% Preference of protists for protists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_P_P=zeros(param.nbr_P,param.nbr_P);

for j=1:param.nbr_P
%error function    
theta_P_P(j,:) = sqrt(pi/2)*param.sigma_P* (...
  erf( (log(param.V_dw)-log(param.V(j)/param.beta_P))/(sqrt(2)*param.sigma_P) ) ...
  - erf( (log(param.V_up)-log(param.V(j)/param.beta_P))/(sqrt(2)*param.sigma_P) ));

end
theta_P_P =theta_P_P./(log(param.V_dw(:))-log(param.V_up(:)));




%%%%%% Preference of copepods for protists %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%it starts with generalists, then species 1 all stages, species 2 all stages etc
%stages are ordered from smallest(top) to adults(bottom)                                                          
pred_vec=param.Wvec;
theta_cop_P=zeros(length(pred_vec),param.nbr_P);
for j=1:length(param.Wvec)
%     if param.C_sp_pass>0
        if j<param.ind_pass(1)
            sigma_c=param.sigma_act;
            beta_c=param.beta_act;
        else
                sigma_c=param.sigma_pass;
                beta_c=param.beta_pass;  
        end
%     else
%         sigma_c=param.sigma_act;
%         beta_c=param.beta_act;
%     end
    
theta_cop_P(j,:) = sqrt(pi/2)*sigma_c* (...
  erf( (log(param.V_dw)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ) ...
  - erf( (log(param.V_up)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ));

end
theta_cop_P =theta_cop_P./(log(param.V_dw(:)')-log(param.V_up(:)'));




%%%%%% Preference of copepods for copepods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_cop_cop=zeros(length(param.W(:)),length(param.W(:)));
for j=1:length(param.W(:))
    if j<param.ind_pass(1)
        sigma_c=param.sigma_act;
        beta_c=param.beta_act;
    else
                sigma_c=param.sigma_pass;
                beta_c=param.beta_pass;              
    end
theta_cop_cop(j,:) = sqrt(pi/2)*sigma_c* (...
  erf( (log(param.C_dw(:))-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ) ...
  - erf( (log(param.C_up(:))-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ));


end
%finish error function for all juveniles on juveniles
theta_cop_cop(param.ind_j,param.ind_j) = theta_cop_cop(param.ind_j,param.ind_j)./(log(param.C_dw(param.ind_j))-log(param.C_up(param.ind_j)));
%finish error function for adults on juveniles
theta_cop_cop(param.ind_a,param.ind_j) = theta_cop_cop(param.ind_a,param.ind_j)./(log(param.C_dw(param.ind_j))-log(param.C_up(param.ind_j)));

%Now we calculate the normal function for the adult preys
for i=1:length(param.W(:))
                if i<param.ind_pass(1)
                    sigma_c=param.sigma_act;
                    beta_c=param.beta_act;
                else
                    sigma_c=param.sigma_pass;
                    beta_c=param.beta_pass;       
                end


theta_cop_cop(i,param.ind_a)=exp(-((log(beta_c*param.W(param.ind_a)./param.W(i))).^2)./(2*sigma_c^2));

end

% lower preference for passive feeding copepods
theta_cop_cop(:,param.ind_pass)=param.p'.*theta_cop_cop(:,param.ind_pass); 




%%%%%% Preference of copepods for fecal pellets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_cop_F=zeros(length(param.W(:)),param.nbr_fp);

for j=1:length(param.W(:))

    if j<param.ind_pass(1)
        sigma_c=param.sigma_act;
        beta_c=param.beta_act;
    else
        sigma_c=param.sigma_pass;
        beta_c=param.beta_pass;              
    end

theta_cop_F(j,:) = sqrt(pi/2)*sigma_c* (...
  erf( (log(param.F_dw)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ) ...
  - erf( (log(param.F_up)-log(pred_vec(j)/beta_c))/(sqrt(2)*sigma_c) ));

end
theta_cop_F=theta_cop_F./(log(param.F_dw)-log(param.F_up));


theta_cop=[theta_cop_P, theta_cop_cop, theta_cop_F];
end