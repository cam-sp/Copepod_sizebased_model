function [fncs_analytical]=z_function_analytical_solutions()

        fncs_analytical.community_spectrum       = @function_analytical_spectrum;
        fncs_analytical.predation_mortality      = @function_analytical_mortality;

end

function Nc=function_analytical_spectrum(f0,beta,sigma,m,h,v,q,n)

lambda=2;

alpha=sqrt(2*pi)*beta^(lambda-2)*exp((lambda-2)^2*sigma^2/2);

Nc=h*f0/(alpha*v*(1-f0))*m.^(-2-q+n);

end

function mu=function_analytical_mortality(f0,beta,sigma,m,h,n)

lambda=2;

alpha=sqrt(2*pi)*beta^(lambda-2)*exp((lambda-2)^2*sigma^2/2);

mu=f0*h*alpha^(-1)*m.^n;

end