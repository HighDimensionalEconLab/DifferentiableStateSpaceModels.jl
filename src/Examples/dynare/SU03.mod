%=========================================================================%
% Schmitt-Grohe and Uribe (2003, JIE)
% Ambrogio Cesa-Bianchi, July 2012
%=========================================================================%
% If you find any bugs when using this file or want to give me comments 
% and suggestions you can email me at ambrogio.cesabianchi@gmail.com


var  d, c, h, y, i, k, a, lambda,  tb, ca, riskpremium, r ;  

varexo e;                                    
                                             
parameters  gamma, omega, rho, sigmae, delta, psi, alpha, phi, beta, r_w, d_bar;
		alpha  = 0.32;
		rho    = 0.42;
		phi    = 0.028;
		r_w    = 0.04;		
        gamma  = 2;
		omega  = 1.455;
		psi    = 0.000742;
		delta  = 0.1;
		sigmae = 0.0129;
		beta   = 1/(1+r_w);
		h_ss   = ((1-alpha)*(alpha/(r_w+delta))^(alpha/(1-alpha)))^(1/(omega-1)); 
		k_ss   = h_ss/(((r_w+delta)/alpha)^(1/(1-alpha)));
		i_ss   = delta*k_ss;                                                     
		y_ss   = (k_ss^alpha)*(h_ss^(1-alpha));                                   
		d_bar  = 0.7442;
		d_ss   = d_bar;                                                        
		c_ss   = y_ss-i_ss-r_w*d_ss;
		tb_ss  = y_ss-c_ss-i_ss;

model;
    d = (1+exp(r(-1)))*d(-1)- exp(y)+exp(c)+exp(i)+(phi/2)*(exp(k)-exp(k(-1)))^2;
    exp(y) = exp(a)*(exp(k(-1))^alpha)*(exp(h)^(1-alpha));
    exp(k) = exp(i)+(1-delta)*exp(k(-1)); 
    exp(lambda)= beta*(1+exp(r))*exp(lambda(+1)); 
    (exp(c)-((exp(h)^omega)/omega))^(-gamma)   = exp(lambda);  
    ((exp(c)-((exp(h)^omega)/omega))^(-gamma))*(exp(h)^omega)  = exp(lambda)*(1-alpha)*exp(y); 
    exp(lambda)*(1+phi*(exp(k)-exp(k(-1)))) = beta*exp(lambda(+1))*(alpha*exp(y(+1))/exp(k)+1-delta+phi*(exp(i(+1))-delta*exp(k))); 
    a = rho*a(-1)+e; 
    tb = 1-((exp(c)+exp(i))/exp(y));
    ca = (1/exp(y))*(d-d(-1));                                   
    riskpremium = psi*(exp(d-d_bar)-1);
    exp(r) = r_w+riskpremium;
end;


initval;
    r     = log((1-beta)/beta);
    d     = d_ss;
    h     = log(h_ss);
    k     = log(k_ss);
    y     = log(y_ss);
    c     = log(c_ss);
    i     = log(i_ss);
    tb    = 1-((exp(c)+exp(i))/exp(y));
    lambda= log((exp(c)-((exp(h)^omega)/omega))^(-gamma));
end;


check;
steady; 


shocks;
    var e; stderr sigmae;
end;


stoch_simul(order=2, irf=24, nograph);