function [ T ] = T_exact( cc, Theta_in, Theta_out, R_in, R_out, alpha )
    
    R = (cc(:,1).^2+cc(:,2).^2).^(1/2);
    
    c1 = (Theta_in - Theta_out) / ( ...
         ( 1/R_out + alpha*log(R_out) ) + ( 1/R_in - alpha*log(R_in) ) ...
         );
    
    c2 = (Theta_in - Theta_out) * (1/R_out + alpha*log(R_out))/ (...
         ( 1/R_out + alpha*log(R_out) ) + ( 1/R_in - alpha*log(R_in) )...
         ) + Theta_out;
     
    T = -c1 * log(R) + c2;
    
end

