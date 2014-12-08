function [ sx, sy ] = s_exact( cc, Theta_in, Theta_out, R_in, R_out, alpha )

    cx = cc(:,1);
    cy = cc(:,2);
    
    R = (cx.^2+cy.^2).^(1/2);
    
    nx = 1./R .* cx;
    ny = 1./R .* cy;
    
    c1 = (Theta_in - Theta_out) / ( ...
         ( 1/R_out + alpha*log(R_out) ) + ( 1/R_in - alpha*log(R_in) ) ...
         );
     
    sx = c1 * 1./R .* nx;
    sy = c1 * 1./R .* ny;
    
end

