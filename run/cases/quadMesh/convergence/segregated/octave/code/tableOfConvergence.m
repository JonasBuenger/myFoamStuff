function [  ] = tableOfConvergence(L2_errors, num_cells)

str = sprintf(['Mesh level' '\t' 'num cells' '\t' 's_err' '\t\t' 'theta_err' '\t' 'ORDER s' '\t\t' 'ORDER Theta' ]);
disp(str);

for i=1:6
    level = i;
    nc = num_cells(i);
    s_err = L2_errors(i,1);
    T_err = L2_errors(i,3);

    s_order = 0;
    T_order = 0;
    if (i>1)
        s_order = log2( L2_errors(i-1,1) / L2_errors(i,1) );
        T_order = log2( L2_errors(i-1,3) / L2_errors(i,3) );
    end

    str = sprintf([ int2str(level) '\t\t'...
                    int2str(nc) '\t\t' ...
                    num2str(s_err,'%e') '\t'...
                    num2str(T_err,'%e') '\t'...
                    num2str(s_order) '\t\t'...
                    num2str(T_order)]);
    disp(str);
end

end

