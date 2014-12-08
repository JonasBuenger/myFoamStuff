function [  ] = plotErrors( L2_errors, localErrors, data, num_cells)

nc = num_cells(:,1);
n = size(nc);

h = sqrt( min(nc) ./ nc(:) );

fstOrder = zeros(n,1);
secOrder = zeros(n,1);
err_init_s = L2_errors(1,1);
err_init_theta = L2_errors(1,3);
err_init_ave = (err_init_s + err_init_theta) * 0.5;
for i=1:n
    fstOrder(i,1) = err_init_ave * h(i)/h(1);
    secOrder(i,1) = err_init_ave * (h(i)/h(1))^2;
end

hold off;
grid on;
loglog(h, L2_errors(:,1),'+-r','LineWidth',5); hold on;
loglog(h, L2_errors(:,2),'o-.g','LineWidth',5);
loglog(h, L2_errors(:,3),'+-b','LineWidth',5);
loglog(h, fstOrder(:,1), ':k','LineWidth',3);
loglog(h, secOrder(:,1), '--k','LineWidth',3);
title('Error L2-norm','FontSize',22);
hl = legend({'s_x', 's_y','theta','1^{st} order','2^{nd} order'},'location','east');
set(hl,'FontSize',18);
chl = get(hl,'child');
set(chl,'LineWidth',5);
axis([min(h) max(h) 0.0001 0.3])
%set(findall(gcf,'-property','FontSize'),'FontSize',16);
set(gca,'XTick',h(end:-1:1,1));
set(gca,'YTick',[.0001 .001 .01 .1 1.]);
set(gca,'FontSize',18);
hx = xlabel('h_{max}');
hy = ylabel('error');%ylabel('e_s');
set(hx,'FontSize',22);
set(hy,'FontSize',22);

end

