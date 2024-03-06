%3 dims
[X,Y,Z]=meshgrid(a_list,k_list,cy_list);
figure(4);
clf;
% slice(X, Y, Z, synctime, k_list, 0.5, cy_list); 
% slice(X, Y, Z, synctime, 0.5,  k_list, cy_list);

% rate-coupling strength
slice(X, Y, Z, error_ave*sqrt(2/3), [],  [], 0.1); 



colorbar;
%同期領域の色を指定
caxis([0 1e-6])
shading interp
shading flat
xlim([0,a_max]);
ylim([0,k_max]);
zlim([0 cy_max]);
% xlim([0.5,0.6]);
% ylim([1,1.1]);
% zlim([0 0.1]);
xlabel('x:on-off rate $$\alpha$$','Interpreter','Latex','FontSize',15)
ylabel('y:Coupling strength  $$\sigma$$','Interpreter','Latex','FontSize',15)
zlabel('z:cycle $$\tau$$','Interpreter','Latex','FontSize',15)
view([90 -90])
% legend('2sys')
% legend('3sys(mixed) error12')