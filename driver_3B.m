%%%% driver for bouncing drop with bessel decomposition model in bath and
%%%% spherical decomposition model for droplet
%%%code verification.
mult=1; % muktiplier (for sweeps)
U_in = linspace(0.05,1.0,40); %impact velocity
%%%%%%%input parameters
tf = 0.02; % simulation end time
final_data = zeros(length(U_in),7); %pre-allocation
for j = 1:length(mult)
% %%% a couple physical parameterss
rho=  998.2; % %in kg/m^3 %%%%forwater bath
sig = 72.2e-3;%21e-3; % bath surface tension % N/m
mu = 0.9793e-3*mult(j); % bath viscosity in Pa-s
rho_s = rho; %drop density
sig_s = sig; % drop surface t
mu_s=mu; %drop viscosity
nu = mu/rho; %kinematic viscosity
nu_s = mu_s/rho_s;
g = 9.81;%gravity
%%%%%%%%%%%%
R = 0.035e-2; %droplet radius
phys(1) = R; % compile input parameters into array
phys(2) = sig;
phys(3) = rho;
phys(4) = nu;
phys(5) = sig_s;
phys(6) = rho_s;
phys(7) = nu_s;
%simulation parameters
b = 25*R; %size of container
l = 0; %bessel function order
% get eigenvalues
u = besselzero(0,1000,1)'; % code from: Jason Nicholson (2022). Bessel Zero Solver (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver)
k = besselzero(-1,1000,1)'; % code from: Jason Nicholson (2022). Bessel Zero Solver (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver)
u = u/b;
k = k/b;
NM = [151]*ones(length(U_in),1); % number of bath modes
MM=[55]*ones(length(U_in),1); % number of drop modes
%%%%%pressure function choice; options are parabola, poly4,poly6
% ptype = 'parabola';
ptype = 'poly6';
% ptype = 'poly4';
%%%%%%%%%start loop over U_in
for i=1:length(U_in)
    printstr = ['Running case for Radius R = ' num2str(R) ' with velocity U_in = ' num2str(U_in(i)) '_mult' num2str(mult(j)) ];
    sprintf('%s',printstr)
    U = -1*U_in(i); %make sign correct
    nm = NM(i);
    mm = MM(i);
    [t4a,Ya,out4,rca] = bounce_alventosa_bessel_implicit(U,tf,u,k,nm,b,l,phys,mm,ptype,g); %make solver function call
    F4 = Ya(:,end);
    %%%%%%%%%%%%%%%%%%%%%%%%%plotting after running
    figure(8) % plot trajectory
    hold on
    plot(t4a*sqrt(sig/(rho*R^3)),(Ya(1:end,1)./R),'b','LineWidth',1.5)
    xlabel('$t/t_{\sigma}$','Interpreter','latex','FontSize',18);ylabel('$z_{cm}/R$','Interpreter','latex','FontSize',18);
    figure(10) %plot velocity
    hold on
    plot(t4a*sqrt(sig/(rho*R^3)),(Ya(1:end,2)./abs(U)),'r','LineWidth',1.5)
    xlabel('$t/t_{sigma}$','Interpreter','latex','FontSize',18);ylabel('$z_{cm}/R$','Interpreter','latex','FontSize',18);
    %%%%%%%%fill output array for plotting impact parameters
    final_data(i,1) = out4(1); %c_r measured when force turns oof (exact loss of contact)
    final_data(i,2) = out4(2); %tc measured when force turns oof (exact loss of contact)
    final_data(i,3) = out4(3); %max_p 
    final_data(i,4) = out4(4); %We
    final_data(i,5) = out4(5); %Bo
    final_data(i,6) = out4(6); %Oh
    final_data(i,7) = out4(7); %c_r method 2 measure when z=R *this is what is used in paper due to experimental limitations*
    final_data(i,8) = out4(8); %tc method 2 measure when z=R *this is what is used in paper due to experimental limitations*
    final_data(i,9) = min(sum(Ya(:,2*mm+1:nm+2*mm),2)); %penetration dpeth
    fnm = ['./sweep_Vi' num2str(U_in(i)*100) 'M_' num2str(nm) '_L' num2str(mm) '.mat']; %input filename to save output
    save(fnm,'t4a','Ya','final_data','rca','k','b','nm','mm'); %save output
%     plot_drop_waves_3(Ya,k,b,nm,mm,R,t4a,1) %make videos! toggle comment on
%    / off
end
%%%%%%%%%%%%more plots
figure(1)
subplot(1,2,1)
cc = colormap(winter(length(mult)));
hold on
idx = final_data(:,7)==0;
if ~isempty(idx)
    plot(final_data(~idx,4),final_data(~idx,7),'Color',cc(j,:),'LineStyle', '--', 'LineWidth',2)
else
    plot(final_data(:,4),final_data(:,7),'Color',cc(j,:),'LineStyle', '--', 'LineWidth',2)
end
xlabel('$We$','Interpreter','latex','FontSize',18);ylabel('$\alpha$','Interpreter','latex','FontSize',18); %title('Coefficient of Restitution vs We')

subplot(1,2,2)
idx2 = final_data(:,8)==0;
hold on
if ~isempty(idx2)
    plot(final_data(~idx2,4),final_data(~idx2,8)/((rho*R^3)/sig)^(1/2),'Color',cc(j,:),'LineStyle', '--', 'LineWidth',2)
else
    plot(final_data(:,4),final_data(:,8)/((rho*R^3)/sig)^(1/2),'Color',cc(j,:),'LineStyle', '--', 'LineWidth',2)
end
xlabel('$We$','Interpreter','latex','FontSize',18);ylabel('$t_c/t_{\sigma}$','Interpreter','latex','FontSize',18);%title('Dimensionless Contact Time vs We')
end
