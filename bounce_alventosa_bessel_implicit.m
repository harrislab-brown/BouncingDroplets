function [t,Y,out,rc] = bounce_alventosa_bessel_implicit(U,tf,u,k,nm,b,l,phys,mm,ptype,g)
% input
% U = impact velo, tf = end time, u = eigenvalues of J0, k =
% eigenvalues of J0', nm = number of modes, b = domain size, % l = order of bessel function, 
% phys = array of physical params phys(1) = R;phys(2) = sig;phys(3) =
% rho;phys(4) = nu; etc
% mm = number of modes in drop, ptype = pressure function, g=gravity

% output, 
% t = time, Y = solution array, contains center of mass trajectory,
% velocity, bath mode amplitudes and velocities, drop mode amplitudes and
% velociteis,a nd force all as functions of time
% out = array of output parameters (index 1 = c_r, 2=tc,3 = max p,
% 4=We,5=Bo,6=c_r different method, 7 = tc different method, rc =contact
% radius 

%%%%%% drop bounce solution with bessel functions
%%%built in physical parameters
nu_a = 1.516e-5; %air viscosity (Not really necessary)
rho_a = 1.204; %air density
% g = 20; %9.81; %m/s^2
H = 0.1; % depth of bath in m
%%%%simulation physical parameters
R = phys(1);
sig = phys(2);
rho = phys(3);
nu  = phys(4);
sig_s = phys(5);
rho_s = phys(6);
nu_s  = phys(7);
%%%%%dimensionless numbers and coefficients
Re = abs(U)*R/nu_a; % reynolds number 
We = (rho_s*U^2*R)/sig_s;
Bo = (rho_s*R^2*g)/sig_s;
Oh = nu_s*sqrt(rho_s/(sig_s*R));
t_cscale = 1/sqrt(sig_s/(rho_s*R^3));
Nw = 8;
w_d = ((Nw*sig_s)/(rho_s*R^3))^(1/2); % natural freq for drop
m = rho_s*(4/3*pi*R^3); % in kg drop mass
c_d = 10*m*nu_s/(R^2); % viscous damping coefficient
%%%%%%%%coefficients for mode amplitudes droplet
ll = 0:mm;
wl2 = ll.*(ll-1).*(ll+2)*(sig_s)/(rho_s*R^3);
alp = 2*(ll-1).*(2.*ll+1)*(nu_s)/(R^2);
Nl = sqrt((2.*ll+1)/(4*pi));
theta = linspace(0,pi,4000); 
sph = nan(mm+1,length(theta));
for i =1:mm+1
    PLM=legendre(ll(i),cos(theta));
    sph(i,:) = Nl(i).*PLM(1,:);
end
% integration of assumed pressure distribution to get d's and c's
r = linspace(0,b,75/R); %domain of integration
%%%%%number of modes wanted
u=u(1:nm);
k=k(1:nm);
%%%%%%some setting up to do
deltat = 3.75e-5;
dtt = tf/deltat;
t = linspace(0,tf,dtt); % pre allocate
a = zeros(length(t)-1,length(u));
alpha = zeros(length(t)-1,length(u));
zcm = zeros(length(t)-1,1);
zeta = zeros(length(t)-1,1);
rv = zeros(length(t)-1,mm-1);
eta = zeros(length(t)-1,mm-1);
F = zeros(length(t)-1,1);
time = nan(length(t)-1,1);

%%% initial conditions
h0 = R;%+0.99*abs(U)*deltat;%+0.05*R;%+1.00*abs(U)*dt;%R-0.01*R;
zcm(1) = h0;
zeta(1) = U; % impact speed
rv(1,:) = 0;
eta(1,:) =  U/2/(mm-1);
a(1,:) = 0;
alpha(1,:) = U/2/(nm); %
r0 = 0; %impact point (the point where we should evaluate the sums that appear in the forcing term)
Y = [zcm,zeta,rv,eta,a,alpha,F]; %fill solution array
clearvars zcm zeta a alpha rv eta F %delete old placeholdrs

dt_init = deltat;
%strings for filenames saving
rstring = num2str(R*100);
idn1 = rstring=='.';
rstring(idn1) = 'p';
vstring = num2str(abs(U)*100);
idn1 = vstring=='.';
vstring(idn1) = 'p';
rhostring = num2str(rho);
idn1 = rhostring=='.';
rhostring(idn1) = 'p';
fname = ['d0Mat_Rho' rhostring 'R' rstring 'NM' num2str(nm) ptype 'b' num2str(b/R) '.mat']; 
fname2 = ['c0Mat_Rho' rhostring 'R' rstring 'MM' num2str(mm) ptype 'b' num2str(b/R) '.mat']; 
%
Rmult = 1;
% [PSI,ZCM,re] = zcmVareaEllipse(rho,R,Rmult,sig);
PSI = linspace(0.01,89,200); % sweep over various values psi 
pa= R^2/2; % integration constant to set C in pressure shape function
if ~exist(fname,'file') || ~exist(fname2,'file')
    aa(1) = 200; %pre allocate
    aa(2) = 1;
    d_0Mat = zeros(length(u),aa(1),aa(2));
    c_0Mat = zeros(length(ll),aa(1));
    for kl = 1:aa(2) %new
        for i = 1:aa(1) % new psi ( aka rc)
            if strcmpi(ptype,'tophat') % determine constant C for pressure function
                 pc = 1*ones(1,length(r));
                 idx = r<=R*Rmult(kl)*sind(PSI(i));
                pc(~idx) = 0;
                press_c(i) = pa/trapz(r,r.*pc); 
            elseif strcmpi(ptype,'parabola')
                pc = 1*(1-(r./(R*sind(PSI(i)))).^2);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                pc(~idx) = 0;
                press_c(i) = pa/trapz(r,r.*pc); 
            elseif strcmpi(ptype,'poly4')
                pc = 1*(1-(r./(R*sind(PSI(i)))).^4);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                pc(~idx) = 0;
                press_c(i) = pa/trapz(r,r.*pc); 
            elseif strcmpi(ptype,'poly6')
                pc = 1*(1-(r./(R*sind(PSI(i)))).^6);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                pc(~idx) = 0;
                press_c(i) = pa/trapz(r,r.*pc); 
            elseif strcmpi(ptype,'poly8')
                pc = 1*(1-(r./(R*sind(PSI(i)))).^8);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                pc(~idx) = 0;
                press_c(i) = pa/trapz(r,r.*pc); 
            elseif strcmpi(ptype,'bell')
                pc = ones(1,length(r));
                lambda=0.75*R*sind(PSI(i));
                idx = find(r<R*sind(PSI(i)) & r>lambda);
                rl = R*sind(PSI(i));
                pc(idx) = 1./(exp(1./(lambda/rl-abs(r(idx)./rl))+1./(1-abs(r(idx)./rl)))+1);
                idx2 = find(r>R*sind(PSI(i)));
                pc(idx2) = 0;
                press_c(i) = pa/trapz(r,r.*pc); %*2*pi
            end
            if strcmpi(ptype,'tophat') % actually make pressure function H_r
                h_r = press_c(i)*ones(1,length(r));
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                h_r(~idx) = 0;
            elseif strcmpi(ptype,'parabola')
                h_r = press_c(i).*(1-(r./(R*sind(PSI(i)))).^2);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                h_r(~idx) = 0;
            elseif strcmpi(ptype,'poly4')
                h_r = press_c(i).*(1-(r./(R*sind(PSI(i)))).^4);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                h_r(~idx) = 0;
            elseif strcmpi(ptype,'poly6')
                h_r = press_c(i).*(1-(r./(R*sind(PSI(i)))).^6);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                h_r(~idx) = 0;
            elseif strcmpi(ptype,'poly8')
                h_r = press_c(i).*(1-(r./(R*sind(PSI(i)))).^8);
                idx = r<=R*Rmult(kl)*sind(PSI(i));
                h_r(~idx) = 0;
            end
            smt = sinc((k)./k(end)); % function includes pi in the argument
            pre = 2./((b^2)*(besselj(l+1,u.*b)).^2); %prefactor
            int = r.*h_r.*besselj(l,r.*u); % this is what u integrate
            if strcmpi(ptype,'tophat')
                d_0Mat(:,i,kl) = smt.*pre.*trapz(r,int,2); % these are the fourier-bessel coefficients with smoothing for top hat by Storey (1968)
            else
                d_0Mat(:,i,kl) = pre.*trapz(r,int,2); % these are the fourier-bessel coefficients
            end
            for li = 1:mm+1 % now do decompostions for droplet!
                PLM = legendre(ll(li),cos(theta));
                pred = 2*pi;
                if strcmpi(ptype,'tophat')
                     hh_r = press_c(i)*ones(1,length(theta));
                elseif strcmpi(ptype,'parabola')
                    hh_r = press_c(i).*(1-(sin(theta)./(sind(PSI(i)))).^2);
                elseif strcmpi(ptype,'poly4')
                    hh_r = press_c(i).*(1-(sin(theta)./(sind(PSI(i)))).^4);
                elseif strcmpi(ptype,'poly6')
                    hh_r = press_c(i).*(1-(sin(theta)./(sind(PSI(i)))).^6);
                elseif strcmpi(ptype,'poly8')
                    hh_r = press_c(i).*(1-(sin(theta)./(sind(PSI(i)))).^8);
                elseif strcmpi(ptype,'bell')
                    hh_r = press_c(i)*ones(1,length(theta));
                    t1 = (pi/180)*PSI(i);
                    t2 = (asin(lambda/R));
                    idx = find(theta<t1 & theta>real(t2));
                    hh_r(idx) = press_c(i)./(exp(1./(t2/t1-abs(theta(idx))./t1)+1./(1-abs(theta(idx))./t1))+1);
                    idx2 = find(theta>t1);
                    hh_r(idx2) = 0;
                end
                idx = 180/pi*theta<=asind(Rmult(kl)*sind(PSI(i)));
                hh_r(~idx) = 0;          
                intd = hh_r.*PLM(1,:).*sin(theta);
                c_0Mat(li,i) = trapz(theta,intd); %pred.*
            end
            smt2 = sinc(ll./ll(end));
            if strcmpi(ptype,'tophat')
                c_0Mat(:,i) = smt2'.*c_0Mat(:,i);
            else
                c_0Mat(:,i) = c_0Mat(:,i);
            end
        end
        pause(1) 
    end
    save(fname,'d_0Mat');
    save(fname2,'c_0Mat');
else
    load(fname)
    load(fname2)
end
clear d_0Mat
load(fname);
load(fname2);
id = find(Rmult == 1);
d_0 = d_0Mat(:,1,id); % initialize with smallest contact radius first
flag=0; %force on or off flag
nn = 1;
time(1) = 0;
%%%%%%%%%%%%%%some final preallocation
iidd = zeros(length(t),1);
iidd2 = zeros(length(t),1);
rc = zeros(length(t),1);
ZCM=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while time(nn) < tf 
    dt = deltat; % set timestep
    Re = abs(Y(nn,2))*R/nu_a; % reynolds number
    c_dd = (1+0.15*(Re^(0.687)/Re));
    c_f = 48*rho_a*nu_a*pi*R*c_dd;
    % get force f !!
    if Y(nn,1) <= R  && ~flag %initializing sim close, therefore we don't relly need this if statement until after
       % the drop rebounds and comes back down 
        time(nn);     
        iid2 = id;       %for debugging
        iidd2(nn) = iid2;
        iid = id;
        iidd(nn)=iid; % save for debugging                              
        %%%%%%%%%%%%%%%%011122 . just check overlap and see if that's
        %%%%%%%%%%%%%%%%good
        io = find(r<1.5*R);
        shape = R+sum(Y(nn,3:mm+1).*sph(3:end,1:end-1000)'./Nl(3:end),2); % get instantaneous drop interface shape
        inte = sum(Y(nn,2*mm+1:nm+2*mm)'.*besselj(0,k.*r(1:io(end))),1);  % get instantaneous bath interface shape   
        xx = intersections(shape'.*cos(theta(1:end-1000)-pi/2),shape'.*sin(theta(1:end-1000)-pi/2)+Y(nn,1),r(1:io(end)),inte,0); % get r_c(nn)
        if ~isempty(xx)
            ra = xx(end);   %grab largest
        else
            ra = rc(nn-1);
        end
        if ra > R*0.9998 %saturate if rc is greater than R (avoid R exactly 1 because of proejction issues)
            ra = R*0.9998;
        elseif ra < sind(PSI(2))*R % set lower limit (really small ~ 0.01R)
            ra = R*sind(PSI(2));
        end
        psi = asind(ra./R); %get angle between vertical and contact raidus
        if nn>=1
            rc(nn) = ra; % save rc
        end

        if nn>=1
            d_0 = interp2(PSI,u,d_0Mat,psi,u); %interpolate between pre computed d_0's and c_0's ( speeds up codE)
            c_0 = interp2(PSI,ll,c_0Mat,psi,ll);
            if isnan(d_0) % make sure it is not outside array bounds( just a ccheck)
                d_0 = d_0Mat(:,iid,iid2);
                c_0 = c_0Mat(:,iid);
            end
        end
        if ~exist('iid','var')
            iid=1;
            d_0 = d_0Mat(:,iid,iid2); %choose the d_0's
            iidd(nn)=iid; %save for debugging
            c_0 = c_0Mat(:,iid); % choose the c_0's 
            rc(nn) = re(iid,iid2)*sind(PSI(iid));
        end
        if Y(nn,end) >= 0 %Y(nn,1)<=Y(1,1)%
            y = Y(nn,:)';        
            %%%assemble matrix Q
            %%%%%%% step
            c_f=0;
            dth = dt;
            Qhalf = [1,-dth,zeros(1,mm-1),zeros(1,mm-1),zeros(1,nm),zeros(1,nm),0;...
            0,1+dth*c_f/m,zeros(1,mm-1),zeros(1,mm-1),zeros(1,nm),zeros(1,nm), -(dth/m);...
            zeros(mm-1,1),zeros(mm-1,1),1*eye(mm-1),-dth*eye(mm-1),zeros(mm-1,nm),zeros(mm-1,nm),zeros(mm-1,1);...
            zeros(mm-1,1),zeros(mm-1,1),wl2(3:end).*dth.*eye(mm-1), eye(mm-1)+alp(3:end).*dth.*eye(mm-1),zeros(mm-1,nm),zeros(mm-1,nm),(dth/(2*pi*rho*R^3)).*c_0(3:end).*ll(3:end)'.*(2*ll(3:end)'+1);...%Nl(3:end)'.* % Nl'.*(dth/(rho*R^4)).*ones(mm,1);... %%
            zeros(nm,1),zeros(nm,1),zeros(nm,mm-1),zeros(nm,mm-1),eye(nm),-dth*eye(nm),zeros(nm,1);...
            zeros(nm,1),zeros(nm,1),zeros(nm,mm-1),zeros(nm,mm-1),eye(nm).*dth.*(k'.^2*(sig/rho)+g).*k'.*tanh(k'.*H),eye(nm).*(1+dth.*4*nu.*k'.^2),dth.*(1/(rho*pi*R^2)).*k.*d_0.*tanh(k.*H);...
            1,0,-sph(3:end,1)'./Nl(3:end),zeros(1,mm-1),-1*ones(1,nm).*besselj(l,k'.*r0),zeros(1,nm),0];
            FYhalf = [y(1);y(2)-g*dth;y(3:mm+1);y(mm+2:2*mm);y(2*mm+1:nm+2*mm);y(nm+2*mm+1:end-1);R];
            Yhalf =Qhalf\FYhalf;
            Y(nn+1,:) = Yhalf;
          
        else %means this is last timestep in current contact
            c_f=0;
            nn=nn-1; % go back 1 step and reset force =0 since last step gave negative value
            y=Y(nn,:)';
            dth = dt;%propgate in free flight
            Qhalf = [1,-dth,zeros(1,mm-1),zeros(1,mm-1),zeros(1,nm),zeros(1,nm);...
            0,1+dth*c_f/m,zeros(1,mm-1),zeros(1,mm-1),zeros(1,nm),zeros(1,nm);...
            zeros(mm-1,1),zeros(mm-1,1),1*eye(mm-1),-dth*eye(mm-1),zeros(mm-1,nm),zeros(mm-1,nm);...
            zeros(mm-1,1),zeros(mm-1,1),wl2(3:end).*dth.*eye(mm-1), eye(mm-1)+alp(3:end).*dth.*eye(mm-1),zeros(mm-1,nm),zeros(mm-1,nm);...
            zeros(nm,1),zeros(nm,1),zeros(nm,mm-1),zeros(nm,mm-1),eye(nm),-dth*eye(nm);...
            zeros(nm,1),zeros(nm,1),zeros(nm,mm-1),zeros(nm,mm-1),eye(nm).*dth.*(k'.^2*(sig/rho)+g).*k'.*tanh(k'.*H),eye(nm).*(1+dth.*4*nu.*k'.^2);];
            FYhalf = [y(1);y(2)-g*dth;y(3:mm+1);y(mm+2:2*mm);y(2*mm+1:nm+2*mm);y(nm+2*mm+1:end-1)];
            f = 0; %force
            flag=1; %contact flag high (not in contact)
            Y(nn+1,1:end-1) =Qhalf\FYhalf;
            Y(nn+1,end) = f;
        end
    else %free flight
    %%%solve all ODE's!!!!!!
        y = Y(nn,:)';        
        Re = abs(Y(nn,2))*R/nu_a; % reynolds number 
        c_dd = (1+0.15*(Re^(0.687)/Re));
        c_f = 48*rho_a*nu_a*pi*R*c_dd; %can turn on skin friciton if you want, does not affect results in our parameter range
        c_f=0;
        %%%assemble matrix Q
        dth = dt;%./2;
        Qhalf = [1,-dth,zeros(1,mm-1),zeros(1,mm-1),zeros(1,nm),zeros(1,nm);...
        0,1+dth*c_f/m,zeros(1,mm-1),zeros(1,mm-1),zeros(1,nm),zeros(1,nm);...
        zeros(mm-1,1),zeros(mm-1,1),1*eye(mm-1),-dth*eye(mm-1),zeros(mm-1,nm),zeros(mm-1,nm);...
        zeros(mm-1,1),zeros(mm-1,1),wl2(3:end).*dth.*eye(mm-1), eye(mm-1)+alp(3:end).*dth.*eye(mm-1),zeros(mm-1,nm),zeros(mm-1,nm);...
        zeros(nm,1),zeros(nm,1),zeros(nm,mm-1),zeros(nm,mm-1),eye(nm),-dth*eye(nm);...
        zeros(nm,1),zeros(nm,1),zeros(nm,mm-1),zeros(nm,mm-1),eye(nm).*dth.*(k'.^2*(sig/rho)+g).*k'.*tanh(k'.*H),eye(nm).*(1+dth.*4*nu.*k'.^2);];
        FYhalf = [y(1);y(2)-g*dth;y(3:mm+1);y(mm+2:2*mm);y(2*mm+1:nm+2*mm);y(nm+2*mm+1:end-1)];
   
        f = 0;
        flag=1;
        Y(nn+1,1:end-1) =Qhalf\FYhalf;
        Y(nn+1,end) = f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% condition for turning force back on. only applicable for
        %%%% multiple bounces
        if max(Y(:,1)) > Y(1,1) && find(Y(:,1)==max(Y(:,1))) < nn-10 %%% aka if it has bounced at least once
            if Y(nn,1) <=R % if we are close, turn flag back on so that we can start checking force
                flag = 0;
            end
        end
        %%%%%%%%
    end
    time(nn+1) = time(nn)+dt; %step forward
    nn=nn+1;
end
%compute output parameters
idx = find(Y(:,end)==0,2); 
if length(idx) < 2 %force never turns off => no bounce
    c_r = 0;
else
    c_r = abs(Y(idx(end),2)/Y(1,2));
end
t=time;
out(1) = c_r;
tc = t(idx(end));
out(2) = tc;
z = abs(min(Y(:,1)));
out(3) = z;
out(4) = We;
out(5) = Bo;
out(6) = Oh;
idx2 = find(Y(:,1) >= Y(1,1),2);
if length(idx2) < 2
    out(7) = 0;
else
    out(7) = abs(Y(idx2(end),2)/Y(idx2(1),2));
end
out(8) = t(idx2(end));
end
