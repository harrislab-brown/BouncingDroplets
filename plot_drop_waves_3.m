function plot_drop_waves_3(YY,k,b,nm,mm,R,t_1P,d2d3) 
    %close al
    TT = t_1P;
    interfaceT = 0.0:0.1:19.5;%49.5;
    t_DNS = interfaceT.*(R./abs(YY(1,2)));
    r =linspace(0,b,5000); %size of domain
    zcm = YY(:,1);
    theta=linspace(0,pi,4000);
    sph = nan(mm,length(theta));
    ll = [0:mm];
    Nl = sqrt((2.*ll+1)/(4*pi));
    for i=1:mm+1
        PLM=legendre(ll(i),cos(theta));
        sph(i,:) = Nl(i).*PLM(1,:);
    end
    a = YY(:,2*mm+1:nm+2*mm); % amplitude functions
    k = k(1:nm);
    n = length(a); %number of time steps to visualize
    N = 1; %how often to compute surface shape (don't want all time steps takes too long)    
    htt = zeros(length(a),length(r)); %preallocation 
    %making histogram plots
    rstring = num2str(R*100);
    idn1 = rstring=='.';
    rstring(idn1) = 'p';
    vstring = num2str(abs(YY(1,2)*100));
    idn1 = vstring=='.';
    vstring(idn1) = 'p';
    vidName = ['AxiSymSphere_R' rstring 'Vi' vstring 'DropModes' num2str(mm) '_waveamplitudevswavenumber.avi'];
    writerObj = VideoWriter(vidName,'Motion JPEG AVI');
    writerObj.FrameRate = 30;
    writerObj.Quality = 50;
    open(writerObj);
    figure('Position',[0 0 512 756])
    for i=1:N:n % make a loop to get surface height htt(t)
        sums = 0; 
        ht = a(i,:);
        for j = 1:length(k)
            ss = ht(j)*besselj(0,k(j)*r);
            sums = ss+sums;
        end
        htt(i,:) = sums;
        bar(k,abs(ht));
        xlabel('$|k|$','Interpreter','latex');
        ylabel('$|a_i|$','Interpreter','latex');
        title(['R = ' num2str(R) ' V_i = ' vstring '  Oil'])
        annotation('textbox',[.75 .75 .1 .1],'String',['t= ' num2str(TT(i))])
        writeVideo(writerObj, getframe(gcf));
        delete(findall(gcf,'type','annotation'))
    end
    close(writerObj);
    rstring = num2str(R*100);
    idn1 = rstring=='.';
    rstring(idn1) = 'p';
    vstring = num2str(abs(YY(1,2)*100));
    idn1 = vstring=='.';
    vstring(idn1) = 'p';
    vidName = ['AxiSymSphere_R' rstring 'Vi' vstring 'DropModes' num2str(mm) 'crossing.avi'];
    writerObj = VideoWriter(vidName,'MPEG-4');
    writerObj.FrameRate = 5;
    writerObj.Quality = 100;
    open(writerObj);
    t2 = [r, fliplr(r)];
    tshift = 0.0005;
    count=1;
    t_DNS = t_DNS - tshift;
    figure(10)
    colormap(winter)
    for tt=1:N:n
         clf('reset')
        set(gca, 'color', [17 17 17]/255)
        if d2d3==0
            inBetween = [htt(tt,:);fliplr(htt(tt,:)-1e-3*ones(1,length(htt(tt,:))))];
            h2=area(r,[htt(tt,:); -1*ones(1,length(r))]','EdgeColor',[0 0.45 0.74]); % plot surface
            h2(1).FaceColor = [0.95 0.95 0.95];
            h2(2).FaceColor= [0 0.45 0.74];
            hold on
            h1=area(-r,[htt(tt,:); -1*ones(1,length(r))]','EdgeColor',[0 0.45 0.74]); % plot surface
            h1(1).FaceColor = [0.95 0.95 0.95];
            h1(2).FaceColor= [0 0.45 0.74];

            reff = R+sum(YY(tt,3:mm+1).*sph(3:end,:)'./Nl(3:end),2);
            count=count+1;
            plot(reff'.*cos(theta-pi/2),+reff'.*sin(theta-pi/2)+zcm(tt),'Color',[0 0.35 0.74],'LineWidth',3.0)
            plot(-reff'.*cos(theta-pi/2),+reff'.*sin(theta-pi/2)+zcm(tt),'Color',[0 0.35 0.75],'LineWidth',3.0)
            plot(0,zcm(tt)-reff(1),'rx') %marker for contact point
        
            io = find(r<1.5*R);
            shape = R+sum(YY(tt,3:mm+1).*sph(3:end,1:end-1000)'./Nl(3:end),2);
            inte = sum(YY(tt,2*mm+1:nm+2*mm)'.*besselj(0,k.*r(1:io(end))),1);     
            [xx,yy] = intersections(shape'.*cos(theta(1:end-1000)-pi/2),shape'.*sin(theta(1:end-1000)-pi/2)+YY(tt,1),r(1:io(end)),inte,0);
            if ~isempty(xx)
                idx =find(r<=xx(end));
                croppedshape = shape(1:idx(end)+1)'.*sin(theta(1:idx(end)+1)-pi/2)+YY(tt,1);
                area_overlap(tt) = 2*trapz(r(1:idx(end)+1),(inte(1:idx(end)+1)-croppedshape))./(pi*R^2);
            end
            axis([-b/6 b/6 -1e-3 2e-3]) %fix axes
            daspect([1 1 1])
            drawnow();
            set(gca,'visible','off')
            set(gca,'xtick',[])
            set(gcf,'units','pixels','position',[0 0 1080 756])
            set(gca,'units','pixels','position',[0 0 1080 756])
            writeVideo(writerObj, getframe(gca));
        else d2d3==1  
            ta = linspace(0,2*pi,50);
            reff = R+sum(YY(tt,3:mm+1).*sph(3:end,:)'./Nl(3:end),2);
            rr=reff'.*cos(theta-pi/2)./R;
            zr=reff'.*sin(theta-pi/2)+zcm(tt);
            xx1=bsxfun(@times,rr',cos(ta));
            yy1=bsxfun(@times,rr',sin(ta));
            zz1=repmat(zr',1,length(ta));
            xx=bsxfun(@times,r(1:2500)',cos(ta));
            yy=bsxfun(@times,r(1:2500)',sin(ta));
            zz=repmat(htt(tt,(1:2500))',1,length(ta));
            sd=surf(xx1*R,yy1*R,zz1);
            set(sd,'FaceColor',[0 .25 1],'FaceAlpha',1,'FaceLighting','gouraud','EdgeColor','none')
            hold on
            sb=surf(xx,yy,zz);
            shading flat
            axis equal            
            if tt>1
                sb.CData = (zz-zzo)./(t_1P(tt)-t_1P(tt-N));
                colormap(winter)
            end
            if tt>1
                sd.CData = (zz1-zz1o)./(t_1P(tt)-t_1P(tt-N))-YY(tt,2);
            else
                set(sd,'FaceColor',[0 .25 1],'FaceAlpha',1,'FaceLighting','gouraud','EdgeColor','none')
            end
            
            view(45,10)
            xlim([-b/2 b/2]) 
            ylim([-b/2 b/2])
            zlim([-1.5*R 6*R])%fix axes
            daspect([1 1 1])
            dim = [.2 .25 .3 .3];
            tsig = 1./sqrt(72.2e-3/(998.2*(0.35e-3)^3));
            str = num2str(TT(tt)/tsig);
            annotation('textbox',dim,'String',['$t/t_{\sigma} =  $' str],'FitBoxToText','on','Interpreter','latex');
            drawnow();
            set(gca,'visible','off')
            set(gca,'xtick',[])
            set(gcf,'units','pixels','position',[0 0 1080 756])
            set(gca,'units','pixels','position',[0 0 1080 756])
            light('style','local')
            writeVideo(writerObj, getframe(gca));
            hold off
            zz1o=zz1;
            zzo=zz;
        end
    end
    close(writerObj);
end