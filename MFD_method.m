%--------------------------------------------------------------------------
%  Code         :  Mode function decomposition (MFD) method is applied to extract
%                  the dominant equatorial Pacific thermocline modes.   
%  Author       :  KANG XU (South China Sea Institute of Oceanology, CAS)
%  E-mail       £∫ xukang0527@hotmail.com
%  Fristwritten :  Sep. 30, 2015
%  Ref.         :  Xu, K., R.-X. Huang, W. Wang, C. Zhu, and R. Lu, 2017: 
%                  Thermocline fluctuations in the equatorial Pacific related 
%                  to the two types of El Ni?o events,submitted to J. Clmite, 
%                  minor revision; 
%  Data source  :  SODA 2.1.6 dataset
%--------------------------------------------------------------------------
clear all; close all;clc;
%-- Input monthly equatorial Pacific thermocline depth (thcd) data from SODA 2.1.6
%-- thcd(nx,nt), nx = 720, nt = 612;
load Eq_thcd_monthly_soda_58-08.mat;
thcd = thcd';
%--- Longtitude resolution
lon = [0.25 : 0.5 : 359.75];
imt = length(lon);
nt  = size(thcd,1);
%--- Choose the region of Pacific Ocean  (130°„E-180°„-80°„W);
for i  = 263 : 558;
    ii = i - 262;
    thcP(:,ii) = thcd(:,i);
end
clear i ii imt;
%--  Scope of longtitude for the Pacific Ocean
llon = lon(263:558);

%-- Cal. thcP anomalies (removed annual cycle)
mon = 12;
nyr = nt/mon;
nx  = size(thcP,2);
thcPm = reshape(thcP',nx,mon,nyr);
thcPm_month = nanmean(thcPm,3);
for ix  = 1 : nx;
    for iyr = 1 : nyr;
        xx = squeeze(thcPm(ix,:,iyr));
        yy = xx - thcPm_month(ix,:);
        slePa(ix,:,iyr) = yy; 
    clear xx yy;
    end
end
sleP = reshape(slePa,nx,nt);
sleP = sleP';
clear thcPm* ix iyr imt slePa;

%-- Set up the related parameters 
ipt   = length(thcP(1,:));      %---  longitudinal grid-point number;
dxip  = 2/(ipt-1);              %---  the intervals of [-1,1];
xip   = [-1:dxip:1];            %---  xip;
nmt   = 70;                     %---  nmt means the frist 70 modes;

%-- Construct a series of functions f(x) satisfying the Eq. A1 in the Referrence;
for n  = 1 : nmt;
    rn = double(n);
    for k  = 1 : ipt;
        xk = xip(k);
        i=double(floor(rn*(1-xk)/2)+1);                 %--- Take the smallest integer;
        fnk=(-1)^(i-1)*(1+rn*(xk-(1-2.0*(i-1)/rn)));    %--- From (+1,+1) to (-1,-1);
        fnk=(-1)^(n-1)*fnk;                             %--- Transfer; From (-1,-1) to (+1,+1); 
        f(n,k)=fnk;
    end
end
clear fnk xk n i k rn thcP; 

%-- Draw Figure [First six mode functions for the equatorial thermocline]
figure;
set(gcf,'paperunit','centimeters');
set(gcf,'paperposition',[1 10 20 12]);
plot(xip,f(1,:),'k','linewidth',3);hold on;
plot(xip,f(2,:),'r','linewidth',3);hold on;
plot(xip,f(3,:),'b','linewidth',2);hold on;
plot(xip,f(4,:),'g','linewidth',2);hold on;
plot(xip,f(5,:),'--m','linewidth',1);hold on
plot(xip,f(6,:),'--c','linewidth',1);hold on;
L1=['M 1'];L2=['M 2'];L3=['M 3'];L4=['M 4'];L5=['M 5'];L6=['M 6'];
hl = legend (L1,L2,L3,L4,L5,L6,3);
    set(hl,'box','off','fontsize',14,'FontWeight','bold');
    axis ([-1 1 -1 1]);
    set(gca,'fontsize',14)
    set(gca,'ytick',[-1.0,-0.5,0,0.5,1.0])
    set(gca,'xtick',[-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0])
ylabel ('Depth (m)');
    set(gca,'yDir','reverse');
title ('Mode Functions'); 
print -depsc Mode_functions.eps;
clear L1 L2 L3 L4 L5 L6;

%---  Cal. Normalization of a series of functions f(x); 
for n   = 1:nmt;
    ff  = squeeze(f(n,1:ipt));
    xms = 0.0;
    for i = 1:ipt;
        xms = xms + ff(i)*ff(i)*dxip;
    end
    ff = ff/xms^0.5;
    fx(n,:) = ff;
end
%  amp=abs(fx(1,1));    %-- the amplitude of f(x);
clear n ff xms i;

%-- Draw Fig. A1 [First six normalized mode functions for the equatorial thermocline]
figure
set(gcf,'paperunit','centimeters');
set(gcf,'paperposition',[1 10 20 12]);
plot(xip,fx(1,:),'k','linewidth',3);hold on;
plot(xip,fx(2,:),'r','linewidth',3);hold on;
plot(xip,fx(3,:),'b','linewidth',2);hold on;
plot(xip,fx(4,:),'g','linewidth',2);hold on;
plot(xip,fx(5,:),'--m','linewidth',1);hold on
plot(xip,fx(6,:),'--c','linewidth',1);hold on;
L1=['M 1'];L2=['M 2'];L3=['M 3'];L4=['M 4'];L5=['M 5'];L6=['M 6'];
hl = legend (L1,L2,L3,L4,L5,L6,3);
    set(hl,'box','off','fontsize',14,'FontWeight','bold');
    axis ([-1 1 -1.2186 1.2186]);
    set(gca,'fontsize',14);
    set(gca,'ytick',[-1.0,-0.5,0,0.5,1.0]);
    set(gca,'xtick',[-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0]);
ylabel ('Depth (m)');
    set(gca,'yDir','reverse');
title ('Normalized Mode Functions'); 
print -depsc Mode_functions_Normalization.eps;
clear L1 L2 L3 L4 L5 L6;

%--  Based on Eq.A2 in ref., the equatorial thermocline depth anomalies D(x, t) can 
%--  be projected onto these mode function series for any given time t;
%--  Subsequently, s(t,i)time evolution for each mode i is obtained. 
for it = 1:nt;
    for n = 1:nmt
        d_m = nanmean(sleP,2);
        d(it,:) = sleP(it,:) - d_m(it);
        s(it,n) = nansum(d(it,:).*fx(n,:))*dxip;
    end
end
clear it n;

%--- Cal. the anomalies of s(t,i);
ky=nt/12.0;
for n  = 1:nmt;
    sx = squeeze(s(1:nt,n));
    ss = [0,0,0,0,0,0,0,0,0,0,0,0];
    for kr  = 1:ky;
        for imon = 1:12;
            kmon = imon+12*(kr-1);           
            ss(imon) = ss(imon) + sx(kmon);
        end    
    end
    for kr = 1:ky;
        for imon = 1:12;
            kmon = imon+12*(kr-1);
            saa(kmon)  = sx(kmon)-ss(imon)/ky;
            sa(kmon,n) = saa(kmon);
        end
    end
end
clear kr imon kmon;

%-- Cal. Normalization of sa(t,i) and Fractional variance
for n  = 1:nmt;
    sy = squeeze(sa(1:nt,n));
    sy_std = std(sy);
    for it = 1:nt;
        sss(it,n) = (sy(it)-mean(sy))/sy_std;
    end
    sd_s(n) = sy_std; 
end

for n = 1 : nmt;
    varexp(n) = sd_s(n)/nansum(sd_s(:))*100;  
end
%--  The rule of thumb of North et al. (1982).
delt_var = varexp * sqrt(2/(nmt));
clear n it;

%-- Draw Figure [Normalized time evolution of the first mode (M1) and Second mode (M2) index] 
%-- Seen in Fig. 3 of the ref.;
figure;
xplot=[1,nt];yplot=[0,0];
    set(gcf,'paperunit','centimeters');
    set(gcf,'paperposition',[1.5 0.5 18.0 28.0])
h1=subplot(2,1,1);
plot(sss(:,1),'b','linewidth',0.8);hold on;
plot(xplot,yplot,'k','linewidth',0.5);
% plot(xplot,yplot1,'--k','linewidth',0.8);
% plot(xplot,yplot2,'--k','linewidth',0.8);
text(0.01,1.07,'(a) Normalized M1  ','sc','FontSize',14);
axis ([1 nt -5 5]);
hl=legend({'M1'},2);set(hl,'box','off','fontsize',12);%,'FontWeight','bold');
    set(gca,'FontSize',14);
    set(gca,'ytick',[-4:2:4]);
    set(gca,'xtick',[25:60:612]);
    set(gca,'xticklabel',['1960';'1965';'1970';'1975';'1980';'1985';'1990';'1995';'2000';'2005']);
%---
h2=subplot(2,1,2);
plot(sss(:,2),'r','linewidth',0.8);hold on;
plot(xplot,yplot,'k','linewidth',0.5);
text(0.01,1.07,'(b) Normalized M2','sc','FontSize',14);
axis ([1 nt -5 5]);
hl=legend({'M2'},2);set(hl,'box','off','fontsize',12);
    set(gca,'FontSize',14);
    set(gca,'ytick',[-4:2:4]);
    set(gca,'xtick',[25:60:612]);
    set(gca,'xticklabel',['1960';'1965';'1970';'1975';'1980';'1985';'1990';'1995';'2000';'2005']);
%---
set(h1,'position',[.10 .77 .8 .2]);
set(h2,'position',[.10 .51 .8 .2]);
print -depsc M1-M2-index_Normalized.eps; 

%--- Composites of the winter (DJF) accumulated sum of modes assocaited 
%    with EP and CP El Nino events;
%--- EP type of El Nino years: 1972, 1976, 1982, 1997, 2006;
%--- CP type of El Nino years: 1958, 1968, 1977, 1994, 2004;

%-- epm means the No. of the month from Jan. 1958 for EP El Nino winters;
epm = [180,181,182,228,229,230,300,301,302,480,481,482,588,589,590]; 
%-- cpm means the No. of the month from Jan. 1958 for CP El Nino winters;
cpm = [12,13,14,132,133,134,240,241,242,444,445,446,564,565,566];

%-- For EP type
for it  = 1 : size(epm,2);
    km  = epm(it);
        dsx0(it,:) = sleP(km,:); 
        dsx1(it,:) = s(km,1)*fx(1,:);
        dsx2(it,:) = dsx1(it,:) + s(km,2)*fx(2,:);
        dsx3(it,:) = dsx2(it,:) + s(km,3)*fx(3,:);
        dsx4(it,:) = dsx3(it,:) + s(km,4)*fx(4,:);
        dsx5(it,:) = dsx4(it,:) + s(km,5)*fx(5,:);
        dsx6(it,:) = dsx5(it,:) + s(km,6)*fx(6,:);
    clear km;
end

d0ep = nanmean(dsx0,1);  
d1ep = nanmean(dsx1,1);d2ep = nanmean(dsx2,1);d3ep = nanmean(dsx3,1);
d4ep = nanmean(dsx4,1);d5ep = nanmean(dsx5,1);d6ep = nanmean(dsx6,1);
clear dsx*;

%-- For CP type
for it  = 1 : size(cpm,2);
    km  = cpm(it);
        dsx0(it,:) = sleP(km,:);
        dsx1(it,:) = s(km,1)*fx(1,:);
        dsx2(it,:) = dsx1(it,:) + s(km,2)*fx(2,:);
        dsx3(it,:) = dsx2(it,:) + s(km,3)*fx(3,:);
        dsx4(it,:) = dsx3(it,:) + s(km,4)*fx(4,:);
        dsx5(it,:) = dsx4(it,:) + s(km,5)*fx(5,:);
        dsx6(it,:) = dsx5(it,:) + s(km,6)*fx(6,:);
    clear km;
end

d0cp = nanmean(dsx0,1);
d1cp = nanmean(dsx1,1);d2cp = nanmean(dsx2,1);d3cp = nanmean(dsx3,1);
d4cp = nanmean(dsx4,1);d5cp = nanmean(dsx5,1);d6cp = nanmean(dsx6,1);
clear dsx*;


%--  Draw [Fig.2 in the ref.] Fractional variances of the modes  
%    and Composited EP/CP cases for the accumulated sum of modes 

figure;
    set(gcf,'paperunit','centimeters');
    set(gcf,'paperposition',[0.5 0.5 24.0 28.0])
x1=[0.3,8.8];y1 = [9.10,9.10];
x2=[0.3,8.8];y2 = [9.54,9.54];
%-- For fractional variances of the modes 
h1=subplot(3,1,1);      %--  draw Fig.2a
    bar(varexp(1:8),'FaceColor',[0 153 255]/255,'BarWidth',0.6,'LineStyle','none');hold on;
    errorbar(varexp(1:8),delt_var(1:8),'k-','linewidth',1.5,'LineStyle','none');hold on;
    axis ([0.3 8.8 0 30]);grid on;
ylabel ('Fractional variance (%)','FontSize',14);
xlabel ('Modes','FontSize',14);
title ('(a) Fractional variance explained by the first eight modes', 'FontSize',14); 
    set(gca,'FontSize',14);
    set(gca,'ytick',[0:5:30]);
    set(gca,'xtick',[1:1:8]);
    set(gca,'tickdir','out','LineWidth',0.8);  

%-- For EP cases 
h2=subplot(3,1,2);    %--  draw Fig.2b
xP = (1:1:size(d0ep,2));
    plot(xP,d0ep,'k-','linewidth',1.0); hold on;
    plot(xP,d1ep,'r-','linewidth',1.0);hold on;
    plot(xP,d2ep,'b-','linewidth',2.0);hold on;
    plot(xP,d3ep,'g-','linewidth',2.0);hold on;
    plot(xP,d4ep,'m--','linewidth',1.5);hold on;
    plot(xP,d5ep,'c--','linewidth',1.5);hold on;
    plot(xP,d6ep,'y--','linewidth',1.2);hold on;
    plot(xP,d0ep,'k-','linewidth',1.5); hold on;
    plot(xP,d1ep,'r-','linewidth',1.0);hold on;
axis ([1 296 -40 40]);grid on;
ylabel ('Depth anomaly (m)', 'FontSize',14);
xlabel ('Longtitude','FontSize',14);  
title ('(b) EP El Nino, D20 anomaly and its reconstructions based on modes', 'FontSize',14); 
    set(gca,'ytick',[-40:20:40],'FontSize',14);
    set(gca,'xtick',[38:60:296],'FontSize',14);
    set(gca,'xticklabel',['150E';'180 ';'150W';'120W';' 90W'],'FontSize',14);
    set(gca,'tickdir','out','LineWidth',0.8); 
    set(gca,'yDir','reverse'); 
L1=['D20 anomaly'];L2=['Mode 1'];L3=['Modes 1-2'];L4=['Modes 1-3'];L5=['Modes 1-4'];L6=['Modes 1-5'];L7=['Modes 1-6'];
hl=legend (L1,L2,L3,L4,L5,L6,L7,3);
    set(hl,'box','off','fontsize',10);

%-- For CP cases 
h3=subplot(3,1,3);     %---     draw Fig.2c
    plot(xP,d0cp,'k-','linewidth',1.0); hold on;
    plot(xP,d1cp,'r-','linewidth',1.0);hold on;
    plot(xP,d2cp,'b-','linewidth',2.0);hold on;
    plot(xP,d3cp,'g-','linewidth',2.0);hold on;
    plot(xP,d4cp,'m--','linewidth',1.5);hold on;
    plot(xP,d5cp,'c--','linewidth',1.5);hold on;
    plot(xP,d6cp,'y--','linewidth',1.2);hold on;
    plot(xP,d0cp,'k-','linewidth',1.5); hold on;
    plot(xP,d2cp,'b-','linewidth',1.0);hold on;
hl=legend (L1,L2,L3,L4,L5,L6,L7,3);
    set(hl,'box','off','fontsize',10);
axis ([1 296 -30 30]);grid on;
ylabel ('Depth anomaly (m)', 'FontSize',14);
xlabel ('Longtitude','FontSize',14);
title ('(c) CP El Nino, D20 anomaly and its reconstructions based on modes', 'FontSize',14); 
    set(gca,'ytick',[-30:10:30],'FontSize',14);
    set(gca,'xtick',[38:60:296],'FontSize',14);
    set(gca,'xticklabel',['150E';'180 ';'150W';'120W';' 90W'],'FontSize',14);
    set(gca,'tickdir','out','LineWidth',0.8);  
    set(gca,'yDir','reverse');
    
set(h1,'position',[.18 .77 .6 .2]);
set(h2,'position',[.18 .48 .6 .2]);
set(h3,'position',[.18 .19 .6 .2]);
print -depsc Fractional_Varance-Composite_MFD_Modes_EP-CP_Elnino.eps; 


