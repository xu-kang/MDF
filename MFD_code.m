clear all; close all;clc;
load E:\Work\thermocline_md\data\Eq_thc_depth_weight_mean_mon_58-08.mat;
thcd=thcd';
%--- lon经度范围
lon=[0.25:0.5:359.75];
imt=length(lon);
nt = size(thcd,1);
%% %--- 截取三大洋格点信息  Atlantic包括两段，需要拼接
%---  draw Basin-mean depth of the 20^oC isothermal
%--- Atlantic (45°W--10°E)
for i=631:imt
    ii=i-630;
thcA(:,ii)=thcd(:,i);
end
for i=1:22
    ii=i+(imt-630);
thcA(:,ii)=thcd(:,i);
end
%--- Indian (45°E--100°E) （x=91-201）
for i=91:201
    ii=i-90;
thcI(:,ii)=thcd(:,i);
end
%--- Pacific (130°E--80°W)
for i=263:558
    ii=i-262;
thcP(:,ii)=thcd(:,i);
end
clear i ii;

llon = lon(263:558);

%--- 三大洋thermocline深度
for it=1:nt
   depth_P(it)=nanmean(thcP(it,:))*(-1.0);
   depth_I(it)=nanmean(thcI(it,:))*(-1.0);
   depth_A(it)=nanmean(thcA(it,:))*(-1.0);
end
% 
% figure;
%     set(gcf,'paperunit','centimeters');
%     set(gcf,'paperposition',[1.5 3.0 18.0 24.0])
% %     set(gcf,'PaperPositionMode','manual');
% h1=subplot(3,1,1);
% plot(depth_P(1:612),'k','linewidth',1);hold on;
%     set(gca,'FontSize',14);
%     set(gca,'xtick',[25:60:612]);
%     set(gca,'xticklabel',['1960';'1965';'1970';'1975';'1980';'1985';'1990';'1995';'2000';'2005']);
%     set(gca,'ytick',[-150:10:-100]);
% text(0.02,0.88,'a) Pacific','sc','FontSize',14);
% ylabel ('Depth (m)','FontSize',14);
% title ('Basin-mean depth of the 20^oC isothermal','FontSize',16);
% axis ([1 612 -155 -95]);
% %------
% h2=subplot(3,1,2);
% plot(depth_A(1:612),'k','linewidth',1);hold on;
%     set(gca,'FontSize',14);
%     set(gca,'xtick',[25:60:612]);
%     set(gca,'xticklabel',['1960';'1965';'1970';'1975';'1980';'1985';'1990';'1995';'2000';'2005']);
%     set(gca,'ytick',[-100:10:-60]);
% text(0.02,0.88,'b) Atlantic','sc','FontSize',14);
% ylabel ('Depth (m)','FontSize',14);
% axis ([1 612 -105 -60]);
% %----
% h3=subplot(3,1,3);
% plot(depth_I(1:612),'k','linewidth',1);hold on;
%     set(gca,'FontSize',14);
%     set(gca,'xtick',[25:60:612]);
%     set(gca,'xticklabel',['1960';'1965';'1970';'1975';'1980';'1985';'1990';'1995';'2000';'2005']);
%     set(gca,'ytick',[-130:10:-100]);
% text(0.02,0.88,'c) Indian','sc','FontSize',14);
% ylabel ('Depth (m)','FontSize',14);
% xlabel ('Year','FontSize',14)
% axis ([1 612 -135 -95]);
% %---
% set(h1,'position',[.12 .60 .85 .2]);
% set(h2,'position',[.12 .35 .85 .2]);
% set(h3,'position',[.12 .10 .85 .2]);
% print -depsc E:\Work\thermocline_md\pic\Eq_thcd_basin_mean_matlab.eps;


%% %--- 计算Equatorial Pacific Thermocline Modes
%--- 参数设置
% x = [0.25:0.5:359.75];
% xP=lon(263:559);  
sleP=thcP;

%---计算异常值;
[slePa,sleP_mean]=anomaly(sleP,12);
sleP = slePa; 
clear slePa sleP_mean;

ipt=length(thcP(1,:));   %--- Pacific格点数
dxip=2/(ipt-1);          %--- -1到1之间格局间隔
xip=[-1:dxip:1];         %--- xip的矩阵
nmt=70;                   %--- nmt指计算前nmt模态

for n=1:nmt
    rn=double(n);
for k=1:ipt
    xk=xip(k);
%     i=double(int16(rn*(1-xk)/2-0.5)+1);           %--- int16取整数,注意需要-0.5,只能变成最小整数
      i=double(floor(rn*(1-xk)/2)+1);               %--- 取最小整数
      fnk=(-1)^(i-1)*(1+rn*(xk-(1-2.0*(i-1)/rn)));  %--- 首先计算（1,1）点出开始积分，但这样不太好理解； 
      fnk=(-1)^(n-1)*fnk;                           %--- 小技巧，将首先在（-1，-1）处开始积分
    f(n,k)=fnk;
end
end
clear fnk xk n i k rn thcP; 

% %%----draw 理想化theromcline前8个模态分解；
% figure
% set(gcf,'paperunit','centimeters');
% set(gcf,'paperposition',[1 10 20 12]);
% plot(xip,f(1,:),'k','linewidth',3);hold on;
% plot(xip,f(2,:),'r','linewidth',3);hold on;
% plot(xip,f(3,:),'b','linewidth',2);hold on;
% plot(xip,f(4,:),'g','linewidth',2);hold on;
% plot(xip,f(5,:),'--m','linewidth',1);hold on
% plot(xip,f(6,:),'--c','linewidth',1);hold on;
% L1=['M 1'];L2=['M 2'];L3=['M 3'];L4=['M 4'];L5=['M 5'];L6=['M 6'];
% hl = legend (L1,L2,L3,L4,L5,L6,3);
%  set(hl,'box','off','fontsize',14,'FontWeight','bold');
% axis ([-1 1 -1 1]);
% set(gca,'fontsize',14)
% set(gca,'ytick',[-1.0,-0.5,0,0.5,1.0])
% set(gca,'xtick',[-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0])
% ylabel ('Depth (m)');
% set(gca,'yDir','reverse');
% title ('Mode Functions'); 
% print -depsc E:\Work\thermocline_md\pic\Mod_functions-y-reverse.eps;
% clear L1 L2 L3 L4 L5 L6;

%---对f（x）函数进行标准化 
for n=1:nmt
    ff=squeeze(f(n,1:ipt));
    xms=0.0;
    for i=1:ipt
        xms=xms+ff(i)*ff(i)*dxip;
    end
    ff=ff/xms^0.5;
    fx(n,:)=ff;
end
amp=abs(fx(1,1))*2;
clear n ff xms i;

% % % for i=1:nmt
% % %     f_std(i)=nanstd(f(i,:));
% % %     fxx(i,:)=f(i,:)/f_std(i);
% % % end

% % % % draw 均一化的理想theromcline前8个模态分解；
% figure
% set(gcf,'paperunit','centimeters');
% set(gcf,'paperposition',[1 10 20 12]);
% plot(xip,fx(1,:),'k','linewidth',3);hold on;
% plot(xip,fx(2,:),'r','linewidth',3);hold on;
% plot(xip,fx(3,:),'b','linewidth',2);hold on;
% plot(xip,fx(4,:),'g','linewidth',2);hold on;
% plot(xip,fx(5,:),'--m','linewidth',1);hold on
% plot(xip,fx(6,:),'--c','linewidth',1);hold on;
% L1=['M 1'];L2=['M 2'];L3=['M 3'];L4=['M 4'];L5=['M 5'];L6=['M 6'];
% hl = legend (L1,L2,L3,L4,L5,L6,3);
%  set(hl,'box','off','fontsize',14,'FontWeight','bold');
% axis ([-1 1 -1.2186 1.2186]);
% set(gca,'fontsize',14);
% set(gca,'ytick',[-1.0,-0.5,0,0.5,1.0]);
% set(gca,'xtick',[-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0]);
% ylabel ('Depth (m)');
% set(gca,'yDir','reverse');
% title ('Normalized Mode Functions'); 
% print -depsc E:\Work\thermocline_md\pic\Mod_functions_Normalization-y-reverse.eps;
% clear L1 L2 L3 L4 L5 L6;

%% %--- 对d(x)*f（x）*delt x 进行积分，求得模态振幅s 
for it=1:nt
    for n=1:nmt
        d_m=nanmean(sleP,2);
        d(it,:)=sleP(it,:)-d_m(it);
        s(it,n)=nansum(d(it,:).*fx(n,:))*dxip;
    end
end

% figure;
% plot(d_m);

%  figure;
%  plot(s(:,1),'r','linewidth',1);
% save s s
% 
% %---------  
% clear all;close all; clc;
% load s;
% nmt=50;
nt=612;
%---  计算前8个模态的anomalise
ky=nt/12.0;
% ss=zeros(12,8);

for n=1:nmt
    sx=squeeze(s(1:nt,n));
    ss=[0,0,0,0,0,0,0,0,0,0,0,0];
    for kr=1:ky
        for imon=1:12
            kmon=imon+12*(kr-1);           
            ss(imon)=ss(imon)+sx(kmon);
        end    
    end
    for kr=1:ky
        for imon=1:12
            kmon=imon+12*(kr-1);
            saa(kmon)=sx(kmon)-ss(imon)/ky;
            sa(kmon,n)=saa(kmon);
        end
    end
end
% figure;
% plot(sa(:,1),'r','linewidth',1);

% % %---- 挑选年月合成，首选1997年11月
% yearst = 1958;
% year   = 1972;
% month  =  12;
% km     = (year - yearst)*12 + month;
% % km     = 600;
% dsx0 = sleP(km,:); %- d_m(km) ;
% dsx1 = d_m(km) + s(km,1)*fx(1,:);
% dsx2 = d_m(km) + s(km,1)*fx(1,:) + s(km,2)*fx(2,:);
% dsx3 = d_m(km) + s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:);
% dsx4 = d_m(km) + s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:) + s(km,4)*fx(4,:);
% dsx5 = d_m(km) + s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:) + s(km,4)*fx(4,:) + s(km,5)*fx(5,:);
% dsx6 = d_m(km) + s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:) + s(km,4)*fx(4,:) + s(km,5)*fx(5,:) + s(km,6)*fx(6,:);
 
% % dsx0 = sleP(km,:) - d_m(km) ;
% % dsx1 = s(km,1)*fx(1,:);
% % dsx2 = s(km,1)*fx(1,:) + s(km,2)*fx(2,:);
% % dsx3 = s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:);kmon
% % dsx4 = s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:) + s(km,4)*fx(4,:);
% % dsx5 = s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:) + s(km,4)*fx(4,:) + s(km,5)*fx(5,:);
% % dsx6 = s(km,1)*fx(1,:) + s(km,2)*fx(2,:) + s(km,3)*fx(3,:) + s(km,4)*fx(4,:) + s(km,5)*fx(5,:) + s(km,6)*fx(6,:);
% % 
% % 
% % b = s(485,2)
% 
% 
% 
% xP = (1:1:size(dsx0,2));
% figure
%    set(gcf,'paperunit','centimeters');
%    set(gcf,'paperposition',[3 10 12 12]);
% plot(xP,dsx0,'r','linewidth',2); hold on;
% plot(xP,dsx1,'k');hold on;
% plot(xP,dsx2,'b');hold on;
% plot(xP,dsx3,'m');hold on;
% % plot(xP,-dsx4,'k--');hold on;
% % plot(xP,-dsx5,'b--');hold on;
% % plot(xP,-dsx6,'g--');hold on;
% L1=['20^o Isothermal'];L2=['Mode 1'];L3=['Modes 1-2'];L4=['Modes 1-3'];L5=['Mode 1-4'];L6=['Modes 1-5'];L7=['Modes 1-6'];
% % legend (L1,L2,L3,L4,L5,L6,L7,4);
% % axis ([130 280 -220 20]);
% ylabel ('depth (m)');
% title ('20C isothermal mode decompostion, Nov. 1989'); 
% set(gca,'xtick',[38:60:296]);
% set(gca,'xticklabel',['150E';'180 ';'150W';'120W';' 90W']);
% % print -depsc E:\Work\thermocline_md\pic\Thermocline+Modes_Pac_Nov_1968.eps;


%--- EP和CP年份冬季DJF合成，我们选择EP：1972, 1976, 1982,1997,2006
%--- CP: 1958,1968,1977,1994,2004
epm = [180,181,182,228,229,230,300,301,302,480,481,482,588,589,590];                            %--- EP年的冬季月份
cpm = [12,13,14,132,133,134,240,241,242,444,445,446,564,565,566];    %--- CP年的冬季月份
% it = (1994 - 1958)*12+12;

aa = nanmean(d_m,1);
bb = nanmean(sleP,1);

for it  = 1:size(epm,2);
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

d0ep = nanmean(dsx0,1);  %- nanmean(d_m,1);
d1ep = nanmean(dsx1,1);d2ep = nanmean(dsx2,1);d3ep = nanmean(dsx3,1);
d4ep = nanmean(dsx4,1);d5ep = nanmean(dsx5,1);d6ep = nanmean(dsx6,1);
clear dsx*;

for it  = 1:size(cpm,2);
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

d0cp = nanmean(dsx0,1);% - nanmean(d_m,1);
d1cp = nanmean(dsx1,1);d2cp = nanmean(dsx2,1);d3cp = nanmean(dsx3,1);
d4cp = nanmean(dsx4,1);d5cp = nanmean(dsx5,1);d6cp = nanmean(dsx6,1);
clear dsx*;

% figure;
% plot(bb,'r');hold on;

% figure;
% plot(d0cp,'r');hold on;
% plot(d1cp,'b--');hold on;
% plot(d2cp,'k--');hold on;
% plot(d3cp,'m--');hold on;
% plot(d4cp,'c--');hold on;
% plot(d5cp,'g--');hold on;
% 
% figure;
% plot(d0ep,'r');hold on;
% plot(d1ep,'b--');hold on;
% plot(d2ep,'k--');hold on;
% plot(d5ep,'g--');hold on;


%% -%-  模态分解进行标准化和每个模态的方差贡献计算
for n=1:nmt
    sy=squeeze(sa(1:nt,n));
    sy_std=std(sy);
    for it=1:nt
        sss(it,n)=(sy(it)-mean(sy))/sy_std;
    end
    sd_s(n)=sy_std;
%    varexp(n)=sd_s(n)/nansum(sd_s(:))*100;    
end

for n=1:nmt
    varexp(n)=sd_s(n)/nansum(sd_s(:))*100;  
end
%--- north 检验
delt_var = varexp * sqrt(2/(nmt));

%% -- draw Standard deviation of the modes and Composited EP/CP cases for mode comdepostion
figure;
    set(gcf,'paperunit','centimeters');
    set(gcf,'paperposition',[0.5 0.5 24.0 28.0])
x1=[0.3,8.8];y1 = [9.10,9.10];
x2=[0.3,8.8];y2 = [9.54,9.54];
h1=subplot(3,1,1);
bar(varexp(1:8),'FaceColor',[0 153 255]/255,'BarWidth',0.6,'LineStyle','none');hold on;
errorbar(varexp(1:8),delt_var(1:8),'k-','linewidth',1.5,'LineStyle','none');hold on;
% plot(x1,y1,'r--','linewidth',0.8); hold on;
% plot(x2,y2,'g--','linewidth',0.8); hold on;
axis ([0.3 8.8 0 30]);grid on;
ylabel ('Fractional variance (%)','FontSize',14);
xlabel ('Modes','FontSize',14);
title ('(a) Fractional variance explained by the first eight modes', 'FontSize',14); 
    set(gca,'FontSize',14);
    set(gca,'ytick',[0:5:30]);
    set(gca,'xtick',[1:1:8]);
    set(gca,'tickdir','out','LineWidth',0.8);  
    
h2=subplot(3,1,2);    
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
% axis ([1 296 0 200]);grid on;
ylabel ('Depth anomaly (m)', 'FontSize',14);
xlabel ('Longtitude','FontSize',14);
title ('(b) EP El Nino, D20 anomaly and its reconstructions based on modes', 'FontSize',14); 
set(gca,'ytick',[-40:20:40],'FontSize',14);
set(gca,'xtick',[38:60:296],'FontSize',14);
set(gca,'xticklabel',['150E';'180 ';'150W';'120W';' 90W'],'FontSize',14);
    set(gca,'tickdir','out','LineWidth',0.8);  
L1=['D20 anomaly'];L2=['Mode 1'];L3=['Modes 1-2'];L4=['Modes 1-3'];L5=['Modes 1-4'];L6=['Modes 1-5'];L7=['Modes 1-6'];
hl=legend (L1,L2,L3,L4,L5,L6,L7,3);
    set(hl,'box','off','fontsize',10);%,'FontWeight','bold'); 
set(gca,'yDir','reverse');    
    
h3=subplot(3,1,3);    
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
    set(hl,'box','off','fontsize',10);%,'FontWeight','bold');
% axis ([1 296 -90 90]);grid on;
axis ([1 296 -30 30]);grid on;
ylabel ('Depth anomaly (m)', 'FontSize',14);
xlabel ('Longtitude','FontSize',14);
title ('(c) CP El Nino, D20 anomaly and its reconstructions based on modes', 'FontSize',14); 
set(gca,'ytick',[-30:10:30],'FontSize',14);
% set(gca,'ytick',[0:30:200],'FontSize',14);
set(gca,'xtick',[38:60:296],'FontSize',14);
set(gca,'xticklabel',['150E';'180 ';'150W';'120W';' 90W'],'FontSize',14);
    set(gca,'tickdir','out','LineWidth',0.8);  
set(gca,'yDir','reverse');
    
set(h1,'position',[.18 .77 .6 .2]);
set(h2,'position',[.18 .48 .6 .2]);
set(h3,'position',[.18 .19 .6 .2]);
print -depsc E:\Work\thermocline_md\pic\Modes_Fuc_Pac_std-Comp_mode-decomposition-revison3-yrev-3.eps; 


%--- only D20 anomaly of EP and CP
figure;
    set(gcf,'paperunit','centimeters');
    set(gcf,'paperposition',[0.5 0.5 24.0 28.0])
    
h1=subplot(2,1,1);    
xP = (1:1:size(d0ep,2));
plot(xP,d0ep,'k-','linewidth',1.0); hold on;
plot(xP,d1ep,'r-','linewidth',1.0);hold on;
plot(xP,d2ep,'b-','linewidth',2.0);hold on;
plot(xP,d3ep,'g-','linewidth',2.0);hold on;
plot(xP,d4ep,'m--','linewidth',1.5);hold on;
plot(xP,d5ep,'c--','linewidth',1.5);hold on;
plot(xP,d6ep,'y--','linewidth',1.2);hold on;
plot(xP,d0ep,'k-','linewidth',1.5); hold on;
plot(xP,d1ep-bb,'r-','linewidth',1.0);hold on;
axis ([1 296 -40 40]);grid on;
% axis ([1 296 0 200]);grid on;
ylabel ('Depth anomaly(m)', 'FontSize',14);
xlabel ('Longtitude','FontSize',14);
title ('(a) EP El Nino, D20 anomaly and its reconstructions based on modes', 'FontSize',12); 
set(gca,'ytick',[-40:20:40],'FontSize',14);
set(gca,'xtick',[38:60:296],'FontSize',14);
set(gca,'xticklabel',['150E';'180 ';'150W';'120W';' 90W'],'FontSize',14);
    set(gca,'tickdir','out','LineWidth',0.8);  
L1=['D20 anomaly'];L2=['Mode 1'];L3=['Modes 1-2'];L4=['Modes 1-3'];L5=['Modes 1-4'];L6=['Modes 1-5'];L7=['Modes 1-6'];
hl=legend (L1,L2,L3,L4,L5,L6,L7,3);
    set(hl,'box','off','fontsize',10);%'FontWeight','bold');
set(gca,'yDir','reverse');    
    
h2=subplot(2,1,2);    
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
    set(hl,'box','off','fontsize',10);%'FontWeight','bold');
% axis ([1 296 -90 90]);grid on;
axis ([1 296 -30 30]);grid on;
ylabel ('Depth anomaly(m)', 'FontSize',14);
xlabel ('Longtitude','FontSize',14);
title ('(b) CP El Nino, D20 anomaly and its reconstructions based on modes', 'FontSize',12); 
set(gca,'ytick',[-30:10:30],'FontSize',14);
% set(gca,'ytick',[0:30:200],'FontSize',14);
set(gca,'xtick',[38:60:296],'FontSize',14);
set(gca,'xticklabel',['150E';'180 ';'150W';'120W';' 90W'],'FontSize',14);
    set(gca,'tickdir','out','LineWidth',0.8);  
set(gca,'yDir','reverse');
    
set(h1,'position',[.18 .77 .55 .2]);
set(h2,'position',[.18 .48 .55 .2]);
print -depsc E:\Work\thermocline_md\pic\Composite_D20_anomaly-revison3-yrev.eps; 


