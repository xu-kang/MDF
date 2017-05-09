function w = omega(u,v,lon,lat,depth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[a b c]=size(u);
ux0=zeros(a,b,c);
vy0=zeros(a,b,c);
ux=zeros(a-2,b-2,c-1);
vy=zeros(a-2,b-2,c-1);
aa=zeros(a-2,b-2,c);

for k=1:c
for j=1:b
for i=2:(a-1)
ux0(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(lon(i+1)-lon(i-1));
end
end
end

ux=ux0(2:(a-1),2:(b-1),1:(c-1));
 clear ux0

for k=1:c
for j=2:(b-1)
for i=1:a
vy0(i,j,k)=(v(i,j+1,k)-v(i,j-1,k))/(lat(j+1)-lat(j-1));
end
end
end

vy=vy0(2:(a-1),2:(b-1),1:(c-1));
 clear vx0
aa(:,:,2:c)=-(ux+vy);
aa(:,:,1)=0;
% aa=mean(aa,3);
% aa=squeeze(aa);
w(:,:,1)=aa(:,:,2)*(depth(1));
w(:,:,2)=aa(:,:,2)*(depth(2));
for j=2:(c-2)
w(:,:,j+1)=aa(:,:,j+1)*(depth(j+1)-depth(j-1))+w(:,:,j-1);
end
