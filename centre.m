function [resultx resulty resultz] = centre(x,lon,lat,depth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[a b c]=size(x);
resultx0=zeros(a,b,c);
resulty0=zeros(a,b,c);
resultz0=zeros(a,b,c);

for k=1:c
for j=1:b
for i=2:(a-1)
resultx0(i,j,k)=(x(i+1,j,k)-x(i-1,j,k))/(lon(i+1)-lon(i-1));
end
end
end
resultx(:,:,:)=resultx0(2:(a-1),2:(b-1),1:(c-1));

for k=1:c
for j=2:(b-1)
for i=1:a
resulty0(i,j,k)=(x(i,j+1,k)-x(i,j-1,k))/(lat(j+1)-lat(j-1));
end
end
end
resulty(:,:,:)=resulty0(2:(a-1),2:(b-1),1:(c-1));

for k=2:(c-1)
for j=1:b
for i=1:a
resultz0(i,j,k)=(x(i,j,k+1)-x(i,j,k-1))/(depth(k+1)-depth(k-1));
end
end
end
resultz0(:,:,1)=(x(:,:,2)-x(:,:,1))/(depth(2)-depth(1));
resultz(:,:,:)=resultz0(2:(a-1),2:(b-1),1:(c-1));

