function [resultx resulty] = centre(x,lon,lat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[a b]=size(x);
resultx0=zeros(a,b);
resulty0=zeros(a,b);

for j=1:b
for i=2:(a-1)
resultx0(i,j)=(x(i+1,j)-x(i-1,j))/(lon(i+1)-lon(i-1));
end
end
resultx(:,:)=resultx0(2:(a-1),2:(b-1));

for j=2:(b-1)
for i=1:a
resulty0(i,j)=(x(i,j+1)-x(i,j-1))/(lat(j+1)-lat(j-1));
end
end
resulty(:,:)=resulty0(2:(a-1),2:(b-1));

