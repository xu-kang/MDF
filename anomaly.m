function [a,y]=anomaly(x,f);
% ANOMALY       constructs a series of anomalies from a timeseries array
%
% [a,y] = anomaly(x,f)
%
% x     = vector or array with columns of time series. Array with any
%	  number of dimensions can be input, as long as the first
%	  dimension (along columns) is time. All columns are treated
%	  separately.  
% f     = the number of samples in one empirical cycle
%	  (default = 12 for seasonal cycle in monthly data).
%
% a	= the anomalies in an array of same size as input x (M-by-N-by-O-...)
% y	= the "normal" period (year) in a similar array but with size 
%         (f-by-N-by-O-...)
%
% EXAMPLES: 
% Input f=1 results in output of the perturbations from the
% overall mean, ie.	    _
%			u = u + u'   <=>   [u_pert,u_mean]=ANOMALY(u,1)
% 
% Removing the empirical seasonal cycle from a timeseries x of monthly data
% can be done by
%			a = ANOMALY(x)
%
% An array of timeseries from different locations must have time running
% along the 1st dimension (columns). It might have longitude along rows and
% latitude inward in the third dimension, and maybe even vertical levels in
% the 4th dimension. In any case, ANOMALY will treat each spatial point
% separately, making it's calculations only along the time-dimension.
%
% INVERSE TRANSFORMATION: 
% After creating the anomalies (a) and the normal period (y), there is
% really no need to let the original dataset take up space. Backtracking can
% be done by input of a and y, ie.
% 
%			x = ANOMALY(a,y)
%
% See also MEAN NANMEAN

%Time-stamp:<Last updated on 02/05/21 at 10:59:29 by even@gfi.uib.no>
%File:<d:/home/matlab/anomaly.m>

error(nargchk(1,2,nargin));
if nargin<2 | isempty(f),	f=12;	end
if isvector(x)==2,		x=x(:);	end

D=size(x);			% Basic size
M=D(1); D=D(2:end);		% Separate length of first dim 
                                % from the rest 
				
if issingle(f)
  Mf=floor(M/f);		% Whole number of jumps (years)
  y=cutmat(x,1:Mf*f,1);		% Use only part with whole "years"
  y=reshape(y,[f,Mf,D]);	% Put years beside eachother and push the
				% other dimensions one right
  y=mean(y,2);			% Find the mean along rows
  y=reshape(y,[f,D]);		% Shape back to "original" form with 
				% N and O as 2nd and 3rd dimensions
  Y=repmat(y,[Mf+1,ones(1,length(D))]);	% Stack years to fit
  Y=cutmat(Y,1:M,1);			% data-matrix for calculations
  a=x-Y;			% Anomalies = remove average
else
  y=f;				% INVERSE:
  d=size(y);
  f=d(1); d=d(2:end);
  if d~=D, error('The length of higher dimensions must match!'); end
  Mf=floor(M/f);			%
  Y=repmat(y,[Mf+1,ones(1,length(D))]);	%
  Y=cutmat(Y,1:M,1);			%
  a=Y+x; % x=Y+a (output is x)
end
