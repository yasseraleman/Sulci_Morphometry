function varargout = Intercept_Plane_with_LinesSegment(varargin);
%
% Syntax :
%      varargout = Intercept_Plane_with_LinesSegment(planEq, P1, P2);
%
% This script computes the interception between planes and different line segments.
%
% Input Parameters:
%       planEq                  : Sulci Surface in matlab format
%       P1 and P2               : Line segments points
%
% Output Parameters:
%       intPoints               : Interception points
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

if nargin <3
    error('Wrong inputs');
    return;
end

planEq = varargin{1};
P1 = varargin{2};
P2 = varargin{3};

n1 = planEq(:,1:3);
D1 = planEq(:,4);

cont = 0;
N = size(P1,1);
intPoints = [0 0 0 0 0];
%% =============================== Main Program ==========================%
for i = 1:size(n1,1) % For each plane

     Num = dot(repmat(n1(i,:),[N 1]),P1,2)+repmat(D1(i,:),[N 1]);
     Den = dot(repmat(n1(i,:),[N 1]),P2-P1,2);
     t = -1*Num./Den;clear Num Den;
     intersc = P1+repmat(t,[1 3]).*(P2-P1);
     dots = dot(P1-intersc,intersc-P2,2)./(sqrt(dot(P1-intersc,P1-intersc,2)).*sqrt(dot(intersc-P2,intersc-P2,2))); 
     ind = find(sign((round(dots*10))/10) ==1);
     if ~isempty(ind)
        intPoints = [intPoints;intersc(ind,:) ind i*ones(length(ind),1)];
     end
     clear t;
end 
intPoints(1,:) = [];
%% ====================== Identifying sulci boundaries ===================%
varargout{1} = intPoints;
return;