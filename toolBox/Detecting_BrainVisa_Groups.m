function varargout = Detecting_BrainVisa_Groups(varargin);
%
% Syntax :
%  [groupsIds, groupsNames, allsulciIds, labelsId] = Detecting_BrainVisa_Groups(ChangeTxtFile, hemi);
%
% This script creates a stats table using the BrainVisa outputs.
%
% Input Parameters:
%       ChangeTxtFile               : Text File for grouping BrainVISA
%                                     sulci
%         hemi                      : hemisphere ('left' or 'right').
%
% Output Parameters:
%       groupsIds                    : Groups IDs (BrainVISA Nomenclature due 
%                                      to mophometry statistics do not work 
%                                      with custom sulcus Ids).
%       groupsNames                  : Custom Names for each groups.
%       allsulciIds                  : All valid sulci Ids.
%        labelsId                    : Group numerical Ids, all sulcal
%                                      nodes that presents the same labelId
%                                      belongs to the same group
%
% See also:
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% September 13th 2014
% Version $1.0

%% =========================== Input parameters  =========================%
% if nargin <2
%     error('Two inputs are mandatory');
%     return;
% end

ChangeTxtFile = varargin{1};
hemi = varargin{2};
%% ======================= End of Input parameters  ======================%

% ============= Reading Nomenclature text file  ==========================%
hemi = lower(hemi);
hemi(1) = upper(hemi(1));
fio = fopen(ChangeTxtFile,'rt');lines = '';cont = 0;
groupsNames = '';
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    if isempty(strfind(line,'#'))
        lines = strvcat(lines,line);
    else
        temp = strread(line,'%s','delimiter','-');  
        strname = temp{2};
        if ~strcmp(lower(strname),'unknown')
        groupsNames = strvcat(groupsNames,[hemi '-' strname]);
        end
    end
end
fclose(fio);
% ========================================================================%
allsulciIds = '';
groupsIds = '';
labelsId = 0;
for j= 1:size(lines,1)
    sts = strread(deblank(lines(j,:)),'%s','delimiter',' ');
    Sulclist = '';
    Nstruct = size(sts,1);
    for k=1:Nstruct 
        if k==1
            Nsulcname = char(sts{k});
            if strcmp(Nsulcname(1:end-1),'unknown')
                Nsulcname = 'unknown';
            else
                Nsulcname = [Nsulcname(1:end-1) lower(hemi)];
                groupsIds = strvcat(groupsIds, Nsulcname);
                labelsId = [labelsId;j*ones(Nstruct-1,1)];
            end
        else
            Sulclist = strvcat(Sulclist,[char(sts{k}) lower(hemi)]);
            if ~strcmp(Nsulcname(1:end-1),'unknown')
                allsulciIds = strvcat(allsulciIds,[char(sts{k}) lower(hemi)]);
            end
        end
    end
    stsl{j} = strvcat(Nsulcname,Sulclist);
    
end
labelsId(1) = [];
labels = [ones(19,1);2*ones(11,1);3*ones(12,1);4*ones(7,1);5*ones(2,1);6*ones(7,1);7*ones(4,1);8;9];
if nargout ==1
    varargout{1} = groupsIds;
end
if nargout ==2
    varargout{1} = groupsIds;
    varargout{2} = groupsNames;
end
if nargout ==4
    varargout{1} = groupsIds;
    varargout{2} = groupsNames;
    varargout{3} = allsulciIds;
end
if nargout ==4
    varargout{1} = groupsIds;
    varargout{2} = groupsNames;
    varargout{3} = allsulciIds;
    varargout{4} = labelsId;
end
return;