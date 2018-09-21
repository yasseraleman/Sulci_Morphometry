function [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(ArgFile);
%
% Syntax :
% [NodeIds, SulcNames,NodeNpoint, SulcLabels, TmtktriIds] = Read_SulcalArgFiles(ArgFile);
%
% This function replace the sulcis names, contained in Sulclist, with a new
% sulc name in the Arg File (Brainvisa Sulci organization format)
%
% Input Parameters:
%   InArg             : Input Arg File
%   Sulclist          : List names for sulcis that will be renamed.
%   Nsulcname         : New Sulcis name
%   Outdir            : Output Directory.
%
% Output Parameters:
%   OutArg            : New/Output Arg file with the sulcis names changed
%
%
% Related references:
%
%
% See also: Multi_Replace_Sulc_name_Argfile
%
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% March 20th 2012
% Version $1.0


%%
[pth, nm, ext] = fileparts(ArgFile);

%=====================Checking Input Parameters===========================%
% if nargin<4
%     Outdir = pth;
% end
%=========================================================================%
%==========================  Reading Arg File  ===========================%
fio = fopen(ArgFile);
contid = 0;
contpo = 0;
conttm = 0;
cont = 0;
SulcNames = '';
SulcLabels = '';
while 1
    cont = cont + 1;
    line = fgetl(fio);
    if ~ischar(line),   break,   end
    % Finding Sulcus Node Id
    if ~isempty(strfind(lower(line),'*begin node fold'))
        contid = contid + 1;
        temp = strread(line,'%s');
        NodeIds(contid) = str2num(temp{4});
    end
    % Finding Sulcus name
    if ~isempty(strfind(lower(line),'name'))
        if strcmp(lower(line(1:4)),'name')
            temp = strread(line,'%s');
            sname = temp{2};
            SulcNames = strvcat(SulcNames,sname);
        end
    end
    % Finding Sulcus number of points
    if ~isempty(strfind(lower(line),'point_number'))
        if strcmp(lower(line(1:12)),'point_number')
            contpo = contpo + 1;
            temp = strread(line,'%s');
            NodeNpoint(contpo) = str2num(temp{2});
        end
    end
    % Finding Sulcus label
    if ~isempty(strfind(lower(line),'label'))
        if strcmp(lower(line(1:6)),'label ')
            temp = strread(line,'%s');
            slabel = temp{2};
            SulcLabels = strvcat(SulcLabels,slabel);
        end
    end
    % Finding Sulcus Tmtktri Id
    if ~isempty(strfind(lower(line),'tmtktri_label'))
        if strcmp(lower(line(1:13)),'tmtktri_label')
            conttm = conttm + 1;
            temp = strread(line,'%s');
            TmtktriIds(conttm) = str2num(temp{2});
        end
    end
    
end
fclose(fio);
return