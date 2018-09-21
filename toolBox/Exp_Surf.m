function [OutFiles, Surfa] = Exp_Surf(SurfF, IMFile, Out_dir, form, ei, sa, join);
%
% Syntax :
% [OutFiles, Surfa] = Exp_Surf(SurfF, IMFile, form, ei, sa);
%
% Export/Import diferent surfaces format files from/to matlab
%
% Input Parameters:
%   SurfF       : Surfaces filenames.
%   IMFile      : Images Filenames(just in case of importing  to matlab format).
%   Out_dir     : Output Directory .
%   form        : Output formats( ie: mat(matlab format),txt or srx (Imagic format),
%                 off(FSL Betall format),obj(MNI surface format),dfs(Brainsuite
%                  surface format), mesh(Brainvisa format),asc(FreeSurfer ASCII format)
%                  ,none(No export or import process)).
%   ei          : Variable indicating export(exp) or import(imp) surface.
%   sa          : Variable to save or not the matlab surfaces.
%   join        : If join varible is empty the surface centers will be corrected.
%                 If join varible is 'joining' do not correct surface centers.
%
% Output Parameters:
%  OutFiles     : List with the exported/imported surfaces.
%  Surfa        : Cell Array with surfaces in matlab variables format.
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 30th 2006
% Version $1.0

warning off;
%=====================Checking Input Parameters===========================%
if nargin==0
    [SurfF,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
    ei = input('Do you want to export or import surfaces(exp/imp):   ','s');
    if (~strcmp(ei,'exp')&~strcmp(ei,'imp'))
        errordlg('Please enter a correct format( ie: exp/imp) ');
        return;
    end
    if strcmp(lower(ei),'exp')
        form = input('Please enter the format( ie: mat,srx,txt,off,obj,dfs,mesh,vtk,frb,non):   ','s');
        if (strcmp(lower(form),'mat')&strcmp(lower(form),'srx')&strcmp(lower(form),'txt')&strcmp(lower(form),'off')&strcmp(lower(form),'obj')&strcmp(lower(form),'dfs')&strcmp(lower(form),'mesh')&strcmp(lower(form),'asc')&strcmp(lower(form),'vtk')&strcmp(lower(form),'non'))
            errordlg('Please enter a correct format( ie: mat,srx,off,obj,dfs,mesh,asc,frb,vtk,non) ');
            return;
        end
        sa = input('Do you want to save the surfaces in a new folder(y/n):   ','s');
        if strcmp(lower(sa),'y')&strcmp(lower(sa),'n')
            errordlg('Please enter a correct answer( ie: y/n) ');
            return;
        end
        if strcmp(sa,'y')
            [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
        elseif strcmp(sa,'n')
            Out_dir = '';
        end
    else
        sa = input('Do you want to save the surfaces(y/n):   ','s');
        if strcmp(lower(sa),'y')&strcmp(lower(sa),'n')
            errordlg('Please enter a correct answer( ie: y/n) ');
            return;
        end
        if strcmp(sa,'y')
            [Out_dir,sts] = spm_select([1],'dir','Selecting Output Directory','',cd);
        end
    end
elseif nargin < 7
    join = '';
end
if ~exist('SurfF','var')|isempty(SurfF)
    [SurfF,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
end
if ~exist('ei','var')|isempty(ei)
    ei = input('Do you want to export or import surfaces(exp/imp):   ','s');
    if strcmp(lower(ei),'exp')&strcmp(lower(ei),'imp')
        errordlg('Please enter a correct format( ie: exp/imp) ');
        return;
    end
    if strcmp(lower(ei),'exp')
        form = input('Please enter the format( ie: mat,srx,txt,off,obj,dfs,mesh,vtk,frb,non):   ','s');
        if (strcmp(lower(form),'mat')&strcmp(lower(form),'srx')&strcmp(lower(form),'txt')&strcmp(lower(form),'off')&strcmp(lower(form),'obj')&strcmp(lower(form),'dfs')&strcmp(lower(form),'mesh')&strcmp(lower(form),'asc')&strcmp(lower(form),'vtk')&strcmp(lower(form),'non'))
            errordlg('Please enter a correct format( ie: mat,srx,off,obj,dfs,mesh,asc,frb,vtk,non) ');
            return;
        end
    end
end
if (~exist('sa','var'))|(isempty(sa))
    sa = input('Do you want to save the surfaces(y/n):   ','s');
    if strcmp(lower(sa),'y')&strcmp(lower(sa),'n')
        errordlg('Please enter a correct answer( ie: y/n) ');
        return;
    end
end
%=========================================================================%
%=========================Main program====================================%
sa = lower(sa);
wh = whos('SurfF');
if (strcmp(wh.class,'struct'))|(strcmp(wh.class,'cell'))
    Ns = length(SurfF);
elseif ischar(SurfF(1,:));
    Ns = size(SurfF,1);
end
OutFiles = '';
if ~isempty(IMFile)
    if strcmp(IMFile,'0')
        IMFile = repmat(IMFile,[Ns 1]);
    end
end
for i = 1:Ns
    if strcmp(wh.class,'struct')
        Surf = SurfF(1,i);
        if (~exist('Out_dir','var'))|(isempty(Out_dir))
            pth = cd;
        else
            pth = Out_dir(i,:);
        end
    elseif strcmp(wh.class,'cell')
        Surf = SurfF{i,1};
        if (~exist('Out_dir','var'))|(isempty(Out_dir))
            pth = cd;
        else
            pth = Out_dir(i,:);
        end
    elseif ischar(SurfF(1,:));
        [pthm,nmm,extm] = fileparts(SurfF(i,:));
        if strcmp(extm(1:4),'.mat')
            Surf = load('-mat',[pthm filesep nmm extm(1:4)]);
            try Surf = Surf.Surf;
            catch
                warndlg([pthm filesep nmm extm(1:4) '  is does not contain surface information '])
            end
        end
    else
        disp('Please select a correct file');
        return;
    end
    switch lower(ei)
        case 'exp'
            Surfa = '';
            if ~exist('Out_dir','var')|isempty(Out_dir)
                pth =pthm;
            else
                pth = Out_dir(i,:);
            end
            if ~isfield(Surf,'SurfData');
                errordlg('This file does not contain surface information');
                return;
            end
            for j = 1:length(Surf)
                %                 switch
                if ~isempty(strfind(form,'mat'))
                    name = Surf(j).Name;
                    if strcmp(sa,'y')
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        temp = strfind(name,'Surf_');
                        if isempty(temp)
                            name = ['Surf_' name];
                        end
                        [Out] = savematfile([npth filesep name '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep name '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                elseif ~isempty(strfind(form,'srx'))
                    name = Surf(j).Name;
                    vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    face = Surf(j).SurfData.faces;
                    Nfaces = size(Surf(j).SurfData.faces,1);
                    Npoints = size(Surf(j).SurfData.vertices,1);
                    mkdir(pth,'SRX_Surfaces'); npth = [pth filesep 'SRX_Surfaces'];
                    dim = Surf(j).Dim;
                    if sum(Surf(j).Orig)==0
                        vert = Surf(j).SurfData.vertices+repmat([91 127 73],[size(Surf(j).SurfData.vertices,1),1]);
                    end
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = ['Surf_' name];
                    end
                    TxtFile = [npth filesep name '.srx'];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    AT = inv([1 0 0 1; 0 0 -1 dim(2); 0 1 0 1; 0 0 0 1]);
                    xyzn = [vert ones(size(vert,1),1)]*AT'; xyzn(:,4) = [];
                    fid = fopen(TxtFile, 'w');
                    fwrite(fid, size(xyzn, 1), 'int32');
                    fwrite(fid, size(face, 1), 'int32');
                    if isempty(join)
                        temp = xyzn'; temp = temp(:);
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        temp = Surf.SurfData.vertices'; temp = temp(:);
                    else
                        temp = xyzn'; temp = temp(:);
                    end
                    fwrite(fid, temp, 'float32');
                    temp = face'-1; temp = temp(:);
                    fwrite(fid, temp, 'int32');
                    fclose(fid);
                    Surfa = '';
                elseif ~isempty(strfind(form,'txt'))
                    name = Surf(j).Name;
                    vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    face = Surf(j).SurfData.faces;
                    Nfaces = size(Surf(j).SurfData.faces,1);
                    tri = Surf(j).Tri;
                    Npoints = size(Surf(j).SurfData.vertices,1);
                    mkdir(pth,'TXT_Surfaces'); npth = [pth filesep 'TXT_Surfaces'];
                    dim = Surf(j).Dim;
                    if sum(Surf(j).Orig)==0
                        vert = Surf(j).SurfData.vertices+repmat([91 127 73],[size(Surf(j).SurfData.vertices,1),1]);
                    end
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = ['Surf_' name];
                    end
                    TxtFile = [npth filesep name '.txt'];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    AT = inv([1 0 0 1; 0 0 -1 dim(2); 0 1 0 1; 0 0 0 1]);
                    xyzn = [vert ones(size(vert,1),1)]*AT'; xyzn(:,4) = [];
                    if isempty(join)
                        xyzn = [(1:size(xyzn,1))' xyzn];
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        xyzn = [(1:size(xyzn,1))' Surf(j).SurfData.vertices];
                    else
                        xyzn = [(1:size(xyzn,1))' xyzn];
                    end
                    face = [(1:size(face,1))' face-1];
                    fid = fopen(TxtFile,'w');
                    fprintf(fid,'%8.2f\n',size(xyzn,1));
                    fprintf(fid,'%8.2f %8.2f %8.2f %8.2f\n',xyzn');
                    fprintf(fid,'%8.2f\n',size(face,1));
                    fprintf(fid,'%8.2f %8.2f %8.2f %8.2f\n',face');
                    fprintf(fid,'%8.2f\n',size(xyzn,1));
                    for z = 1:size(tri,1)
                        ind = find(tri(z,:)~=0);
                        temp =tri(z,ind);
                        if size(temp,2)~=1
                            s = '%8.2f '; s = repmat(s,[1,temp(1,2)+2]);s(end:end+1)='\n';
                        else
                            s = '%8.2f\n';
                        end
                        fprintf(fid,s,temp);
                    end
                    fclose(fid);
                    Surfa = '';
                elseif ~isempty(strfind(form,'off'))
                    name = Surf(j).Name;
                    if isempty(join)
                        vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        vert = Surf(j).SurfData.vertices;
                    else
                        vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    end
                    face = Surf(j).SurfData.faces;
                    Nfaces = size(Surf(j).SurfData.faces,1);
                    Npoints = size(Surf(j).SurfData.vertices,1);
                    mkdir(pth,'OFF_Surfaces'); npth = [pth filesep 'OFF_Surfaces'];
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = ['Surf_' name];
                    end
                    TxtFile = [npth filesep name '.off'];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    xyzn = [Npoints Nfaces 0; vert];
                    faces = [3*ones(Nfaces,1) face-1];
                    fid = fopen(TxtFile,'w');
                    fprintf(fid,'%s\n','OFF');
                    fprintf(fid,'%f %f %f\n',xyzn');
                    fprintf(fid,'%u %u %u %u\n',faces');
                    fclose(fid);
                    Surfa = '';
                elseif ~isempty(strfind(form,'obj'))
                    name = Surf(j).Name;
                    vert = Surf(j).SurfData.vertices;
                    face = Surf(j).SurfData.faces;
                    if isfield(Surf(j).SurfData,'VertexNormals')
                        normals = Surf(j).SurfData.VertexNormals;
                    elseif (~isfield(Surf(j).SurfData,'VertexNormals'))&((strcmp(form,'dfs')|strcmp(form,'obj')|strcmp(form,'mesh')))
                        fv.vertices = vert;
                        fv.faces = face;
                        h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
                        norma = sqrt(sum((normals').^2));
                        normals = normals./repmat(norma',[1 3]);
                    end
                    Nfaces = size(Surf(j).SurfData.faces,1);
                    Npoints = size(Surf(j).SurfData.vertices,1);
                    mkdir(pth,'OBJ_Surfaces'); npth = [pth filesep 'OBJ_Surfaces'];
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = ['Surf_' name];
                    end
                    TxtFile = [npth filesep name '.obj'];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    fid = fopen(TxtFile,'wt');
                    fprintf(fid,'%s %1.1f %1.1f %1.1f %u %1.1f %u\n','P',0.3,0.3,0.4,10,1,Npoints);
                    fprintf(fid,' %6.4f %6.4f %6.4f\n',vert');
                    fprintf(fid,'\n');
                    fprintf(fid,' %1.6f %1.6f %1.6f\n',normals');
                    fprintf(fid,'\n');
                    fprintf(fid,' %u\n',Nfaces);
                    if isfield(Surf,'Is');
                        [Colors] = Surf_Color(Surf,'jet');
                        fprintf(fid,' %u %1.1f %1.1f %1.1f %u\n',[2 Colors(1,:) 1]);
                        fprintf(fid,' %1.1f %1.1f %1.1f %u\n',[Colors(2:end,:) ones(Npoints-1,1)]');
                        fprintf(fid,'\n');
                    else
                        fprintf(fid,' %u %u %u %u %u\n',[0 1 1 1 1]);
                        fprintf(fid,'\n');
                    end
                    temp = [3:3:3*Nfaces]';
                    temps = repmat(Nfaces,[1 8]);modp = mod(temps,[2:9]);
                    ind = find(modp==0);
                    if ~isempty(ind)
                        Nl = max(ind)+1;
                    else
                        pr = primes(Nfaces);
                        temps = repmat(Nfaces,[1 size(pr,1)]);
                        modp = mod(temps,pr);
                        ind = find(modp==0);Nl = pr(min(ind));
                    end
                    tempt = reshape(temp,[Nfaces/Nl Nl]);
                    cad = [repmat(' %u',[1 Nl-1])  ' %u\n'];
                    fprintf(fid,cad,tempt);
                    fprintf(fid,'\n');
                    tempf = face';tempf = tempf(:)-1;
                    temps = repmat(Nfaces*3,[1 8]);modp = mod(temps,[2:9]);
                    ind = find(modp==0);
                    if ~isempty(ind)
                        Nl = max(ind)+1;
                    else
                        pr = primes(Nfaces*3);
                        temps = repmat(Nfaces*3,[1 size(pr,1)]);
                        modp = mod(temps,pr);
                        ind = find(modp==0);Nl = pr(min(ind));
                    end
                    tempf = reshape(tempf,[3*Nfaces/Nl Nl]);
                    cad = [repmat(' %u',[1 Nl-1])  ' %u'];
                    fprintf(fid,cad,tempf);
                    fclose(fid);
                    Surfa = '';
                elseif ~isempty(strfind(form,'dfs'))
                    name = Surf(j).Name;
                    if isempty(join)
                        vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        vert = Surf(j).SurfData.vertices;
                    else
                        vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    end
                    face = Surf(j).SurfData.faces;
                    if isfield(Surf(j).SurfData,'VertexNormals')
                        normals = Surf(j).SurfData.VertexNormals;
                    elseif (~isfield(Surf(j).SurfData,'VertexNormals'))&((strcmp(form,'dfs')|strcmp(form,'obj')|strcmp(form,'mesh')))
                        fv.vertices = vert;
                        fv.faces = face;
                        h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
                        norma = sqrt(sum((normals').^2));
                        normals = normals./repmat(norma',[1 3]);
                    end
                    Nfaces = size(Surf(j).SurfData.faces,1);
                    Npoints = size(Surf(j).SurfData.vertices,1);
                    mkdir(pth,'DFS_Surfaces'); npth = [pth filesep 'DFS_Surfaces'];
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = ['Surf_' name];
                    end
                    TxtFile = [npth filesep name '.dfs'];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    fid=fopen(TxtFile,'wb','ieee-le');
                    magic = ['D' 'U' 'F' 'F' 'S' 'U' 'R' 'F']';
                    version = [1 0 0 0];
                    hdrsize = 184;
                    mdoffset = 0;
                    pdoffset = 0;
                    nTriangles = size(face,1);
                    nVertices  = size(vert,1);
                    nStrips = 0;
                    stripSize = 0;
                    normals = 0;
                    uvStart = 0;
                    vcoffset = 0;
                    precision = 0;
                    pad=[0 0 0];
                    orientation=eye(4);
                    fwrite(fid,magic,'char');
                    fwrite(fid,version,'char');
                    fwrite(fid,hdrsize,'int32');
                    fwrite(fid,mdoffset,'int32');
                    fwrite(fid,pdoffset,'int32');
                    fwrite(fid,nTriangles,'int32');
                    fwrite(fid,nVertices,'int32');
                    fwrite(fid,nStrips,'int32');
                    fwrite(fid,stripSize,'int32');
                    fwrite(fid,normals,'int32');
                    fwrite(fid,uvStart,'int32');
                    fwrite(fid,vcoffset,'int32');
                    fwrite(fid,precision,'int32');
                    fwrite(fid,orientation,'float64');
                    fwrite(fid,(face-1)','int32');
                    fwrite(fid,vert','float32');
                    fclose(fid);
                    Surfa = '';
                    
                elseif ~isempty(strfind(form,'mesh'))
                    name = Surf(j).Name;
                    vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    face = Surf(j).SurfData.faces;
                    if isfield(Surf(j).SurfData,'VertexNormals')
                        normals = Surf(j).SurfData.VertexNormals;
                    elseif (~isfield(Surf(j).SurfData,'VertexNormals'))&((strcmp(form,'dfs')|strcmp(form,'obj')|strcmp(form,'mesh')))
                        fv.vertices = vert;
                        fv.faces = face;
                        h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
                        norma = sqrt(sum((normals').^2));
                        normals = normals./repmat(norma',[1 3]);
                    end
                    Nfaces = size(Surf(j).SurfData.faces,1);
                    tri = Surf(j).Tri;
                    Npoints = size(Surf(j).SurfData.vertices,1);
                    if j == 1
                        mkdir(pth,'MESH_Surfaces'); npth = [pth filesep 'MESH_Surfaces'];
                        temp = strfind(name,'Surf_');
                        if isempty(temp)
                            name = deblank(['Surf_' name]);
                        end
                        TxtFile = [npth filesep name '.mesh'];
                        OutFiles = strvcat(OutFiles,TxtFile);
                        fid = fopen(TxtFile,'wb');
                        fwrite(fid, 'binar', 'uchar');
                        fwrite(fid, hex2dec('41424344'), 'uint32');
                        fwrite(fid, 4, 'uint32');
                        fwrite(fid, 'VOID', 'uchar');
                        fwrite(fid, size(Surf(1).SurfData.faces,2), 'uint32');
                        fwrite(fid, size(Surf,2), 'uint32');
                        for k=1:size(Surf,2)
                            Mat = inv([1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]);
                            vert = Surf(k).SurfData.vertices+repmat(Surf(k).Orig,[size(Surf(k).SurfData.vertices,1),1]);
                            face = Surf(k).SurfData.faces;
                            if isfield(Surf(k).SurfData,'VertexNormals')
                                normals = Surf(k).SurfData.VertexNormals;
                            elseif (~isfield(Surf(k).SurfData,'VertexNormals'))
                                fv.vertices = vert;
                                fv.faces = face;
                                h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
                                norma = sqrt(sum((normals').^2));
                                normals = normals./repmat(norma',[1 3]);
                            end
                            if isempty(join)
                                vert = Mat*[vert ones(size(vert,1),1)]'; vert(4,:) = [];vert = vert(:);
                            elseif ~isempty(join)&strcmp(lower(join),'joining')
                                vert = Surf(j).SurfData.vertices'; vert = vert(:);
                            else
                                vert = Mat*[vert ones(size(vert,1),1)]'; vert(4,:) = [];vert = vert(:);
                            end
                            face = face';face = face(:)-1;
                            normals = normals';normals = normals(:);
                            fwrite(fid, k-1, 'uint32');
                            fwrite(fid, size(vert,1)/3, 'uint32');
                            fwrite(fid, vert', 'float32');
                            fwrite(fid, size(normals,1)/3, 'uint32');
                            fwrite(fid, normals', 'float32');
                            fwrite(fid, 0, 'uint32');
                            fwrite(fid, size(face,1)/3, 'uint32');
                            fwrite(fid, face', 'uint32');
                        end
                        fclose(fid);
                        Surfa = '';
                    end
                elseif ~isempty(strfind(form,'asc'))
                    name = Surf(j).Name;
                    vert = Surf(j).SurfData.vertices;
                    face = Surf(j).SurfData.faces;
                    mkdir(pth,'ASC_Surfaces'); npth = [pth filesep 'ASC_Surfaces'];
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = ['Surf_' name];
                    end
                    TxtFile = [npth filesep name '.asc'];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    fid = fopen(TxtFile,'wt');
                    vert =[vert zeros(size(vert,1),1)];
                    face =[face-1 zeros(size(face,1),1)];
                    fprintf(fid,'%s\n',['#!ascii version of ' name]);
                    fprintf(fid,'%u %u\n',size(vert,1),size(face,1));
                    fprintf(fid,'%8.6f %8.6f %8.6f %u\n',vert');
                    fprintf(fid,'%u %u %u %u\n',face');
                    fclose(fid);
                    Surfa = '';
                elseif ~isempty(strfind(form,'frb'))|~isempty(strfind(form,'inflated'))|~isempty(strfind(form,'pial'))|~isempty(strfind(form,'white'))|~isempty(strfind(form,'sphere'))|~isempty(strfind(form,'orig'))
                    %case ('frb')|('inf')|('pia')|('sph')|('whi')
                    name = Surf(j).Name;
                    mkdir(pth,'FRB_Surfaces'); npth = [pth filesep 'FRB_Surfaces'];
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = deblank(['Surf_' name]);
                    end
                    TxtFile = [npth filesep name '.' deblank(form)];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    if size(Surf(j).SurfData.faces,2) == 3
                        MNumber = 16777214;
                    elseif size(Surf(j).SurfData.faces,2) == 4
                        MNumber  = 16777215;
                        
                    end
                    fid = fopen(TxtFile, 'wb', 'b') ;
                    fwrite(fid, bitand(bitshift(MNumber, -16), 255), 'uchar') ;
                    fwrite(fid, bitand(bitshift(MNumber, -8), 255), 'uchar') ;
                    fwrite(fid, bitand(MNumber, 255), 'uchar') ;
                    fwrite(fid, sprintf('created from IBASPM on %s\n',datestr(now)),'char');
                    fwrite(fid, sprintf('created from IBASPM on %s\n',datestr(now)),'char');
                    Npoints = size(Surf(j).SurfData.vertices,1) ;  % number of vertices
                    Nfaces = size(Surf(j).SurfData.faces,1) ;  % number of faces
                    fwrite(fid, Npoints,'int32');
                    fwrite(fid, Nfaces,'int32');
                    vertcol = reshape(Surf(j).SurfData.vertices',size(Surf(j).SurfData.vertices,1)*size(Surf(j).SurfData.vertices,2),1);
                    fwrite(fid, vertcol,'float32');
                    facecol = reshape(Surf(j).SurfData.faces',size(Surf(j).SurfData.faces,1)*size(Surf(j).SurfData.faces,2),1)-1;
                    fwrite(fid, facecol,'int32');
                    fclose(fid) ;
                    Surfa = '';
                elseif ~isempty(strfind(form,'vtk'))
                    % case 'vtk'
                    name = Surf(j).Name;
                    if isempty(join)
                        Surf(j).SurfData.vertices = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        Surf(j).SurfData.vertices = Surf(j).SurfData.vertices;
                    else
                        Surf(j).SurfData.vertices = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                    end
                    mkdir(pth,'VTK_Surfaces'); npth = [pth filesep 'VTK_Surfaces'];
                    temp = strfind(name,'Surf_');
                    if isempty(temp)
                        name = deblank(['Surf_' name]);
                    end
                    TxtFile = [npth filesep name '.vtk'];
                    OutFiles = strvcat(OutFiles,TxtFile);
                    fid = fopen(TxtFile, 'wt') ;
                    Npoints = size(Surf(j).SurfData.vertices,1) ;  % number of vertices
                    fprintf(fid,'%s\n','# vtk DataFile Version 3.0');
                    fprintf(fid,'%s\n',['created from IBASPM on ' datestr(now)]);
                    fprintf(fid,'%s\n','ASCII');
                    fprintf(fid,'%s\n','DATASET POLYDATA');
                    fprintf(fid,'%s %u %s\n','POINTS',Npoints,'float');
                    fprintf(fid,'%8.6f %8.6f %8.6f\n',Surf(j).SurfData.vertices');
                    Nfaces = size(Surf(j).SurfData.faces,1) ;  % number of faces
                    Pol = size(Surf(j).SurfData.faces,2);
                    face = [repmat(Pol,[Nfaces 1]) Surf(j).SurfData.faces-1];
                    fprintf(fid,'%s %u %u\n','POLYGONS',Nfaces,Nfaces*(Pol+1));
                    fprintf(fid,'%u %u %u %u\n',face');
                    if isfield(Surf(j),'Is')
                        fprintf(fid,'%s %u\n','POINT_DATA',Npoints);
                        fprintf(fid,'%s\n','SCALARS Scalars float');
                        fprintf(fid,'%s\n','LOOKUP_TABLE default');
                        fprintf(fid,'%8.6f\n',Surf(j).Is');
                        fprintf(fid,'%s\n','VECTORS Vectors float');
                        fprintf(fid,'%8.6f %8.6f %8.6f\n',Surf(j).SurfData.VertexNormals');
                    end
                    fclose(fid);
                    Surfa = '';
                end
            end
        case 'imp'
            [pth,nm,ext] = fileparts(deblank(SurfF(i,:)));
            %             if ~isempty(pth)
            %                 Sfile = [pth filesep nm ext(1:4)];
            %             else
            %                 Sfile = [nm ext(1:4)];
            %             end
            %form = lower(ext(2:4));
            Sfile = deblank(SurfF(i,:));
            form = deblank(lower(ext(2:end)));
            switch form
                case 'mat'
                    disp(['The surface ' [nm ext(1:4)] ' is already in matlab format']);
                    if strcmp(sa,'y')
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                case 'srx'
                    Surf.Imp = 'srx';
                    Surf.Name = nm;
                    Surf.Area = 0;
                    fid = fopen(Sfile, 'r');
                    Npoints = fread(fid, 1, 'int32');
                    Nfaces = fread(fid, 1, 'int32');
                    
                    % Reading Vertices
                    xyzn = fread(fid, 3*Npoints, 'float32');
                    
                    % Reading Faces
                    face = fread(fid, 3*Nfaces, 'int32');
                    fclose(fid);
                    vert = reshape(xyzn, [3 Npoints])';xyzn=vert;
                    face = reshape(face, [3 Nfaces])';face = face+1;
                    Surf.SurfData.faces = face;clear face;
                    if ~isempty(IMFile(i,:))&~strcmp(IMFile(i,:),'0')
                        V = spm_vol(IMFile(i,:));
                        Surf.Dim = V.dim(1:3);
                        AT = [1 0 0 1; 0 0 -1 V.dim(2); 0 1 0 1; 0 0 0 1];
                        vert = [vert ones(Npoints,1)]*AT'; vert(:,4) = [];
                        Surf.Orig = abs(V.mat(1:3,4)');
                        Surf.VoxSize = abs(diag(V.mat(1:3,1:3))');
                    else
                        dx = ceil(max(vert(:,1))-min(vert(:,1)))+2;
                        dy = ceil(max(vert(:,2))-min(vert(:,2)))+2;
                        dz = ceil(max(vert(:,3))-min(vert(:,3)))+2;
                        Surf.Orig = [abs(min(vert(:,1)))+ abs(dx/2)+1 abs(min(vert(:,2)))+ abs(dy/2)+1 abs(min(vert(:,3)))+ abs(dz/2)+1];
                        Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                        AT = [1 0 0 1; 0 0 -1 dy; 0 1 0 1; 0 0 0 1];
                        vert = [vert ones(Npoints,1)]*AT'; vert(:,4) = [];
                        Surf.VoxSize = [1 1 1];
                    end
                    if isempty(join)
                        Surf.SurfData.vertices = vert-repmat(Surf.Orig,[Npoints,1]);
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        Surf.SurfData.vertices = xyzn;
                    else
                        Surf.SurfData.vertices = vert-repmat(Surf.Orig,[Npoints,1]);
                    end
                    Surf.SurfData.vertices(:,1) = [];
                    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf.Tri = Tri; clear Tri;
                    h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
                    norma = sqrt(sum((Normals').^2));
                    Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
                    Surf.Type = 'Mask';
                    Surf.Area = Area_Comp(Surf.SurfData);
                    if strcmp(sa,'y')
                        if exist('Out_dir','var')
                            pth = Out_dir;
                        end
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                case 'txt'
                    Surf.Imp = 'txt';
                    Surf.Name = nm;
                    Surf.Area = 0;
                    fid = fopen(Sfile, 'r');
                    
                    % Reading Vertices
                    Npoints = textread(Sfile,'%f',1,'headerlines',0);
                    vert = reshape(textread(Sfile,'%f',Npoints*4,'headerlines',1),4,Npoints)';
                    xyzn=vert;
                    % Reading Faces
                    Nfaces = textread(Sfile,'%f',1,'headerlines',Npoints+1);
                    face = reshape(textread(Sfile,'%f',Nfaces*4,'headerlines',Npoints+2),4,Nfaces)'+1;
                    face(:,1) = [];vert(:,1) = [];
                    Surf.SurfData.faces = face;clear face;
                    if ~isempty(IMFile(i,:))&~strcmp(IMFile(i,:),'0')
                        V = spm_vol(IMFile(i,:));
                        Surf.Dim = V.dim(1:3);
                        AT = [1 0 0 1; 0 0 -1 V.dim(2); 0 1 0 1; 0 0 0 1];
                        vert = [vert ones(Npoints,1)]*AT'; vert(:,4) = [];
                        Surf.Orig = abs(V.mat(1:3,4)');
                        Surf.VoxSize = abs(diag(V.mat(1:3,1:3))');
                    else
                        dx = ceil(max(vert(:,1))-min(vert(:,1)))+2;
                        dy = ceil(max(vert(:,2))-min(vert(:,2)))+2;
                        dz = ceil(max(vert(:,3))-min(vert(:,3)))+2;
                        Surf.Orig = [abs(min(vert(:,1)))+ abs(dx/2)+1 abs(min(vert(:,2)))+ abs(dy/2)+1 abs(min(vert(:,3)))+ abs(dz/2)+1];
                        Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                        AT = [1 0 0 1; 0 0 -1 dy; 0 1 0 1; 0 0 0 1];
                        vert = [vert ones(Npoints,1)]*AT'; vert(:,4) = [];
                        Surf.VoxSize = [1 1 1];
                    end
                    if isempty(join)
                        Surf.SurfData.vertices = vert-repmat(Surf.Orig,[Npoints,1]);
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        Surf.SurfData.vertices = xyzn;
                    else
                        Surf.SurfData.vertices = vert-repmat(Surf.Orig,[Npoints,1]);
                    end
                    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                    Surf.SurfData.vertices(:,1) = [];
                    % Computing Normals
                    h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
                    norma = sqrt(sum((Normals').^2));
                    Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf.Tri = Tri; clear Tri;
                    Surf.Type = 'Mask';
                    Surf.Area = Area_Comp(Surf.SurfData);
                    if strcmp(sa,'y')
                        if exist('Out_dir','var')
                            pth = Out_dir;
                        end
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                case 'off'
                    Surf.Imp = 'off';
                    fid = fopen(deblank(Sfile),'r');
                    lin = fgetl(fid);
                    lin = fgetl(fid);
                    [Npoints,Nfaces,Car] = strread(lin,'%f%f%f','delimiter',' ');
                    
                    % Reading Vertices
                    Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',2),3,Npoints)';
                    xyzn=Surf.SurfData.vertices;
                    % Reading faces
                    Surf.SurfData.faces = uint32(reshape(textread(Sfile,'%f',Nfaces*4,'headerlines',Npoints+2),4,Nfaces)')+1;
                    Surf.SurfData.faces(:,1) = [];
                    if ~isempty(IMFile(i,:))&~strcmp(IMFile(i,:),'0')
                        V = spm_vol(IMFile(i,:));
                        Surf.Dim = V.dim(1:3);
                        Surf.Orig = abs(V.mat(1:3,4)');
                        Surf.VoxSize = abs(diag(V.mat(1:3,1:3))');
                    else
                        dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
                        dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
                        dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
                        Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
                        Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                        Surf.VoxSize = [1 1 1];
                    end
                    if isempty(join)
                        Surf.SurfData.vertices = Surf.SurfData.vertices-repmat(Surf.Orig,[Npoints,1]);
                    elseif ~isempty(join)&strcmp(lower(join),'joining')
                        Surf.SurfData.vertices = xyzn;
                    else
                        Surf.SurfData.vertices = Surf.SurfData.vertices-repmat(Surf.Orig,[Npoints,1]);
                    end
                    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                    
                    % Computing Normals
                    h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
                    norma = sqrt(sum((Normals').^2));
                    Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf.Tri = Tri;
                    clear Tri;
                    Surf.Name = nm;
                    Surf.Area = Area_Comp(Surf.SurfData);
                    Surf.VoxSize = [1 1 1];
                    Surf.Type = 'Mask';
                    if strcmp(sa,'y')
                        if exist('Out_dir','var')
                            pth = Out_dir;
                        end
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                case 'obj'
                    Surf.Imp = 'obj';
                    fid = fopen(deblank(SurfF(i,:)),'r');
                    lin = fgetl(fid);
                    [typ,ac,dif,spr,spc,tr,Npoints] = strread(lin,'%s%n%n%n%u%n%u','delimiter',' ');
                    if lower(char(typ)) == 'p';
                        
                        % Reading Vertices
                        Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',1),3,Npoints)';
                        
                        % Reading Normals
                        Surf.SurfData.VertexNormals = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',Npoints+2),3,Npoints)';
                        Surf.SurfData.FaceColor = 'interp';
                        
                        % Reading Faces
                        Nfaces = textread(Sfile,'%f',1,'headerlines',Npoints*2+3);
                        cl_flag=int32(textread(Sfile,'%f',1,'headerlines',Npoints*2+4));
                        if cl_flag == 0
                            cl=int32(textread(Sfile,'%f',5,'headerlines',Npoints*2+4));
                            end_ind=textread(Sfile,'%f',Nfaces,'headerlines',Npoints*2+6);
                            temps = repmat(Nfaces,[1 8]);modp = mod(temps,[2:9]);
                            ind = find(modp==0);
                            if ~isempty(ind)
                                Nl = max(ind)+1;
                            else
                                pr = primes(Nfaces);
                                temps = repmat(Nfaces,[1 size(pr,1)]);
                                modp = mod(temps,pr);
                                ind = find(modp==0);Nl = pr(max(ind));
                            end
                            leng = max(end_ind);
                            ind=textread(Sfile,'%f',double(leng),'headerlines',Npoints*2+6+Nfaces/Nl);
                            ind = ind+1;
                            Surf.SurfData.FaceColor = [1 1 1];
                        elseif cl_flag == 1
                        elseif cl_flag == 2
                            temp = textread(Sfile,'%f',Npoints*4+1,'headerlines',Npoints*2+4); temp(1) = [];temp(4:4:end) = [];
                            Surf.SurfData.FaceVertexCData = reshape(temp,3,Npoints)';clear temp;
                            end_ind=textread(Sfile,'%f',Nfaces,'headerlines',Npoints*3+5);
                            temps = repmat(Nfaces,[1 8]);modp = mod(temps,[2:9]);
                            ind = find(modp==0);
                            if ~isempty(ind)
                                Nl = max(ind)+1;
                            else
                                pr = primes(Nfaces);
                                temps = repmat(Nfaces,[1 size(pr,1)]);
                                modp = mod(temps,pr);
                                ind = find(modp==0);Nl = pr(min(ind));
                            end
                            leng = max(end_ind);
                            ind=textread(Sfile,'%f',double(leng),'headerlines',Npoints*3+5+Nfaces/Nl);
                            ind = ind+1;
                        end
                        faces =  [ind(end_ind-2) ind(end_ind-1) ind(end_ind)];
                        if ~isempty(IMFile(i,:))&~strcmp(IMFile(i,:),'0')
                            V = spm_vol(IMFile(i,:));
                            Surf.Dim = V.dim(1:3);
                            Surf.Orig = abs(V.mat(1:3,4)');
                            Surf.VoxSize = abs(diag(V.mat(1:3,1:3))');
                        else
                            dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
                            dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
                            dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
                            Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
                            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                            Surf.VoxSize = [1 1 1];
                        end
                        Surf.SurfData.faces = faces; clear faces;
                        clear ind end_ind;
                        [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                        Temp = sum(Tri);
                        Tri(:,Temp==0) = [];
                        Surf.Tri = Tri;
                        clear Tri;
                        Surf.Name = nm;
                        Surf.Area = Area_Comp(Surf.SurfData);
                        Surf.Type = 'Mask';
                        if strcmp(sa,'y')
                            mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                            [Out] = savematfile([npth filesep nm '.mat'], Surf);
                            OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                            %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                            Surfa = '';
                        elseif strcmp(sa,'n')
                            OutFiles = '';
                            Surfa{i,1} = Surf;
                        end
                    end
                case 'dfs'
                    Surf.Imp = 'dfs';
                    Sfile = [pth filesep nm ext(1:4)];
                    fid=fopen(deblank(Sfile),'rb','ieee-le');
                    magic=char(fread(fid,8,'char'));
                    version=fread(fid,4,'char');
                    hdrsize=fread(fid,1,'int32');
                    mdoffset=fread(fid,1,'int32');
                    pdoffset=fread(fid,1,'int32');
                    Nfaces=fread(fid,1,'int32');
                    Npoints=fread(fid,1,'int32');
                    nStrips=fread(fid,1,'int32');
                    stripSize=fread(fid,1,'int32');
                    Normals=fread(fid,1,'int32');
                    uvStart=fread(fid,1,'int32');
                    vcoffset=fread(fid,1,'int32');
                    precision=fread(fid,1,'int32');
                    orientation=fread(fid,[4 4],'float64');
                    fseek(fid,hdrsize,-1);
                    
                    % Reading Faces
                    faces = fread(fid,[3 Nfaces],'int32')+1;
                    
                    % Reading Vertices
                    vert=fread(fid,[3 Npoints],'float32');
                    Surf.SurfData.faces = faces' ;clear faces;
                    Surf.SurfData.vertices = vert' ;xyzn=vert';clear vert;
                    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf.Tri = Tri;
                    Surf.Name = nm;
                    Surf.Area = Area_Comp(Surf.SurfData);
                    Surf.VoxSize = [1 1 1];
                    Surf.Type = 'Mask';
                    
                    % Reading Normals
                    if (Normals>0)
                        fseek(fid,Normals,-1);
                        Surf.SurfData.VertexNormals = fread(fid,[3 Npoints],'float32')';
                    end;
                case 'gii'
                    Sfile = [pth filesep nm ext(1:4)];
                    S = gifti(deblank(Sfile));
                    vertices = S.vertices;
                    faces    = S.faces;
                    mat    = S.mat;
                    if ~iscell(vertices)
                        vertices = {vertices};
                        faces = {faces};
                        mat = {mat};
                    end
                    Tsteps = length(vertices);
                    for t=1:Tsteps
                        Surf(t).Imp = 'gii';
                        Surf(t).SurfData.vertices = double(vertices{t});
                        Surf(t).SurfData.faces = faces{t};
                        if ~isempty(mat{t})
                            Surf(t).Mat = mat{t};
                        end
                        Npoints = size(Surf(t).SurfData.vertices,1);
                        Nfaces = size(Surf(t).SurfData.faces,1);
                        
                        [Tri] = Vert_Neib(double(Surf(t).SurfData.faces),Npoints,Nfaces);
                        Temp = sum(Tri);
                        Tri(:,Temp==0) = [];
                        Surf(t).Tri = Tri;
                        Surf(t).Name = nm;
                        Surf(t).Area = Area_Comp(Surf(t).SurfData);
                    end
                    Surfa{i,1} = Surf;
                case 'mesh'
                    fid = fopen(deblank(Sfile),'rb');
                    Type = char(fread(fid, 5, 'uchar'));  %- 'ascii' or 'binar'
                    if strcmp(Type','binar');
                        [byte_swapping, COUNT]     = fread(fid, 1, 'uint32'); %- 'ABCD' or 'DCBA'
                        ff = strcmp(dec2hex(byte_swapping),'41424344');
                        if ~ff
                            [fn, pm, mf] = fopen(1); %- machine format
                            fclose(fid);
                            if strmatch(mf,'ieee-le');
                                fid = fopen(deblank(Sfile),'r','ieee-be');
                            else
                                fid = fopen(deblank(Sfile),'r','ieee-le');
                            end
                            [file_format, COUNT]   = fread(fid, 5, 'uchar');
                            [byte_swapping, COUNT] = fread(fid, 1, 'uint32');
                        end
                        Val = fread(fid,1,'uint32');
                        Text = char(fread(fid,4,'char'));
                        Pol = fread(fid,1,'uint32');
                        Tsteps = fread(fid,1,'uint32');
                        for t=1:Tsteps
                            Tinst = fread(fid,1,'uint32');
                            
                            % Reading vertices
                            Npoints = fread(fid,1,'uint32');
                            vert = fread(fid,Npoints*3,'float32')';
                            vert = reshape(vert,[3,Npoints])';xyzn=vert;
                            Surf(t).SurfData.vertices =vert;
                            % Reading normals
                            T = fread(fid,1,'uint32');
                            normals = fread(fid,Npoints*3,'float32');
                            normals = reshape(normals,[3,Npoints])';
                            
                            % Reading Faces
                            T = fread(fid,1,'uint32');
                            Nfaces = fread(fid,1,'uint32');
                            faces = fread(fid,Nfaces*Pol,'uint32');
                            faces = reshape(faces,[Pol,Nfaces])'+1;
                            
                            Surf(t).SurfData.VertexNormals = normals;clear normals
                            Surf(t).SurfData.faces = faces;clear faces vert;
                            if size(Surf(t).SurfData.faces,2) ==3
                                [Tri] = Vert_Neib(double(Surf(t).SurfData.faces),Npoints,Nfaces);
                                Temp = sum(Tri);
                                Tri(:,Temp==0) = [];
                                Surf(t).Tri = Tri;
                            end
                            if Tsteps~=1
                                Surf(t).Name = [nm '_' sprintf('%.3d',t)];
                            else
                                Surf(t).Name = nm;
                            end
                            Surf(t).Area = 0;
                            Surf(t).Type = 'Mask';
                        end
                    elseif strcmp(Type(1:5,1)','ascii');
                        Text = fscanf(fid,'%s',1);
                        Pol = fscanf(fid,'%d',1);
                        Tsteps = fscanf(fid,'%d',1);
                        for t=1:Tsteps
                            ms = fscanf(fid,'\n%d',1);
                            
                            % Reading vertices
                            Npoints = fscanf(fid,'\n%d\n',1);
                            vert = fscanf(fid,'(%f ,%f ,%f) ',3*Npoints);
                            vert = reshape(vert,[3,Npoints])';
                            
                            % Reading Normals
                            if Pol==3
                                T = fscanf(fid,'\n%d\n',1);
                                normals = fscanf(fid,'(%f ,%f ,%f) ',3*Npoints);
                                normals = reshape(normals,[3,Npoints])';
                                T = fscanf(fid,'\n%d\n',1);
                            end
                            
                            % Reading Faces
                            Nfaces = fscanf(fid,'\n%d\n',1);
                            faces = fscanf(fid,'(%d ,%d ,%d) ',Pol*Nfaces);
                            faces = reshape(faces,[Pol,Nfaces])'+1;
                            Surf(t).Name = nm;
                            Surf(t).SurfData.vertices = vert;
                            Surf(t).SurfData.VertexNormals = normals;clear normals
                            Surf(t).SurfData.faces = faces;clear faces vert;
                            if size(Surf(t).SurfData.faces,2) ==3
                                [Tri] = Vert_Neib(double(Surf(t).SurfData.faces),Npoints,Nfaces);
                                Temp = sum(Tri);
                                Tri(:,Temp==0) = [];
                                Surf(t).Tri = Tri;
                            end
                            Surf(t).Area = 0;
                            Surf(t).Type = 'Mask';
                        end
                    end
                    fclose(fid);
                    if strcmp(sa,'y')
                        if exist('Out_dir','var')
                            pth = Out_dir;
                        end
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                case 'asc'
                    Surf.Name = nm;
                    fid = fopen(Sfile, 'rt');
                    line = fgetl(fid);
                    line = fgetl(fid);
                    [Npoints,Nfaces] = strread(line,'%u%u','delimiter',' ');
                    
                    % Reading Vertices
                    Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*4,'headerlines',2),4,Npoints)';
                    Surf.SurfData.vertices(:,4) = [];
                    
                    % Reading Faces
                    Surf.SurfData.faces = reshape(textread(Sfile,'%u',Nfaces*4,'headerlines',2+Npoints),4,Nfaces)'+1;
                    Surf.SurfData.faces(:,4) = [];
                    fclose(fid);
                    if ~isempty(IMFile(i,:))&~strcmp(IMFile(i,:),'0')
                        V = spm_vol(IMFile(i,:));
                        Surf.Dim = V.dim(1:3);
                        Surf.Orig = abs(V.mat(1:3,4)');
                        Surf.VoxSize = abs(diag(V.mat(1:3,1:3))');
                    else
                        dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
                        dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
                        dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
                        Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
                        Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                        Surf.VoxSize = [1 1 1];
                    end
                    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf.Tri = Tri; clear Tri;
                    h = figure;h1 = patch(Surf.SurfData);  Normals = get(h1,'VertexNormals');close(h);
                    norma = sqrt(sum((Normals').^2));
                    Surf.SurfData.VertexNormals = Normals./repmat(norma',[1 3]);
                    Surf.Area = Area_Comp(Surf.SurfData);
                    Surf.Imp = 'asc';
                    Surf.Type = 'Mask';
                    if strcmp(sa,'y')
                        if exist('Out_dir','var')
                            pth = Out_dir;
                        end
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                case 'meta'
                    [str,div,num] = textread(Sfile,'%s%s%s',1,'headerlines',15);
                    Npoints = cell2mat(num);
                    Npoints=str2num(Npoints);
                    Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*4,'headerlines',17),4,Npoints)';
                    Surf.SurfData.vertices(:,1) =[];
                    [str,div,num] = textread(Sfile,'%s%s%s',1,'headerlines',18+Npoints);
                    Nfaces = cell2mat(num);
                    Nfaces=str2num(Nfaces);
                    Surf.SurfData.faces = reshape(textread(Sfile,'%f',Nfaces*4,'headerlines',20+Npoints),4,Nfaces)';
                    Surf.SurfData.faces(:,1) =[];
                    Surf.SurfData.faces = Surf.SurfData.faces+1;
                    [pth,nm,ext] = fileparts(Sfile);
                    Surf.Name = nm;
                    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf.Tri = Tri; clear Tri;
                    Surf.Area = Area_Comp(Surf.SurfData);
                    Surf.Imp = 'meta';
                    if strcmp(sa,'y')
                        if exist('Out_dir','var')
                            pth = Out_dir;
                        end
                        mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                        
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                case 'vtk'
                    Surf.Name = nm;
                    fid = fopen(Sfile, 'rt');
                    line = fgetl(fid);line = fgetl(fid);line = fgetl(fid);
                    line = fgetl(fid);
                    [Info,typ] = strread(line,'%s%s','delimiter',' ');
                    if ~strcmp(lower(typ),'polydata')
                        errordlg('Please select a correct surface format');
                        return;
                    end
                    line = fgetl(fid);
                    
                    % Reading Vertices
                    [txt,Npoints,typ] = strread(line,'%s%n%s','delimiter',' ');
                    Surf.SurfData.vertices = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',5),3,Npoints)';
                    xyzn =Surf.SurfData.vertices;
                                        
                    % Reading Faces
                    [line,Nfaces] = textread(Sfile,'%s %u',1,'headerlines',5+Npoints);
                    pol = textread(Sfile,'%u',1,'headerlines',6+Npoints);pol = pol+1;
                    Surf.SurfData.faces = reshape(textread(Sfile,'%f',Nfaces*pol,'headerlines',6+Npoints),pol,Nfaces)';
                    Surf.SurfData.faces(:,1)=[]; Surf.SurfData.faces= Surf.SurfData.faces+1;
                    tempChar = textread(Sfile,'%s',1,'headerlines',7+Npoints+Nfaces);
                    switch lower(char(tempChar))
                        case 'scalars'
                            Surf.Is = textread(Sfile,'%f',Npoints,'headerlines',9+Npoints+Nfaces);
                        case 'color_scalars'
                            Surf.SurfData.FaceVertexCData = reshape(textread(Sfile,'%f',Npoints*4,'headerlines',8+Npoints+Nfaces),4,Npoints)';
                        otherwise
                            Surf.SurfData.VertexNormals = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',10+2*Npoints+Nfaces),3,Npoints)';
                            
                    end
                    %                     if ~isempty(Surf.Is)
                    %                         Surf.SurfData.VertexNormals = reshape(textread(Sfile,'%f',Npoints*3,'headerlines',10+2*Npoints+Nfaces),3,Npoints)';
                    %                     else
                    %                         Surf =rmfield(Surf,'Is');
                    %                     end
                    [Tri] = Vert_Neib(double(Surf.SurfData.faces),Npoints,Nfaces);
                    Temp = sum(Tri);
                    Tri(:,Temp==0) = [];
                    Surf.Tri = Tri; clear Tri;
                    Surf.Area = Area_Comp(Surf.SurfData);
                    Surf.Type = 'Mask';
                    Surf.Imp = 'vtk';
                    if strcmp(sa,'y')
                        if exist('Out_dir','var')
                            pth = Out_dir;
                        end
                        mkdir(pth,'VTK_Surfaces'); npth = [pth filesep 'VTK_Surfaces'];
                        [Out] = savematfile([npth filesep nm '.mat'], Surf);
                        OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                        Surfa = '';
                    elseif strcmp(sa,'n')
                        OutFiles = '';
                        Surfa{i,1} = Surf;
                    end
                otherwise
                    Tmn =  16777214 ;
                    Qmn =  16777215 ;
                    [pth,nm,ext] = fileparts(deblank(SurfF(i,:)));
                    fid = fopen(deblank(SurfF(i,:)),'rb','b');
                    if fid >0
                        MNumber = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
                        if(MNumber == Qmn)
                            Surf.Imp = 'frb';
                            Npoints = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
                            Nfaces = bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar');
                            Surf.SurfData.vertices = fread(fid, Npoints*3, 'int16') ./ 100 ;
                            for t=1:Npoints
                                for k=1:4
                                    Surf.SurfData.faces(t,k) =  bitshift(fread(fid,1,'uchar'), 16) + bitshift(fread(fid,1,'uchar'),8) + fread(fid,1,'uchar'); ;
                                end
                            end
                            Surf.SurfData.faces = Surf.SurfData.faces+1;
                        elseif (MNumber == Tmn)
                            Surf.Imp = 'frb';
                            tline = fgets(fid);
                            tline = fgets(fid);
                            Np = fread(fid,1,'int32');
                            Nf = fread(fid,1,'int32');
                            vertices = fread(fid,Np*3,'float32');
                            Surf.Name = [nm ext];
                            Surf.SurfData.vertices = reshape(vertices,3,Np)';
                            faces = fread(fid,Nf*3,'int32');
                            Surf.SurfData.faces = reshape(faces,3,Nf)'+1;
                            fclose(fid);
                        else
                            errordlg('This file is not a known Surface File. Please try again');
                            return
                        end
                        if ~isempty(IMFile(i,:))&~strcmp(IMFile(i,:),'0')
                            V = spm_vol(IMFile(i,:));
                            Surf.Dim = V.dim(1:3);
                            Surf.Orig = abs(V.mat(1:3,4)');
                            Surf.VoxSize = abs(diag(V.mat(1:3,1:3))');
                        else
                            dx = ceil(max(Surf.SurfData.vertices(:,1))-min(Surf.SurfData.vertices(:,1)))+2;
                            dy = ceil(max(Surf.SurfData.vertices(:,2))-min(Surf.SurfData.vertices(:,2)))+2;
                            dz = ceil(max(Surf.SurfData.vertices(:,3))-min(Surf.SurfData.vertices(:,3)))+2;
                            Surf.Orig = [abs(min(Surf.SurfData.vertices(:,1)))+ abs(dx/2)+1 abs(min(Surf.SurfData.vertices(:,2)))+ abs(dy/2)+1 abs(min(Surf.SurfData.vertices(:,3)))+ abs(dz/2)+1];
                            Surf.Dim = [abs(dx) abs(dy) abs(dz)];
                            Surf.VoxSize = [1 1 1];
                        end
                        [Tri] = Vert_Neib(double(Surf.SurfData.faces),Np,Nf);
                        Temp = sum(Tri);
                        Tri(:,Temp==0) = [];
                        Surf.Tri = Tri;
                        Surf.Name = nm;
                        Surf.Area = Area_Comp(Surf.SurfData);
                        Surf.Type = 'Mask';
                        if strcmp(sa,'y')
                            if exist('Out_dir','var')
                                pth = Out_dir;
                            end
                            mkdir(pth,'MAT_Surfaces'); npth = [pth filesep 'MAT_Surfaces'];
                            [Out] = savematfile([npth filesep nm '.mat'], Surf);
                            OutFiles = strvcat(OutFiles,[npth filesep nm '.mat']);
                            %eval(['save ' deblank(OutFiles(i,:)) ' Surf']);
                            Surfa = '';
                        elseif strcmp(sa,'n')
                            OutFiles = '';
                            Surfa{i,1} = Surf;
                        end
                    else
                        errordlg('This file is not a known Surface File. Please try again');
                        return
                    end
            end
    end
end
%=========================End of main program=============================%
return;







