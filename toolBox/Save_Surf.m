function [OutFiles] = Save_Surf(SurfF, SurfFiles);
%
% Syntax :
% [OutFiles] = Save_Surf(SurfF, SurfFile);
%
% Save the matlab surface SurfF in the file SurfFile
%
% Input Parameters:
%   SurfF       : Surfaces filenames or matlab variables.
%   SurfFiles   : New surface name.
%
% Output Parameters:
%  OutFiles     : List with the saved surfaces.

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
    [SurfF,sts] = spm_select([1],'any','Selecting Surface Files','',cd);
    [FileName, Out_dir] = uiputfile('Saving Surface, please, specify surface format...');
    SurfFiles = [Out_dir FileName];
end
if ~exist('SurfF','var')|isempty(SurfF)
    [SurfF,sts] = spm_select([1 Inf],'any','Selecting Surface Files','',cd);
end
if ~exist('SurfFiles','var')|isempty(SurfFiles)
    [FileName, Out_dir] = uiputfile('Saving Surface, please, specify surface format...');
    SurfFiles = [Out_dir FileName];
end
%=========================================================================%
%=========================Main program====================================%
wh = whos('SurfF');
if (strcmp(wh.class,'struct'))
    Ns = length(SurfF);
elseif ischar(SurfF(1,:));
    Ns = size(SurfF,1);
elseif (strcmp(wh.class,'cell'))
    Ns = length(SurfF);
    for k = 1:Ns
        Surf(k) = SurfF(1,k);
        clear SurfF; SurfF = Surf;
    end
    Ns = size(Surf,1);
end
OutFiles = '';
wh = whos('SurfF');
for i = 1:Ns
    [pthm,nmm,extm] = fileparts(SurfFiles(i,:));
    if strcmp(wh.class,'struct')
        Surf = SurfF(i);
    elseif ischar(SurfF(1,:));
        [OutputFiles, Surfa] = Exp_Surf(SurfF(i,:), '0', '','', 'imp','n');
    end
    form = extm(2:4);
    for j = 1:length(Surf)
        if ~isfield(Surf(j),'Name')
            Surf(j).Name =  ['Surface_' num2str(j)];
        end
        if length(Surf)==1
            TxtFile = SurfFiles(i,:);
        else
            TxtFile = [SurfFiles(i,:) '_' num2str(j) '_' Surf(j).Name];
        end
        switch lower(form)
            case 'mat'
                [Out] = savematfile(TxtFile, Surf);
                OutFiles = strvcat(OutFiles,TxtFile);
            case 'srx'
                vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                face = Surf(j).SurfData.faces;
                Nfaces = size(Surf(j).SurfData.faces,1);
                Npoints = size(Surf(j).SurfData.vertices,1);
                dim = Surf(j).Dim;
                if sum(Surf(j).Orig)==0
                    vert = Surf(j).SurfData.vertices+repmat([91 127 73],[size(Surf(j).SurfData.vertices,1),1]);
                end
                OutFiles = strvcat(OutFiles,TxtFile);
                AT = inv([1 0 0 1; 0 0 -1 dim(2); 0 1 0 1; 0 0 0 1]);
                xyzn = [vert ones(size(vert,1),1)]*AT'; xyzn(:,4) = [];
                fid = fopen(TxtFile, 'w');
                fwrite(fid, size(xyzn, 1), 'int32');
                fwrite(fid, size(face, 1), 'int32');
                temp = xyzn'; temp = temp(:);
                fwrite(fid, temp, 'float32');
                temp = face'-1; temp = temp(:);
                fwrite(fid, temp, 'int32');
                fclose(fid);
            case 'txt'
                vert = Surf(j).SurfData.vertices+repmat(Surf(j).Orig,[size(Surf(j).SurfData.vertices,1),1]);
                face = Surf(j).SurfData.faces;
                Nfaces = size(Surf(j).SurfData.faces,1);
                tri = Surf(j).Tri;
                Npoints = size(Surf(j).SurfData.vertices,1);
                dim = Surf(j).Dim;
                if sum(Surf(j).Orig)==0
                    vert = Surf(j).SurfData.vertices+repmat([91 127 73],[size(Surf(j).SurfData.vertices,1),1]);
                end
                OutFiles = strvcat(OutFiles,TxtFile);
                AT = inv([1 0 0 1; 0 0 -1 dim(2); 0 1 0 1; 0 0 0 1]);
                xyzn = [vert ones(size(vert,1),1)]*AT'; xyzn(:,4) = [];
                xyzn = [(1:size(xyzn,1))' xyzn];
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
            case 'off'
                vert = Surf(j).SurfData.vertices;
                face = Surf(j).SurfData.faces;
                Nfaces = size(Surf(j).SurfData.faces,1);
                Npoints = size(Surf(j).SurfData.vertices,1);
                OutFiles = strvcat(OutFiles,TxtFile);
                xyzn = [Npoints Nfaces 0; vert];
                faces = [3*ones(Nfaces,1) face-1];
                fid = fopen(TxtFile,'w');
                fprintf(fid,'%s\n','OFF');
                fprintf(fid,'%f %f %f\n',xyzn');
                fprintf(fid,'%u %u %u %u\n',faces');
                fclose(fid);
            case 'obj'
                
                fid = fopen(TxtFile,'wt');
                fprintf(fid,'%s\n','# Max2Obj Version 4.0 Mar 10th, 2001');
                fprintf(fid,'%s\n','#');
                [pth,nm,ext] = fileparts(TxtFile);
                ind = strfind(nm,'.');nm(ind) = [];
                Ns = length(Surf);
                Nfaces = 0;
                Surft = Surf(j);
                
                %% Saving Vertices
                fprintf(fid, '\n');
                fprintf(fid, '%s\n','g');
                fprintf(fid, '%s\n',['# object Surf_' num2str(i) ' to come ...']);
                %     fprintf(fid, '%s\n',['# object ' Matname ' to come ...']);
                fprintf(fid, '%s\n','#');
                
                if isfield(Surft.SurfData,'FaceVertexCData' )
                    Mat = [Surft.SurfData.vertices Surft.SurfData.FaceVertexCData]';
                    fprintf(fid,'v %.6f %.6f %.6f %.6f %.6f %.6f\n', Mat(:));
                else
                    Mat = Surft.SurfData.vertices';
                    fprintf(fid,'v %.6f %.6f %.6f\n', Mat(:));
                end
                fprintf(fid, '%s\n',['# ' num2str(size(Surft.SurfData.vertices,1)) ' vertices']);
                fprintf(fid, '\n');
                fprintf(fid, '%s\n',['g Surf_' num2str(i)]);
                %% Saving Structure Faces
                
                Mat = Surft.SurfData.faces'+Nfaces;
                Newfaces = size(Surft.SurfData.vertices,1);
                Nfaces = Nfaces+max(Surft.SurfData.faces(:));
                fprintf(fid,'%s\n', ['s 2']);
                fprintf(fid,'f %u %u %u\n', Mat(:));
                fprintf(fid, '%s\n',['# ' num2str(size(Surft.SurfData.faces,1)) ' faces']);
                % -----------------------------------------
                
                
                fprintf(fid, '\n');
                fprintf(fid, '%s','g');
                fclose(fid);
                
%                 
%                 
%                 
%                 
%                 
%                 
%                 vert = Surf(j).SurfData.vertices;
%                 face = Surf(j).SurfData.faces;
%                 if isfield(Surf(j).SurfData,'VertexNormals')
%                     normals = Surf(j).SurfData.VertexNormals;
%                 elseif (~isfield(Surf(j).SurfData,'VertexNormals'))&((strcmp(form,'dfs')|strcmp(form,'obj')|strcmp(form,'mesh')))
%                     fv.vertices = vert;
%                     fv.faces = face;
%                     h = figure;h1 = patch(fv);  normals = get(h1,'VertexNormals');close(h);clear fv;
%                     norma = sqrt(sum((normals').^2));
%                     normals = normals./repmat(norma',[1 3]);
%                 end
%                 Nfaces = size(Surf(j).SurfData.faces,1);
%                 Npoints = size(Surf(j).SurfData.vertices,1);
%                 OutFiles = strvcat(OutFiles,TxtFile);
%                 fid = fopen(TxtFile,'wt');
%                 fprintf(fid,'%s %1.1f %1.1f %1.1f %u %1.1f %u\n','P',0.3,0.3,0.4,10,1,Npoints);
%                 fprintf(fid,' %6.4f %6.4f %6.4f\n',vert');
%                 fprintf(fid,'\n');
%                 fprintf(fid,' %1.6f %1.6f %1.6f\n',normals');
%                 fprintf(fid,'\n');
%                 fprintf(fid,' %u\n',Nfaces);
%                 if isfield(Surf,'Is');
%                     [Colors] = Surf_Color(Surf,'jet');
%                     fprintf(fid,' %u %1.1f %1.1f %1.1f %u\n',[2 Colors(1,:) 1]);
%                     fprintf(fid,' %1.1f %1.1f %1.1f %u\n',[Colors(2:end,:) ones(Npoints-1,1)]');
%                     fprintf(fid,'\n');
%                 else
%                     fprintf(fid,' %u %u %u %u %u\n',[0 1 1 1 1]);
%                     fprintf(fid,'\n');
%                 end
%                 temp = [3:3:3*Nfaces]';
%                 temps = repmat(Nfaces,[1 8]);modp = mod(temps,[2:9]);
%                 ind = find(modp==0);
%                 if ~isempty(ind)
%                     Nl = max(ind)+1;
%                 else
%                     pr = primes(Nfaces);
%                     temps = repmat(Nfaces,[1 size(pr,1)]);
%                     modp = mod(temps,pr);
%                     ind = find(modp==0);Nl = pr(min(ind));
%                 end
%                 tempt = reshape(temp,[Nfaces/Nl Nl]);
%                 cad = [repmat(' %u',[1 Nl-1])  ' %u\n'];
%                 fprintf(fid,cad,tempt);
%                 fprintf(fid,'\n');
%                 tempf = face';tempf = tempf(:)-1;
%                 temps = repmat(Nfaces*3,[1 8]);modp = mod(temps,[2:9]);
%                 ind = find(modp==0);
%                 if ~isempty(ind)
%                     Nl = max(ind)+1;
%                 else
%                     pr = primes(Nfaces*3);
%                     temps = repmat(Nfaces*3,[1 size(pr,1)]);
%                     modp = mod(temps,pr);
%                     ind = find(modp==0);Nl = pr(min(ind));
%                 end
%                 tempf = reshape(tempf,[3*Nfaces/Nl Nl]);
%                 cad = [repmat(' %u',[1 Nl-1])  ' %u'];
%                 fprintf(fid,cad,tempf);
%                 fclose(fid);
            case 'dfs'
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
            case 'mes'
                OutFiles = strvcat(OutFiles,TxtFile);
                fid = fopen(TxtFile,'wb');
                fwrite(fid, 'binar', 'uchar');
                fwrite(fid, hex2dec('41424344'), 'uint32');
                fwrite(fid, 4, 'uint32');
                fwrite(fid, 'VOID', 'uchar');
                fwrite(fid, size(Surf(1).SurfData.faces,2), 'uint32');
                fwrite(fid, size(Surf,2), 'uint32');
                for k=1:size(Surf,2)
                    vert = Surf(k).SurfData.vertices;
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
                        Nfaces = size(Surf(k).SurfData.faces,1);
                        Npoints = size(Surf(k).SurfData.vertices,1);
                        face = face';face = face(:)-1;
                        normals = normals';normals = normals(:);
                        vert = vert';vert = vert(:);
                        if size(Surf,2) == 1
                            fwrite(fid, k-1, 'uint32');
                        else
                            fwrite(fid, k, 'uint32');
                        end
                        fwrite(fid, size(vert,1)/3, 'uint32');
                        fwrite(fid, vert', 'float32');
                        fwrite(fid, size(normals,1)/3, 'uint32');
                        fwrite(fid, normals', 'float32');
                        fwrite(fid, 0, 'uint32');
                        fwrite(fid, size(face,1)/3, 'uint32');
                        fwrite(fid, face', 'uint32');
                    end
                    fclose(fid);
            case 'asc'
                name = Surf(j).Name;
                vert = Surf(j).SurfData.vertices;
                face = Surf(j).SurfData.faces;
                OutFiles = strvcat(OutFiles,TxtFile);
                fid = fopen(TxtFile,'wt');
                vert =[vert zeros(size(vert,1),1)];
                face =[face-1 zeros(size(face,1),1)];
                fprintf(fid,'%s\n',['#!ascii version of ' name]);
                fprintf(fid,'%u %u\n',size(vert,1),size(face,1));
                fprintf(fid,'%8.6f %8.6f %8.6f %u\n',vert');
                fprintf(fid,'%u %u %u %u\n',face');
                fclose(fid);
            case 'frb'
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
            case 'vtk'
                Surftemp = Surf(j);
                name = Surftemp.Name;
                OutFiles = strvcat(OutFiles,TxtFile);
                fid = fopen(TxtFile, 'wt') ;
                Npoints = size(Surftemp.SurfData.vertices,1) ;  % number of vertices
                fprintf(fid,'%s\n','# vtk DataFile Version 3.0');
                fprintf(fid,'%s\n',['created from IBASPM on ' datestr(now)]);
                fprintf(fid,'%s\n','ASCII');
                fprintf(fid,'%s\n','DATASET POLYDATA');
                fprintf(fid,'%s %u %s\n','POINTS',Npoints,'float');
                fprintf(fid,'%8.6f %8.6f %8.6f\n',Surftemp.SurfData.vertices');
                Nfaces = size(Surftemp.SurfData.faces,1) ;  % number of faces
                Pol = size(Surftemp.SurfData.faces,2);
                face = [repmat(Pol,[Nfaces 1]) Surftemp.SurfData.faces-1];
                fprintf(fid,'%s %u %u\n','POLYGONS',Nfaces,Nfaces*(Pol+1));
                fprintf(fid,'%u %u %u %u\n',face');
                
                
                if ~isfield(Surf(j).SurfData,'VertexNormals')
                    Surftemp = Compute_Surface_Normals(Surftemp);
                end
                if ~isfield(Surftemp.SurfData,'FaceVertexCData')
                    if isfield(Surftemp,'Is')
% %                         fprintf(fid,'%s %u\n','POINT_DATA',Npoints);
% %                         fprintf(fid,'%s\n','SCALARS Scalars float');
% %                         fprintf(fid,'%s\n','LOOKUP_TABLE default');
% %                         fprintf(fid,'%8.6f\n',Surftemp.Is');
                        Surftemp.SurfData.FaceVertexCData = Surf_Color(Surftemp,'jet');
                    end
                end
                fprintf(fid,'%s %u\n','POINT_DATA',Npoints);
                fprintf(fid,'%s\n','COLOR_SCALARS lut 4');
                fprintf(fid,'%8.6f %8.6f %8.6f %8.6f\n',[Surftemp.SurfData.FaceVertexCData ones(Npoints,1)]');
                fprintf(fid,'%s\n','VECTORS Vectors float');
                fprintf(fid,'%8.6f %8.6f %8.6f\n',Surftemp.SurfData.VertexNormals');
                
                
                fclose(fid);
            otherwise
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
        end
    end
end
    %=========================End of main program=============================%
    return;

