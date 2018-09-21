function Out_mesh = save_mesh(Surf,Out_mesh);
%
% Syntax :
% Out_mesh = save_mesh(Surf,Out_mesh);
%
% Script file to save brainvisa mesh files
%
% Input Parameters:
%   Surf              :  Surface variable
%   Out_mesh          :  Output mesh filename
%
% Output Parameters:
%   Out_mesh          :  Output mesh filename
%
%
% Related references:
%
%
% See also: load_mesh
% 
%
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% April 9th 2012
% Version $1.0

%=================== Main Program ========================================%

fid = fopen(Out_mesh,'wb');
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
        Surft.SurfData = fv;
        Surft = Compute_Surface_Normals(Surft);
        normals = Surft.SurfData.VertexNormals;
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
%========================================================================%
return
