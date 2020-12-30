function [H] = alist2sparse(fname)
% reads binary parity check matrix in "alist" format from file FNAME and
% converts it to sparse matrix used in MATLAB routines.
% This is an interface to matrices at https://github.com/aff3ct/configuration_files/tree/master/dec/LDPC
%
% Example
%        [H] = alist2sparse('A');   % A is the ascii file in alist format


%   Adapted from Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
%   $Revision: 1.1 $  $Date: 2000/03/23 $ Bug fixed by Hatim Behairy

fid = fopen(fname);

line = fgetl(fid);
while (strcmp(line(1), '#'))
    disp(['Reading : ', line]);
    line = fgetl(fid);
end
nbrs = sscanf(line,'%d');
n = nbrs(1);
m = nbrs(2);

line = fgetl(fid);
while (strcmp(line(1), '#'))
    disp(['Reading : ', line]);
    line = fgetl(fid);
end
nbrs = sscanf(line,'%d');
dmax_VN = nbrs(1);
dmax_CN = nbrs(2);

line = fgetl(fid);
while (strcmp(line(1), '#'))
    disp(['Reading : ', line]);
    line = fgetl(fid);
end
d_VN = sscanf(line,'%d',[1 n]); % number of elements in rows
num2(1:n)=dmax_VN;

line = fgetl(fid);
while (strcmp(line(1), '#'))
    disp(['Reading : ', line]);
    line = fgetl(fid);
end
d_CN = sscanf(line,'%d',[1 m]); % no need

line = fgetl(fid);
while (strcmp(line(1), '#'))
    disp(['Reading : ', line]);
    line = fgetl(fid);
end
position = zeros(n,dmax_VN);
nbrs = sscanf(line,'%d');
for j=1:length(nbrs)
    position(1,j) = nbrs(j);
end

for i=2:n
    line = fgetl(fid);
    nbrs = sscanf(line,'%d');
    for j=1:length(nbrs)
        position(i,j) = nbrs(j);
    end
end

ii = zeros(1,sum(d_VN));
jj = ii;
k = 1;
for i=1:n
    for j=1:d_VN(i)
        jj(k) = i;
        ii(k) = position(i,j);
        ss = 1;
        k = k+1 ;
    end
end
H = sparse(ii,jj,ss,m,n);
fclose(fid);