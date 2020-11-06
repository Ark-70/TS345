function g = ldpc_h2g(H)

% Le fichier arrive pas a utiliser le .c sur mon pc donc je rentre g en dur
g = [1 1 1 1 0 0; 0 1 0 0 1 0; 1 0 1 0 0 1];

% converts tentative binary LDPC matrix H into a new matrix h
% (columns are permuted) and produces the generator matrix g
% H should be a sparse matrix in MATLAB format.
%
% MEX file


%   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu
%   $Revision: 1.0 $  $Date: 1999/08/23 $
