function [T,p] = statchol(sigma)
%STATCHOL  Do Cholesky-like decomposition, allowing zero eigenvalues
%   SIGMA must be symmetric.  In general T is not square or triangular.
%   P is the number of negative eigenvalues, and T is empty if P>0.

%   $Revision: 1.1 $  $Date: 2001/08/01 20:57:49 $

[T p] = chol(sigma);

if p > 0
   % Can get factors of the form sigma==T'*T using the eigenvalue
   % decomposition of a symmetric matrix, so long as the matrix
   % is positive semi-definite.
   [U,D] = eig((sigma+sigma')/2);
   D = diag(D);
   
   tol = max(D) * length(D) * eps;
   t = (abs(D) > tol);
   D = D(t);
   p = sum(D<0);
   
   if (p==0)
      T = diag(sqrt(D)) * U(:,t)';
   else
      T = [];
   end
end
