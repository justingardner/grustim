function y = qadd(a,b,base)
  
%    qadd(a,b,base)
%
% addition in a Galois Field for mod (power of prime)
% Reference: K. Godfrey, Perturbation Signals for System Identificaton,
% 1993
%
% send comments and questions to ttliu@ucsd.edu  
  

if (a >= base) | (b >= base)
  error(sprintf('qadd(a,b), a and b must be < %d',base));
end
  
switch base
  case 4
   amat = [0 1 2 3;
	   1 0 3 2;
	   2 3 0 1;
	   3 2 1 0];
 case 8
   amat = [0 1 2 3 4 5 6 7;
	   1 0 3 2 5 4 7 6;
	   2 3 0 1 6 7 4 5;
	   3 2 1 0 7 6 5 4;
	   4 5 6 7 0 1 2 3;
	   5 4 7 6 1 0 3 2;
	   6 7 4 5 2 3 0 1;
	   7 6 5 4 3 2 1 0;];
 case 9
  amat = [ 0 1 2 3 4 5 6 7 8;
	   1 2 0 4 5 3 7 8 6;
	   2 0 1 5 3 4 8 6 7;
	   3 4 5 6 7 8 0 1 2;
	   4 5 3 7 8 6 1 2 0;
	   5 3 4 8 6 7 2 0 1;
	   6 7 8 0 1 2 3 4 5;
	   7 8 6 1 2 0 4 5 3;
	   8 6 7 2 0 1 5 3 4;];
   
 otherwise;
  error(sprintf('qadd base %d not supported yet',base));
end

y = amat(a+1,b+1);



