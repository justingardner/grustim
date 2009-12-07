function y = qmult(a,b,base)

% qmult(a,b,base)
%
% multiplication in a Galois Field when base is a power of a prime number
% Reference: K. Godfrey, Perturbation Signals for System Identificaton,
% 1993
%
% send comments and questions to ttliu@ucsd.edu  

if (a >= base) | (b >= base)
  error(sprintf('qadd(a,b), a and b must be < %d',base));
end
  
switch base
  case 4
   amult = [ 0 0 0 0;
	     0 1 2 3;
	     0 2 3 1;
	     0 3 1 2;];
	   
	   
	   
 case 8
  amult = [ 0 0 0 0 0 0 0 0;
	    0 1 2 3 4 5 6 7;
	    0 2 4 6 5 7 1 3;
	    0 3 6 5 1 2 7 4;
	    0 4 5 1 7 3 2 6;
	    0 5 7 2 3 6 4 1;
	    0 6 1 7 2 4 3 5;
	    0 7 3 4 6 1 5 2;];

 case 9
    amult = [0 0 0 0 0 0 0 0 0;
	     0 1 2 3 4 5 6 7 8;
	     0 2 1 6 8 7 3 5 4;
	     0 3 6 4 7 1 8 2 5;
	     0 4 8 7 2 3 5 6 1;
	     0 5 7 1 3 8 2 4 6;
	     0 6 3 8 5 2 4 1 7;
	     0 7 5 2 6 4 1 8 3;
	     0 8 4 5 1 6 7 3 2;];
   
 otherwise;
  error(sprintf('qmult base %d not supported yet',base));
end

y = amult(a+1,b+1);



