function out = randint(n, m, range, seed);
%RANDINT Random integer matrix generator.
%       OUT = RANDINT(N) generates an N-by-N matrix of random binary numbers. 
%       The appearance of "0" and "1" have even probability.
%
%       OUT = RANDINT(N, M) generates an N-by-M matrix of random binary 
%       numbers. The appearance of "0" and "1" have even probability.
%
%       OUT = RANDINT(N, M, RANGE) generates an N-by-M matrix of random 
%       integer numbers.
%
%       The minimum and maximum output integer numbers are specified in 
%       RANGE(1) and RANGE(2). All integers in range [RANGE(1), RANGE(2)] 
%       have even probability and appear in the output variable OUT. When
%       RANGE is a scalar positive integer, the output integer range is 
%       [0, RANGE-1].
%       
%       OUT = RANDINT(N, M, RANGE, SEED) specifies the random number seed 
%       in SEED.
%
%       See also RAND, DE2BI, BI2DE.

%       Wes Wang 6/13/94, 10/11/95.
%       Copyright (c) 1996-98 by The MathWorks, Inc.
%       $Revision: 1.8 $

%rand('uniform')

if nargin > 3
    rand('seed', seed');
end;
if nargin < 2
    m = n;
end;
if nargin < 3
    range = [0, 1];
end;

len_range = length(range);
if len_range < 2
    if len_range < 1
        range = [0, 1];
    elseif range < 0
        error('Incorrect range assignment in RANDINT');
    else
        range = [0, range-1];
    end
end;
range = sort(range);

range(1) = ceil(range(1));
range(2) = floor(range(2));
if range(1) == range(2) 
    out = ones(n, m) * range(1);
    return;
end;

distance = range(2) - range(1);

r = rand(n, m);
out = ones(n,m)*range(1);
for i = 1:distance
    ind = find(r >= i/(distance+1));
    out(ind) = (range(1) + i) * ind./ind;
end;

%--end randint.m--
