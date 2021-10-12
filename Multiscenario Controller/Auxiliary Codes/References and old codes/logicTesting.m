clear 
close all
clc

% number of ramifications
Nr = 2;
% number of levels
Nl = 4;
% number of scenarios
Ns = Nl^Nr;
% control horizon
Np = 10;

values = 1:Nl;
temp = dec2base(0:Ns-1,Nl)
temp2 = temp -'0'+1
combs = values(temp2)

% values = 1:Nl;                              %//  data
% k = Nr;                                      %//  data
% n = length(values);                          %//  number of values
% combs = values(dec2base(0:n^k-1,n)-'0'+1)  %//  generate all tuples


%combs = combs(all(diff(combs.')>=0, 1),:)   %'// keep only those that are sorted

