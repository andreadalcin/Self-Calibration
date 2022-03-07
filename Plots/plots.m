clc, clear, close all

rng default;

A = [normrnd(17502, 6695, 100, 1), normrnd(12436, 3884, 100, 1), normrnd(14546, 2771, 100, 1)];
B = [normrnd(12436, 3884, 100, 1), normrnd(12436, 3884, 100, 1), normrnd(14546, 2771, 100, 1)];
C = [normrnd(12436, 3884, 100, 1), normrnd(12436, 3884, 100, 1), normrnd(14546, 2771, 100, 1)];
D = [normrnd(12436, 3884, 100, 1), normrnd(12436, 3884, 100, 1), normrnd(14546, 2771, 100, 1)];

D = squeeze(mat2cell(permute(cat(3,A,B,C,D),[1,3,2]),size(A,1),4,ones(1,size(A,2))));
%                                        Number of matricies ^
clf
boxplotGroup(D','PrimaryLabels', {'S' 'SKV' 'Ours'}, ...
  'SecondaryLabels',{'# F = 15', '# F = 25' '# F = 50', '# F = 75'}, 'InterGroupSpace', 1)