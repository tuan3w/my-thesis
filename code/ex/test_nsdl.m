addpath(genpath('../ksvdbox/'));
lambda = 0.012;
subdictsize = 20;
load('featurevectors.mat');
%load ./data_usps.mat
training_feats = normcols(training_feats);
testing_feats = normcols(testing_feats);
nclass = size(H_train, 1);


% init params
dic_params.iternum=50;
dic_params.data = training_feats;
dic_params.H = H_train;
dic_params.lambda2 = lambda;
dic_params.subdictsize= subdictsize;
dictsize = subdictsize * nclass;
dic_params.dictsize = dictsize;


% init H
H = zeros(nclass, dictsize);
offset = 1;
[l,~] = find(H_train);
idx = {};
for i = 1:nclass
    H(i,offset:offset + subdictsize -1) = 1;
    idx{i} = find(l == i);
    offset = offset + subdictsize;
end
dic_params.idx = idx;

% main
acc = 0.0;
c = 4;
t = 0.0;
for i=1:c
	%[D,~] = nsdl(dic_params, 'i');
    [D,~] = nsdl(dic_params, '');
   	[acci, dt, erri] = CRC_RLS(D, H, testing_feats, H_test, lambda);
	acc = acc + acci;
	t = t + dt;
	fprintf('\n[%d] Lambda=%f, acc = %f', i, lambda, acci);
end
av = acc/c;
avt = t/c;

fprintf('\nAverage acc = %f, average testing time %f', av, avt);
