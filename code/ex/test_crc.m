addpath(genpath('../ksvdbox/'));
%load('featurevectors.mat');
load ./data_usps.mat
training_feats = normcols(training_feats);
testing_feats = normcols(testing_feats);


acc = 0.0;
t = 0.0;
lambda = 0.1;
l = 0.1:0.1:10
for lambda =l
subdictsize = 90;   % use full samples
for i=1:1
	%[~,~, D, H] = obtaintraingtestingsamples(training_feats, H_train,subdictsize);  % init dictionary by sampling from data
    [acci, dt] = CRC_RLS(training_feats, H_train, testing_feats, H_test, lambda);
    t = t + dt;
    acc = acc + acci;
    fprintf('\n[%d] lambda=%f, acc = %f', i, lambda, acci);
end
end
av = acc/5;
t = t/5;
fprintf('\nAverage acc = %f, average testing time %f', av, t);
