addpath(genpath('../ksvdbox/'));
load('featurevectors.mat');
training_feats = normcols(training_feats);
testing_feats = normcols(testing_feats);


acc = 0.0;
t = 0.0;
T = 70;
for i=1:5
    [acci, dt] = SRC(training_feats, H_train, testing_feats, H_test, T);
    t = t + dt;
    acc = acc + acci;
    fprintf('\n[%d] T=%f, acc = %f', i, T, acci);
end
av = acc/5;
t = t/5;
fprintf('\nAverage acc = %f, average testing time %f', av, t);
