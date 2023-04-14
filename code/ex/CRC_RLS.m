function [acc,dt, err2] = SRC_LS(training_feats, H_train, testing_feats, H_test, lambda)

D = normcols(training_feats);

[l,~] = find(H_train);
idx={};
nclass = size(H_test,1);
total = size(H_test,2);
P = (D'*D  + lambda * eye(size(D,2))) \ D'; 
%imagesc(Gamma);
%Gamma(Gamma<0) = 0.0;
Di = {};
for i = 1:nclass
    idx{i} = find(l==i);
	Di{i} = D(:,idx{i});
end
predict = [];
truePredict = 0;
[testlabel, ~]= find(H_test == 1);

%id = find(testlabel = 11)(1);
err = zeros(nclass, total);
tic
Gamma = P * testing_feats;
for i=1:nclass
    id = idx{i};
	Gi = Gamma(id,:);
    err(i,:) = sum((testing_feats - Di{i} * Gi).^2)./sum(Gi.^2);
end
[~, predict] = min(err);
dt = toc;
truePredict = sum(predict == testlabel');
acc = truePredict * 1.0/total;
dt = dt/total;
err2 = sum(sum((testing_feats - D * Gamma).^2));
%fprintf('\nAverage time : %f', dt);
