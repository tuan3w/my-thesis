load ./featurevectors

dic_params.dictsize = 20;
dic_params.iternum = 30;
dic_params.Tdata = 4;
[l,~] = find(H);
D = [];
nclass = size(H, 1);
for i = 1:nclass
    idx = find(l == i);
    data_i = training_feats(:,idx);
    dic_params.data = data_i;


end

