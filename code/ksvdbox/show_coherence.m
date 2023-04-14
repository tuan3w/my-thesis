function [] = show_coherence(dic)

h = [];
n = size(dic,2);
for i=1:n-1
    for j=i+1:n
        h = [h, abs(dic(:,i)' * dic(:,j))];
    end
end
[n,x] = hist(h,0:0.1:1);
bar(x, n./sum(n),.25, 'hist');

end
