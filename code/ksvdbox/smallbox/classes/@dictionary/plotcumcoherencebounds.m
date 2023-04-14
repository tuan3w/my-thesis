function plotcumcoherencebounds(obj)

[d N] = size(obj.phi);

mu1 = cumcoherence(obj,1);
mu2 = cumcoherence(obj,2);
mu1_min = (1:d)*sqrt((N-d)/(d*(N-1)));

figure, subplot(2,1,1)
hold on, grid on
plot(1:d,mu1,'k-')
plot(1:d,mu2,'k--')
plot(1:d,mu1_min,'b-.')
set(gca,'XScale','log','YScale','log');
axis tight
ylabel('p-cumulative coherence')
xlabel('k')
legend('\mu_1(k)','\mu_2(k)','\mu_{1_{min}}(k)')

% temp = conv(mu1,[1 1]);
% kMax = find(temp<1,1,'last');
% title(['Worst case bound: k_{max} = ' num2str(kMax)])
% 
% % subplot(2,1,2)
% p_BP_fails = exp(-1./(8*mu1(1)^2*(1:d)));
% p_TH_fails = exp(-1./(128*mu2.^2));
% plot(1:d,p_BP_fails,'k-+');
% plot(1:d,p_TH_fails,'k-o');
% legend('P(BP fails)','P(TH fails)');

% % Plot cumulative coherence with lower and upper bounds
% mumin = (1:d)*sqrt((N-d)/(d*(N-1)));
% mumax = (1:d);
% figure,
% subplot(1,6,1:2)
% hold on, grid on
% plot(1,mumax(1),'k-s');
% plot(1,mu(1),'ko');
% plot(1,mumin(1),'k-d')
% 
% 
% subplot(1,6,3:6)
% hold on, grid on
% plot(2:d,mumax(2:end),'k-s');
% plot(2:d,mu(2:end),'k-o');
% plot(2:d,mumin(2:end),'k-d');
% set(gca,'XScale','log','YScale','log');
% axis tight
% xlabel('k');
% ylabel('\mu(k)');
% ylim([mumin(1) 10])
% 
% % Plot recovery bounds
% plot(2:d,1-mu(1:d-1),'r-o')
% plot(2:d,1-mumin(1:d-1)','r-d')
% plot([2 d],[1/3 1/3],'b');
% legend('\mu_{max}(k)','\mu(k)','\mu_{min}(k)','Exact-Sparse \mu','Exact-Sparse \mu_{min}','Sparse');


% v = conv(mu,[1 1]);
% ind = find(v<1, 1, 'last');
% 
% line([ind ind], [min(mu) max(mu)],'Color','red');
% title(['Minimum allowed sparsity (Tanner):' num2str(ind/obj.len)]);
end
