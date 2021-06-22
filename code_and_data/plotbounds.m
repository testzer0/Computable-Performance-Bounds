load('plotbounds.mat');

vals = 1:40;
mbyns = ['m/n=0.3'; 'm/n=0.4'; 'm/n=0.5'; 'm/n=0.6'; 'm/n=0.7'; 'm/n=0.8'; 'm/n=0.9'];
colors2 = ['--b','--g','--r','--c','--m','--y','--k'];
colors1 = ['b','g','r','c','m','y','k'];

figure(1);
clf(1);
set(gca, 'YScale', 'log')
hold on;

plots = zeros(1,7);

for i = 1:7
    idx = ~any(isnan(omegabounds(i,:)),1);
    plots(1,i) = plot(vals(idx),omegabounds(i,idx),colors1(i));
    plot(vals(idx),mcbounds(i,idx),colors2(i));
end

legend(plots,mbyns);
title('Omega vs MC bounds');