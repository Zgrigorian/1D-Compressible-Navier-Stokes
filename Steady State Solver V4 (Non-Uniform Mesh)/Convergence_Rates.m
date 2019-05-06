clc;
clear;
close all;
tfs=[0;.1;.2;.3;.4;.5;1;5];
dynamics=[0;548;556;563;571;578;616;917];
dogleg=[18468;15904;12825;7695;5130;4618;2025;2025];
trust=[28215;14364;18981;7695;5643;4104;2052;1539];
leven=[59050;44158;33888;23618;13348;6159;2565;1026];
figure(1);
plot(tfs,dynamics);
xlabel('Integration Time','interpreter','latex');
ylabel('Function Evaluations','interpreter','latex');
title('Dynamic Function Evaluations','interpreter','latex');
figure(2);
semilogy(tfs,dogleg,tfs,trust,tfs,leven);
legend('trust-region-dogleg','trust-region','levenberg-marquet');
xlabel('Integration Time');
ylabel('Function Evaluations');
title('Steady State Solver Function Evaluations');