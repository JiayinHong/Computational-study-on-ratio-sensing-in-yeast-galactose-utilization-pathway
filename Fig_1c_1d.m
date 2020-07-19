% this script generate diagrams to illustrate the biological implications
% of decision front shift and slope variation in Fig 1c & 1d
% note that the diagrams are only for illustration purpose,
% not real data/simulation

cmap = parula;
a = logspace(5,9,100);  % gal titration
b = logspace(5,9,100);  % glu titration
[A,B] = meshgrid(a,b);

figure
set(gcf,'position',[360 77 714 621]);

param_values = [10^4, 10^4, .5*10^5, 10^7; 10^4, 10^6, .5*10^5, 10^7; 2, 2, 1, .1];
nA = 2;
i = 1;
for param_value = param_values
    Kglu = param_value(1);
    Kgal = param_value(2);
    nB = param_value(3);
    ff = 1./(1+(Kgal^nA)./(Kglu^nB).*(B.^nB)./(A.^nA)+(Kgal^nA)./(A.^nA));
    
    subplot(2,2,i)
    contourf(a,b,ff,5);hold on;
    set(gca,'xscale','log','yscale','log','zscale','linear');
    set(gca,'fontsize',15,'fontname','Times New Roman','fontweight','normal')
    set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5))
    set(gca,'xticklabel',{'2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}'})
    set(gca,'yticklabel',{'2^{-5}','2^{-4}','2^{-3}','2^{-2}','2^{-1}'})   
    caxis manual
    caxis([0 1]);
    if i==2
        hc = colorbar('Ticks',[0, 0.2, 0.4, 0.6, 0.8, 1]);
        colormap(parula(6))
        hc.Position(1) = hc.Position(1)+0.07;
    end
    i = i+1;  
end

export_fig('../figures/Fig-1c-1d', '-pdf','-transparent','-c[NaN NaN NaN NaN]')

