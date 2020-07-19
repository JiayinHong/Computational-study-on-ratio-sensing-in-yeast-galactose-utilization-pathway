% competitive binding at transcriptional level can enable ratio sensing
% different parameter regimes give rise to distinct signal integration
% pattern, such as glu(/gal) threshold sensing
% simulation results for Fig 3b - 3e

cmap=cbrewer('seq', 'PuBuGn', 9);
a = logspace(5,9,100);  % gal titration
b = logspace(5,9,100);  % glu titration
[A,B] = meshgrid(a,b);

figure
set(gcf,'position',[360 190 542 508]);
ha = tight_subplot(2,2,[.15 .12],[.1 .1],[.1 .12]);
param_values = [10^9, 10^7, 10^9, 10^8; 10^9, 10^8, 10^9, 10^7; 10^4, 10^3,100, 100; 10^4, 100, 10^3, 10^3];
i = 1;
for param_value = param_values
    KG = param_value(1);
    KM = param_value(2);
    phiR = param_value(3);
    phiA = param_value(4);
    ind_level = 1./(1+1./phiA.*(KG./A+1)+phiR./phiA.*(KG./A+1)./(KM./B+1));
    axes(ha(i));
    colormap(cmap);
    contourf(a,b,ind_level,5);hold on;
    set(gca,'xscale','log','yscale','log','zscale','linear');
    set(gca,'fontsize',15,'fontname','Times New Roman')
    set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5))
    
    if i==1
        title('$$\bf\frac{K_G}{K_M}=1\hspace{3mm}\bf\frac{\phi_R}{\phi_A}=1$$'...
            ,'interpreter','latex','FontSize',15)
    end
    if i==2
        title('$$\bf\frac{K_G}{K_M}=0.1\hspace{3mm}\bf\frac{\phi_R}{\phi_A}=10$$'...
            ,'interpreter','latex','FontSize',15)
        cb = colorbar('location','east');   cb.Limits(1)=0;
        cb.Ticks=[0:.2:.8];
        cb.Position(1)=cb.Position(1)+.09;
    end
    if i==3
        title('$$\bf\frac{K_G}{K_M}=1\hspace{3mm}\bf\frac{\phi_R}{\phi_A}=0.1$$'...
            ,'interpreter','latex','FontSize',15)
    end
    if i==4
        title('$$\bf\frac{K_G}{K_M}=10\hspace{3mm}\bf\frac{\phi_R}{\phi_A}=0.1$$'...
            ,'interpreter','latex','FontSize',15)
    end
    i = i+1;
end

export_fig('../figures/Fig-3b-3e', '-pdf','-transparent','-c[NaN NaN NaN NaN]')


