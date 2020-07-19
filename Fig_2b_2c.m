% competitive binding at transporter level can enable ratio sensing
% simulation results for Fig 2b & 2c

cmap=cbrewer('seq', 'YlGnBu', 9);
a = logspace(5,9,100);  % gal titration
b = logspace(5,9,100);  % glu titration
[A,B] = meshgrid(a,b);

%% Fig 2b, show decision front shift & threshold regime 
figure
set(gcf,'position',[360 224 780 474]);
ha = tight_subplot(2,3,[.1 .07],[.1 .1],[.1 .1]);
i = 1;
nA = 1;
nB = 1;
for Kglu = [10^5 10^7]
    for r = [0.1 1 10]
        Kgal = r*Kglu;
        ind_level = 1./(1+(Kgal^nA)./(Kglu^nB).*(B.^nB)./(A.^nA)+(Kgal^nA)./(A.^nA));
        axes(ha(i));
        colormap(cmap)
        contourf(a,b,ind_level,5);hold on;
        set(gca,'xscale','log','yscale','log','zscale','linear');
        set(gca,'fontsize',15,'fontname','Times New Roman')
        set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5))
        
        if i==1
            title('$$\bf\frac{K_{gal}}{K_{gluc}} = 0.1$$'...
                ,'interpreter','latex','FontSize',16)
            % \bf for bold font, set fontsize inside latex interpreter
            yl1 = ylabel('K_{gluc} = 10^5','FontName','Times New Roman'...
                ,'FontSize',16, 'FontWeight','bold');
        end
        if i==2
            title('$$\bf\frac{K_{gal}}{K_{gluc}} = 1$$'...
                ,'interpreter','latex','FontSize',16)
        end
        if i==3
            title('$$\bf\frac{K_{gal}}{K_{gluc}} = 10$$'...
                ,'interpreter','latex','FontSize',16)
            cb = colorbar('location', 'east');
            cb.Limits(1)=0; cb.Ticks=[0, 0.2, 0.4, 0.6, 0.8];
            cb.Position(1)=cb.Position(1)+0.08; % manually shift the position of the colorbar
        end
        if i==4
            yl4 = ylabel('K_{gluc} = 10^7','FontName','Times New Roman'...
                ,'FontSize',16, 'FontWeight','bold');
            yl1.Position(1) = yl4.Position(1);
        end
        i = i+1;
    end
end

export_fig('../figures/Fig-2b', '-pdf','-transparent','-c[NaN NaN NaN NaN]')

%% Fig 2c, show decision front slope variation
figure
set(gcf,'position',[360 185 594 513]);
ha = tight_subplot(2,2,[.17 .2],[.1 .1],[.1 .12]);
param_values = [10^5, 10^7, 10^5, 10^7; 2, 2, 1, 0.1; 1, 0.1, 2, 2];
i = 1;
for param_value = param_values
        Kglu = param_value(1);
        Kgal = Kglu;
        nA = param_value(2);
        nB = param_value(3);
        ind_level = 1./(1+(Kgal^nA)./(Kglu^nB).*(B.^nB)./(A.^nA)+(Kgal^nA)./(A.^nA));
        axes(ha(i));
        colormap(cmap);
        contourf(a,b,ind_level,5);hold on;  
        set(gca,'xscale','log','yscale','log','zscale','linear');
        set(gca,'fontsize',15,'fontname','Times New Roman')
        set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5))
        cb = colorbar('location', 'east');  cb.Limits(1)=0;
        cb.Ticks=[0:0.2:0.8];
        cb.Position(1)=cb.Position(1)+0.09;
        
        if i==1
            title('$$\bf\frac{n_{gal}}{n_{gluc}} = 2$$'...
                ,'interpreter','latex','FontSize',15)
        end      
        if i==2
            title('$$\bf\frac{n_{gal}}{n_{gluc}} = 20$$'...
                ,'interpreter','latex','FontSize',15)
        end
        if i==3
            title('$$\bf\frac{n_{gal}}{n_{gluc}} = \frac{1}{2}$$'...
                ,'interpreter','latex','FontSize',15)
        end
        if i==4
             title('$$\bf\frac{n_{gal}}{n_{gluc}} = \frac{1}{20}$$'...
                ,'interpreter','latex','FontSize',15)
            cb.Ticks=[0:0.1:0.5];
        end
        i = i+1;
end

export_fig('../figures/Fig-2c', '-pdf','-transparent','-c[NaN NaN NaN NaN]')
