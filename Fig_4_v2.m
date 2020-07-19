% concatenation of transporter and transcriptional circuits
% figure showing expanded dynamical range of ratio-sensing
% simulation results for Fig 4, version 2

%% typical signal integration patterns in transporter circuit
cmap=cbrewer('seq', 'YlGnBu', 9);
a = logspace(5,9,100);  % gal titration
b = logspace(5,9,100);  % glu titration
[A,B] = meshgrid(a,b);

figure
set(gcf,'position',[339 231 508 215]);
ha = tight_subplot(1,2,[.03 .1],[.13 .13],[.1 .13]);
nA = 1;
nB = 1;
i = 1;
for Kglu = [10^5,10^7]
    Kgal = Kglu;
    ind_level = 1./(1+(Kgal^nA)./(Kglu^nB).*(B.^nB)./(A.^nA)+(Kgal^nA)./(A.^nA));
    axes(ha(i));
    colormap(cmap);
    contourf(a,b,ind_level,5);hold on;
    set(gca,'xscale','log','yscale','log','zscale','linear');
    set(gca,'fontsize',15,'fontname','Times New Roman')
    set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5))
    
    if i==1
        title('K_{glu}=K_{gal}=10^5','FontName','Times New Roman'...
            ,'FontSize',14, 'FontWeight','bold');
    end
    if i==2
        title('K_{glu}=K_{gal}=10^7','FontName','Times New Roman'...
            ,'FontSize',14, 'FontWeight','bold');
        cb = colorbar('location','east');   cb.Limits(1)=0;
        cb.Ticks=[0:.2:.8]; cb.FontName = 'Times New Roman';
        cb.Position(1)=cb.Position(1)+.12;
        cb.FontSize = 16;
    end
    i = i+1;
end

export_fig('../figures/Fig-4a-v2', '-pdf','-transparent','-c[NaN NaN NaN NaN]')

%% typical signal integration patterns in transcription circuit
cmap=cbrewer('seq', 'PuBuGn', 9);
a = logspace(5,9,100);  % gal titration
b = logspace(5,9,100);  % glu titration
[A,B] = meshgrid(a,b);

figure
set(gcf,'position',[133 250 291 939]);
ha = tight_subplot(4,1,[.06 .05],[.05 .05],[.3 .2]);
param_values = [10^9, 10^9, 10^8, 10^7; 10^9, 10^9, 10^7, 10^8; 10^4, 100, 100, 10^3; 10^4, 10^3, 10^3, 100];
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
        ylabel({'$\bf K_G=10^9\hspace{2mm} K_M=10^9$'...
            ,'$\bf \phi_R=10^4\hspace{3mm} \phi_A=10^4$'}...
            ,'interpreter','latex','FontSize',14)
        cb = colorbar('location','east');   cb.Limits(1)=0;
        cb.Ticks=[0:.2:.8]; cb.FontName = 'Times New Roman';
        cb.Position(1)=cb.Position(1)+.18;
        cb.FontSize = 16;
    end  
    if i==2
        ylabel({'$\bf K_G=10^9\hspace{2mm} K_M=10^9$'...
            ,'$\bf \phi_R=100\hspace{3mm} \phi_A=10^3$'}...
            ,'interpreter','latex','FontSize',14)
    end   
    if i==3
        ylabel({'$\bf K_G=10^8\hspace{2mm} K_M=10^7$'...
            ,'$\bf \phi_R=100\hspace{3mm} \phi_A=10^3$'}...
            ,'interpreter','latex','FontSize',14)
    end  
    if i==4
        ylabel({'$\bf K_G=10^7\hspace{2mm} K_M=10^8$'...
            ,'$\bf \phi_R=10^3\hspace{3mm} \phi_A=100$'}...
            ,'interpreter','latex','FontSize',14)
    end
    i = i+1;
end

export_fig('../figures/Fig-4b-v2', '-pdf','-transparent','-c[NaN NaN NaN NaN]')

%% concatenation of the two circuits
cmap=cbrewer('seq', 'YlOrBr', 9);
a = logspace(5,9,100);  % gal titration
b = logspace(5,9,100);  % glu titration
[A,B] = meshgrid(a,b);

figure
set(gcf,'position',[210 256 442 979]);
ha = tight_subplot(4,2,[.06 .1],[.05 .08],[.18 .12]);
param_values = [10^9, 10^9, 10^8, 10^7; 10^9, 10^9, 10^7, 10^8; 10^4, 100, 100, 10^3; 10^4, 10^3, 10^3, 100];
Psi = 10^8;    % sugar transportation capacity
i = 1;
for param_value = param_values
    for Kgal = [10^5,10^7]
        Kglu = Kgal;
        KG = param_value(1);
        KM = param_value(2);
        phiR = param_value(3);
        phiA = param_value(4);
        GALin = Psi ./(1+Kgal./Kglu.*B./A+Kgal./A);
        GLUin = Psi ./(1+Kglu./Kgal.*A./B+Kglu./B);
        ind_level = 1./(1+1./phiA.*(KG./GALin+1)+phiR./phiA.*(KG./GALin+1)./(KM./GLUin+1));
        axes(ha(i));
        colormap(cmap);
        contourf(a,b,ind_level,5);hold on;
        set(gca,'xscale','log','yscale','log','zscale','linear');
        set(gca,'fontsize',15, 'FontName','Times New Roman')
        set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5))
        
        if i==1
            ylabel({'$\bf K_G=10^9\hspace{2mm} K_M=10^9$'...
                ,'$\bf \phi_R=10^4\hspace{3mm} \phi_A=10^4$'}...
                ,'interpreter','latex','FontSize',14)
            title({'$\bf K_{gal}=10^5$'...
                '$\bf K_{glu}=10^5$'},'interpreter','latex'...
                ,'FontName','Times New Roman'...
                ,'FontSize',15, 'FontWeight','bold');
        end
        if i==3
            ylabel({'$\bf K_G=10^9\hspace{2mm} K_M=10^9$'...
                ,'$\bf \phi_R=100\hspace{3mm} \phi_A=10^3$'}...
                ,'interpreter','latex','FontSize',14)
        end
        if i==5
            ylabel({'$\bf K_G=10^8\hspace{2mm} K_M=10^7$'...
                ,'$\bf \phi_R=100\hspace{3mm} \phi_A=10^3$'}...
                ,'interpreter','latex','FontSize',14)
        end
        if i==7
            ylabel({'$\bf K_G=10^7\hspace{2mm} K_M=10^8$'...
                ,'$\bf \phi_R=10^3\hspace{3mm} \phi_A=100$'}...
                ,'interpreter','latex','FontSize',14)
        end
        if i==2
            set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5)...
                ,'FontName','Times New Roman','FontSize',14)
            title({'$\bf K_{gal}=10^7$'...
                '$\bf K_{glu}=10^7$'},'interpreter','latex'...
                ,'FontName','Times New Roman'...
                ,'FontSize',15, 'FontWeight','bold');
            cb = colorbar('location','east');   cb.Limits(1)=0;
            cb.Ticks=[0:.2:.8]; cb.FontName = 'Times New Roman';
            cb.FontSize = 16;
            cb.Position(1)=cb.Position(1)+.12;
        end
        i = i+1;
    end
end

export_fig('../figures/Fig-4c-v2', '-pdf','-transparent','-c[NaN NaN NaN NaN]')


