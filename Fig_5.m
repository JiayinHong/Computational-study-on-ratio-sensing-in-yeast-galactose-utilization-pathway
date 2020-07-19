% topology modifications can alter the slope of the decision front
% simulation results to supplement Fig 5

cmap=cbrewer('seq', 'PuBuGn', 9);
a = logspace(5,9,50);  % gal titration
b = logspace(5,9,50);  % glu titration
[A,B] = meshgrid(a,b);
tic
figure; set(gcf,'position',[289 232 696 448]);
ii = 1;
% topologies that halve the slope of the decision front
for topos={'NAR_null','IE_null','NAR_AE','NAR_PAR','IE_PAR','IE_AE'}
    topo = topos{1};
    ind_level = get_ind_level(topo,A,B);
    subplot(2,3,ii)
    contourf(a,b,ind_level,5);hold on;
    title(changeunderscore(topo))
    colormap(cmap)
    set(gca,'xscale','log','yscale','log','zscale','linear');
    set(gca,'xtick','','ytick','')
    if ii==3
        cb = colorbar('location','east');   cb.Limits(1)=0;
        cb.Ticks=[0:.2:.8]; cb.FontName = 'Times New Roman';
        cb.Position(1)=cb.Position(1)+.09;
    end 
    set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5)...
        ,'FontName','Times New Roman','FontSize',14)
    ii=ii+1;
end
export_fig('../figures/Fig-5b', '-pdf','-transparent','-c[NaN NaN NaN NaN]')

figure; set(gcf,'position',[289 232 696 448]);
ii = 1;
% topologies that double the slope of the decision front
for topos={'null_NAR','null_IE','AE_NAR','PAR_NAR','PAR_IE','AE_IE'}
    topo = topos{1};
    ind_level = get_ind_level(topo,A,B);
    subplot(2,3,ii)
    contourf(a,b,ind_level,5);hold on;
    title(changeunderscore(topo))
    colormap(cmap)
    set(gca,'xscale','log','yscale','log','zscale','linear');
    set(gca,'xtick','','ytick','')
    if ii==3
        cb = colorbar('location','east');   cb.Limits(1)=0;
        cb.Ticks=[0:.2:.8]; cb.FontName = 'Times New Roman';
        cb.Position(1)=cb.Position(1)+.09;
    end
    set(gca,'xtick',logspace(5,9,5),'ytick',logspace(5,9,5)...
        ,'FontName','Times New Roman','FontSize',14)
    ii=ii+1;
end
export_fig('../figures/Fig-5c', '-pdf','-transparent','-c[NaN NaN NaN NaN]')


function ind_level = get_ind_level(topo,A,B)
% set parameter values
beta1 = 1000; % production rate for activator or repressor
beta2 = 10; % production rate for repressor or activator
alpha = 0.001;  % degradation rate for activator and repressor
KA = 1; % dissociation constant for activator binding to UAS
KR = 1; % dissociation constant for repressor binding to URS
Kx = 1; % auto-regulation Michaelis constant for activator
Ky = 1; % auto-regulation Michaelis constant for repressor
KM = 10^9;  % half-saturation constant for glu binding to Mig1p
KG = 10^9;  % half-saturation constant for gal binding to Gal3p
switch topo
    case 'NAR_null'
        Atotal = sqrt(beta1./alpha.*Kx.*(KG+A)./A);
        Rtotal = beta2/alpha;
    case 'IE_null'
        Rtotal = beta2/alpha;
        Atotal = beta1./alpha./(1+Rtotal./Kx.*B./(KM+B));
    case 'NAR_AE'
        Atotal = sqrt(beta1./alpha.*Kx.*(KG+A)./A);
        Rtotal = beta2./alpha./(1+Ky.*(KG+A)./A./Atotal);
    case 'NAR_PAR'
        Atotal = sqrt(beta1./alpha.*Kx.*(KG+A)./A);
        Rtotal = beta2./alpha - Ky.*(KM./B+1);
    case 'IE_PAR'
        Rtotal = beta2./alpha - Ky.*(KM./B+1);
        Atotal = beta1./alpha./(1+Rtotal./Kx.*B./(KM+B));
    case 'IE_AE'
        fun = @IE_AE;
    case 'null_NAR'
        Atotal = beta2/alpha;
        Rtotal = sqrt(beta1./alpha.*Ky.*(KM+B)./B);
    case 'null_IE'
        Atotal = beta2/alpha;
        Rtotal = beta1./alpha./(1+Atotal./Ky.*A./(KG+A));
    case 'AE_NAR'
        Rtotal = sqrt(beta1./alpha.*Ky.*(KM+B)./B);
        Atotal = beta2./alpha./(1+Kx.*(KM+B)./B./Rtotal);
    case 'PAR_NAR'
        Atotal = beta2./alpha - Kx.*(KG./A+1);
        Rtotal = sqrt(beta1./alpha.*Ky.*(KM+B)./B);
    case 'PAR_IE'
        Atotal = beta2./alpha - Kx.*(KG./A+1);
        Rtotal = beta1./alpha./(1+Atotal./Ky.*A./(KG+A));
    case 'AE_IE'
        fun = @AE_IE;
end

if ismember(topo,{'IE_AE','AE_IE'})
    x0 = [0,0];
    for i=1:50
        for j=1:50
            f = @(x) fun(x,A(i,j),B(i,j));
            [out,~] = fsolve(f,x0);
            Atotal(i,j) = out(1);
            Rtotal(i,j) = out(2);
        end
    end
end

ind_level = 1./(1+KA./Atotal.*(KG./A+1).*(1+Rtotal./KR.*B./(KM+B)));
end

function F=IE_AE(x,SA,SR)
beta1 = 1000;
beta2 = 10;
alpha = 0.001;
Kx = 1;
Ky = 30;
KM = 10^9;
KG = 10^9;
F(1) = beta1./(1+SR.*x(2)./(Kx.*(KM+SR)))-alpha.*x(1);  % dynamics of activator
F(2) = beta2./(1+Ky./x(1).*(1+KG./SA))-alpha.*x(2); % dynamics of repressor
end

function F=AE_IE(x,SA,SR)
beta1 = 10;
beta2 = 1000;
alpha = 0.001;
Kx = 30;
Ky = 1;
KM = 10^9;
KG = 10^9;
F(1) = beta1./(1+Kx./x(2).*(1+KM./SR))-alpha.*x(1); % dynamics of activator
F(2) = beta2./(1+SA.*x(1)./(Ky.*(KG+SA)))-alpha.*x(2);  % dynamics of repressor
end
