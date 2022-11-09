%% Spinel structural formula

function [APFU]=spinel_fe3(data,headers)

[m,n]=size(data); %finds the x and y size of the input data matrix

I(1,1)=find(strcmp(headers,'TiO2'));
I(1,2)=find(strcmp(headers,'Al2O3'));
I(1,3)=find(strcmp(headers,'Cr2O3'));

%makes NiO optional
if strcmp(headers,'NiO')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'NiO'));
end

%makes ZnO optional
if strcmp(headers,'ZnO')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'ZnO'));
end

I(1,6)=find(strcmp(headers,'FeO'));
I(1,7)=find(strcmp(headers,'MnO'));
I(1,8)=find(strcmp(headers,'MgO'));

cat=3.0; %cations per formula unit
Opfu=4.0; %oxygens per formula unit

%% Molecular weights

SiO2=60.083;
TiO2=79.865;
Al2O3=101.961;
Cr2O3=151.989;
Fe2O3=159.6874;
Y2O3=225.809;
NiO=74.692;
ZnO=81.381;
FeO=71.8442;
MnO=70.937;
MgO=40.304;
CaO=56.0774;
Na2O=61.979;
K2O=94.195;
BaO=153.329;
F=18.998;
Cl=35.45;

W=[TiO2,Al2O3,Cr2O3,NiO,ZnO,FeO,MnO,MgO];

%% Calculate cations units
MC(:,1)=data(:,I(1,1))./W(:,1); %for TiO2
MC(:,2)=(data(:,I(1,2))./W(:,2)).*2; %for Al2O3
MC(:,3)=(data(:,I(1,3))./W(:,3)).*2; %for Cr2O3

%adds a column of zeros if NiO is not included in the calculation
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=data(:,I(1,4))./W(:,4); %for NiO
end

%adds a column of zeros if ZnO is not included in the calculation
if I(1,5)==0
    MC(:,5)=zeros(m,1);
else
    MC(:,5)=data(:,I(1,5))./W(:,5); %for ZnO
end

MC(:,6)=data(:,I(1,6))./W(:,6); %for FeO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MnO
MC(:,8)=data(:,I(1,8))./W(:,8); %for MgO

MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cation

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for TiO2
O2(:,2)=MCnorm(:,2).*(3/2); %for Al2O3
O2(:,3)=MCnorm(:,3).*(3/2); %for Cr2O3
O2(:,4)=MCnorm(:,4); %for NiO
O2(:,5)=MCnorm(:,5); %for ZnO
O2(:,6)=MCnorm(:,6); %for FeO
O2(:,7)=MCnorm(:,7); %for MnO
O2(:,8)=MCnorm(:,8); %for MgO

O2total=sum(O2,2); %O2 totals

%% Atoms PFU

APFU(:,1)=MCnorm(:,1); %for Ti
APFU(:,2)=MCnorm(:,2); %for Al
APFU(:,3)=MCnorm(:,3); %for Cr
APFU(:,5)=MCnorm(:,4); %for Ni
APFU(:,6)=MCnorm(:,5); %for Zn
APFU(:,8)=MCnorm(:,7); %for Mn
APFU(:,9)=MCnorm(:,8); %for Mg

%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 4
%if so, then there is no Fe3+
%if totalO2 < 4, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(4-totalO2) then the amount
%of Fe3+ = 2*(4-totalO2), if false then, all Fe is Fe3+

for c=1:m
    if (Opfu-O2total(c,1)) >= 0
        if MCnorm(c,6) > 2*(Opfu-O2total(c,1))
            APFU(c,4)=2*(Opfu-O2total(c,1)); 
        else
            APFU(c,4)=MCnorm(c,6);
        end
    else
        APFU(c,4)=0;
    end
end

APFU(:,7)=MCnorm(:,6)-APFU(:,4); %the APFU of Fe2+ equals totalFe-Fe3+

APFU(:,10)=sum(APFU,2); %calculations the total, which should be 3

% Oxygen deficiency 
APFU(:,11)=Opfu-O2total; %must be greater than zero

%% plots

prompt1='Do you wish to plot ternary diagrams? (y|n): ';
wantplot=input(prompt1, 's');

if strcmp(wantplot, 'y')

    %prompts the user to determine which symbols to use
    prompt2='What symbols do you want to use?:';
    disp('Options are (CASE SENSITIVE): circle, square, diamond, and triangle.') %for simplicity only 4 options are available
    wantsymbols=input(prompt2,'s');

    %assigns the variable the appropriate symbol marker
    if strcmp(wantsymbols,'circle')
        symb='o';
    end

    if strcmp(wantsymbols,'square')
        symb='s';
    end

    if strcmp(wantsymbols,'diamond')
        symb='d';
    end

    if strcmp(wantsymbols,'triangle')
        symb='^';
    end

    %prompts the user to determine which symbol fill color to use
    prompt3='Specify the fill color:';
    disp('Options are (CASE SENSITIVE): blue, orange, yellow, purple, green, cyan, & red.')
    wantfil=input(prompt3,'s');

    %assigns the variable the appropriate fill color
    if strcmp(wantfil,'blue')
        fil=[0 0.4470 0.7410];
    end

    if strcmp(wantfil,'orange')
        fil=[0.8500 0.3250 0.0980];
    end

    if strcmp(wantfil,'yellow')
        fil=[0.9290 0.6940 0.1250];
    end

    if strcmp(wantfil,'purple')
        fil=[0.4940 0.1840 0.5560];
    end

    if strcmp(wantfil,'green')
        fil=[0.4660 0.6740 0.1880];
    end

    if strcmp(wantfil,'cyan')
        fil=[0.3010 0.7450 0.9330];
    end

    if strcmp(wantfil,'red')
        fil=[0.6350 0.0780 0.1840];
    end

    %prompts the user to determine which symbol fill color to use
    prompt4='Specify symbol size (numeric scalar):';
    disp('Note: Between 50 & 200 is good for most applications.')
    symbsize=input(prompt4);

    %plots a ternary for Al-Fe3-Cr
    figure('Name','Al-Fe3-Cr Plot')

    %plot a grid intervals of 0.2 for different endmembers
    pgon=polyshape([0 0.5 1],[0 sqrt(3)/2 0]);
    plot(pgon,'FaceColor','w')
    hold on

    %plot grid (spacing of 10 %)
    plot([0.45 0.55], [0.779422863 0.779422863],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.9
    hold on
    plot([0.4 0.6], [0.692820323 0.692820323],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.8
    hold on
    plot([0.35 0.65], [0.606217783 0.606217783],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.7
    hold on
    plot([0.3 0.7], [0.519615242 0.519615242],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.6
    hold on
    plot([0.25 0.75], [0.433012702 0.433012702],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.5
    hold on
    plot([0.2 0.8], [0.346410162 0.346410162],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.4
    hold on
    plot([0.15 0.85], [0.259807621 0.259807621],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.3
    hold on
    plot([0.1 0.9], [0.173205081 0.173205081],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.2
    hold on
    plot([0.05 0.95], [0.08660254 0.08660254],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.1
    hold on
    plot([0.05 0.1], [0.08660254 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.9
    hold on
    plot([0.1 0.2], [0.173205081 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.8
    hold on
    plot([0.15 0.3], [0.259807621 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.7
    hold on
    plot([0.2 0.4], [0.346410162 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.6
    hold on
    plot([0.25 0.5], [0.433012702 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.5
    hold on
    plot([0.3 0.6], [0.519615242 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.4
    hold on
    plot([0.35 0.7], [0.606217783 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.3
    hold on
    plot([0.4 0.8], [0.692820323 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.2
    hold on
    plot([0.45 0.9], [0.779422863 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.1
    hold on
    plot([0.95 0.9], [0.08660254 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.9
    hold on
    plot([0.9 0.8], [0.173205081 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.8
    hold on
    plot([0.85 0.7], [0.259807621 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.7
    hold on
    plot([0.8 0.6], [0.346410162 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.6
    hold on
    plot([0.75 0.5], [0.433012702 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.5
    hold on
    plot([0.7 0.4], [0.519615242 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.4
    hold on
    plot([0.65 0.3], [0.606217783 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.3
    hold on
    plot([0.6 0.2], [0.692820323 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.2
    hold on
    plot([0.55 0.1], [0.779422863 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.1
    hold on

    %plot boundaries of the triangle
    plot([0 1],[0 0],'k','linewidth',1.5)
    hold on
    plot([0 0.5],[0 sqrt(3)/2],'k','linewidth',1.5)
    hold on
    plot([0.5 1],[sqrt(3)/2 0],'k','linewidth',1.5)
    hold on

    %tick labels
    text(-0.02,-0.039,'0.0','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.18,-0.039,'0.2','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.38,-0.039,'0.4','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.58,-0.039,'0.6','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.78,-0.039,'0.8','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.98,-0.039,'1.0','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(1.02,0.01,'0.0','FontSize',12)
    text(0.92,0.18,'0.2','FontSize',12)
    text(0.82,0.35,'0.4','FontSize',12)
    text(0.72,0.52,'0.6','FontSize',12)
    text(0.62,0.70,'0.8','FontSize',12)
    text(0.52,0.87,'1.0','FontSize',12)
    text(-0.025,0.040,'1.0','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.075,0.215,'0.8','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.178,0.386,'0.6','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.278,0.558,'0.4','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.382,0.735,'0.2','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.482,0.901,'0.0','FontSize',12,'HorizontalAlignment','center','Rotation',300)

    %labels
    text(-0.10,0.0,'Cr','FontSize',14)
    text(0.51,0.92,'Fe^{3+}+2Ti','FontSize',14)
    text(1.02,-0.04,'Al','FontSize',14)

    %Plot
    Fe3Ti=(2.*APFU(:,1)+APFU(:,4))./(2.*APFU(:,1)+APFU(:,4)+APFU(:,2)+APFU(:,3)); %fraction of Fe3+ and Ti
    Al=(APFU(:,2))./(2.*APFU(:,1)+APFU(:,4)+APFU(:,2)+APFU(:,3)); %fraction of Al

    %transforms the data to ternary space
    X1=0.5.*(Fe3Ti)+(Al);
    Y1=(Fe3Ti)*(cos(30*pi()/180));

    scatter(X1(:),Y1(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
    hold off

    axis image
    axis off

    %plots a ternary for Fe2-Mg-Mn
    figure('Name','Fe2-Mg-Mn Plot')

    %plot a grid intervals of 0.2 for different endmembers
    pgon=polyshape([0 0.5 1],[0 sqrt(3)/2 0]);
    plot(pgon,'FaceColor','w')
    hold on

    %plot grid (spacing of 10 %)
    plot([0.45 0.55], [0.779422863 0.779422863],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.9
    hold on
    plot([0.4 0.6], [0.692820323 0.692820323],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.8
    hold on
    plot([0.35 0.65], [0.606217783 0.606217783],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.7
    hold on
    plot([0.3 0.7], [0.519615242 0.519615242],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.6
    hold on
    plot([0.25 0.75], [0.433012702 0.433012702],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.5
    hold on
    plot([0.2 0.8], [0.346410162 0.346410162],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.4
    hold on
    plot([0.15 0.85], [0.259807621 0.259807621],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.3
    hold on
    plot([0.1 0.9], [0.173205081 0.173205081],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.2
    hold on
    plot([0.05 0.95], [0.08660254 0.08660254],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.1
    hold on
    plot([0.05 0.1], [0.08660254 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.9
    hold on
    plot([0.1 0.2], [0.173205081 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.8
    hold on
    plot([0.15 0.3], [0.259807621 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.7
    hold on
    plot([0.2 0.4], [0.346410162 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.6
    hold on
    plot([0.25 0.5], [0.433012702 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.5
    hold on
    plot([0.3 0.6], [0.519615242 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.4
    hold on
    plot([0.35 0.7], [0.606217783 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.3
    hold on
    plot([0.4 0.8], [0.692820323 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.2
    hold on
    plot([0.45 0.9], [0.779422863 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.1
    hold on
    plot([0.95 0.9], [0.08660254 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.9
    hold on
    plot([0.9 0.8], [0.173205081 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.8
    hold on
    plot([0.85 0.7], [0.259807621 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.7
    hold on
    plot([0.8 0.6], [0.346410162 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.6
    hold on
    plot([0.75 0.5], [0.433012702 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.5
    hold on
    plot([0.7 0.4], [0.519615242 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.4
    hold on
    plot([0.65 0.3], [0.606217783 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.3
    hold on
    plot([0.6 0.2], [0.692820323 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.2
    hold on
    plot([0.55 0.1], [0.779422863 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.1
    hold on

    %plot boundaries of the triangle
    plot([0 1],[0 0],'k','linewidth',1.5)
    hold on
    plot([0 0.5],[0 sqrt(3)/2],'k','linewidth',1.5)
    hold on
    plot([0.5 1],[sqrt(3)/2 0],'k','linewidth',1.5)
    hold on

    %tick labels
    text(-0.02,-0.039,'0.0','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.18,-0.039,'0.2','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.38,-0.039,'0.4','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.58,-0.039,'0.6','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.78,-0.039,'0.8','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(0.98,-0.039,'1.0','FontSize',12,'HorizontalAlignment','center','Rotation',60)
    text(1.02,0.01,'0.0','FontSize',12)
    text(0.92,0.18,'0.2','FontSize',12)
    text(0.82,0.35,'0.4','FontSize',12)
    text(0.72,0.52,'0.6','FontSize',12)
    text(0.62,0.70,'0.8','FontSize',12)
    text(0.52,0.87,'1.0','FontSize',12)
    text(-0.025,0.040,'1.0','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.075,0.215,'0.8','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.178,0.386,'0.6','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.278,0.558,'0.4','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.382,0.735,'0.2','FontSize',12,'HorizontalAlignment','center','Rotation',300)
    text(0.482,0.901,'0.0','FontSize',12,'HorizontalAlignment','center','Rotation',300)

    %labels
    text(-0.25,0.0,'Mn+Zn+Ni','FontSize',14)
    text(0.51,0.92,'Fe^{2+}','FontSize',14)
    text(1.02,-0.04,'Mg','FontSize',14)

    %Plot
    Fe2=APFU(:,7)./(APFU(:,5)+APFU(:,6)+APFU(:,7)+APFU(:,8)+APFU(:,9));
    Mg=APFU(:,9)./(APFU(:,5)+APFU(:,6)+APFU(:,7)+APFU(:,8)+APFU(:,9));

    %transforms the data to ternary space
    X2=0.5.*(Fe2)+(Mg);
    Y2=(Fe2)*(cos(30*pi()/180));

    scatter(X2(:),Y2(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
    hold off

    axis image
    axis off
end

APFU=array2table(APFU,'VariableNames',{'Ti_B','Al_B','Cr_B','Fe3_B','Ni_A','Zn_A','Fe2_A','Mn_A','Mg_A','Sum','O2_deficiency'});

end
            
