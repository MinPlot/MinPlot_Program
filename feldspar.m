%%feldspar 

function [StrctFrm]=feldspar(data,headers)

[m,n]=size(data); %finds the x and y size of the input data matrix

I(1,1)=find(strcmp(headers,'SiO2'));
I(1,2)=find(strcmp(headers,'Al2O3'));

%Makes FeO optional
if strcmp(headers,'FeO')==zeros(1,length(headers))
    I(1,3)=0;
else
    I(1,3)=find(strcmp(headers,'FeO'));
end

%Makes MnO optional
if strcmp(headers,'MnO')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'MnO'));
end

%Makes MgO optional
if strcmp(headers,'MgO')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'MgO'));
end

I(1,6)=find(strcmp(headers,'CaO'));
I(1,7)=find(strcmp(headers,'Na2O'));
I(1,8)=find(strcmp(headers,'K2O'));

%Makes BaO optional
if strcmp(headers,'BaO')==zeros(1,length(headers))
    I(1,9)=0;
else
    I(1,9)=find(strcmp(headers,'BaO'));
end

Opfu=8.0; %oxygens per formula unit


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

%% moles of cations


MC(:,1)=data(:,I(1,1))./SiO2; %for SiO2
MC(:,2)=(data(:,I(1,2))./Al2O3)*2; %for Al2O3

%calculates for FeO if it is included in the analysis 
if I(1,3)==0
    MC(:,3)=zeros(m,1); 
else
    MC(:,3)=data(:,I(1,3))./FeO; %for FeO
end 

%calculates for MnO if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1); 
else
    MC(:,4)=data(:,I(1,4))./MnO; %for MnO
end 

%calculates for MgO if it is included in the analysis 
if I(1,5)==0
    MC(:,5)=zeros(m,1); 
else
    MC(:,5)=data(:,I(1,5))./MgO; %for MgO
end 

MC(:,6)=data(:,I(1,6))./CaO; %for CaO
MC(:,7)=(data(:,I(1,7))./Na2O).*2; %for Na2O
MC(:,8)=(data(:,I(1,8))./K2O).*2; %for K2O

%calculates for BaO if it is included in the analysis 
if I(1,9)==0
    MC(:,9)=zeros(m,1); 
else
    MC(:,9)=data(:,I(1,9))./BaO; %for BaO
end 

%% Oxygen Units
O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*(3/2); %for Al2O3
O2(:,3)=MC(:,3); %for FeO
O2(:,4)=MC(:,4); %for MnO
O2(:,5)=MC(:,5); %for MgO
O2(:,6)=MC(:,6); %for CaO
O2(:,7)=MC(:,7)*(1/2); %for Na2O
O2(:,8)=MC(:,8)*(1/2); %for K2O
O2(:,9)=MC(:,9); %for BaO

O2total=sum(O2,2); %O2 totals

O2total=sum(O2,2); %O2 totals
MCnormfact=Opfu./sum(O2,2); %normalization factor

%% Structural Formula
StrctFrm=MCnormfact.*MC; %creates a matrix of normalized cations
StrctFrm(:,10)=sum(StrctFrm,2); %calculations the total, which should be close to 5

%% endmember calculation

Endmembers(:,1)=StrctFrm(:,6)./(StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)); %anorthite 
Endmembers(:,2)=StrctFrm(:,7)./(StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)); %albite
Endmembers(:,3)=StrctFrm(:,8)./(StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)); %orthoclase


%% Ternary plot 

%prompts the user if they wish to plot the CPX data
prompt1='Do you wish to plot the feldspar ternary? (y|n): ';
wantplots=input(prompt1, 's');

if strcmp(wantplots, 'y')
    
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

    figure('Name','Feldspar Ternary');
    
    %makes white background
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
    
    %labels
    text(-0.10,0.0,'ab','FontSize',14)
    text(0.51,0.92,'an','FontSize',14)
    text(1.02,-0.04,'or','FontSize',14)
    
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
    
    %prompts the user to plot the subdivisions
    prompt5='Do you want to plot the feldspar compositional subdivisionns? (y|n):';
    feld=input(prompt5,'s');
    
    if strcmp(feld,'y')
    
        %boundaries of the feldspar subdivisions
        plot([0.05 0.145],[0.08660254 0.077942286],'k','linewidth',1.5)
        hold on
        plot([0.1 0.145],[0 0.077942286],'k','linewidth',1.5)
        hold on
        plot([0.145 0.225],[0.077942286 0.129903811],'k','linewidth',1.5)
        hold on
        plot([0.15 0.235],[0.259807621 0.233826859],'k','linewidth',1.5)
        hold on
        plot([0.225 0.235],[0.129903811 0.233826859],'k','linewidth',1.5)
        hold on
        plot([0.25 0.325],[0.433012702 0.389711432],'k','linewidth',1.5)
        hold on
        plot([0.235 0.325],[0.233826859 0.389711432],'k','linewidth',1.5)
        hold on
        plot([0.35 0.415],[0.606217783 0.545596004],'k','linewidth',1.5)
        hold on
        plot([0.325 0.415],[0.389711432 0.545596004],'k','linewidth',1.5)
        hold on
        plot([0.45 0.505],[0.779422863 0.701480577],'k','linewidth',1.5)
        hold on
        plot([0.415 0.505],[0.545596004 0.701480577],'k','linewidth',1.5)
        hold on
        plot([0.55 0.505],[0.779422863 0.701480577],'k','linewidth',1.5)
        hold on
        plot([0.37 0.383],[0 0.08660254],'k','linewidth',1.5)
        hold on
        plot([0.383 0.225],[0.08660254 0.129903811],'k','linewidth',1.5)
        hold on
        plot([0.383 0.95],[0.08660254 0.08660254],'k','linewidth',1.5)
        hold on
        text(0.24,0.18,'olig','FontSize',12)
        text(0.29,0.30,'ands','FontSize',12)
        text(0.38,0.46,'labr','FontSize',12)
        text(0.48,0.62,'bytw','FontSize',12)
        text(0.30,0.14,'ano','FontSize',12)
        text(0.62,0.11,'mc/or/sa','FontSize',12)
    end

    %transforms the data to ternary space
    X1=0.5.*(Endmembers(:,1))+(Endmembers(:,3));
    Y1=(Endmembers(:,1))*(cos(30*pi()/180));
    
    
    scatter(X1(:),Y1(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
    hold off
            
    axis image
    axis off
end

all=[StrctFrm Endmembers];
StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe2_M1','Mn_M1','Mg_M1','Ca_M1','Na_M1','K_M1','Ba_M1','Sum','Xan','Xab','Xor'}); 
end
