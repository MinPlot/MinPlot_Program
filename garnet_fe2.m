%%Garnet Structural Formula

function [StrctFrm, APFU]=garnet_fe2(data,headers,wantstrctfrm)

[m,n]=size(data); %finds the x and y size of the input data matrix

%finds the column position oxide headers in any order to assign data to the correct array
%positions 
I(1,1)=find(strcmp(headers,'SiO2'));

%makes TiO2 optional
if strcmp(headers,'TiO2')==zeros(1,length(headers))
    I(1,2)=0;
else
    I(1,2)=find(strcmp(headers,'TiO2'));
end

I(1,3)=find(strcmp(headers,'Al2O3'));

%Makes Cr2O3 optional
if strcmp(headers,'Cr2O3')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'Cr2O3'));
end

%Makes Y2O3 optional
if strcmp(headers,'Y2O3')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'Y2O3'));
end

I(1,6)=find(strcmp(headers,'FeO'));
I(1,7)=find(strcmp(headers,'MnO'));
I(1,8)=find(strcmp(headers,'MgO'));
I(1,9)=find(strcmp(headers,'CaO'));

%Makes Na2O optional
if strcmp(headers,'Na2O')==zeros(1,length(headers))
    I(1,10)=0;
else
    I(1,10)=find(strcmp(headers,'Na2O'));
end

Opfu=12.0; %oxygens per formula unit

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

W=[SiO2,TiO2,Al2O3,Cr2O3,Y2O3,FeO,MnO,MgO,CaO,Na2O];

%% Calculate cations units


MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%adds a column of zeros if Ti is not included in the calculation
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

MC(:,3)=(data(:,I(1,3))./W(:,3)).*2; %for Al2O3

%adds a column of zeros if Cr is not included in the calculation
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=(data(:,I(1,4))./W(:,4)).*2; %for Cr2O3
end
%adds a column of zeros if Y is not included in the calculation
if I(1,5)==0
    MC(:,5)=zeros(m,1);
else
    MC(:,5)=(data(:,I(1,5))./W(:,5)).*2; %for Y2O3
end

MC(:,6)=data(:,I(1,6))./W(:,6); %for FeO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MnO
MC(:,8)=data(:,I(1,8))./W(:,8); %for MgO
MC(:,9)=data(:,I(1,9))./W(:,9); %for CaO

%adds a column of zeros if Na is not included in the calculation
if I(1,10)==0
    MC(:,10)=zeros(m,1);
else
    MC(:,10)=(data(:,I(1,10))./W(:,10)).*2; %for Na2O
end

%% Calculate Oxygen Units

O2=zeros(size(data));
O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4).*(3/2); %for Cr2O3
O2(:,5)=MC(:,5).*(3/2); %for Y2O3
O2(:,6)=MC(:,6); %for FeO
O2(:,7)=MC(:,7); %for MnO
O2(:,8)=MC(:,8); %for MgO
O2(:,9)=MC(:,9); %for CaO
O2(:,10)=MC(:,10)./2; %for Na2O


O2total=sum(O2,2); %O2 totals

MCnormfact=Opfu./sum(O2,2); %normalization factor

%% Atoms pfu

APFU=MCnormfact.*MC; %creates a matrix of normalized cations
APFU(:,11)=sum(APFU,2); %calculations the total, which should be close to 8

%% structural formula calculation

%T SITE
%Si 
for c=1:m
    if APFU(c,1)<3.000
        StrctFrm(c,1)=APFU(c,1); %If Si < 3, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=3; %If Si is in excess, then Si(T) = 3
    end
end

%Al(T)
for c=1:m
    if 3-StrctFrm(c,1)>0 %Is 3-Si > 0? If y, then some Al goes into T
        if 3-StrctFrm(c,1)>APFU(c,3) %For low Al Grt, 3-Si may be > Al
            StrctFrm(c,2)=APFU(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=3-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to Y
        end
    else
        StrctFrm(c,2)=0; %if Si=3, then no Al goes into T
    end
end


%Sum of T site
StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2);

%Y SITE

%Al (Y)
StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2);

%Ti (Y)
StrctFrm(:,5)=APFU(:,2);

%Cr (Y) 
StrctFrm(:,6)=APFU(:,4);

%Y sum
StrctFrm(:,7)=StrctFrm(:,6)+StrctFrm(:,5)+StrctFrm(:,4);

%X SITE

%Y (X)
StrctFrm(:,8)=APFU(:,5);

%Fe2+ (X)
StrctFrm(:,9)=APFU(:,6);

%Mn (X)
StrctFrm(:,10)=APFU(:,7);

%Mg (X)
StrctFrm(:,11)=APFU(:,8);

%Ca (X)
StrctFrm(:,12)=APFU(:,9);

%Na (X)
StrctFrm(:,13)=APFU(:,10);

% X Sum
StrctFrm(:,14)=StrctFrm(:,13)+StrctFrm(:,12)+StrctFrm(:,11)+StrctFrm(:,10)+StrctFrm(:,9)+StrctFrm(:,8);

%% end member calculations



Endmembers(:,1)=APFU(:,6)./(APFU(:,6)+APFU(:,7)+APFU(:,8)+APFU(:,9)); % XAlm
Endmembers(:,2)=APFU(:,8)./(APFU(:,6)+APFU(:,7)+APFU(:,8)+APFU(:,9)); % XPrp
Endmembers(:,3)=APFU(:,7)./(APFU(:,6)+APFU(:,7)+APFU(:,8)+APFU(:,9)); % XSps
Endmembers(:,4)=APFU(:,9)./(APFU(:,6)+APFU(:,7)+APFU(:,8)+APFU(:,9)); % XGrs

if strcmp(wantstrctfrm, 'y')

    %prompts the user if they wish to plot the garnet data
    prompt2='Do you wish to plot garnet compositions? (y|n): ';
    wantplots=input(prompt2, 's');

    if strcmp(wantplots, 'y')

        %prompts the user to determine which symbols to use
        prompt3='What symbols do you want to use?:';
        disp('Options are (CASE SENSITIVE): circle, square, diamond, and triangle.') %for simplicity only 4 options are available
        wantsymbols=input(prompt3,'s');

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
        prompt4='Specify the fill color:';
        disp('Options are (CASE SENSITIVE): blue, orange, yellow, purple, green, cyan, & red.')
        wantfil=input(prompt4,'s');

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
        prompt5='Specify symbol size (numeric scalar):';
        disp('Note: Between 50 & 200 is good for most applications.')
        symbsize=input(prompt5);

        %plots a ternary for Alm+Sps-Prp-Grs
        figure('Name','Alm+Sps-Prp-Grs garnet')

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

        %labels
        text(-0.14,0.02,'Alm','FontSize',14)
        text(-0.16,-0.02,'+Sps','FontSize',14)
        text(0.51,0.92,'Grs','FontSize',14)
        text(1.02,-0.04,'Prp','FontSize',14) 

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

        %transforms the data to ternary space
        X2=0.5.*(Endmembers(:,4))+(Endmembers(:,2));
        Y2=(Endmembers(:,4))*(cos(30*pi()/180));

        scatter(X2(:),Y2(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
        hold off

        axis image
        axis off

    end
end

all=[StrctFrm Endmembers];
StrctFrm=array2table(all,'VariableNames',{'Si','Aliv','Sum_Tet','Alvi','Ti','Cr','Sum_Oct','Y','Fe2','Mn','Mg','Ca','Na','Sum_Dod','Xalm','Xprp','Xsps','Xgrs'});

APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Cr','Y','Fe2','Mn','Mg','Ca','Na','Total'});
end


