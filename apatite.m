%apatite structural formula 
function [StrctFrm,APFU]=apatite(data,headers,wantstrctfrm)

%choose is sulfur included?
prompt1 = 'Is sulfur analyzed? (y|n): ';
wantS = input(prompt1, 's');

%choose the speciation of sulfur:
if strcmp(wantS, 'y') 
    prompt2 = 'Enter S6+/S, S4+/S, S1-/S, and S2-/S (e.g, [0.6; 0.2; 0.2; 0.0]): ';
    Sratio = input(prompt2);
else
    Sratio = [0;0;0;0];
end

[m,n]=size(data); %finds the x and y size of the input data matrix

%finds the column position oxide headers in any order to assign data to the correct array
%positions 

%Makes SiO2 optional
if strcmp(headers,'SiO2')==zeros(1,length(headers))
    I(1,1)=0;
else
    I(1,1)=find(strcmp(headers,'SiO2'));
end

I(1,2)=find(strcmp(headers,'P2O5'));

%Makes TiO2 optional
if strcmp(headers,'TiO2')==zeros(1,length(headers))
    I(1,3)=0;
else
    I(1,3)=find(strcmp(headers,'TiO2'));
end

%Makes Al2O3 optional
if strcmp(headers,'Al2O3')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'Al2O3'));
end

%Makes FeO optional
if strcmp(headers,'FeO')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'FeO'));
end

%Makes MnO optional
if strcmp(headers,'MnO')==zeros(1,length(headers))
    I(1,6)=0;
else
    I(1,6)=find(strcmp(headers,'MnO'));
end

%Makes MgO optional
if strcmp(headers,'MgO')==zeros(1,length(headers))
    I(1,7)=0;
else
    I(1,7)=find(strcmp(headers,'MgO'));
end

I(1,8)=find(strcmp(headers,'CaO'));

%Makes BaO optional
if strcmp(headers,'BaO')==zeros(1,length(headers))
    I(1,9)=0;
else
    I(1,9)=find(strcmp(headers,'BaO'));
end

%Makes SrO optional
if strcmp(headers,'SrO')==zeros(1,length(headers))
    I(1,10)=0;
else
    I(1,10)=find(strcmp(headers,'SrO'));
end

%Makes K2O optional
if strcmp(headers,'K2O')==zeros(1,length(headers))
    I(1,11)=0;
else
    I(1,11)=find(strcmp(headers,'K2O'));
end

%Makes Na2O optional
if strcmp(headers,'Na2O')==zeros(1,length(headers))
    I(1,12)=0;
else
    I(1,12)=find(strcmp(headers,'Na2O'));
end

%Makes Ce2O3 optional
if strcmp(headers,'Ce2O3')==zeros(1,length(headers))
    I(1,13)=0;
else
    I(1,13)=find(strcmp(headers,'Ce2O3'));
end

%Makes La2O3 optional
if strcmp(headers,'La2O3')==zeros(1,length(headers))
    I(1,14)=0;
else
    I(1,14)=find(strcmp(headers,'La2O3'));
end

%Makes SO3 optional
if strcmp(headers,'SO3')==zeros(1,length(headers))
    I(1,15)=0;
else
    I(1,15)=find(strcmp(headers,'SO3'));
end

I(1,16)=find(strcmp(headers,'F'));
I(1,17)=find(strcmp(headers,'Cl'));

%% Molecular weights

SiO2=60.083;
TiO2=79.865;
Al2O3=101.961;
Cr2O3=151.989;
Fe2O3=159.6874;
Y2O3=225.809;
Ce2O3=328.229;
La2O3=325.80794;
NiO=74.692;
ZnO=81.381;
FeO=71.8442;
MnO=70.937;
MgO=40.304;
CaO=56.0774;
Na2O=61.979;
K2O=94.195;
BaO=153.329;
SrO=103.619;
SO3=80.0594;
F=18.998;
Cl=35.45;
P2O5=141.942524;

W=[SiO2,P2O5,TiO2,Al2O3,FeO,MnO,MgO,CaO,BaO,SrO,K2O,Na2O,Ce2O3,La2O3,SO3,F,Cl];

%% Moles of oxides

%calculates for SiO2 if it is included in the analysis 
if I(1,1)==0
    MC(:,1)=zeros(m,1);
else
    MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2
end 

MC(:,2)=data(:,I(1,2))./W(:,2); %for P2O5

%calculates for TiO2 if it is included in the analysis 
if I(1,3)==0
    MC(:,3)=zeros(m,1);
else
    MC(:,3)=data(:,I(1,3))./W(:,3); %for TiO2
end 

%calculates for Al2O3 if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=data(:,I(1,4))./W(:,4); %for Al2O3
end 

%calculates for FeO if it is included in the analysis 
if I(1,5)==0
    MC(:,5)=zeros(m,1);
else
    MC(:,5)=data(:,I(1,5))./W(:,5); %for FeO
end 

%calculates for MnO if it is included in the analysis 
if I(1,6)==0
    MC(:,6)=zeros(m,1);
else
    MC(:,6)=data(:,I(1,6))./W(:,6); %for MnO
end 

%calculates for MgO if it is included in the analysis 
if I(1,7)==0
    MC(:,7)=zeros(m,1);
else
    MC(:,7)=data(:,I(1,7))./W(:,7); %for MgO
end 


MC(:,8)=data(:,I(1,8))./W(:,8); %for CaO


%calculates for BaO if it is included in the analysis 
if I(1,9)==0
    MC(:,9)=zeros(m,1);
else
    MC(:,9)=data(:,I(1,9))./W(:,9); %for BaO
end 

%calculates for SrO if it is included in the analysis 
if I(1,10)==0
    MC(:,10)=zeros(m,1);
else
    MC(:,10)=data(:,I(1,10))./W(:,10); %for SrO
end 

%calculates for K2O if it is included in the analysis 
if I(1,11)==0
    MC(:,11)=zeros(m,1);
else
    MC(:,11)=data(:,I(1,11))./W(:,11); %for K2O
end 

%calculates for Na2O if it is included in the analysis 
if I(1,12)==0
    MC(:,12)=zeros(m,1);
else
    MC(:,12)=data(:,I(1,12))./W(:,12); %for Na2O
end 

%calculates for Ce2O3 if it is included in the analysis 
if I(1,13)==0
    MC(:,13)=zeros(m,1);
else
    MC(:,13)=data(:,I(1,13))./W(:,13); %for Ce2O3
end 

%calculates for La2O3 if it is included in the analysis 
if I(1,14)==0
    MC(:,14)=zeros(m,1);
else
    MC(:,14)=data(:,I(1,14))./W(:,14); %for La2O3
end 

%calculates for SO3 if it is included in the analysis 
if I(1,15)==0
    MC(:,15)=zeros(m,1);
else
    MC(:,15)=data(:,I(1,15))./W(:,15); %for SO3
end 

MC(:,16)=data(:,I(1,16))./W(:,16); %for F
MC(:,17)=data(:,I(1,17))./W(:,17); %for Cl

%% Moles of O2 units
O2(:,1)=MC(:,1).*2; %SiO2
O2(:,2)=MC(:,2).*5; %P2O5
O2(:,3)=MC(:,3).*2; %TiO2
O2(:,4)=MC(:,4).*3; %Al2O3
O2(:,5)=MC(:,5); %FeO
O2(:,6)=MC(:,6); %MnO
O2(:,7)=MC(:,7); %MgO
O2(:,8)=MC(:,8); %CaO
O2(:,9)=MC(:,9); %BaO
O2(:,10)=MC(:,10); %SrO
O2(:,11)=MC(:,11); %K2O
O2(:,12)=MC(:,12); %Na2O
O2(:,13)=MC(:,13).*3; %Ce2O3
O2(:,14)=MC(:,14).*3; %La2O3
O2(:,15)=(MC(:,15).*Sratio(1,1)).*3; %S6+
O2(:,16)=(MC(:,15).*Sratio(2,1)).*2; %S4+
O2(:,17)=MC(:,15).*Sratio(3,1); %S1-
O2(:,18)=MC(:,15).*Sratio(4,1); %S2-
O2(:,19)=MC(:,16); %F
O2(:,20)=MC(:,17); %Cl

O2_Sum=sum(O2(:,1:16),2); %sum of O2, excluding F, S1-, S2-, and Cl

O2_N=(25)./O2_Sum; %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

%% atoms pfu

APFU(:,1)=N_Ox(:,1)./2; %Si
APFU(:,2)=N_Ox(:,2).*(2/5); %P
APFU(:,3)=N_Ox(:,3)./2; %Ti
APFU(:,4)=N_Ox(:,4).*(2/3); %Al
APFU(:,5)=N_Ox(:,5); %Fe
APFU(:,6)=N_Ox(:,6); %Mn
APFU(:,7)=N_Ox(:,7); %Mg
APFU(:,8)=N_Ox(:,8); %Ca
APFU(:,9)=N_Ox(:,9); %Ba
APFU(:,10)=N_Ox(:,10); %Sr
APFU(:,11)=N_Ox(:,11).*2; %K
APFU(:,12)=N_Ox(:,12).*2; %Na
APFU(:,13)=N_Ox(:,13).*(2/3); %Ce
APFU(:,14)=N_Ox(:,14).*(2/3); %La
APFU(:,15)=N_Ox(:,15)./3; %S6+
APFU(:,16)=N_Ox(:,16)./2; %S4+
APFU(:,17)=N_Ox(:,17); %S1-
APFU(:,18)=N_Ox(:,18); %S2-
APFU(:,19)=N_Ox(:,19); %F
APFU(:,20)=N_Ox(:,20); %Cl
APFU(:,21)=2-(N_Ox(:,20)+N_Ox(:,19)+N_Ox(:,17)+(2.*N_Ox(:,18))); %H = 2 - (F + Cl + S1- + 2 * S2-)

%% Structural Formula

%T site
StrctFrm(:,1)=APFU(:,2); %P (T)
StrctFrm(:,2)=APFU(:,1); %Si (T)
StrctFrm(:,3)=APFU(:,15); %S6+ (T)
StrctFrm(:,4)=APFU(:,16); %S4+ (T)

%M sites
StrctFrm(:,5)=APFU(:,3); %Ti (M)
StrctFrm(:,6)=APFU(:,4); %Al (M)
StrctFrm(:,7)=APFU(:,5); %Fe (M)
StrctFrm(:,8)=APFU(:,6); %Mn (M)
StrctFrm(:,9)=APFU(:,7); %Mg (M)
StrctFrm(:,10)=APFU(:,8); %Ca (M)
StrctFrm(:,11)=APFU(:,9); %Ba (M)
StrctFrm(:,12)=APFU(:,10); %Sr (M)
StrctFrm(:,13)=APFU(:,11); %K (M)
StrctFrm(:,14)=APFU(:,12); %Na (M)
StrctFrm(:,15)=APFU(:,13); %Ce (M)
StrctFrm(:,16)=APFU(:,14); %La (M)
StrctFrm(:,17)=sum(StrctFrm(:,1:16),2); %cation sum

%Z site
StrctFrm(:,18)=APFU(:,21); %OH
StrctFrm(:,19)=APFU(:,19); %F
StrctFrm(:,20)=APFU(:,20); %Cl
StrctFrm(:,21)=APFU(:,17); %S1-
StrctFrm(:,22)=APFU(:,18); %S2-
StrctFrm(:,23)=sum(StrctFrm(:,18:22),2); %anion sum

%variables for the ternary plot
OH = (APFU(:,21) + APFU(:,17) + APFU(:,18))./StrctFrm(:,23); %OH + S 
F = APFU(:,19)./StrctFrm(:,23); %F

StrctFrm=array2table(StrctFrm,'VariableNames',{'P_T','Si_T','S6_T','S4_T','Ti_M','Al_M','Fe_M' ...
    ,'Mn_M','Mg_M','Ca_M','Ba_M','Sr_M','K_M','Na_M','Ce_M','La_M','Cation_Sum','OH_Z', ...
    'F_Z','Cl_Z','S1-_Z','S2-_Z','Anion Sum'});

APFU=array2table(APFU,'VariableNames',{'Si','P','Ti','Al','Fe','Mn','Mg','Ca','Ba' ...
    ,'Sr','K','Na','Ce','La','S6','S4','S1-','S2-','F','Cl','OH'});

if strcmp(wantstrctfrm, 'y')
    
    %prompts the user if they wish to plot the mica data
    prompt1='Do you wish to plot the OH-F-Cl ternary? (y|n): ';
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

        figure('Name','Cl-F-OH Ternary');
        pgon=polyshape([0 0.5 1],[0 sqrt(3)/2 0]);
        plot(pgon,'FaceColor','w')
        hold on
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
        plot([0 1],[0 0],'k','linewidth',1.5)
        hold on
        plot([0 0.5],[0 sqrt(3)/2],'k','linewidth',1.5)
        hold on
        plot([0.5 1],[sqrt(3)/2 0],'k','linewidth',1.5)
        hold on

        %labels
        text(-0.10,0.0,'Cl','FontSize',14)
        text(0.50,0.92,'F','FontSize',14)
        text(1.02,-0.04,'OH + S','FontSize',14)

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
        X3=0.5.*(F)+(OH);
        Y3=(F)*(cos(30*pi()/180));

        scatter(X3(:),Y3(:),symbsize,symb,'filled','MarkerFaceAlpha',2/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
        hold off

        axis image
        axis off
    end
end



end
















