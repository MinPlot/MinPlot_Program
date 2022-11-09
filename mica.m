%mica structural formula 
function [StrctFrm, APFU]=mica(data,headers,wantstrctfrm)

[m,n]=size(data); %finds the x and y size of the input data matrix

%finds the column position oxide headers in any order to assign data to the correct array
%positions 
I(1,1)=find(strcmp(headers,'SiO2'));

%Makes TiO2 optional
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

I(1,5)=find(strcmp(headers,'FeO'));
I(1,6)=find(strcmp(headers,'MnO'));
I(1,7)=find(strcmp(headers,'MgO'));

%Makes CaO optional
if strcmp(headers,'CaO')==zeros(1,length(headers))
    I(1,8)=0;
else
    I(1,8)=find(strcmp(headers,'CaO'));
end


I(1,9)=find(strcmp(headers,'Na2O'));
I(1,10)=find(strcmp(headers,'K2O'));

%Makes BaO optional
if strcmp(headers,'BaO')==zeros(1,length(headers))
    I(1,11)=0;
else
    I(1,11)=find(strcmp(headers,'BaO'));
end

%Makes F optional
if strcmp(headers,'F')==zeros(1,length(headers))
    I(1,12)=0;
else
    I(1,12)=find(strcmp(headers,'F'));
end

%Makes Cl optional
if strcmp(headers,'Cl')==zeros(1,length(headers))
    I(1,13)=0;
else
    I(1,13)=find(strcmp(headers,'Cl'));
end

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

W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O,BaO,F,Cl];

%% Moles of oxides

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%calculates for TiO2 if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

MC(:,3)=(data(:,I(1,3))./W(:,3)); %for Al2O3

%calculates for Cr2O3 if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=data(:,I(1,4))./W(:,4); %for Cr2O3
end 

MC(:,5)=(data(:,I(1,5))./W(:,5)); %for FeO
MC(:,6)=(data(:,I(1,6))./W(:,6)); %for MnO
MC(:,7)=(data(:,I(1,7))./W(:,7)); %for MgO

%calculates for CaO if it is included in the analysis 
if I(1,8)==0
    MC(:,8)=zeros(m,1);
else
    MC(:,8)=data(:,I(1,8))./W(:,8); %for CaO
end 

MC(:,9)=(data(:,I(1,9))./W(:,9)); %for Na2O
MC(:,10)=(data(:,I(1,10))./W(:,10)); %for K2O

%calculates for BaO if it is included in the analysis 
if I(1,11)==0
    MC(:,11)=zeros(m,1);
else
    MC(:,11)=data(:,I(1,11))./W(:,11); %for BaO
end 

%calculates for F if it is included in the analysis 
if I(1,12)==0
    MC(:,12)=zeros(m,1);
else
    MC(:,12)=data(:,I(1,12))./W(:,12); %for F
end 

%calculates for Cl if it is included in the analysis 
if I(1,13)==0
    MC(:,13)=zeros(m,1);
else
    MC(:,13)=data(:,I(1,13))./W(:,13); %for Cl
end 


%% Moles of O2 units
O2(:,1)=MC(:,1).*2; %SiO2
O2(:,2)=MC(:,2).*2; %TiO2
O2(:,3)=MC(:,3).*3; %Al2O3
O2(:,4)=MC(:,4).*3; %Cr2O3
O2(:,5)=MC(:,5); %FeO
O2(:,6)=MC(:,6); %MnO
O2(:,7)=MC(:,7); %MgO
O2(:,8)=MC(:,8); %CaO
O2(:,9)=MC(:,9); %Na2O
O2(:,10)=MC(:,10); %K2O
O2(:,11)=MC(:,11); %BaO
O2(:,12)=MC(:,12); %F
O2(:,13)=MC(:,13); %Cl

O2_Sum=sum(O2(:,1:13),2)-0.5.*(O2(:,12)+O2(:,13)); %sum of O2, including F and Cl

O2_N=(11)./O2_Sum; %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

%% atoms pfu

APFU(:,1)=N_Ox(:,1)./2; %Si
APFU(:,2)=N_Ox(:,2)./2; %Ti
APFU(:,3)=N_Ox(:,3).*(2/3); %Al
APFU(:,4)=N_Ox(:,4).*(2/3); %Cr
APFU(:,5)=N_Ox(:,5); %Fe
APFU(:,6)=N_Ox(:,6); %Mn
APFU(:,7)=N_Ox(:,7); %Mg
APFU(:,8)=N_Ox(:,8); %Ca
APFU(:,9)=N_Ox(:,9).*2; %Na
APFU(:,10)=N_Ox(:,10).*2; %K
APFU(:,11)=N_Ox(:,11); %Ba
APFU(:,12)=N_Ox(:,12); %F
APFU(:,13)=N_Ox(:,13); %Cl
APFU(:,14)=sum(APFU,2); %calculations the total

%% Structural Formula

%T site
StrctFrm(:,1)=APFU(:,1); %Si (T)

%Al (T)
for c=1:m
    if 4-StrctFrm(c,1) > APFU(c,3)
        StrctFrm(c,2)=APFU(c,3); 
    else
        StrctFrm(c,2)=4-StrctFrm(c,1);
    end
end

StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2); %Sum of T

StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2); %Altotal-Al(T)=Al(M)
StrctFrm(:,5)=APFU(:,2); %Ti (M)
StrctFrm(:,6)=APFU(:,4); %Cr (M)
StrctFrm(:,7)=APFU(:,5); %Fe2+ (M)
StrctFrm(:,8)=APFU(:,6); %Mn (M)
StrctFrm(:,9)=APFU(:,7); %Mg (M)
StrctFrm(:,10)=sum(StrctFrm(:,4:1:9),2); %sum (M)

StrctFrm(:,11)=APFU(:,8); %Ca (I)
StrctFrm(:,12)=APFU(:,9); %Na (I)
StrctFrm(:,13)=APFU(:,10); %K (I)
StrctFrm(:,14)=APFU(:,11); %Ba (I)

StrctFrm(:,15)=sum(StrctFrm(:,11:1:14),2)+StrctFrm(:,10)+StrctFrm(:,3); %cation sum

StrctFrm(:,16)=APFU(:,12); % F (A)
StrctFrm(:,17)=APFU(:,13); % Cl (A)
StrctFrm(:,18)=2-(APFU(:,12)+APFU(:,13)); % OH (A)

%% Endmembers

%Dioctahedral Micas
%the sum of M=3 cations in TriOct micas, and M=2 cations in DiOct, 

for c=1:m
   if StrctFrm(c,10)>2 %If M>2, there is a TriOct component
       XTriOct(c,:)=StrctFrm(c,10)-2; %Sum of M-2, scales from 0 to 1
   else
       if StrctFrm(c,10)>3
           XTriOct(c,:)=1;  %If M>3, there is a no DiOct Component
       else
           XTriOct(c,:)=0; %If M<3 and M<2 then there is no TriOct component
       end
   end
end

XDiOct=1-XTriOct; %fraction of dioctahedral mica

%calculate fractionn of Ms + Pg + Mrg + Pyl
%looks at octahedrally coordinated Al
for c=1:m
   if StrctFrm(c,4)>2 
       XMsum(c,:)=1; %If viAl > 2, fix it as 1 (upper limit for muscovite, paragonite, margarite, pyrophyllite)
   else
       if StrctFrm(c,4)<1
           XMsum(c,:)=0;  %If viAl < 1, there may be a celadonite (Fe3+) component, which is not considered here
       else
           XMsum(c,:)=StrctFrm(c,4)-1; %scales viAl between 0 and 1
       end
   end
end

%Al-Celadonite component 
XCel=1-XMsum; %total amout (Fe Al-Celadonite + Mg Al-Celadonite)

%Calculate XMg
for c=1:m
   if StrctFrm(c,9)<=0
       XMg(c,:)=0; %Avoids devision by 0 if no Mg is present
   else
       XMg(c,:)=StrctFrm(c,9)./(StrctFrm(c,9)+StrctFrm(c,7)); %calculates XMg
   end
end

XMgCel=XMg.*XCel; %fraction of Al-Celadonite (Mg-bearing endmember)
XFeCel=XCel-XMgCel; %fraction of Ferro Al-Celadonite

%fraction Ms, Pg, and Mrg
XMPM=(StrctFrm(:,11)+StrctFrm(:,12)+StrctFrm(:,13)).*XMsum; %Sum of Ca + Na + K multiplied by the total fraction of Ms, Pg, Mrg, and Prl

XPrl=XMsum-XMPM; %fraction of pyrophyllite

%Calculate XMrg: XCa * total proportion of Ms + Pg + Mrg
for c=1:m
   if (StrctFrm(c,11)+StrctFrm(c,12)+StrctFrm(c,13))<=0
       XMrg(c,:)=0; %If there is no Ca, Na, and K, then no XCa, avoids division by 0
   else
       XMrg(c,:)=(StrctFrm(c,11)./(StrctFrm(c,11)+StrctFrm(c,12)+StrctFrm(c,13))).*XMPM(c,:);
   end
end

%Calculate XPg: XNa * total proportion of Ms + Pg + Mrg
for c=1:m
   if (StrctFrm(c,11)+StrctFrm(c,12)+StrctFrm(c,13))<=0
       XPg(c,:)=0; %If there is no Ca, Na, and K, then no XNa, avoids division by 0
   else
       XPg(c,:)=(StrctFrm(c,12)./(StrctFrm(c,11)+StrctFrm(c,12)+StrctFrm(c,13))).*XMPM(c,:);
   end
end

%Calculate XMs: XK * total proportion of Ms + Pg + Mrg
for c=1:m
   if (StrctFrm(c,11)+StrctFrm(c,12)+StrctFrm(c,13))<=0
       XMs(c,:)=0; %If there is no Ca, Na, and K, then no XNa, avoids division by 0
   else
       XMs(c,:)=(StrctFrm(c,13)./(StrctFrm(c,11)+StrctFrm(c,12)+StrctFrm(c,13))).*XMPM(c,:);
   end
end

DiOct_Endmembers(:,1)=XMgCel.*XDiOct; %Al-Celadonite
DiOct_Endmembers(:,2)=XFeCel.*XDiOct; %Fe Al-Celadonite
DiOct_Endmembers(:,3)=XPrl.*XDiOct; %Pyrophyllite
DiOct_Endmembers(:,4)=XMrg.*XDiOct; %Margarite
DiOct_Endmembers(:,5)=XPg.*XDiOct; %Paragonite
DiOct_Endmembers(:,6)=XMs.*XDiOct; %Muscovite
DiOct_Endmembers(:,7)=XTriOct; %trioctahedral component

%Trioctahedral Micas

%pholgopite-annite solid solution
XPhlAnn=APFU(:,1)-2; %fraction of the Phl+Ann endmembers
                     %Phl and Ann has 3 Si, whereas Sid and East has 2
XPhl=XPhlAnn.*XMg; %unnormalized fraction of Phl
XAnn=XPhlAnn-XPhl; %unnormalized fraction of Ann

%siderophyllite-eastonite solid solution
XSidEast=1-XPhlAnn; %fraction of the Sid+Eastendmembers
XEast=XSidEast.*XMg;
XSid=XSidEast-XEast;

TriOct_Endmembers(:,1)=XPhl.*XTriOct; %phlogopite
TriOct_Endmembers(:,2)=XAnn.*XTriOct; %annite
TriOct_Endmembers(:,3)=XEast.*XTriOct; %eastonite 
TriOct_Endmembers(:,4)=XSid.*XTriOct; %siderophyllite
TriOct_Endmembers(:,5)=XDiOct; %dioctohedral component

%% Plots 

if strcmp(wantstrctfrm, 'y')
    
    %prompts the user if they wish to plot the mica data
    prompt1='Do you wish to plot Mica compositions? (y|n): ';
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
        
        prompt5='Do you wish to plot TriOctahedral Micas? (y|n): ';
        wanttri=input(prompt5, 's');
        
        if strcmp(wanttri, 'y')
            
            figure('Name','TriOctahedral Mica Plot')
            %annite to eastonite grid lines 
            plot([0.2 0.0], [0.0 0.2],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5)
            hold on
            plot([0.4 0.0], [0.0 0.4],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([0.6 0.0], [0.0 0.6],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5)
            hold on
            plot([0.8 0.0], [0.0 0.8],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([1.0 0.0], [0.0 1.0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5)
            hold on
            plot([0.2 1.0], [1.0 0.2],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([0.4 1.0], [1.0 0.4],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5)
            hold on
            plot([0.6 1.0], [1.0 0.6],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([0.8 1.0], [1.0 0.8],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5)
            hold on

            
            %phlogopite to siderophyllite grid lines
            plot([0.8 1.0], [0.0 0.2],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) 
            hold on
            plot([0.6 1.0], [0.0 0.4],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([0.4 1.0], [0.0 0.6],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) 
            hold on
            plot([0.2 1.0], [0.0 0.8],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([0.0 1.0], [0.0 1.0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) 
            hold on
            plot([0.0 0.8], [0.2 1.0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([0.0 0.6], [0.4 1.0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) 
            hold on
            plot([0.0 0.4], [0.6 1.0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5)
            hold on
            plot([0.0 0.2], [0.8 1.0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) 
            hold on
            
            scatter(XMg(:,1),StrctFrm(:,4),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            xlabel('X_{Mg}')
            ylabel('Al_{M}')
            axis equal %makes the spacing of the axes intervals equal
            %set axis font size to 12
            ax=gca;
            fon=ax.FontSize;
            ax.FontSize=12;
            xlim([0 1])
            ylim([0 1])
            text(0,-0.06,'Annite','FontSize',12,'HorizontalAlignment','center')
            text(0,1.05,'Siderophyllite','FontSize',12,'HorizontalAlignment','center')
            text(1,1.05,'Eastonite','FontSize',12,'HorizontalAlignment','center')
            text(1,-0.06,'Phlogopite','FontSize',12,'HorizontalAlignment','center')
        end
        
        prompt6='Do you wish to plot the Ms-Cel-Prl ternary? (y|n): ';
        wanttern=input(prompt6, 's');
        
        if strcmp(wanttern, 'y')
            
            figure('Name','Ms-Cel-Prl Ternary');
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
            text(-0.10,0.0,'ms','FontSize',14)
            text(0.50,0.92,'Alcel','FontSize',14)
            text(1.02,-0.04,'prl','FontSize',14)
            
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

            XAlCel=DiOct_Endmembers(:,1)+DiOct_Endmembers(:,2); %sum Fe and Mg endmembers

            %normalize the ternary components to 1
            XAlCelN=XAlCel./(XAlCel+XPrl+XMs); 
            XPrlN=XPrl./(XAlCel+XPrl+XMs); 
            
            %transforms the data to ternary space
            X2=0.5.*(XAlCelN)+(XPrlN);
            Y2=(XAlCelN)*(cos(30*pi()/180));
            
            scatter(X2(:),Y2(:),symbsize,symb,'filled','MarkerFaceAlpha',2/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            hold off
            
            axis image
            axis off
            
        end
        
        prompt7='Do you wish to plot Muscovite-Celadonite solid solution (Al vs Si)? (y|n): ';
        wantSiAl=input(prompt7, 's');
        
        if strcmp(wantSiAl, 'y')
            figure('Name','Muscovite-Celadonite Plot')
            
            plot([3 1],[3 4],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([2.78 2.82],[3.06 3.14],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([2.58 2.62],[3.16 3.24],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([2.38 2.42],[3.26 3.34],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([2.18 2.22],[3.36 3.44],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([1.98 2.02],[3.46 3.54],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([1.78 1.82],[3.56 3.64],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([1.58 1.62],[3.66 3.74],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([1.38 1.42],[3.76 3.84],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            plot([1.18 1.22],[3.86 3.94],'color',[0.3 0.3 0.3],'linewidth',1.0)
            hold on
            scatter(APFU(:,3),APFU(:,1),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            xlabel('Al (APFU)','FontSize',12)
            ylabel('Si (APFU)','FontSize',12)
            

            %set axis font size to 12
            ax=gca;
            ax.FontSize=12;
            axis equal %makes the spacing of the axes intervals equal
            xlim([1 3])
            ylim([3 4])
            text(1,4.07,'Celadonite','FontSize',14,'HorizontalAlignment','center')
            text(3,2.89,'Muscovite','FontSize',14,'HorizontalAlignment','center')
            text(3,2.82,'+ Paragonite','FontSize',14,'HorizontalAlignment','center')
            text(2.85,3.2,'0.1','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(2.65,3.3,'0.2','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(2.45,3.4,'0.3','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(2.25,3.5,'0.4','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(2.05,3.6,'0.5','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(1.85,3.7,'0.6','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(1.65,3.8,'0.7','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(1.45,3.9,'0.8','FontSize',12,'HorizontalAlignment','center','Rotation',60)
            text(1.15,3.8,'0.9','FontSize',12,'HorizontalAlignment','center','Rotation',60)           
        end
        
        prompt8='Do you wish to plot Na vs Si? (y|n): ';
        wantSiNa=input(prompt8, 's');
        
        if strcmp(wantSiNa, 'y')
            figure('Name','Na (APFU) vs. Si (APFU) Plot')

            scatter(APFU(:,9),APFU(:,1),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            
            xlabel('Na (APFU)','FontSize',12)
            ylabel('Si (APFU)','FontSize',12)
            
            %set axis font size to 12
            ax=gca;
            ax.FontSize=12;
            axis equal %makes the spacing of the axes intervals equal
            ax.Box = 'on';
            xlim([0 1])
            ylim([3 4])

        end

        prompt9='Do you wish to plot the Cl-F-OH ternary? (y|n): ';
        wantOH=input(prompt9, 's');
        
        if strcmp(wantOH, 'y')
            
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
            text(1.02,-0.04,'OH','FontSize',14)
            
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
            
            F=StrctFrm(:,16)./(StrctFrm(:,16)+StrctFrm(:,17)+StrctFrm(:,18));
            OH=StrctFrm(:,18)./(StrctFrm(:,16)+StrctFrm(:,17)+StrctFrm(:,18));
            
            %transforms the data to ternary space
            X3=0.5.*(F)+(OH);
            Y3=(F)*(cos(30*pi()/180));

            scatter(X3(:),Y3(:),symbsize,symb,'filled','MarkerFaceAlpha',2/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            hold off
            
            axis image
            axis off
            
        end
    end
    
all=[StrctFrm DiOct_Endmembers TriOct_Endmembers];
StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Sum_T','Al_M','Ti_M','Cr_M','Fe_M','Mn_M','Mg_M','M_Sum','Ca_I','Na_I','K_I','Ba_I','Total_Cations','F_A','Cl_A','OH_A','XAlcel','XFeAlcel','Xprl','Xmrg','Xpg','Xms','XTriOct','Xphl','Xann','Xeas','Xsid','XDiOct'});
APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Cr','Fe','Mn','Mg','Ca','Na','K','Ba','F','Cl','Sum'});

end


