%Clinopyroxene structural formula with Fe3+ estimation
function [StrctFrm, APFU]=pyroxene_fe3(data,headers,wantstrctfrm)
         
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
I(1,8)=find(strcmp(headers,'CaO'));

%Makes Na2O optional
if strcmp(headers,'Na2O')==zeros(1,length(headers))
    I(1,9)=0;
else
    I(1,9)=find(strcmp(headers,'Na2O'));
end

%Makes K2O optional
if strcmp(headers,'K2O')==zeros(1,length(headers))
    I(1,10)=0;
else
    I(1,10)=find(strcmp(headers,'K2O'));
end

cat=4.0; %cations per formula unit
Opfu=6.0; %oxygens per formula unit


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

W=[SiO2,TiO2,Al2O3,Cr2O3,FeO,MnO,MgO,CaO,Na2O,K2O];

%% Calculate cations units

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%calculates for TiO2 if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

MC(:,3)=(data(:,I(1,3))./W(:,3)).*2; %for Al2O3

%calculates for Cr2O3 if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=(data(:,I(1,4))./W(:,4)).*2; %for Cr2O3
end

MC(:,5)=data(:,I(1,5))./W(:,5); %for FeO
MC(:,6)=data(:,I(1,6))./W(:,6); %for MnO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MgO
MC(:,8)=data(:,I(1,8))./W(:,8); %for CaO

%calculates for K2O if it is included in the analysis 
if I(1,9)==0
    MC(:,9)=zeros(m,1);
else
    MC(:,9)=(data(:,I(1,9))./W(:,9)).*2; %for Na2O
end


%calculates for K2O if it is included in the analysis 
if I(1,10)==0
    MC(:,10)=zeros(m,1);
else
    MC(:,10)=(data(:,I(1,10))./W(:,10)).*2; %for K2O
end

MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5); %for FeO
O2(:,6)=MCnorm(:,6); %for MnO
O2(:,7)=MCnorm(:,7); %for MgO
O2(:,8)=MCnorm(:,8); %for CaO
O2(:,9)=MCnorm(:,9)./2; %for Na2O
O2(:,10)=MCnorm(:,10)./2; %for K2O

O2total=sum(O2,2); %O2 totals

%% Atoms pfu

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,5)=MCnorm(:,4); %for Cr
APFU(:,7)=MCnorm(:,6); %for Mn
APFU(:,8)=MCnorm(:,7); %for Mg
APFU(:,9)=MCnorm(:,8); %for Ca
APFU(:,10)=MCnorm(:,9); %for Na
APFU(:,11)=MCnorm(:,10); %for K

%calculation of Fe3+ from stoichiometry and charge balance
%the following if statement firsts checks if totalO2 = 6
%if so, then there is no Fe3+
%if totalO2 < 6, then we assume that the deficiency is caused by the
%assumption Fetotal = Fe2+
%in the nested if statement, if FeTotal > 2*(6-totalO2) then the amount
%of Fe3+ = 2*(6-totalO2), if false then, all Fe is Fe3+
for c=1:m
    if (Opfu-O2total(c,1)) > 0
        if MCnorm(c,5) > 2.*(Opfu-O2total(c,1))
            APFU(c,4)=2.*(Opfu-O2total(c,1)); 
        else
            APFU(c,4)=MCnorm(c,5);
        end
    else
        APFU(c,4)=0;
    end
end

APFU(:,6)=MCnorm(:,5)-APFU(:,4); %the APFU of Fe2+ = totalFe-Fe3+

APFU(:,12)=sum(APFU,2); %calculations the total, which should be 4

% Oxygen deficiency 
APFU(:,13)=Opfu-O2total; %must be greater than zero

%XMg
XMg=APFU(:,8)./(APFU(:,8)+APFU(:,6)); 

%% structural formula calculation

%T SITE
%Si 
for c=1:m
    if APFU(c,1)<2.000
        StrctFrm(c,1)=APFU(c,1); %If Si < 2, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=2; %If Si is in excess, then Si(T) = 2
    end
end

%Al(T)
for c=1:m
    if 2-StrctFrm(c,1)>0 %Is 2-Si > 0? If y, then some Al goes into T
        if 2-StrctFrm(c,1)>APFU(c,3) %For low Al cpx, 2-Si may be > Al
            StrctFrm(c,2)=APFU(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=2-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to M1
        end
    else
        StrctFrm(c,2)=0; %if Si=2, then no Al goes into T
    end
end

%Fe3+(T)
for c=1:m
    if 2-StrctFrm(c,1)-StrctFrm(c,2)>0 %Is 2-(Si+Al) > 0? If y, then some Fe3+ goes into T
        if 2-StrctFrm(c,1)-StrctFrm(c,2)>APFU(c,4) %For low Fe3+ cpx, 2-(Si+Al) may be > Fe3+
            StrctFrm(c,3)=APFU(c,4); %All Fe3+ goes into T
        else
            StrctFrm(c,3)=2-StrctFrm(c,1)-StrctFrm(c,2); %if there isn't enough space in T for all Fe3+, the rest will go to M1
        end
    else
        StrctFrm(c,3)=0; %if Si+Al=2, then no Fe3+ goes into T
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3); %Si + Al + Fe3+ in T

%M1 SITE

%Al(M1)
StrctFrm(:,5)=APFU(:,3)-StrctFrm(:,2); %Al(M1) = Total Al - Al(T)

%Ti (M1)
StrctFrm(:,6)=APFU(:,2);

%Fe3+ (M1)
StrctFrm(:,8)=APFU(:,4)-StrctFrm(:,3); %Fe3+(M1) = Total Fe3+ - Fe3+(T)

%Cr3+ (M1)
StrctFrm(:,7)=APFU(:,5);

%Mn (M1)
StrctFrm(:,9)=APFU(:,7);

%Mg (M1)
for c=1:m
    if XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5))<APFU(c,8) %if XMg*(1-Sum(Al to Mn) in M1) is < Mg
        StrctFrm(c,10)=XMg(c).*(1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,5)); %Mg(M1)=XMg*(1-Sum(Al to Mn) in M1) and some Mg goes into M2
    else
        StrctFrm(c,10)=APFU(c,8); %if not, all Mg goes into M1
    end
end

%Fe2+ (M1)
for c=1:m
    if 1-(StrctFrm(c,6)+StrctFrm(c,7)+StrctFrm(c,8)+StrctFrm(c,9)+StrctFrm(c,10)+StrctFrm(c,5))>0 %Is 1-Sum(Al to Mg) in M1 >0?, if yes then:
        if 1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,10)-StrctFrm(c,5)>APFU(c,6) %Is 1-Sum(Al to Mg) in M1 > Fe2+?, if yes then:
            StrctFrm(c,11)=APFU(c,6); %all Fe2+ goes into M1
        else
            StrctFrm(c,11)=1-StrctFrm(c,6)-StrctFrm(c,7)-StrctFrm(c,8)-StrctFrm(c,9)-StrctFrm(c,10)-StrctFrm(c,5); %if there isn't enough space in M1 for all Fe2+, some goes into M2
        end
    else
        StrctFrm(c,11)=0; %If M1 is already filled, then no Fe2+ goes into M1
    end
end


%Sum of M1 site
StrctFrm(:,12)=StrctFrm(:,5)+StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8)+StrctFrm(:,9)+StrctFrm(:,10)+StrctFrm(:,11);

%M2 SITE
%Mg (M2)
for c=1:m
    if APFU(c,8)-StrctFrm(c,10)>0 %Is Mgtotal-Mg(M1) > 0?
        StrctFrm(c,13)=APFU(c,8)-StrctFrm(c,10); %if yes, then some Mg goes into M2
    else
        StrctFrm(c,13)=0; %if no, then all Mg is in M1
    end
end

%Fe2+ (M2)
for c=1:m
    if APFU(c,6)-StrctFrm(c,11)>0 %Is Fe2+total-Fe2+(M1) > 0?
        StrctFrm(c,14)=APFU(c,6)-StrctFrm(c,11); %if yes, then some Fe2+ goes into M2
    else
        StrctFrm(c,14)=0; %if no, then all Fe2+ is in M1
    end
end

%Ca (M2)
StrctFrm(:,15)=APFU(:,9);

%Na (M2)
StrctFrm(:,16)=APFU(:,10);

%K (M2)
StrctFrm(:,17)=APFU(:,11);

%Sum of M2 site
StrctFrm(:,18)=StrctFrm(:,13)+StrctFrm(:,14)+StrctFrm(:,15)+StrctFrm(:,16)+StrctFrm(:,17);

%% end member calculations

%only do endmember calculation if structural formula with endembers is selected 

if strcmp(wantstrctfrm, 'y')
    %choose to calculate endmember fractions with Na-Fe3+-Cr 
    prompt1 = 'Do you wish to include Na-Fe3+-Cr endmembers? (y|n): ';
    wantNaFe3 = input(prompt1, 's');
    
    if strcmp(wantNaFe3, 'y')
        %Na-Ca endmembers 
        A(:,1)=APFU(:,9); %Ca
        A(:,2)=APFU(:,3); %Al
        A(:,3)=APFU(:,4); %Fe3+
        A(:,4)=APFU(:,8); %Mg
        A(:,5)=APFU(:,6); %Fe2+
        A(:,6)=APFU(:,5); %Cr
        AT=transpose(A); %transpose of A
        
        M=[2 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 2 0 0 0; 0 2 0 0 0 0; 0 0 0 0 0 1];
        X=zeros(6,m);
        
        for c=1:m
            X(:,c)=inv(M)*AT(:,c); %calculates endmembers
        end
        
        Xtot=sum(X);%sum of endmembers
        Xnorm=X./sum(X); %normalizes the endmembers to 1
        XnormT=transpose(Xnorm); %transposes back
        
        Endmembers(:,1)=XnormT(:,1); % XWo
        Endmembers(:,2)=XnormT(:,2); % XFs
        Endmembers(:,3)=XnormT(:,3); % XEn
        Endmembers(:,4)=XnormT(:,4); % XJd
        Endmembers(:,5)=XnormT(:,5); % XAeg
        Endmembers(:,6)=XnormT(:,6); % XKos
        Endmembers(:,7)=(XnormT(:,1)+XnormT(:,2)+XnormT(:,3))./(XnormT(:,1)+XnormT(:,2)+XnormT(:,3)+XnormT(:,4)+XnormT(:,5)+XnormT(:,6)); % XQuad
        
        EnWoFs(:,1)=(APFU(:,8)./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XEn
        EnWoFs(:,2)=(APFU(:,9)./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XWo
        EnWoFs(:,3)=((APFU(:,6)+APFU(:,7)+APFU(:,4))./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XFs
        
        %prompts the user if they wish to plot the CPX data
        prompt2='Do you wish to plot clinopyroxene compositions? (y|n): ';
        wantplots=input(prompt2, 's');
        
        %plots for Na-Fe3+-Cr choice
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
            
            %Q-J diagram
            figure('Name','Q-J Diagram')
            plot([0 2],[2 0],'k','linewidth',1.5)
            hold on
            plot([0 1.5],[1.5 0],'k','linewidth',1.5)
            hold on
            plot([0.3 0.4],[1.2 1.6],'k','linewidth',1.5)
            hold on
            plot([1.2 1.6],[0.3 0.4],'k','linewidth',1.5)
            hold on
            
            Y3=APFU(:,8)+APFU(:,6)+APFU(:,9);
            X3=APFU(:,10).*2;
            
            scatter(X3(:),Y3(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            hold off
            text(0.25,1.8,'Quad','FontSize',12)
            text(1.05,1,'Ca-Na','FontSize',12)
            text(1.85,0.2,'Na','FontSize',12)
            xlabel('J=2Na')
            ylabel('Q=(Ca+Mg+Fe^{2+})')
            
            %plots a ternary for Ca-Mg-Fe pyroxenes
            figure('Name','Ca-Mg-Fe Pyroxenes')
            
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
            %plot mineral compositional fields
            plot([0.25 0.75],[0.4330127 0.4330127],'k','linewidth',1.5)
            hold on
            plot([0.225 0.775],[0.38971143 0.38971143],'k','linewidth',1.5)
            hold on
            plot([0.1 0.9],[0.17320508 0.17320508],'k','linewidth',1.5)
            hold on
            plot([0.025 0.975],[0.04330127 0.04330127],'k','linewidth',1.5)
            hold on
            plot([0.5 0.5],[0.38971143 0.4330127],'k','linewidth',1.5)
            hold on
            plot([0.5 0.5],[0 0.04330127],'k','linewidth',1.5)
            hold on
            
            %labels
            text(-0.10,0.0,'en','FontSize',14)
            text(0.51,0.92,'wo','FontSize',14)
            text(1.02,-0.04,'fs','FontSize',14) 
            text(0.20,0.02,'Enstatite','FontSize',12)
            text(0.65,0.02,'Ferrosilite','FontSize',12)
            text(0.425,0.1,'Pigeonite','FontSize',12)
            text(0.45,0.28,'Augite','FontSize',12)
            text(0.29,0.412,'Diopside','FontSize',12)
            text(0.52,0.412,'Hedenbergite','FontSize',12)
            
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
            X2=0.5.*(EnWoFs(:,2))+(EnWoFs(:,3));
            Y2=(EnWoFs(:,2))*(cos(30*pi()/180));
            
            scatter(X2(:),Y2(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            hold off
            
            axis image
            axis off
            
            %plots a ternary figure for sodic-calcic pyroxenes
            figure('Name','Sodic-Calcic Pyroxenes');
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
            plot([0.4 0.6],[0.69282032 0.69282032],'k','linewidth',1.5)
            hold on
            plot([0.1 0.9],[0.17320508 0.17320508],'k','linewidth',1.5)
            hold on
            plot([0.5 0.5],[0 0.69282032],'k','linewidth',1.5)
            hold on
            
            %labels
            text(-0.10,0.0,'jd','FontSize',14)
            text(0.52,0.92,'Quad','FontSize',14)
            text(1.02,-0.04,'aeg','FontSize',14)
            text(0.20,0.09,'Jadeite','FontSize',12)
            text(0.66,0.09,'Aegirine','FontSize',12)
            text(0.26,0.35,'Omphacite','FontSize',12)
            text(0.58,0.37,'Aegirine-','FontSize',12)
            text(0.60,0.33,'Augite','FontSize',12)
            
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
            X1=0.5.*(Endmembers(:,7))+(Endmembers(:,5));
            Y1=(Endmembers(:,7))*(cos(30*pi()/180));
            
            scatter(X1(:),Y1(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            hold off
            
            axis image
            axis off
            
        end
            
        all=[StrctFrm Endmembers APFU(:,13)];
        StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M1','Ti_M1','Cr_M1','Fe3_M1','Mn_M1','Mg_M1','Fe2_M1','Sum_M1','Mg_M2','Fe2_M2','Ca_M2','Na_M2','K_M2','Sum_M2','Xwo','Xfs','Xen','Xjd','Xaeg','Xkos','XQuad','O2_deficiency'});            
        
    else
        %ortho and calcic pyroxene
        %Normalization procedure follows Morimoto et al. (1988)
        Endmembers(:,3)=(APFU(:,8)./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XEn
        Endmembers(:,1)=(APFU(:,9)./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XWo
        Endmembers(:,2)=((APFU(:,6)+APFU(:,5)+APFU(:,4))./(APFU(:,9)+APFU(:,6)+APFU(:,8)+APFU(:,7)+APFU(:,4))); % XFs
        prompt2 = 'Do you wish to plot clinopyroxene compositions? (y|n): ';
        wantplots = input(prompt2, 's');
        
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
            
            %Q-J diagram
            figure('Name','Q-J Diagram')
            plot([0 2],[2 0],'k','linewidth',1.5)
            hold on
            plot([0 1.5],[1.5 0],'k','linewidth',1.5)
            hold on
            plot([0.3 0.4],[1.2 1.6],'k','linewidth',1.5)
            hold on
            plot([1.2 1.6],[0.3 0.4],'k','linewidth',1.5)
            hold on
            
            Y3=APFU(:,8)+APFU(:,6)+APFU(:,9);
            X3=APFU(:,10).*2;
            
            scatter(X3(:),Y3(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            hold off
            text(0.25,1.8,'Quad','FontSize',12)
            text(1.05,1,'Ca-Na','FontSize',12)
            text(1.85,0.2,'Na','FontSize',12)
            xlabel('J=2Na')
            ylabel('Q=(Ca+Mg+Fe^{2+})')
            
            %plots a ternary for Ca-Mg-Fe pyroxenes
            figure('Name','Ca-Mg-Fe Pyroxenes');
            
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
            plot([0.25 0.75],[0.4330127 0.4330127],'k','linewidth',1.5)
            hold on
            plot([0.225 0.775],[0.38971143 0.38971143],'k','linewidth',1.5)
            hold on
            plot([0.1 0.9],[0.17320508 0.17320508],'k','linewidth',1.5)
            hold on
            plot([0.025 0.975],[0.04330127 0.04330127],'k','linewidth',1.5)
            hold on
            plot([0.5 0.5],[0.38971143 0.4330127],'k','linewidth',1.5)
            hold on
            plot([0.5 0.5],[0 0.04330127],'k','linewidth',1.5)
            hold on
            
            %labels
            text(-0.10,0.0,'en','FontSize',14)
            text(0.51,0.92,'wo','FontSize',14)
            text(1.02,-0.04,'fs','FontSize',14)  
            text(0.20,0.02,'Enstatite','FontSize',12)
            text(0.65,0.02,'Ferrosilite','FontSize',12)
            text(0.425,0.1,'Pigeonite','FontSize',12)
            text(0.45,0.28,'Augite','FontSize',12)
            text(0.29,0.412,'Diopside','FontSize',12)
            text(0.52,0.412,'Hedenbergite','FontSize',12)
            
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
            X2=0.5.*(Endmembers(:,1))+(Endmembers(:,2));
            Y2=(Endmembers(:,1))*(cos(30*pi()/180));
      
            scatter(X2(:),Y2(:),symbsize,symb,'filled','MarkerFaceAlpha',2/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            hold off
            
            axis image
            axis off
            
        end

        all=[StrctFrm Endmembers APFU(:,13)];
        StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M1','Ti_M1','Cr_M1','Fe3_M1','Mn_M1','Mg_M1','Fe2_M1','Sum_M1','Mg_M2','Fe2_M2','Ca_M2','Na_M2','K_M2','Sum_M2','Xwo','Xfs','Xen','O2_deficiency'});
      
    end
    
    
else
    APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn','Mg','Ca','Na','K','Sum','O2_deficiency'});
end

end

