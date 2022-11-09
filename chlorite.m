%%chlorite Structural Formula

function [StrctFrm, APFU]=chlorite(data,headers,wantstrctfrm)

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

%Makes NiO optional
if strcmp(headers,'NiO')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'NiO'));
end

I(1,5)=find(strcmp(headers,'FeO'));
I(1,6)=find(strcmp(headers,'MnO'));
I(1,7)=find(strcmp(headers,'MgO'));

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

W=[SiO2,TiO2,Al2O3,NiO,FeO,MnO,MgO];

%% Calculate cations units

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%calculates for TiO2 if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1); 
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

MC(:,3)=(data(:,I(1,3))./W(:,3)).*2; %for Al2O3


%calculates for NiO if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=data(:,I(1,4))./W(:,4); %for NiO
end

MC(:,5)=data(:,I(1,5))./W(:,5); %for FeO
MC(:,6)=data(:,I(1,6))./W(:,6); %for MnO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MgO


%% Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*(3/2); %for Al2O3
O2(:,4)=MC(:,4); %for NiO
O2(:,5)=MC(:,5); %for FeO
O2(:,6)=MC(:,6); %for MnO
O2(:,7)=MC(:,7); %for MgO

O2total=sum(O2,2); %O2 totals
Opfu=14.0; %oxygens per formula unit
MCnormfact=Opfu./O2total; %normalization factor

APFU=MCnormfact.*MC; %creates a matrix of normalized cations
APFU(:,8)=sum(APFU,2); %calculations the total, which should be close to 10

%% Structural formula

% Si (T)
for c=1:m
    if APFU(c,1)<4.000 %if Si is <4.00, then the Si content is the APFU of Si
        StrctFrm(c,1)=APFU(c,1);
    else
        StrctFrm(c,1)=4; %or else, it assumes 4
    end
end

%Al(T)
for c=1:m
    if 4-StrctFrm(c,1)>0 %Is 4-Si > 0? If y, then some Al goes into T
        if 4-StrctFrm(c,1)>APFU(c,3) %For low Al , 4-Si may be > Al
            StrctFrm(c,2)=APFU(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=4-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to M
        end
    else
        StrctFrm(c,2)=0; %if Si=4, then no Al goes into T
    end
end

%T site sum
StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2); %should sum to 4

%Al (M)
StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2); 

StrctFrm(:,5)=APFU(:,2); %Ti (M)
StrctFrm(:,6)=APFU(:,4); %Ni (M)
StrctFrm(:,7)=APFU(:,5); %Fe total (M)
StrctFrm(:,8)=APFU(:,6); %Mn (M)
StrctFrm(:,9)=APFU(:,7); %Mg (M)
StrctFrm(:,10)=(StrctFrm(:,4)-StrctFrm(:,2))./2; %Vacancies on M1 = (Alvi-Aliv)/2
StrctFrm(:,11)=StrctFrm(:,9)+StrctFrm(:,8)+StrctFrm(:,7)+StrctFrm(:,6)+StrctFrm(c,5)+StrctFrm(c,4)+StrctFrm(c,3); %Total Cations

XMg(:,1)=StrctFrm(:,9)./(StrctFrm(:,7)+StrctFrm(:,9)); %XMg = Mg/(Mg+Fe)


%% Plots 

if strcmp(wantstrctfrm, 'y')
    
    %prompts the user if they wish to plot the CPX data
    prompt1='Do you wish to plot chlorite compositions? (y|n): ';
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
        
        figure('Name','Chlorite Plot')
        
        scatter(StrctFrm(:,4),XMg(:,1),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
        axis square %makes the spacing of the axes intervals equal
        xlim([1 3])
        ylim([0 1])
        ylabel('X_{Mg}')
        xlabel('Al_{M} (APFU)')
        
        %labels
        text(0.8,-0.06,'Chamosite','FontSize',14)
        text(0.8,1.05,'Clinochlore','FontSize',14)
        text(2.85,1.05,'Sudoite','FontSize',14)
        %set axis font size to 12
        ax=gca;
        fon=ax.FontSize;
        ax.FontSize=12;
        ax.Box = 'on';


    end
end

all=[StrctFrm XMg];
StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Sum_T','Al_M','Ti_M','Ni_M','Fe_M','Mn_M','Mg_M','vac_M1','Total_Cations','XMg'});
APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Ni','Fe','Mn','Mg','Sum'});

end



            