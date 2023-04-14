%% Ilmenite structural formula

function [APFU]=ilmenite_fe3(data,headers)

[m,n]=size(data); %finds the x and y size of the input data matrix

I(1,1)=find(strcmp(headers,'TiO2'));

%makes SiO2 optional
if strcmp(headers,'SiO2')==zeros(1,length(headers))
    I(1,2)=0;
else
    I(1,2)=find(strcmp(headers,'SiO2'));
end

I(1,3)=find(strcmp(headers,'Al2O3'));
I(1,4)=find(strcmp(headers,'Cr2O3'));

%makes V2O3 optional
if strcmp(headers,'V2O3')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'V2O3'));
end

I(1,6)=find(strcmp(headers,'FeO'));
I(1,7)=find(strcmp(headers,'MnO'));

%makes MgO optional
if strcmp(headers,'MgO')==zeros(1,length(headers))
    I(1,8)=0;
else
    I(1,8)=find(strcmp(headers,'MgO'));
end

%makes CaO optional
if strcmp(headers,'CaO')==zeros(1,length(headers))
    I(1,9)=0;
else
    I(1,9)=find(strcmp(headers,'CaO'));
end

cat=2.0; %cations per formula unit
Opfu=3.0; %oxygens per formula unit

%% Molecular weights

SiO2=60.083;
TiO2=79.865;
Al2O3=101.961;
Cr2O3=151.989;
Fe2O3=159.6874;
V2O3=149.881;
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

W=[TiO2,SiO2,Al2O3,Cr2O3,V2O3,FeO,MnO,MgO,CaO];

%% Calculate cations units

MC(:,1)=data(:,I(1,1))./W(:,1); %for TiO2

%adds a column of zeros if SiO2 is not included in the calculation
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for SiO2
end

MC(:,3)=(data(:,I(1,3))./W(:,3)).*2; %for Al2O3
MC(:,4)=(data(:,I(1,4))./W(:,4)).*2; %for Cr2O3

%adds a column of zeros if V2O3 is not included in the calculation
if I(1,5)==0
    MC(:,5)=zeros(m,1);
else
    MC(:,5)=(data(:,I(1,5))./W(:,5)).*2; %for V2O3
end

MC(:,6)=data(:,I(1,6))./W(:,6); %for FeO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MnO
MC(:,8)=data(:,I(1,8))./W(:,8); %for MgO

%adds a column of zeros if MgO is not included in the calculation
if I(1,8)==0
    MC(:,8)=zeros(m,1);
else
    MC(:,8)=data(:,I(1,8))./W(:,8); %for MgO
end

%adds a column of zeros if CaO is not included in the calculation
if I(1,9)==0
    MC(:,9)=zeros(m,1);
else
    MC(:,9)=data(:,I(1,9))./W(:,9); %for CaO
end

MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cation

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for TiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5).*(3/2); %for V2O3
O2(:,6)=MCnorm(:,6); %for FeO
O2(:,7)=MCnorm(:,7); %for MnO
O2(:,8)=MCnorm(:,8); %for MgO
O2(:,9)=MCnorm(:,9); %for CaO

O2total=sum(O2,2); %O2 totals

%% Atoms PFU

APFU(:,1)=MCnorm(:,1); %for Ti
APFU(:,2)=MCnorm(:,2); %for Si
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,4)=MCnorm(:,4); %for Cr
APFU(:,5)=MCnorm(:,5); %for V
APFU(:,8)=MCnorm(:,7); %for Mn
APFU(:,9)=MCnorm(:,8); %for Mg
APFU(:,10)=MCnorm(:,9); %for Ca

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
            APFU(c,6)=2*(Opfu-O2total(c,1)); 
        else
            APFU(c,6)=MCnorm(c,6);
        end
    else
        APFU(c,6)=0;
    end
end

APFU(:,7)=MCnorm(:,6)-APFU(:,6); %the APFU of Fe2+ equals totalFe-Fe3+

APFU(:,11)=sum(APFU,2); %calculations the total, which should be 3

% Oxygen deficiency 
APFU(:,12)=Opfu-O2total; %must be greater than zero

%% Endmembers 

XHem=APFU(:,6)./(APFU(:,6)+APFU(:,1));
XIlGkPy=1-XHem;
XIlm=XIlGkPy.*(APFU(:,7)./(APFU(:,7)+APFU(:,8)+APFU(:,9)));
XPph=XIlGkPy.*(APFU(:,8)./(APFU(:,7)+APFU(:,8)+APFU(:,9)));
XGk=XIlGkPy.*(APFU(:,9)./(APFU(:,7)+APFU(:,8)+APFU(:,9)));

all=[APFU XIlm XGk XPph XHem];

APFU=array2table(all,'VariableNames',{'Ti','Si','Al','Cr','V','Fe3','Fe2','Mn','Mg','Ca','Sum','O2_deficiency','XIlm','XGk','XPph','XHem'});

end 

