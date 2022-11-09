%%staurolite Structural Formula

function [StrctFrm]=staurolite(data,headers)

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

%Makes ZnO optional
if strcmp(headers,'ZnO')==zeros(1,length(headers))
    I(1,8)=0;
else
    I(1,8)=find(strcmp(headers,'ZnO'));
end

%Makes CaO optional
if strcmp(headers,'CaO')==zeros(1,length(headers))
    I(1,9)=0;
else
    I(1,9)=find(strcmp(headers,'CaO'));
end

%Makes Na2O optional
if strcmp(headers,'Na2O')==zeros(1,length(headers))
    I(1,10)=0;
else
    I(1,10)=find(strcmp(headers,'Na2O'));
end

%Makes K2O optional
if strcmp(headers,'K2O')==zeros(1,length(headers))
    I(1,11)=0;
else
    I(1,11)=find(strcmp(headers,'K2O'));
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

W=[SiO2,TiO2,Al2O3,Fe2O3,Cr2O3,FeO,MnO,MgO,ZnO,CaO,Na2O,K2O,F,Cl];

%% Moles of Oxides

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%calculates for TiO2 if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

MC(:,3)=data(:,I(1,3))./W(:,3); %for Al2O3

%adjust Fe2O3 and FeO to a predetermined ratio
prompt1 = 'Assume an Fe3+/Fetotal ratio: ';
disp('Holdaway et al. (1991) recommends 0.035 for ilm-bearing rocks (XHem<0.1)');
disp('and 0.07 for hem-ilm rocks (XHem>0.1).');
wantFe3ratio = input(prompt1);

FeOT=data(:,I(1,5));
Fe2O3T=data(:,I(1,5)).*(0.5*Fe2O3/FeO);
MC(:,4)=(Fe2O3T-((1-wantFe3ratio).*FeOT*(0.5*Fe2O3/FeO)))./W(:,4);

%calculates for Cr2O3 if it is included in the analysis 
if I(1,4)==0
    MC(:,5)=zeros(m,1);
else
    MC(:,5)=data(:,I(1,4))./W(:,5); %for Cr2O3
end 

%For FeO with adjustment for Fe3+/Fetotal ratio
MC(:,6)=(((1-wantFe3ratio).*Fe2O3T)./(0.5*Fe2O3/FeO))./W(:,6);

MC(:,7)=data(:,I(1,6))./W(:,7); %for MnO
MC(:,8)=data(:,I(1,7))./W(:,8); %for MgO

%calculates for ZnO if it is included in the analysis 
if I(1,8)==0
    MC(:,9)=zeros(m,1);
else
    MC(:,9)=data(:,I(1,8))./W(:,9); %for ZnO
end 

%calculates for CaO if it is included in the analysis 
if I(1,9)==0
    MC(:,10)=zeros(m,1);
else
    MC(:,10)=data(:,I(1,9))./W(:,10); %for CaO
end 

%calculates for Na2O if it is included in the analysis 
if I(1,10)==0
    MC(:,11)=zeros(m,1);
else
    MC(:,11)=data(:,I(1,10))./W(:,11); %for Na2O
end 

%calculates for K2O if it is included in the analysis 
if I(1,11)==0
    MC(:,12)=zeros(m,1);
else
    MC(:,12)=data(:,I(1,11))./W(:,12); %for K2O
end 

%calculates for F if it is included in the analysis 
if I(1,12)==0
    MC(:,13)=zeros(m,1);
else
    MC(:,13)=data(:,I(1,12))./W(:,13); %for F
end 

%calculates for Cl if it is included in the analysis 
if I(1,13)==0
    MC(:,14)=zeros(m,1);
else
    MC(:,14)=data(:,I(1,13))./W(:,14); %for Cl
end 

%% Moles of O2 units

O2(:,1)=MC(:,1).*2; %SiO2
O2(:,2)=MC(:,2).*2; %TiO2
O2(:,3)=MC(:,3).*3; %Al2O3
O2(:,4)=MC(:,4).*3; %Fe2O3
O2(:,5)=MC(:,5).*3; %Cr2O3
O2(:,6)=MC(:,6); %FeO
O2(:,7)=MC(:,7); %MnO
O2(:,8)=MC(:,8); %MgO
O2(:,9)=MC(:,9); %ZnO
O2(:,10)=MC(:,10); %CaO
O2(:,11)=MC(:,11); %Na2O
O2(:,12)=MC(:,12); %K2O
O2(:,13)=MC(:,13); %F
O2(:,14)=MC(:,14); %Cl

O2_Sum=sum(O2(:,1:14),2)-0.5.*(O2(:,13)+O2(:,14)); %sum of O2, including F and Cl

O2_N=(48)./O2_Sum; %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

%% Iterate OH

%initial oxygens of OH
Hin=4-(N_Ox(:,14)+N_Ox(:,13)); %H = 4 - (F+Cl)
O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O
H2Oi=O_OH.*18.015; %initial H2O wt %

for z=1:10
    O2(:,1)=MC(:,1).*2; %SiO2
    O2(:,2)=MC(:,2).*2; %TiO2
    O2(:,3)=MC(:,3).*3; %Al2O3
    O2(:,4)=MC(:,4).*3; %Fe2O3
    O2(:,5)=MC(:,5).*3; %Cr2O3
    O2(:,6)=MC(:,6); %FeO
    O2(:,7)=MC(:,7); %MnO
    O2(:,8)=MC(:,8); %MgO
    O2(:,9)=MC(:,9); %ZnO
    O2(:,10)=MC(:,10); %CaO
    O2(:,11)=MC(:,11); %Na2O
    O2(:,12)=MC(:,12); %K2O
    O2(:,13)=MC(:,13); %F
    O2(:,14)=MC(:,14); %Cl
    O2(:,15)=O_OH; %H2O

    O2_Sum=sum(O2(:,1:15),2)-0.5.*(O2(:,13)+O2(:,14)); %sum of O2, including F and Cl

    O2_N=(48)./O2_Sum; %normalization factor

    %normalized moles of anions
    N_Ox=O2.*O2_N;

    Hin=4-(N_Ox(:,14)+N_Ox(:,13)); %H = 4 - (F+Cl)
    O_OH=(0.5.*Hin)./O2_N; %oxygen moles of H2O
    H2Oi=O_OH.*18.015; %initial H2O wt %
end

%% Moles of cations (Oxygen Normalization)

APFU_O(:,1)=N_Ox(:,1)./2; %Si
APFU_O(:,2)=N_Ox(:,2)./2; %Ti
APFU_O(:,3)=N_Ox(:,3).*(2/3); %Al
APFU_O(:,4)=N_Ox(:,4).*(2/3); %Fe3
APFU_O(:,5)=N_Ox(:,5).*(2/3); %Cr3
APFU_O(:,6)=N_Ox(:,6); %FeO
APFU_O(:,7)=N_Ox(:,7); %MnO
APFU_O(:,8)=N_Ox(:,8); %MgO
APFU_O(:,9)=N_Ox(:,9); %ZnO
APFU_O(:,10)=N_Ox(:,10); %CaO
APFU_O(:,11)=N_Ox(:,11).*2; %Na2O
APFU_O(:,12)=N_Ox(:,12).*2; %K2O
APFU_O(:,13)=N_Ox(:,13); %F
APFU_O(:,14)=N_Ox(:,14); %Cl
APFU_O(:,15)=N_Ox(:,15)*2; %OH

%% normalization to 25.53 Si + Al

%STEP 1
%Normmalization to Si + Al 
Nfact=25.53./(APFU_O(:,1)+APFU_O(:,3)); %normalization

APFU_N=APFU_O.*Nfact; %normalization of the cations

%cation charges
Chrg(:,1)=APFU_N(:,1).*4; %Si4+
Chrg(:,2)=APFU_N(:,2).*4; %Ti4+
Chrg(:,3)=APFU_N(:,3).*3; %Al3+
Chrg(:,4)=APFU_N(:,4).*3; %Fe3+
Chrg(:,5)=APFU_N(:,5).*3; %Cr3+
Chrg(:,6)=APFU_N(:,6).*2; %Fe2+
Chrg(:,7)=APFU_N(:,7).*2; %Mn2+
Chrg(:,8)=APFU_N(:,8).*2; %Mg2+
Chrg(:,9)=APFU_N(:,9).*2; %Zn2+
Chrg(:,10)=APFU_N(:,10).*2; %Ca2+
Chrg(:,11)=APFU_N(:,11); %Na1+
Chrg(:,12)=APFU_N(:,12); %K1+
Chrg(:,13)=APFU_N(:,13); %K1+
Chrg(:,14)=APFU_N(:,14); %K1+

Crg_Ex=sum(Chrg(:,1:14),2)-92; %excess charge (>92)

%reformulates OH, for O substitution following normalization
APFU_N(:,15)=96-sum(Chrg(:,1:14),2); 

APFU_N(:,16)=sum(APFU_N(:,1:12),2); %calculations the total
APFU_N(:,17)=30-APFU_N(:,16); %total number of vacancies

%% organization

%cation normalized data
StrctFrm_N=APFU_N;
StrctFrm_N(:,13)=APFU_N(:,17);
StrctFrm_N(:,14)=APFU_N(:,16);
StrctFrm_N(:,15)=APFU_N(:,13);
StrctFrm_N(:,16)=APFU_N(:,14);
StrctFrm_N(:,17)=APFU_N(:,15);

for c=1:m
    if Crg_Ex(c,1)<0
        StrctFrm_N(c,18)=0;
    else
        StrctFrm_N(c,18)=Crg_Ex(c,1); 
    end
end

StrctFrm_N(:,19)=APFU_N(:,15)+APFU_N(:,14)+APFU_N(:,13)+StrctFrm_N(:,18);

%48 O normalized data
StrctFrm_O=APFU_O;

APFU_O(:,16)=sum(APFU_O(:,1:12),2); %calculations the total
APFU_O(:,17)=30-APFU_O(:,16); %total number of vacancies
StrctFrm_O(:,13)=APFU_O(:,17);
StrctFrm_O(:,14)=APFU_O(:,16);
StrctFrm_O(:,15)=APFU_O(:,13);
StrctFrm_O(:,16)=APFU_O(:,14);
StrctFrm_O(:,17)=APFU_O(:,15);
StrctFrm_O(:,18)=APFU_O(:,15)+APFU_O(:,14)+APFU_O(:,13);

%prompts the user what normalization to output
prompt1='Normalization by 25.53 Al + Si? (y|n): ';
wantnorm=input(prompt1, 's');

if strcmp(wantnorm, 'y')

    StrctFrm=array2table(StrctFrm_N,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn','Mg','Zn','Ca','Na','K','vac','cations_sum','F','Cl','OH','O','X_sum'});

else

    StrctFrm=array2table(StrctFrm_O,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Fe2','Mn','Mg','Zn','Ca','Na','K','vac','cations_sum','F','Cl','OH','X_sum'});

end


