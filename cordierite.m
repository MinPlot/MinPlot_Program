%Cordierite structural formula assumnig Fetotal = FeO

function [StrctFrm, APFU]=cordierite(data,headers,wantstrctfrm)

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
I(1,4)=find(strcmp(headers,'FeO'));
I(1,5)=find(strcmp(headers,'MnO'));
I(1,6)=find(strcmp(headers,'MgO'));

%Makes CaO optional
if strcmp(headers,'CaO')==zeros(1,length(headers))
    I(1,7)=0;
else
    I(1,7)=find(strcmp(headers,'CaO'));
end

I(1,8)=find(strcmp(headers,'Na2O'));
I(1,9)=find(strcmp(headers,'K2O'));

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

W=[SiO2,TiO2,Al2O3,FeO,MnO,MgO,CaO,Na2O,K2O];


%% Moles of oxides

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%calculates for TiO2 if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

MC(:,3)=(data(:,I(1,3))./W(:,3)); %for Al2O3
MC(:,4)=(data(:,I(1,4))./W(:,4)); %for FeO
MC(:,5)=(data(:,I(1,5))./W(:,5)); %for MnO
MC(:,6)=(data(:,I(1,6))./W(:,6)); %for MgO

%calculates for CaO if it is included in the analysis 
if I(1,7)==0
    MC(:,7)=zeros(m,1);
else
    MC(:,7)=data(:,I(1,7))./W(:,7); %for CaO
end 

MC(:,8)=(data(:,I(1,8))./W(:,8)); %for Na2O
MC(:,9)=(data(:,I(1,9))./W(:,9)); %for K2O

%% Moles of O2 units
O2(:,1)=MC(:,1).*2; %SiO2
O2(:,2)=MC(:,2).*2; %TiO2
O2(:,3)=MC(:,3).*3; %Al2O3
O2(:,4)=MC(:,4); %FeO
O2(:,5)=MC(:,5); %MnO
O2(:,6)=MC(:,6); %MgO
O2(:,7)=MC(:,7); %CaO
O2(:,8)=MC(:,8); %Na2O
O2(:,9)=MC(:,9); %K2O

O2_N=(18)./sum(O2,2); %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

%% atoms pfu

APFU(:,1)=N_Ox(:,1)./2; %Si
APFU(:,2)=N_Ox(:,2)./2; %Ti
APFU(:,3)=N_Ox(:,3).*(2/3); %Al
APFU(:,4)=N_Ox(:,4); %Fe
APFU(:,5)=N_Ox(:,5); %Mn
APFU(:,6)=N_Ox(:,6); %Mg
APFU(:,7)=N_Ox(:,7); %Ca
APFU(:,8)=N_Ox(:,8).*2; %Na
APFU(:,9)=N_Ox(:,9).*2; %K
APFU(:,10)=sum(APFU,2); %calculations the total


%% structural formula

%T1 site
StrctFrm(:,1)=APFU(:,1); %Si (T1)

%Al (T1)
for c=1:m
    if 6-StrctFrm(c,1) > APFU(c,3)
        StrctFrm(c,2)=APFU(c,3); 
    else
        StrctFrm(c,2)=6-StrctFrm(c,1);
    end
end

StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2); %Sum of T1

%T2 site
StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2); %Altotal-Al(T1)=Al (T2)
StrctFrm(:,5)=APFU(:,2); %Ti (T2)
StrctFrm(:,6)=StrctFrm(:,5)+StrctFrm(:,4); % sum of T2 site

%Octahedral Site

StrctFrm(:,7)=APFU(:,4); %Fe (Oct)
StrctFrm(:,8)=APFU(:,5); %Mn (Oct)
StrctFrm(:,9)=APFU(:,6); %Mg (Oct)
StrctFrm(:,10)=StrctFrm(:,9)+StrctFrm(:,8)+StrctFrm(:,7); %Sum of Oct site

%A site
StrctFrm(:,11)=APFU(:,7); %Ca (A)
StrctFrm(:,12)=APFU(:,8); %Na (A)
StrctFrm(:,13)=APFU(:,9); %K (A)
StrctFrm(:,14)=StrctFrm(:,13)+StrctFrm(:,12)+StrctFrm(:,11); %A site sum

StrctFrm(:,15)=APFU(:,6)./(APFU(:,6)+APFU(:,4)); %XMg

StrctFrm=array2table(StrctFrm,'VariableNames',{'Si_T1','Al_T1','Sum_T1','Al_T2','Ti_T2','Sum_T2','Fe_Oct','Mn_Oct','Mg_Oct','Sum_Oct','Ca_A','Na_A','K_A','A_Sum','XMg'});
APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Fe','Mn','Mg','Ca','Na','K','Sum'});

end






