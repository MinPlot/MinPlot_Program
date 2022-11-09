%%talc Structural Formula

function [StrctFrm, APFU]=talc(data,headers,wantstrctfrm)

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

%Makes CaO optional
if strcmp(headers,'CaO')==zeros(1,length(headers))
    I(1,8)=0;
else
    I(1,8)=find(strcmp(headers,'CaO'));
end

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

W=[SiO2,TiO2,Al2O3,NiO,FeO,MnO,MgO,CaO,Na2O,K2O];

%% Moles of cations

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%calculates for TiO2 if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

MC(:,3)=(data(:,I(1,3))./W(:,3)); %for Al2O3

%calculates for NiO if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=data(:,I(1,4))./W(:,4); %for NiO
end

MC(:,5)=data(:,I(1,5))./W(:,5); %for FeO
MC(:,6)=data(:,I(1,6))./W(:,6); %for MnO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MgO

%calculates for CaO if it is included in the analysis 
if I(1,8)==0
    MC(:,8)=zeros(m,1);
else
    MC(:,8)=data(:,I(1,8))./W(:,8); %for CaO
end

%calculates for Na2O if it is included in the analysis 
if I(1,9)==0
    MC(:,9)=zeros(m,1);
else
    MC(:,9)=data(:,I(1,9))./W(:,9); %for Na2O
end

%calculates for K2O if it is included in the analysis 
if I(1,10)==0
    MC(:,10)=zeros(m,1);
else
    MC(:,10)=data(:,I(1,10))./W(:,8); %for K2O
end

%% Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*3; %for Al2O3
O2(:,4)=MC(:,4); %for NiO
O2(:,5)=MC(:,5); %for FeO
O2(:,6)=MC(:,6); %for MnO
O2(:,7)=MC(:,7); %for MgO
O2(:,8)=MC(:,8); %CaO
O2(:,9)=MC(:,9); %Na2O
O2(:,10)=MC(:,10); %K2O

O2_N=(11)./sum(O2,2); %normalization factor

%normalized moles of anions
N_Ox=O2.*O2_N;

%% atoms pfu

APFU(:,1)=N_Ox(:,1)./2; %Si
APFU(:,2)=N_Ox(:,2)./2; %Ti
APFU(:,3)=N_Ox(:,3).*(2/3); %Al
APFU(:,4)=N_Ox(:,4); %Ni
APFU(:,5)=N_Ox(:,5); %Fe
APFU(:,6)=N_Ox(:,6); %Mn
APFU(:,7)=N_Ox(:,7); %Mg
APFU(:,8)=N_Ox(:,8); %Ca
APFU(:,9)=N_Ox(:,9).*2; %Na
APFU(:,10)=N_Ox(:,10).*2; %K
APFU(:,11)=sum(APFU,2); %calculations the total

%% structural formula

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

StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2); %Sum of T

StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2); %Altotal-Al(T)=Al(M)
StrctFrm(:,5)=APFU(:,2); %Ti (M)
StrctFrm(:,6)=APFU(:,4); %Ni (M)
StrctFrm(:,7)=APFU(:,5); %Fe2+ (M)
StrctFrm(:,8)=APFU(:,6); %Mn (M)
StrctFrm(:,9)=APFU(:,7); %Mg (M)
StrctFrm(:,10)=APFU(:,8); %Ca (M)
StrctFrm(:,11)=APFU(:,9); %Na (M)
StrctFrm(:,12)=APFU(:,10); %K (M)
StrctFrm(:,13)=APFU(:,11); %Sum

StrctFrm=array2table(StrctFrm,'VariableNames',{'Si_T','Al_T','Sum_T','Al_M','Ti_M','Ni_M','Fe_M','Mn_M','Mg_M','Ca_M','Na_M','K_M','Total_Cations'});
APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Ni','Fe','Mn','Mg','Ca','Na','K','Sum'});

end


