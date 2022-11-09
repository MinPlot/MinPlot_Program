%Calculate titanite structural formula 

function [StrctFrm, APFU]=titanite(data,headers)

[m,n]=size(data); %finds the x and y size of the input data matrix

%finds the column position oxide headers in any order to assign data to the correct array
%positions 
I(1,1)=find(strcmp(headers,'SiO2'));
I(1,2)=find(strcmp(headers,'TiO2'));
I(1,3)=find(strcmp(headers,'Al2O3'));

%Makes Y2O3 optional
if strcmp(headers,'Y2O3')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'Y2O3'));
end

I(1,5)=find(strcmp(headers,'FeO'));
I(1,6)=find(strcmp(headers,'MnO'));
I(1,7)=find(strcmp(headers,'MgO'));
I(1,8)=find(strcmp(headers,'CaO'));
I(1,9)=find(strcmp(headers,'Na2O'));

%Makes K2O optional
if strcmp(headers,'K2O')==zeros(1,length(headers))
    I(1,10)=0;
else
    I(1,10)=find(strcmp(headers,'K2O'));
end

%Makes F optional
if strcmp(headers,'F')==zeros(1,length(headers))
    I(1,11)=0;
else
    I(1,11)=find(strcmp(headers,'F'));
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

W=[SiO2,TiO2,Al2O3,Fe2O3,Y2O3,MnO,MgO,CaO,Na2O,K2O,F];


%% moles of cations

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2
MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
MC(:,3)=(data(:,I(1,3))./W(:,3)).*2; %for Al2O3

%calculates for Y2O3 if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=(data(:,I(1,4))./W(:,4)).*2; %for Y2O3
end 

MC(:,5)=((data(:,I(1,5)).*1.111378)./W(:,5)).*2; %converts FeO to Fe2O3
MC(:,6)=data(:,I(1,6))./W(:,6); %for MnO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MgO
MC(:,8)=data(:,I(1,8))./W(:,8); %for CaO
MC(:,9)=(data(:,I(1,9))./W(:,9)).*2; %for Na2O

%calculates for K2O if it is included in the analysis 
if I(1,10)==0
    MC(:,10)=zeros(m,1);
else
    MC(:,10)=(data(:,I(1,10))./W(:,10)).*2; %for K2O
end 

%calculates for F if it is included in the analysis 
if I(1,11)==0
    MC(:,11)=zeros(m,1);
else
    MC(:,11)=data(:,I(1,11))./W(:,11); %for F
end 

%% Normalized moles of cations 

% Normalization is calculated assuming that the tetrahedral and octahedral
% sites are fully occupied by Si, Ti, Al, Fe3+, Mn, and Mg
N_fact=2./(MC(:,1)+MC(:,2)+MC(:,3)+MC(:,5)+MC(:,6)+MC(:,7)); %normalization factor

APFU=MC.*N_fact;
APFU(:,12)=sum(APFU(:,1:10),2);

%% structural formula

%T site
StrctFrm(:,1)=APFU(:,1); %Si (T)

%Al(T)
for c=1:m
    if 1-StrctFrm(c,1)>0 %Is 1-Si > 0? If y, then some Al goes into T
        if 1-StrctFrm(c,1)>APFU(c,3) %For low Al , 1-Si may be > Al
            StrctFrm(c,2)=APFU(c,3); %All Al goes into T
        else
            StrctFrm(c,2)=1-StrctFrm(c,1); %if there isn't enough space in T for all Al, the rest will go to M
        end
    else
        StrctFrm(c,2)=0; %if Si=1, then no Al goes into T
    end
end

%T site sum
StrctFrm(:,3)=StrctFrm(:,1)+StrctFrm(:,2); %should sum to 1

%Octahedral Site
StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2); %Al 
StrctFrm(:,5)=APFU(:,2); %Ti4+
StrctFrm(:,6)=APFU(:,5); %Fe3+ 
StrctFrm(:,7)=APFU(:,6); %Mn 
StrctFrm(:,8)=APFU(:,7); %Mg 

%Decahedral Site
StrctFrm(:,9)=APFU(:,8); %Ca 
StrctFrm(:,10)=APFU(:,4); %Y3+ 
StrctFrm(:,11)=APFU(:,9); %Na
StrctFrm(:,12)=APFU(:,10); %K

%cation sum
StrctFrm(:,13)=sum(StrctFrm(:,1:2),2)+sum(StrctFrm(:,4:12),2);

%anions
StrctFrm(:,14)=APFU(:,11); %F 
StrctFrm(:,15)=(StrctFrm(:,4)+StrctFrm(:,6))-StrctFrm(:,14); %OH = (Al + Fe) - F
StrctFrm(:,16)=5-(10-(4.*APFU(:,1)+4.*APFU(:,2)+3.*APFU(:,3)+3.*APFU(:,4)+3.*APFU(:,5)+2.*APFU(:,6)+2.*APFU(:,7)+2.*APFU(:,8)+APFU(:,9)+APFU(:,10)-0.5.*StrctFrm(:,14)-0.5.*StrctFrm(:,15))); %O2
%O2 anions are calculated as the sum of the charges of cations minus 1/2
%for F + OH

%anion sum
StrctFrm(:,17)=StrctFrm(:,14)+StrctFrm(:,15)+StrctFrm(:,16);

%XTtn
StrctFrm(:,18)=StrctFrm(:,5)./(StrctFrm(:,5)+StrctFrm(:,4)+StrctFrm(:,6)+StrctFrm(:,7)+StrctFrm(:,8));

StrctFrm=array2table(StrctFrm,'VariableNames',{'Si_T','Al_T','Sum_T','Al_Oct','Ti_Oct','Fe3_Oct','Mn_Oct','Mg_Oct','Ca_Dec','Y_Dec','Na_Dec','K_Dec','Cation_Sum','F','OH','O','Anion_Sum','XTtn'});
APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Y','Fe3','Mn','Mg','Ca','Na','K','F','Cation_Sum'});

end



