function [O2_Ti,APFU_I,MC]=Amph_O2Ti(D,W,wantTiOH)

%Determines normalization for the the O2-Ti adjustment 

[m,n]=size(D); %finds the x and y size of the input data matrix

%moles of cations 
MC(:,1)=D(:,1)./W(:,1); %SiO2
MC(:,2)=D(:,2)./W(:,2); %TiO2
MC(:,3)=(D(:,3)./W(:,3)).*2; %Al2O3
MC(:,4)=(D(:,4)./W(:,4)).*2; %Cr2O3
MC(:,5)=(D(:,5)./W(:,5)).*2; %Fe2O3
MC(:,6)=D(:,6)./W(:,6); %FeO
MC(:,7)=D(:,7)./W(:,7); %MnO
MC(:,8)=D(:,8)./W(:,8); %MgO
MC(:,9)=D(:,9)./W(:,9); %CaO
MC(:,10)=(D(:,10)./W(:,10)).*2; %Na2O
MC(:,11)=(D(:,11)./W(:,11)).*2; %K2O
MC(:,12)=(D(:,12)./W(:,12)).*2; %H2O
MC(:,13)=D(:,13)./W(:,13); %F
MC(:,14)=D(:,14)./W(:,14); %Cl

%Moles of O2
O2(:,1)=MC(:,1)*2; %SiO2
O2(:,2)=MC(:,2)*2; %TiO2
O2(:,3)=MC(:,3)*1.5; %Al2O3
O2(:,4)=MC(:,4)*1.5; %Cr2O3
O2(:,5)=MC(:,5)*1.5; %Fe2O3
O2(:,6)=MC(:,6); %FeO
O2(:,7)=MC(:,7); %MnO
O2(:,8)=MC(:,8); %MgO
O2(:,9)=MC(:,9); %CaO
O2(:,10)=MC(:,10)*0.5; %Na2O
O2(:,11)=MC(:,11)*0.5; %K2O
O2(:,12)=MC(:,12)*0.5; %H2O
O2(:,13)=MC(:,13); %F
O2(:,14)=MC(:,14); %Cl

%correction of O2 sum for F and Cl
Anion_Sum=sum(O2,2)-(0.5*O2(:,13))-(0.5*O2(:,14)); %anion sum, O2 equivalents minus 1/2 F and 1/2 Cl
Anhydrous_Sum=Anion_Sum-O2(:,12)-(0.5*O2(:,13))-(0.5*O2(:,14)); %sum not including F or OH
F_Cl=((O2(:,13)+O2(:,14))./Anion_Sum).*24; %F+Cl per 24 anions

%correction of O2 sum for Ti-O substitution
Ti_adj(:,1)=zeros(m,1);
for c=1:m 
    %if the moles of Ti > 8-(Si + Al) & 8-(Si + Al) > 0, then some Ti goes into the T site 
    if (((MC(c,2).*24)./Anion_Sum(c,1)) > (8-(((MC(c,1)+MC(c,3)).*24)./Anion_Sum(c,1)))) && ((8-(((MC(c,1)+MC(c,3)).*24)./Anion_Sum(c,1))) > 0)
        %The Ti in T is subtracted from the Ti that is used to estimate OH
        Ti_adj(c,1)=((MC(c,2).*24)./Anion_Sum(c,1)) - (8-(((MC(c,1)+MC(c,3)).*24)./Anion_Sum(c,1)));
    else
        Ti_adj(c,1)=(MC(c,2).*24)./Anion_Sum(c,1); %all Ti is used to estimate OH
    end
end

%decides whether to correct for OH=2-2Ti
if strcmp(wantTiOH, 'y')
    %decides whether to use the O2 equivalent from the Ti-O estimate or
    %24-(F+Cl)
    O2_Ti(:,1)=zeros(m,1);
    for c=1:m
        if (23+Ti_adj(c,1))>24
            O2_Ti(c,1)=24-0.5*F_Cl(c,1);
        else
            O2_Ti(c,1)=23+Ti_adj(c,1);
        end
    end
else
    O2_Ti(:,1)=zeros(m,1);
    for c=1:m
        O2_Ti(c,1)=23;
    end
end

O2_N=O2_Ti./Anhydrous_Sum;
APFU_I=MC.*O2_N; %normalized moles of cations + F & Cl

end