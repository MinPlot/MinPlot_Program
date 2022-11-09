function APFU=Amph_KA(APFU_Fe,O2_Ti)
% Normalize to [Si + Al + Ti + Fe + Mn + Mg + Ca + Na] = 15 cations

[m,n]=size(APFU_Fe); %finds the x and y size of the input data matrix

Cat=APFU_Fe(:,1:11);

Norm(:,1)=15./sum(Cat(:,1:1:10),2); %Normalize to 15 cations

N_Cat=Cat.*Norm; %normalized cation units

%Oxygen units
O2(:,1)=N_Cat(:,1).*2; %Si
O2(:,2)=N_Cat(:,2).*2; %Ti
O2(:,3)=N_Cat(:,3).*(3/2); %Al
O2(:,4)=N_Cat(:,4).*(3/2); %Cr
O2(:,5)=N_Cat(:,5).*(3/2); %Fe3+ (blank)
O2(:,6)=N_Cat(:,6); %Fe2+
O2(:,7)=N_Cat(:,7); %Mn
O2(:,8)=N_Cat(:,8); %Mg
O2(:,9)=N_Cat(:,9); %Ca
O2(:,10)=N_Cat(:,10)./2; %Na
O2(:,11)=N_Cat(:,11)./2; %K


O2total=sum(O2,2); %Sum of oxygen units 

APFU=N_Cat; 

%correction for electroneutrality 
for c=1:m
    if (O2_Ti(c,1)-O2total(c,1)) > 0 %if the sum of the O equivalents minus the O2 sum after
        %normalization is >1, then Fe3+ may be calculated 
        if O2(c,6) > 2*(O2_Ti(c,1)-O2total(c,1)) %if Fe2+ anions are > than 2*(sum of the 
            %O equivalents - the O2 sum after normalisation
            APFU(c,5)=2*(O2_Ti(c,1)-O2total(c,1)); %then only some Fe3+ is calculated 
            APFU(c,6)=N_Cat(c,6)-APFU(c,5); %the Fe2+ is whatever is left
        else
            APFU(c,5)=N_Cat(c,6); %if no, all Fe is Fe3+
            APFU(c,6)=0; %Fe2+ zero
        end
    else
        APFU(c,5)=N_Cat(c,5); %No Fe3+ is calculated
        APFU(c,6)=N_Cat(c,6); %Fe2+ = Fetotal
    end
end


APFU(:,12)=8./APFU(:,1); %Criteria 1-1: Si cannot be more than 8 cations

APFU(:,13)=16./sum(APFU(:,1:1:11),2); %Criteria 1-2: Should be <= 16 APFU

APFU(:,14)=15./sum(APFU(:,1:1:9),2); %Criteria 1-3: Should be <=15 APFU

APFU(:,15)=8./(APFU(:,1)+APFU(:,3)); % Criteria 2-1: Si + Al = 8

APFU(:,16)=15./sum(APFU(:,1:1:10),2); %Criteria 2-2: Should be equal to 15 APFU

APFU(:,17)=13./sum(APFU(:,1:1:8),2); %Criteria 2-3: Should be equal to 13 APFU

APFU(:,18)=36./(46-sum(APFU(:,1:4),2)); %Criteria 2-4: tetrahedral sites completely filled by 3+ and 4+ cation

APFU(:,19)=Norm(:); 

end