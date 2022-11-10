%%Calculates amphibole formula following Hawthorne et al. (2012) and
%Leake et al., (1997)

function [Strct_Frm, APFU, APFU_FeT, Fe3_limits, Fe3_class, APFU_SiT, APFU_AfullT, APFU_NaAT, APFU_Fe2O3T, APFU_SiAlT, APFU_KAT, APFU_CaNaBT]=amphibole(data,headers,wantstrctfrm,wantTiOH,wantplot,wantferric)

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
I(1,9)=find(strcmp(headers,'Na2O'));
I(1,10)=find(strcmp(headers,'K2O'));

%Makes F optional
if strcmp(headers,'F')==zeros(1,length(headers))
    I(1,11)=0;
else
    I(1,11)=find(strcmp(headers,'F'));
end

%Makes Cl optional
if strcmp(headers,'Cl')==zeros(1,length(headers))
    I(1,12)=0;
else
    I(1,12)=find(strcmp(headers,'Cl'));
end


%% makes a new array, D, where the oxides are in the correct order and
%filling in zeros where optional oxides are not included 

D(:,1)=data(:,I(1,1)); %SiO2

%adds TiO2 if it is included in the analysis 
if I(1,2)==0
    D(:,2)=zeros(m,1);
else
    D(:,2)=data(:,I(1,2));
end 

D(:,3)=data(:,I(1,3)); %Al2O3


%adds Cr2O3 if it is included in the analysis 
if I(1,4)==0
    D(:,4)=zeros(m,1);
else
    D(:,4)=data(:,I(1,4));
end 

D(:,5)=zeros(m,1); %Fe2O3, all zeros (place holder)
D(:,6)=data(:,I(1,5)); %FeO
D(:,7)=data(:,I(1,6)); %MnO
D(:,8)=data(:,I(1,7)); %MgO
D(:,9)=data(:,I(1,8)); %CaO
D(:,10)=data(:,I(1,9)); %Na2O
D(:,11)=data(:,I(1,10)); %K2O
D(:,12)=zeros(m,1); %H2O, all zeros (place holder)

%adds F if it is included in the analysis 
if I(1,11)==0
    D(:,13)=zeros(m,1);
else
    D(:,13)=data(:,I(1,11));
end 

%adds Cl if it is included in the analysis 
if I(1,12)==0
    D(:,14)=zeros(m,1);
else
    D(:,14)=data(:,I(1,12));
end 

% Molecular weights

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
H2O=18.01528;

W=[SiO2,TiO2,Al2O3,Cr2O3,Fe2O3,FeO,MnO,MgO,CaO,Na2O,K2O,H2O,F,Cl];


%% Initial Fe3+ estimation

%calculate the initial APFU of cations and anions 

[O2_Ti,APFU_I,MC]=Amph_O2Ti(D,W,wantTiOH); %calls the initial amphibole APFU calculation

%Lower Limits on Fe3+:

%Atoms per formula unit (APFU) assuming that all Fe is FeO
APFU_Fe=Amph_Fe(D,W,O2_Ti);

if strcmp(wantferric, 'y')
    %APFU assuming 8 Si cations
    APFU_Si=Amph_Si(APFU_Fe,O2_Ti);
    
    %APFU assuming 16 cations (no vacancies on A)
    APFU_Afull=Amph_Afull(APFU_Fe,O2_Ti);
    
    %APFU assuming Na in A site only
    APFU_NaA=Amph_NaA(APFU_Fe,O2_Ti);
    
    APFU_low=zeros(m,11);
    % Fe3+ minimum estimate: highest of 3 normalization factors
    for c=1:m
        %if the 3 minima criteria have lower normalization factors than the
        %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
        if find(APFU_Fe(c,15:17)>APFU_Fe(c,12:14))
            APFU_low(c,:)=APFU_Fe(c,1:1:11);
            low_check{c,1}='Fe3+ cannot be estimated';
        else %otherwise,a Fe3+ lower estimate can be made
            %if the normalization factors are all greater than 1, then all Fe is assumed to be FeO
            if (APFU_Fe(c,12) > 1) && (APFU_Fe(c,13) > 1) && (APFU_Fe(c,14) > 1)
                APFU_low(c,:)=APFU_Fe(c,1:1:11);
                low_check{c,1}={'All Fe is FeO'};
            else
                %If (8/Si) < (16/sum(all cations)) &  < (15/sum(cations Si to Ca))
                %, then the Si normalized formula is correct
                if (APFU_Fe(c,12) < APFU_Fe(c,13)) && (APFU_Fe(c,12) < APFU_Fe(c,14))
                    APFU_low(c,:)=APFU_Si(c,1:1:11);
                    low_check{c,1}={'Criterion 1-1: 8 Si cations on T'};
                else %otherwise, consider the reamining two options
                    
                    %If 16/sum(all cations) < (8/Si) & < (15/sum(cations Si to Ca))
                    %, then the formula normalized to 16 cations is correct
                    if (APFU_Fe(c,13) < APFU_Fe(c,12)) && (APFU_Fe(c,13) < APFU_Fe(c,14))
                        APFU_low(c,:)=APFU_Afull(c,1:1:11);
                        low_check{c,1}={'Criterion 1-2: 16 total cations (no vac. on A)'};
                    else
                        %otherwise, 15/sum(cations Si to Ca) is the only remaining
                        %option for a normalization factor
                        APFU_low(c,:)=APFU_NaA(c,1:1:11);
                        low_check{c,1}={'Criterion 1-3: Na on A site only'};
                        
                    end
                end
            end
        end
    end
    
    
    % Upper Limits on Fe3+
    %APFU assuming all Fe is Fe2O3
    APFU_Fe2O3=Amph_Fe2O3(D,W,O2_Ti);
    O2_Nfact=APFU_Fe2O3(:,19)./APFU_Fe(:,19); %Normalization factor for all Fe as Fe2O3 relative to all Fe as FeO
    
    %APFU assuming Al + Si in T site only
    APFU_SiAlT=Amph_SiAlT(APFU_Fe,O2_Ti);
    
    %APFU assuming K only on the A site
    APFU_KA=Amph_KA(APFU_Fe,O2_Ti);
    
    %APFU assuming fully occupied T & B sites
    APFU_CaNaB=Amph_CaNaB(APFU_Fe,O2_Ti);
    
    %APFU assuming tetrahedral sites completely filled by 3+ and 4+ cations
    APFU_10Fe3=Amph_10Fe3(APFU_Fe,O2_Ti);
    
    
    %Fe3+ maximum estimate: highest of 4 normalization factors
    APFU_hi=zeros(m,11);
    for c=1:m
        
        %if the 3 minima criteria have lower normalization factors than the
        %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
        if find(APFU_Fe(c,15:17)>APFU_Fe(c,12:14))
            APFU_hi(c,:)=APFU_Fe(c,1:1:11);
            hi_check{c,1}={'Fe3+ cannot be estimated'};
        else
            %If the all Fe3+ normalization factor is grater than the other 
            if (O2_Nfact(c,1) > APFU_Fe(c,15)) && (O2_Nfact(c,1) > APFU_Fe(c,16)) && (O2_Nfact(c,1) > APFU_Fe(c,17)) && (O2_Nfact(c,1) > APFU_Fe(c,18))
                APFU_hi(c,:)=APFU_Fe2O3(c,1:1:11);
                hi_check{c,1}={'Criterion 2-5: All Fe is Fe2O3'};
            else
                if (APFU_Fe(c,15) > O2_Nfact(c,1)) && (APFU_Fe(c,15) > APFU_Fe(c,16)) && (APFU_Fe(c,15) > APFU_Fe(c,17)) && (APFU_Fe(c,15) > APFU_Fe(c,18))
                    APFU_hi(c,:)=APFU_SiAlT(c,1:1:11);
                    hi_check{c,1}={'Criterion 2-1: Al + Si normalized to 8 in T site'};
                else
                    if (APFU_Fe(c,16) > O2_Nfact(c,1)) && (APFU_Fe(c,16) > APFU_Fe(c,15)) && (APFU_Fe(c,16) > APFU_Fe(c,17)) && (APFU_Fe(c,16) > APFU_Fe(c,18))
                        APFU_hi(c,:)=APFU_KA(c,1:1:11);
                        hi_check{c,1}={'Criterion 2-2: K only on the A site, fully occupied T, B, & C sites'};
                    else
                        if (APFU_Fe(c,17) > O2_Nfact(c,1)) && (APFU_Fe(c,17) > APFU_Fe(c,15)) && (APFU_Fe(c,17) > APFU_Fe(c,16)) && (APFU_Fe(c,17) > APFU_Fe(c,18))
                            APFU_hi(c,:)=APFU_CaNaB(c,1:1:11);
                            hi_check{c,1}={'Criterion 2-3: Fully occupied T & B sites'};
                        else
                            APFU_hi(c,:)=APFU_10Fe3(c,1:1:11);
                            hi_check{c,1}={'Criterion 2-4: tetrahedral sites completely filled by 3+ and 4+ cations'};
                        end
                    end
                end
            end
        end
        
    end
    
    % calculate the formula with Fe3+/Fetotal determination from high and low
    % limits
    APFU=(APFU_hi+APFU_low)./2;
    
    %Adds OH, scales for subsequent itterations
    for c=1:m
        if APFU_I(c,12) > 0
            APFU(c,12) = APFU_I(c,12); %OH from the normalized moles cations
        else
            APFU(c,12) = MC(c,12).*(MC(c,1)./APFU(c,1)); %OH is scaled based on the change in Si content
        end
    end
    
    %Adds F
    for c=1:m
        if APFU_I(c,13) > 0
            APFU(c,13) = APFU_I(c,13); %F from the normalized moles cations
        else
            APFU(c,13) = MC(c,13).*(MC(c,1)./APFU(c,1)); %F is scaled based on the change in Si content
        end
    end
    
    %Adds Cl
    for c=1:m
        if APFU_I(c,14) > 0
            APFU(c,14) = APFU_I(c,14); %Cl from the normalized moles cations
        else
            APFU(c,14) = MC(c,14).*(MC(c,1)./APFU(c,1)); %Cl is scaled based on the change in Si content
        end
    end
    
    for c=1:m
        if (APFU(c,12)+APFU(c,13)+APFU(c,14)) > 2 %if (OH_initial + F + Cl) > 2
            if (APFU(c,13)+APFU(c,14)) < 2 % if (F + Cl) < 2
                OH_1(c,1)=2-(APFU(c,13)+APFU(c,14)); %initial OH guess is 2-(F+Cl)
            else
                OH_1(c,1)=0; %W site is filled by Cl and F
            end
        else
            if (APFU(c,12)+APFU(c,13)+APFU(c,14)) < 2 %if (OH_initial + F + Cl) < 2
                if (2-APFU(c,13)-APFU(c,14)) > 0 %if (2-F-Cl) < 2
                    OH_1(c,1)=(2-APFU(c,13)-APFU(c,14)); % OH = 2-F-Cl
                else
                    OH_1(c,1)=0; %W site is filled by Cl and F
                end
            else
                OH_1(c,1)=APFU(c,12);
            end
        end
    end
    
    %OH estimation, Step 2: Adjustment for OH=2-2Ti
    
    %some Ti goes into the T site instead of the M site
    for c=1:m
        %if Ti > (8-(Si+Al)) & (8-(Si+Al))>0
        if APFU(c,2)>(8-(APFU(c,1)+APFU(c,3))) && (8-(APFU(c,1)+APFU(c,3))) > 0
            Ticorr(c,1)=APFU(c,2)-(8-(APFU(c,1)+APFU(c,3))); %Ti in M = Titotal - Ti in T
        else
            Ticorr(c,1)=APFU(c,2); %otherwise all Ti is in M
        end
    end
    
    % Determine appropriate OH correction
    for c=1:m
        if strcmp(wantTiOH, 'y')
            if (APFU(c,12) + APFU(c,13) + APFU(c,14)) > 2 %If (OH + F + Cl) > 2
                if ((2.*APFU(c,12))./((APFU(c,12)+ APFU(c,13) + APFU(c,14))-2.*APFU(c,2)))>0 %if ((2*OH)/((OH+F+Cl)-2*Ti)) >0
                    APFU(c,12)=((2.*APFU(c,12))./((APFU(c,12)+ APFU(c,13) + APFU(c,14))-2.*APFU(c,2)));
                else
                    APFU(c,12)=0; %W is filled by F, Cl, and O
                end
            else
                if (APFU(c,12) + APFU(c,13) + APFU(c,14)) < 2 %If (OH + F + Cl) < 2
                    if ((2-APFU(c,13)-APFU(c,14))-2.*Ticorr(c,1))>0 %if 2-F-Cl-2*Ti > 0
                        APFU(c,12)=((2-APFU(c,13)-APFU(c,14))-2.*Ticorr(c,1)); %OH=(2-F-Cl-2*Ti)
                    else
                        APFU(c,12)=0; %W is filled by F, Cl, and O
                    end
                else
                    APFU(c,12)=APFU(c,12); %OH content is not adjusted for Ti
                end
            end
        else
            APFU(c,12)=OH_1(c,1);
        end
    end
    
    %Moles of anions
    APFU_an(:,1)=APFU(:,1).*2; %SiO2
    APFU_an(:,2)=APFU(:,2).*2; %TiO2
    APFU_an(:,3)=APFU(:,3).*1.5; %Al2O3
    APFU_an(:,4)=APFU(:,4).*1.5; %Cr2O3
    APFU_an(:,5)=APFU(:,5).*1.5; %Fe2O3
    APFU_an(:,6)=APFU(:,6); %FeO
    APFU_an(:,7)=APFU(:,7); %MnO
    APFU_an(:,8)=APFU(:,8); %MgO
    APFU_an(:,9)=APFU(:,9); %CaO
    APFU_an(:,10)=APFU(:,10).*0.5; %Na2O
    APFU_an(:,11)=APFU(:,11).*0.5; %K2O
    APFU_an(:,12)=APFU(:,12).*0.5; %H2O
    
    %for F
    for c=1:m
        if (APFU(c,13)+APFU(c,14))>2 %F + Cl cannot > 2
            APFU_an(c,13)=2.*(APFU(c,13)./(APFU(c,13)+APFU(c,14))); %scales F
        else
            APFU_an(c,13)=APFU(c,13); %Does not scale F
        end
    end
    
    %for Cl
    for c=1:m
        if (APFU(c,13)+APFU(c,14))>2 %F + Cl cannot > 2
            APFU_an(c,14)=2.*(APFU(c,14)./(APFU(c,13)+APFU(c,14))); %scales Cl
        else
            APFU_an(c,14)=APFU(c,14); %Does not scale Cl
        end
    end
    
    
    %the anion sum may be > 24
    anion_sum=sum(APFU_an(:,1:12),2)+0.5*APFU_an(:,13)+0.5*APFU_an(:,14);
    anion_norm=24./anion_sum;
    
    APFU_n=APFU.*anion_norm; %cations normalized to 24 anions
    APFU_ann=APFU_an.*anion_norm; %anions normalized to 24 anions
    Fe3_ratio(:,1)=APFU(:,5)./(APFU(:,5)+APFU(:,6));
    
    molar_mass=APFU_n(:,1).*SiO2+APFU_n(:,2).*TiO2+(APFU_n(:,3)./2).*Al2O3+(APFU_n(:,4)./2).*Cr2O3+(APFU_n(:,5)./2).*Fe2O3+APFU_n(:,6).*FeO+APFU_n(:,7).*MnO+APFU_n(:,8).*MgO+APFU_n(:,9).*CaO+(APFU_n(:,10)./2).*Na2O+(APFU_n(:,11)./2).*K2O+(APFU_n(:,12)./2).*H2O+APFU_ann(:,13).*F+APFU_ann(:,14).*Cl-APFU_ann(:,13).*15.9994.*0.5-APFU_ann(:,14).*15.9994.*0.5;
    
    %new wt% of oxides
    D2=D;
    D2(:,5)=Fe3_ratio(:,1).*(D(:,5).*(2*FeO/Fe2O3)+D(:,6))./(2*FeO/Fe2O3);
    D2(:,6)=(1-Fe3_ratio(:,1)).*(D(:,5).*(2*FeO/Fe2O3)+D(:,6));
    D2(:,12)=APFU_n(:,12).*H2O.*0.5.*(1./molar_mass)*100;
    
    %% Itterates the calculation 10 times
    for z=1:10
        %calculate the initial APFU of cations and anions
        
        [O2_Ti2,APFU_I2,MC2]=Amph_O2Ti(D2,W,wantTiOH); %calls the initial amphibole APFU calculation
        
        %Lower Limits on Fe3+:
        
        %Atoms per formula unit (APFU) assuming that all Fe is FeO
        APFU_Fe2=Amph_Fe(D2,W,O2_Ti2);
        
        %APFU assuming 8 Si cations
        APFU_Si2=Amph_Si(APFU_Fe2,O2_Ti2);
        
        %APFU assuming 16 cations (no vacancies on A)
        APFU_Afull2=Amph_Afull(APFU_Fe2,O2_Ti2);
        
        %APFU assuming Na in A site only
        APFU_NaA2=Amph_NaA(APFU_Fe2,O2_Ti2);
        
        APFU_low2=zeros(m,11);
        % Fe3+ minimum estimate: highest of 3 normalization factors
        for c=1:m
            %if the 3 minima criteria have lower normalization factors than the
            %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
            if find(APFU_Fe2(c,15:17)>APFU_Fe2(c,12:14))
                APFU_low2(c,:)=APFU_Fe2(c,1:1:11);
                low_check2{c,1}='Fe3+ cannot be estimated';
            else %otherwise,a Fe3+ lower estimate can be made
                %if the normalization factors are all greater than 1, then all Fe is assumed to be FeO
                if (APFU_Fe2(c,12) > 1) && (APFU_Fe2(c,13) > 1) && (APFU_Fe2(c,14) > 1)
                    APFU_low2(c,:)=APFU_Fe2(c,1:1:11);
                    low_check2{c,1}={'All Fe is FeO'};
                else
                    %If (8/Si) < (16/sum(all cations)) &  < (15/sum(cations Si to Ca))
                    %, then the Si normalized formula is correct
                    if (APFU_Fe2(c,12) < APFU_Fe2(c,13)) && (APFU_Fe2(c,12) < APFU_Fe2(c,14))
                        APFU_low2(c,:)=APFU_Si2(c,1:1:11);
                        low_check2{c,1}={'Criterion 1-1: 8 Si cations on T'};
                    else %otherwise, consider the reamining two options
                        
                        %If 16/sum(all cations) < (8/Si) & < (15/sum(cations Si to Ca))
                        %, then the formula normalized to 16 cations is correct
                        if (APFU_Fe2(c,13) < APFU_Fe2(c,12)) && (APFU_Fe2(c,13) < APFU_Fe2(c,14))
                            APFU_low2(c,:)=APFU_Afull2(c,1:1:11);
                            low_check2{c,1}={'Criterion 1-2: 16 total cations (no vac. on A)'};
                        else
                            %otherwise, 15/sum(cations Si to Ca) is the only remaining
                            %option for a normalization factor
                            APFU_low2(c,:)=APFU_NaA2(c,1:1:11);
                            low_check2{c,1}={'Criterion 1-3: Na on A site only'};
                            
                        end
                    end
                end
            end
        end
        
        
        % Upper Limits on Fe3+
        %APFU assuming all Fe is Fe2O3
        APFU_Fe2O32=Amph_Fe2O3(D2,W,O2_Ti2);
        O2_Nfact2=APFU_Fe2O32(:,19)./APFU_Fe2(:,19); %Normalization factor for all Fe as Fe2O3 relative to all Fe as FeO
        
        %APFU assuming Al + Si in T site only
        APFU_SiAlT2=Amph_SiAlT(APFU_Fe2,O2_Ti2);
        
        %APFU assuming K only on the A site
        APFU_KA2=Amph_KA(APFU_Fe2,O2_Ti2);
        
        %APFU assuming fully occupied T & B sites
        APFU_CaNaB2=Amph_CaNaB(APFU_Fe2,O2_Ti2);
        
        %APFU assuming tetrahedral sites completely filled by 3+ and 4+ cations
        APFU_10Fe32=Amph_10Fe3(APFU_Fe2,O2_Ti2);
        
        
        %Fe3+ maximum estimate: highest of 4 normalization factors
        APFU_hi2=zeros(m,11);
        for c=1:m
            
            %if the 3 minima criteria have lower normalization factors than the
            %four maximum criteria, then Fe3+ cannot be estimated and FeTotal=Fe2+
            if find(APFU_Fe2(c,15:17)>APFU_Fe2(c,12:14))
                APFU_hi2(c,:)=APFU_Fe2(c,1:1:11);
                hi_check2{c,1}={'Fe3+ cannot be estimated'};
            else
                %If the all Fe3+ normalization factor is grater than the other th
                if (O2_Nfact2(c,1) > APFU_Fe2(c,15)) && (O2_Nfact2(c,1) > APFU_Fe2(c,16)) && (O2_Nfact2(c,1) > APFU_Fe2(c,17)) && (O2_Nfact2(c,1) > APFU_Fe2(c,18))
                    APFU_hi2(c,:)=APFU_Fe2O32(c,1:1:11);
                    hi_check2{c,1}={'Criterion 2-5: All Fe is Fe2O3'};
                else
                    if (APFU_Fe2(c,15) > O2_Nfact2(c,1)) && (APFU_Fe2(c,15) > APFU_Fe2(c,16)) && (APFU_Fe2(c,15) > APFU_Fe2(c,17)) && (APFU_Fe2(c,15) > APFU_Fe2(c,18))
                        APFU_hi2(c,:)=APFU_SiAlT2(c,1:1:11);
                        hi_check2{c,1}={'Criterion 2-1: Al + Si normalized to 8 in T site'};
                    else
                        if (APFU_Fe2(c,16) > O2_Nfact2(c,1)) && (APFU_Fe2(c,16) > APFU_Fe2(c,15)) && (APFU_Fe2(c,16) > APFU_Fe2(c,17)) && (APFU_Fe2(c,16) > APFU_Fe2(c,18))
                            APFU_hi2(c,:)=APFU_KA2(c,1:1:11);
                            hi_check2{c,1}={'Criterion 2-2: K only on the A site, fully occupied T, B, & C sites'};
                        else
                            if (APFU_Fe2(c,17) > O2_Nfact2(c,1)) && (APFU_Fe2(c,17) > APFU_Fe2(c,15)) && (APFU_Fe2(c,17) > APFU_Fe2(c,16)) && (APFU_Fe2(c,17) > APFU_Fe2(c,18))
                                APFU_hi2(c,:)=APFU_CaNaB2(c,1:1:11);
                                hi_check2{c,1}={'Criterion 2-3: Fully occupied T & B sites'};
                            else
                                APFU_hi2(c,:)=APFU_10Fe32(c,1:1:11);
                                hi_check2{c,1}={'Criterion 2-4: tetrahedral sites completely filled by 3+ and 4+ cations'};
                            end
                        end
                    end
                end
            end
        end
        
        % calculate the formula with Fe3+/Fetotal determination from high and low
        % limits
        APFU2=(APFU_hi2+APFU_low2)./2;
        
        %Adds OH
        for c=1:m
            if APFU_I2(c,12) > 0
                APFU2(c,12) = APFU_I2(c,12); %OH from the normalized moles cations
            else
                APFU2(c,12) = MC2(c,12).*(MC2(c,1)./APFU2(c,1)); %OH is scaled based on the change in Si content
            end
        end
        
        %Adds F
        for c=1:m
            if APFU_I2(c,13) > 0
                APFU2(c,13) = APFU_I2(c,13); %F from the normalized moles cations
            else
                APFU2(c,13) = MC2(c,13).*(MC2(c,1)./APFU2(c,1)); %F is scaled based on the change in Si content
            end
        end
        
        %Adds Cl
        for c=1:m
            if APFU_I2(c,14) > 0
                APFU2(c,14) = APFU_I2(c,14); %Cl from the normalized moles cations
            else
                APFU2(c,14) = MC2(c,14).*(MC2(c,1)./APFU2(c,1)); %Cl is scaled based on the change in Si content
            end
        end
        
        %initial OH estimate
        for c=1:m
            if (APFU2(c,12)+APFU2(c,13)+APFU2(c,14)) > 2 %if (OH_initial + F + Cl) > 2
                if (APFU2(c,13)+APFU2(c,14)) < 2 % if (F + Cl) < 2
                    OH_2(c,1)=2-(APFU2(c,13)+APFU2(c,14)); %initial OH guess is 2-(F+Cl)
                else
                    OH_2(c,1)=0; %W site is filled by Cl and F
                end
            else
                if (APFU2(c,12)+APFU2(c,13)+APFU2(c,14)) < 2 %if (OH_initial + F + Cl) < 2
                    if (2-APFU2(c,13)-APFU2(c,14)) > 0 %if (2-F-Cl) < 2
                        OH_2(c,1)=(2-APFU2(c,13)-APFU2(c,14)); % OH = 2-F-Cl
                    else
                        OH_2(c,1)=0; %W site is filled by Cl and F
                    end
                else
                    OH_2(c,1)=APFU2(c,12);
                end
            end
        end
        
        %OH estimation, Step 2: Adjustment for OH=2-2Ti
        
        %some Ti goes into the T site instead of the M site
        for c=1:m
            %if Ti > (8-(Si+Al)) & (8-(Si+Al))>0
            if APFU2(c,2)>(8-(APFU2(c,1)+APFU2(c,3))) && (8-(APFU2(c,1)+APFU2(c,3))) > 0
                Ticorr2(c,1)=APFU2(c,2)-(8-(APFU2(c,1)+APFU2(c,3))); %Ti in M = Titotal - Ti in T
            else
                Ticorr2(c,1)=APFU2(c,2); %otherwise all Ti is in M
            end
        end
        
        %Determine appropriate OH correction
        for c=1:m
            if strcmp(wantTiOH, 'y')
                if (APFU2(c,12) + APFU2(c,13) + APFU2(c,14)) > 2 %If (OH + F + Cl) > 2
                    if ((2.*APFU2(c,12))./((APFU2(c,12)+ APFU2(c,13) + APFU2(c,14))-2.*APFU2(c,2)))>0 %if ((2*OH)/((OH+F+Cl)-2*Ti)) >0
                        APFU2(c,12)=((2.*APFU2(c,12))./((APFU2(c,12)+ APFU2(c,13) + APFU2(c,14))-2.*APFU2(c,2)));
                    else
                        APFU2(c,12)=0; %W is filled by F, Cl, and O
                    end
                else
                    if (APFU2(c,12) + APFU2(c,13) + APFU2(c,14)) < 2 %If (OH + F + Cl) < 2
                        if ((2-APFU2(c,13)-APFU2(c,14))-2.*Ticorr2(c,1))>0 %if 2-F-Cl-2*Ti > 0
                            APFU2(c,12)=((2-APFU2(c,13)-APFU2(c,14))-2.*Ticorr2(c,1)); %OH=(2-F-Cl-2*Ti)
                        else
                            APFU2(c,12)=0; %W is filled by F, Cl, and O
                        end
                    else
                        APFU2(c,12)=APFU2(c,12); %OH content is not adjusted for Ti
                    end
                end
            else
                APFU2(c,12)=OH_2(c,1);
            end
        end
        
        %Moles of anions
        APFU_an2(:,1)=APFU2(:,1)*2; %SiO2
        APFU_an2(:,2)=APFU2(:,2)*2; %TiO2
        APFU_an2(:,3)=APFU2(:,3)*1.5; %Al2O3
        APFU_an2(:,4)=APFU2(:,4)*1.5; %Cr2O3
        APFU_an2(:,5)=APFU2(:,5)*1.5; %Fe2O3
        APFU_an2(:,6)=APFU2(:,6); %FeO
        APFU_an2(:,7)=APFU2(:,7); %MnO
        APFU_an2(:,8)=APFU2(:,8); %MgO
        APFU_an2(:,9)=APFU2(:,9); %CaO
        APFU_an2(:,10)=APFU2(:,10)*0.5; %Na2O
        APFU_an2(:,11)=APFU2(:,11)*0.5; %K2O
        APFU_an2(:,12)=APFU2(:,12)*0.5; %H2O
        
        %for F
        for c=1:m
            if (APFU2(c,13)+APFU2(c,14))>2 %F + Cl cannot > 2
                APFU_an2(c,13)=2.*(APFU2(c,13)./(APFU2(c,13)+APFU2(c,14))); %scales F
            else
                APFU_an2(c,13)=APFU2(c,13); %Does not scale F
            end
        end
        
        %for Cl
        for c=1:m
            if (APFU2(c,13)+APFU2(c,14))>2 %F + Cl cannot > 2
                APFU_an2(c,14)=2.*(APFU2(c,14)./(APFU2(c,13)+APFU2(c,14))); %scales Cl
            else
                APFU_an2(c,14)=APFU2(c,14); %Does not scale Cl
            end
        end
        
        %the anionn sum may be > 24
        anion_sum2=sum(APFU_an2(:,1:12),2)+0.5*APFU_an2(:,13)+0.5*APFU_an2(:,14);
        anion_norm2=24./anion_sum2;
        
        APFU_n2=APFU2.*anion_norm2; %cations normalized to 24 anions
        APFU_ann2=APFU_an2.*anion_norm2; %anions normalized to 24 anions
        Fe3_ratio2(:,1)=APFU_n2(:,5)./(APFU_n2(:,5)+APFU_n2(:,6));
        
        molar_mass2=APFU_n2(:,1).*SiO2+APFU_n2(:,2).*TiO2+(APFU_n2(:,3)./2).*Al2O3+(APFU_n2(:,4)./2).*Cr2O3+(APFU_n2(:,5)./2).*Fe2O3+APFU_n2(:,6).*FeO+APFU_n2(:,7).*MnO+APFU_n2(:,8).*MgO+APFU_n2(:,9).*CaO+(APFU_n2(:,10)./2).*Na2O+(APFU_n2(:,11)./2).*K2O+(APFU_n2(:,12)./2).*H2O+APFU_ann2(:,13).*F+APFU_ann2(:,14).*Cl-APFU_ann2(:,13).*15.9994.*0.5-APFU_ann2(:,14).*15.9994.*0.5;
        
        %new wt% of oxides
        D2=D;
        D2(:,5)=Fe3_ratio(:,1).*(D(:,5).*(2*FeO/Fe2O3)+D(:,6))./(2*FeO/Fe2O3);
        D2(:,6)=(1-Fe3_ratio(:,1)).*(D(:,5).*(2*FeO/Fe2O3)+D(:,6));
        D2(:,12)=APFU_n(:,12).*H2O.*0.5.*(1./molar_mass)*100;
    end
    
    %% Calculate structural formula
    
    Strct_Frm=StrctFrm(APFU_n2);
    
    %% Data Plotting
    
    if strcmp(wantstrctfrm, 'y')
        
        if strcmp(wantplot, 'y')
            
            Amp_Plot(:,1)=Fe3_ratio2(:,1); %Fe3+/Fetotal
            Amp_Plot(:,2)=(APFU_n2(:,5)+APFU_n2(:,6)); %Fetotal
            Amp_Plot(:,3)=Strct_Frm(:,15)./(Strct_Frm(:,15)+Strct_Frm(:,16)); %Ca/(Ca+Na) in B
            Amp_Plot(:,4)=Strct_Frm(:,19)+Strct_Frm(:,20)+2.*Strct_Frm(:,18); %A(Na + K + 2Ca), Ca=0
            Amp_Plot(:,5)=Strct_Frm(:,4)+Strct_Frm(:,7)+2.*Strct_Frm(:,5); %C(Al + Fe3+ +2Ti)
            Amp_Plot(:,6)=Strct_Frm(:,8)./(Strct_Frm(:,8)+Strct_Frm(:,9)+Strct_Frm(:,13)); %XMg
            Amp_Plot(:,7)=Strct_Frm(:,9)./(Strct_Frm(:,9)+Strct_Frm(:,8)+Strct_Frm(:,10)); %Fe2+/(Fe2+ + Mg + Mn) in C
            Amp_Plot(:,8)=Strct_Frm(:,7)./(Strct_Frm(:,7)+Strct_Frm(:,4)+Strct_Frm(:,5)); %Fe3+/(Fe3+ + Al + Ti) in C
            Amp_Plot(:,9)=Strct_Frm(:,1); %Si (T)
            
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
            
            
            prompt5='Do you want to plot Ca, Na-Ca, and Na amphibole classification diagrams? (y|n):';
            want_classplot= input(prompt5, 's');
            
            if strcmp(want_classplot, 'y')
                
                %Plots for Ca amphiboles
                figure('Name','Calcium Amphiboles');
                xlim([0 2])
                hold on
                ylim([0 1])
                hold on
                plot([0 2],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                plot([1.5 1.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                text(0.25,0.25,'Tremolite','FontSize',12,'HorizontalAlignment','center')
                text(0.25,0.75,'Edenite','FontSize',12,'HorizontalAlignment','center')
                text(1,0.25,'Magnesio-Hornblende','FontSize',12,'HorizontalAlignment','center')
                text(1,0.75,'Pargasite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.75,'Sadanagaite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.25,'Tschermakite','FontSize',12,'HorizontalAlignment','center')
                xlabel('^{C}(Al + Fe^{3+} + 2Ti) APFU')
                ylabel('^{A}(Na + K + 2Ca) APFU')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) >= 0.75
                        scatter(Amp_Plot(c,5),Amp_Plot(c,4),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                figure('Name','Calcium Amphiboles 2');
                xlim([5.5 8])
                hold on
                ylim([0 1])
                hold on
                plot([5.5 8],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([6.5 6.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                plot([7.5 7.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                plot([7.5 8],[0.9 0.9],'k','linewidth',0.5)
                hold on
                text(7.0,0.25,'Ferro-hornblende','FontSize',12,'HorizontalAlignment','center')
                text(6.0,0.25,'Ferro-tschermakite','FontSize',12,'HorizontalAlignment','center')
                text(7.0,0.75,'Magnesio-Hornblende','FontSize',12,'HorizontalAlignment','center')
                text(6.0,0.75,'Tschermakite','FontSize',12,'HorizontalAlignment','center')
                text(7.75,0.25,'Ferro-actinolite','FontSize',12,'HorizontalAlignment','center')
                text(7.75,0.75,'Actinolite','FontSize',12,'HorizontalAlignment','center')
                text(7.75,0.95,'Tremolite','FontSize',12,'HorizontalAlignment','center')
                xlabel('Si APFU')
                ylabel('X_{Mg}')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) >= 0.75
                        scatter(Amp_Plot(c,9),Amp_Plot(c,6),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                
                %Na-Ca amphiboles
                figure('Name','Sodium-Calcium Amphiboles');
                xlim([0 2])
                hold on
                ylim([0 1])
                hold on
                plot([0 2],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0],[0 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0.5],[0.5 1.0],'k','linewidth',0.5)
                hold on
                plot([1.5 1.5],[0.0 1.0],'k','linewidth',0.5)
                hold on
                text(0.25,0.75,'Richterite','FontSize',12,'HorizontalAlignment','center')
                text(1,0.25,'Winchite','FontSize',12,'HorizontalAlignment','center')
                text(1,0.75,'Katophorite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.75,'Taramite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.25,'Barroisite','FontSize',12,'HorizontalAlignment','center')
                xlabel('^{C}(Al + Fe^{3+} + 2Ti) APFU')
                ylabel('^{A}(Na + K + 2Ca) APFU')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) < 0.75 && Amp_Plot(c,3) > 0.25
                        scatter(Amp_Plot(c,5),Amp_Plot(c,4),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                %Na amphiboles
                
                figure('Name','Sodium Amphiboles');
                xlim([0 2])
                hold on
                ylim([0 1])
                hold on
                plot([1 2],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([1.5 1.5],[0.5 1.0],'k','linewidth',0.5)
                hold on
                plot([0.5 1.5],[1 0],'k','linewidth',0.5)
                hold on
                text(1.15,0.75,'Eckermannite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.75,'Nyboite','FontSize',12,'HorizontalAlignment','center')
                text(1.75,0.25,'Glaucophane','FontSize',12,'HorizontalAlignment','center')
                xlabel('^{C}(Al + Fe^{3+} + 2Ti) APFU')
                ylabel('^{A}(Na + K + 2Ca) APFU')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) <= 0.25
                        scatter(Amp_Plot(c,5),Amp_Plot(c,4),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
                
                figure('Name','Sodium Amphiboles 2');
                xlim([0 1])
                hold on
                ylim([0 1])
                hold on
                plot([0 1],[0.5 0.5],'k','linewidth',0.5)
                hold on
                plot([0.5 0.5],[0 1],'k','linewidth',0.5)
                hold on
                text(0.25,0.25,'Glaucophane','FontSize',12,'HorizontalAlignment','center')
                text(0.25,0.75,'Ferro-Glaucophane','FontSize',12,'HorizontalAlignment','center')
                text(0.75,0.25,'Magnesio-Riebeckite','FontSize',12,'HorizontalAlignment','center')
                text(0.75,0.75,'Riebeckite','FontSize',12,'HorizontalAlignment','center')
                xlabel('Fe^{3+}/(Fe^{3+} + Al + Ti)')
                ylabel('Fe^{2+}/(Fe^{2+} + Mg + Mn)')
                box on
                hold on
                
                for c=1:m
                    if Amp_Plot(c,3) <= 0.25
                        scatter(Amp_Plot(c,8),Amp_Plot(c,7),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                        hold on
                    end
                end
            end
            
            prompt6='Do you want to plot Fe3+ vs FeTotal? (y|n):';
            want_FePlot= input(prompt6, 's');
            
            if strcmp(want_FePlot, 'y')
                
                figure('Name','Fe3+/Fetotal vs Fetotal');
                xlim([0 5])
                hold on
                ylim([0 1])
                hold on
                xlabel('\SigmaFe (APFU)')
                ylabel('Fe^{3+}/\SigmaFe')
                box on
                hold on
                scatter(Amp_Plot(:,2),Amp_Plot(:,1),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
            end
        end
    end
    
    %% Calculates Outputs
    
    %automatic Fe3+ estimation
    
    %W cations
    W_site(:,1)=APFU_n2(:,12); %OH
    W_site(:,2)=APFU_n2(:,13); %F
    W_site(:,3)=APFU_n2(:,14); %Cl
    W_site(:,4)=APFU_n2(:,12)+APFU_n2(:,13)+APFU_n2(:,14); %W sum
    
    all=[Strct_Frm W_site];
    
    Strct_Frm=array2table(all,'VariableNames',{'Si_T','Al_T','Sum_T','Al_C','Ti_C','Cr_C','Fe3_C','Mg_C','Fe2_C','Mn_C','Sum_C','Mg_B','Fe2_B','Mn_B','Ca_B','Na_B','Sum_B','Ca_A','Na_A','K_A','Cation_Sum','OH_W','F_W','Cl_W','W_sum'});
    APFU=array2table(APFU_n2,'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K','OH','F','Cl'});
    
    limits=[low_check hi_check];
    
    Fe3_limits=array2table(limits,'VariableNames',{'Low_Fe3_limits','High_Fe3_limits'});
    
    amph_class=[APFU_Fe(:,12:18) O2_Nfact];
    
    Fe3_class=array2table(amph_class, 'VariableNames',{'Crit1_1','Crit1_2','Crit1_3','Crit2_1','Crit2_2','Crit2_3','Crit2_4','Crit2_5'});
    
    %Fe3+ low estimates:
    
    %Fe2+ Only formula
    APFU_FeT=array2table(APFU_Fe(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to 8Si
    APFU_SiT=array2table(APFU_Si2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe3+ + Fe2+ + Mn + Mg + Ca + Na + K] = 16
    APFU_AfullT=array2table(APFU_Afull2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe3+ + Fe2+ + Mn + Mg + Ca] = 15
    APFU_NaAT=array2table(APFU_NaA2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Fe3+ high estimates:
    %Fe3+ only
    
    APFU_Fe2O3T=array2table(APFU_Fe2O32(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al] = 8 cations
    APFU_SiAlT=array2table(APFU_SiAlT2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe + Mn + Mg + Ca + Na] = 15 cations
    APFU_KAT=array2table(APFU_KA2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
    
    %Normalized to [Si + Al + Ti + Fe + Mn + Mg + Ca + Na] = 15 cations
    APFU_CaNaBT=array2table(APFU_CaNaB2(:,1:11),'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K'});
else %Fe2+ only calculation

    APFU=APFU_Fe(:,1:11);

    %Adds OH, scales for subsequent itterations
    for c=1:m
        if APFU_I(c,12) > 0
            APFU(c,12) = APFU_I(c,12); %OH from the normalized moles cations
        else
            APFU(c,12) = MC(c,12).*(MC(c,1)./APFU(c,1)); %OH is scaled based on the change in Si content
        end
    end
    
    %Adds F
    for c=1:m
        if APFU_I(c,13) > 0
            APFU(c,13) = APFU_I(c,13); %F from the normalized moles cations
        else
            APFU(c,13) = MC(c,13).*(MC(c,1)./APFU(c,1)); %F is scaled based on the change in Si content
        end
    end
    
    %Adds Cl
    for c=1:m
        if APFU_I(c,14) > 0
            APFU(c,14) = APFU_I(c,14); %Cl from the normalized moles cations
        else
            APFU(c,14) = MC(c,14).*(MC(c,1)./APFU(c,1)); %Cl is scaled based on the change in Si content
        end
    end
    
    for c=1:m
        if (APFU(c,12)+APFU(c,13)+APFU(c,14)) > 2 %if (OH_initial + F + Cl) > 2
            if (APFU(c,13)+APFU(c,14)) < 2 % if (F + Cl) < 2
                OH_1(c,1)=2-(APFU(c,13)+APFU(c,14)); %initial OH guess is 2-(F+Cl)
            else
                OH_1(c,1)=0; %W site is filled by Cl and F
            end
        else
            if (APFU(c,12)+APFU(c,13)+APFU(c,14)) < 2 %if (OH_initial + F + Cl) < 2
                if (2-APFU(c,13)-APFU(c,14)) > 0 %if (2-F-Cl) < 2
                    OH_1(c,1)=(2-APFU(c,13)-APFU(c,14)); % OH = 2-F-Cl
                else
                    OH_1(c,1)=0; %W site is filled by Cl and F
                end
            else
                OH_1(c,1)=APFU(c,12);
            end
        end
    end
    
    %OH estimation, Step 2: Adjustment for OH=2-2Ti
    
    %some Ti goes into the T site instead of the M site
    for c=1:m
        %if Ti > (8-(Si+Al)) & (8-(Si+Al))>0
        if APFU(c,2)>(8-(APFU(c,1)+APFU(c,3))) && (8-(APFU(c,1)+APFU(c,3))) > 0
            Ticorr(c,1)=APFU(c,2)-(8-(APFU(c,1)+APFU(c,3))); %Ti in M = Titotal - Ti in T
        else
            Ticorr(c,1)=APFU(c,2); %otherwise all Ti is in M
        end
    end
    
    % Determine appropriate OH correction
    for c=1:m
        if strcmp(wantTiOH, 'y')
            if (APFU(c,12) + APFU(c,13) + APFU(c,14)) > 2 %If (OH + F + Cl) > 2
                if ((2.*APFU(c,12))./((APFU(c,12)+ APFU(c,13) + APFU(c,14))-2.*APFU(c,2)))>0 %if ((2*OH)/((OH+F+Cl)-2*Ti)) >0
                    APFU(c,12)=((2.*APFU(c,12))./((APFU(c,12)+ APFU(c,13) + APFU(c,14))-2.*APFU(c,2)));
                else
                    APFU(c,12)=0; %W is filled by F, Cl, and O
                end
            else
                if (APFU(c,12) + APFU(c,13) + APFU(c,14)) < 2 %If (OH + F + Cl) < 2
                    if ((2-APFU(c,13)-APFU(c,14))-2.*Ticorr(c,1))>0 %if 2-F-Cl-2*Ti > 0
                        APFU(c,12)=((2-APFU(c,13)-APFU(c,14))-2.*Ticorr(c,1)); %OH=(2-F-Cl-2*Ti)
                    else
                        APFU(c,12)=0; %W is filled by F, Cl, and O
                    end
                else
                    APFU(c,12)=APFU(c,12); %OH content is not adjusted for Ti
                end
            end
        else
            APFU(c,12)=OH_1(c,1);
        end
    end
    
    %Moles of anions
    APFU_an(:,1)=APFU(:,1).*2; %SiO2
    APFU_an(:,2)=APFU(:,2).*2; %TiO2
    APFU_an(:,3)=APFU(:,3).*1.5; %Al2O3
    APFU_an(:,4)=APFU(:,4).*1.5; %Cr2O3
    APFU_an(:,5)=APFU(:,5).*1.5; %Fe2O3
    APFU_an(:,6)=APFU(:,6); %FeO
    APFU_an(:,7)=APFU(:,7); %MnO
    APFU_an(:,8)=APFU(:,8); %MgO
    APFU_an(:,9)=APFU(:,9); %CaO
    APFU_an(:,10)=APFU(:,10).*0.5; %Na2O
    APFU_an(:,11)=APFU(:,11).*0.5; %K2O
    APFU_an(:,12)=APFU(:,12).*0.5; %H2O
    
    %for F
    for c=1:m
        if (APFU(c,13)+APFU(c,14))>2 %F + Cl cannot > 2
            APFU_an(c,13)=2.*(APFU(c,13)./(APFU(c,13)+APFU(c,14))); %scales F
        else
            APFU_an(c,13)=APFU(c,13); %Does not scale F
        end
    end
    
    %for Cl
    for c=1:m
        if (APFU(c,13)+APFU(c,14))>2 %F + Cl cannot > 2
            APFU_an(c,14)=2.*(APFU(c,14)./(APFU(c,13)+APFU(c,14))); %scales Cl
        else
            APFU_an(c,14)=APFU(c,14); %Does not scale Cl
        end
    end
    
    %the anionn sum may be > 24
    anion_sum=sum(APFU_an(:,1:12),2)+0.5*APFU_an(:,13)+0.5*APFU_an(:,14);
    anion_norm=24./anion_sum;
    
    APFU_n=APFU.*anion_norm; %cations normalized to 24 anions

    %only APFU_FeT is output 
    APFU_FeT=array2table(APFU,'VariableNames',{'Si','Ti','Al','Cr','Fe3','Fe2','Mn','Mg','Ca','Na','K','OH','F','Cl'});
    
    %initializes the other required variables, but they are zero
    Strct_Frm=0;
    APFU=0;
    Fe3_limits=0;
    Fe3_class=0;
    APFU_SiT=0;
    APFU_AfullT=0;
    APFU_NaAT=0;
    APFU_Fe2O3T=0;
    APFU_SiAlT=0;
    APFU_KAT=0;
    APFU_CaNaBT=0;
end

end
