%%EPMA data processing script 

clear all
close all

disp('The file must be a .txt file with oxide headings (see template)');
[FileName, PathName, FilterIndex] = uigetfile({'*.txt';'*.csv';'*.xls';'*.xlsx'},'File Selector');
T=readtable(strcat(PathName,FileName));

data=table2array(T); %data converted to array
data(isnan(data))=0;
headers=T.Properties.VariableNames; %saves a list of the oxides in the header

%call functions for data processing 
prompt2='What mineral formula do you want to recalculate?:';
disp('Options are (CASE SENSITIVE): garnet, pyroxene, olivine, amphibole,')
disp('feldspar, mica, staurolite, cordierite, chlorite, chloritoid,')
disp('talc, epidote, titanite, oxyspinel, and sulfide.')
wantformula = input(prompt2, 's');

%garnet 
if strcmp(wantformula, 'garnet')
    prompt3='Calculate Fe3+ from stoichiometry? (y|n):'; %asks if you want the ferric iron calculation
    wantferric=input(prompt3, 's');
    %if you do not want ferric iron, the following procedure occurs 
    if strcmp(wantferric, 'y') 
        prompt4='Do you want cation assignment and endmember calculations? (y|n):';
        wantstrctfrm=input(prompt4, 's');
        if strcmp(wantstrctfrm, 'y')
            StrctFrm=garnet_fe3(data,headers,wantstrctfrm); %calls garnet (w/Fe3+) calculation function, outputs the structural formula with cation assignment
            disp('Mineral formula and endmembers was output.')
        else 
            [~,APFU]=garnet_fe3(data,headers,wantstrctfrm); %calls the garnet (w/Fe3+) calculation, but only outputs APFU and no endmembers
            disp('The calculated APFU was output without cation assignment.')
        end
    else
        %If an FeO only calculation was selected, the following procedure
        %occurs 
        prompt4='Do you want cation assignment and endmember calculations? (y|n):';
        wantstrctfrm=input(prompt4, 's');
        if strcmp(wantstrctfrm, 'y') 
            StrctFrm=garnet_fe2(data,headers,wantstrctfrm); %calls garnet (FeO only) calculation function, outputs the structural formula with cation assignment
            disp('Mineral formula and endmembers was output.')
        else 
            [~,APFU]=garnet_fe2(data,headers,wantstrctfrm); %calls garnet (FeO only) calculation function, outputs the APFU only
            disp('The calculated APFU was output without cation assignment.')
        end           
    end 
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%pyroxene
if strcmp(wantformula, 'pyroxene')
    prompt3='Calculate Fe3+ from stoichiometry? (y|n):'; %asks if you want the ferric iron calculation
    wantferric=input(prompt3, 's');
    %if you do not want ferric iron, the following procedure occurs 
    if strcmp(wantferric, 'y') 
        prompt4='Do you want cation assignment and endmember calculations? (y|n):';
        wantstrctfrm=input(prompt4, 's');
        if strcmp(wantstrctfrm, 'y')

            prompt5='Use the calculation scheme for high-pressure/temperature pyroxene? (y|n):';
            disp('For K-, Ca-Eskola-, & Ca-Tschermaks-rich pyroxenes:')

            wantHP=input(prompt5, 's');

            if strcmp(wantHP, 'y')
                StrctFrm=pyroxene_fe3_HP(data,headers,wantstrctfrm); %calls pyroxene (w/Fe3+) calculation function, outputs the structural formula with cation assignment
                disp('Mineral formula and endmembers was output.')
            else
                StrctFrm=pyroxene_fe3(data,headers,wantstrctfrm); %calls pyroxene (w/Fe3+) calculation function, outputs the structural formula with cation assignment
                disp('Mineral formula and endmembers was output.')
            end
        else 
            [~,APFU]=pyroxene_fe3(data,headers,wantstrctfrm); %calls the pyroxene (w/Fe3+) calculation, but only outputs APFU and no endmembers
            disp('The calculated APFU was output without cation assignment.')
        end
    else
        %If an FeO only calculation was selected, the following procedure
        %occurs 
        prompt4='Do you want cation assignment and endmember calculations? (y|n):';
        wantstrctfrm=input(prompt4, 's');
        if strcmp(wantstrctfrm, 'y') 
            StrctFrm=pyroxene_fe2(data,headers,wantstrctfrm); %calls pyroxene (FeO only) calculation function, outputs the structural formula with cation assignment
            disp('Mineral formula and endmembers was output.')
        else 
            [~,APFU]=pyroxene_fe2(data,headers,wantstrctfrm); %calls pyroxene (FeO only) calculation function, outputs the APFU only
            disp('The calculated APFU was output without cation assignment.')
        end           
    end 
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end  
end

%olivine
if strcmp(wantformula, 'olivine')
    prompt3='Calculate Fe3+ from stoichiometry? (y|n):'; %asks if you want the ferric iron calculation
    wantferric=input(prompt3, 's');
    %if you do not want ferric iron, the following procedure occurs 
    if strcmp(wantferric, 'y') 
        prompt4='Do you want cation assignment and endmember calculations? (y|n):';
        wantstrctfrm=input(prompt4, 's');
        if strcmp(wantstrctfrm, 'y')
            StrctFrm=olivine_fe3(data,headers,wantstrctfrm); %calls olivine (w/Fe3+) calculation function, outputs the structural formula with cation assignment
            disp('Mineral formula and endmembers was output.')
        else 
            [~,APFU]=olivine_fe3(data,headers,wantstrctfrm); %calls the olivine (w/Fe3+) calculation, but only outputs APFU and no endmembers
            disp('The calculated APFU was output without cation assignment.')
        end
    else
        %If an FeO only calculation was selected, the following procedure
        %occurs 
        prompt4='Do you want cation assignment and endmember calculations? (y|n):';
        wantstrctfrm=input(prompt4, 's');
        if strcmp(wantstrctfrm, 'y') 
            StrctFrm=olivine_fe2(data,headers,wantstrctfrm); %calls olivine (FeO only) calculation function, outputs the structural formula with cation assignment
            disp('Mineral formula and endmembers was output.')
        else 
            [~,APFU]=olivine_fe2(data,headers,wantstrctfrm); %calls olivine (FeO only) calculation function, outputs the APFU only
            disp('The calculated APFU was output without cation assignment.')
        end           
    end 
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%feldspar
if strcmp(wantformula, 'feldspar')
    StrctFrm=feldspar(data,headers); %calls feldspar calculation
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
    end
end

%mica
if strcmp(wantformula, 'mica')
    prompt4='Do you want cation assignment and endmember calculations? (y|n):';
    wantstrctfrm=input(prompt4, 's');
    if strcmp(wantstrctfrm, 'y')
        StrctFrm=mica(data,headers,wantstrctfrm); %outputs the structural formula with cation assignment
        disp('Mineral formula and endmembers was output.')
    else
        [~,APFU]=mica(data,headers,wantstrctfrm); %saves only APFU
        disp('The calculated APFU was output without cation assignment.')
    end
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%amphibole
if strcmp(wantformula, 'amphibole')
    
    prompt3='Calculate Fe3+ from stoichiometry? (y|n):'; %asks if you want the ferric iron calculation
    wantferric=input(prompt3, 's');
    
    if strcmp(wantferric, 'y')
        disp('Remember to closely check the 6 criteria for Fe3+ estimation.')
        disp('Not applicable for some compositions, see appendices in')
        disp('Leake et al. (1997) and Hawthorne et al. (2012).')
        
        prompt4='Do you want to correct for OH=2-2Ti? (y|n):';
        wantTiOH= input(prompt4, 's');
       
        prompt5='Do you want cation assignment? (y|n):';
        wantstrctfrm=input(prompt5, 's');
        
        if strcmp(wantstrctfrm, 'y')
            prompt6='Do you want to plot amphibole compositions? (y|n):';
            wantplot= input(prompt6, 's');
            [Strct_Frm, ~, ~, Fe3_limits, Fe3_class, ~, ~, ~, ~, ~, ~, ~]=amphibole(data,headers,wantstrctfrm,wantTiOH,wantplot,wantferric); %outputs the structural formula for the automated Fe3+ estimation only
            disp('Mineral formula was output.')
            
            %prompts you to save the output
            prompt7='Do you wish to save the recalculated data? (y|n):';
            save=input(prompt7,'s');
            if strcmp(save,'y')
                writetable(Strct_Frm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
                writetable(Fe3_limits,strcat(PathName,FileName(1:size(FileName,2)-4),'_Fe3limits.txt'),'Delimiter','\t');
                writetable(Fe3_class,strcat(PathName,FileName(1:size(FileName,2)-4),'_Fe3limitclassification.txt'),'Delimiter','\t');
            end
        else
            wantplot='n';
            [~, APFU, ~, Fe3_limits, Fe3_class, ~, ~, ~, ~, ~, ~, ~]=amphibole(data,headers,wantstrctfrm,wantTiOH,wantplot,wantferric); %outputs the AFPU only
            disp('The calculated APFU was output without cation assignment.')
            
            %prompts you to save the output
            prompt7='Do you wish to save the recalculated data? (y|n):';
            save=input(prompt7,'s');
            if strcmp(save,'y')
                writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
                writetable(Fe3_limits,strcat(PathName,FileName(1:size(FileName,2)-4),'_Fe3limits.txt'),'Delimiter','\t');
                writetable(Fe3_class,strcat(PathName,FileName(1:size(FileName,2)-4),'_Fe3limitclassification.txt'),'Delimiter','\t');
            end
            
        end
    
    else
        wantplot='n';
        wantstrctfrm='n';
        
        prompt4='Do you want to correct for OH=2-2Ti? (y|n):';
        wantTiOH= input(prompt4, 's');

        [~, ~, APFU_FeT, ~, ~, ~, ~, ~, ~, ~, ~, ~]=amphibole(data,headers,wantstrctfrm,wantTiOH,wantplot,wantferric); %outputs the AFPU only
        disp('The calculated FeO only APFU was output without cation assignment.')
        
        %prompts you to save the output
        prompt4='Do you wish to save the recalculated data? (y|n):';
        save=input(prompt4,'s');
        if strcmp(save,'y')
            writetable(APFU_FeT,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
        
    end
end

%epidote
if strcmp(wantformula, 'epidote')
    
    prompt4='Do you want cation assignment and endmember calculations? (y|n):';
    wantstrctfrm=input(prompt4, 's');
    
    if strcmp(wantstrctfrm, 'y')
        [StrctFrm,~]=epidote(data,headers,wantstrctfrm); %calls epidote calculation
    else
        [~,APFU]=epidote(data,headers,wantstrctfrm); %calls epidote calculation
    end
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
end

%chlorite
if strcmp(wantformula, 'chlorite')
    prompt4='Do you want cation assignment and endmember calculations? (y|n):';
    wantstrctfrm=input(prompt4, 's');
    if strcmp(wantstrctfrm, 'y')
        StrctFrm=chlorite(data,headers,wantstrctfrm); %outputs the structural formula with cation assignment
        disp('Mineral formula and endmembers was output.')
    else
        [~,APFU]=chlorite(data,headers,wantstrctfrm); %saves only APFU
        disp('The calculated APFU was output without cation assignment.')
    end
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%titanite
if strcmp(wantformula, 'titanite')
    prompt4='Do you want cation assignment and endmember calculations? (y|n):';
    wantstrctfrm=input(prompt4, 's');
    if strcmp(wantstrctfrm, 'y')
        StrctFrm=titanite(data,headers); %outputs the structural formula with cation assignment
        disp('Mineral formula and endmembers was output.')
    else
        [~,APFU]=titanite(data,headers); %saves only APFU
        disp('The calculated APFU was output without cation assignment.')
    end
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%talc
if strcmp(wantformula, 'talc')
    prompt4='Do you want cation assignment? (y|n):';
    wantstrctfrm=input(prompt4, 's');
    if strcmp(wantstrctfrm, 'y')
        StrctFrm=talc(data,headers,wantstrctfrm); %outputs the structural formula with cation assignment
        disp('Mineral formula and endmembers was output.')
    else
        [~,APFU]=talc(data,headers,wantstrctfrm); %saves only APFU
        disp('The calculated APFU was output without cation assignment.')
    end
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%chloritoid
if strcmp(wantformula, 'chloritoid')
    prompt4='Do you want cation assignment? (y|n):';
    wantstrctfrm=input(prompt4, 's');
    if strcmp(wantstrctfrm, 'y')
        [StrctFrm,~]=chloritoid(data,headers,wantstrctfrm); %outputs the structural formula with cation assignment
        disp('Mineral formula and endmembers was output.')
    else
        [~,APFU]=chloritoid(data,headers,wantstrctfrm); %saves only APFU
        disp('The calculated APFU was output without cation assignment.')
    end
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%staurolite
if strcmp(wantformula, 'staurolite')

    [APFU]=staurolite(data,headers); %saves only APFU
    disp('The calculated APFU was output without cation assignment.')
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')

        writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
    
    end     
end

%spinel
if strcmp(wantformula, 'oxyspinel')
    APFU=spinel_fe3(data,headers); %calls feldspar calculation
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
    end
end

%cordierite
if strcmp(wantformula, 'cordierite')
    prompt4='Do you want cation assignment? (y|n):';
    wantstrctfrm=input(prompt4, 's');
    if strcmp(wantstrctfrm, 'y')
        [StrctFrm,~]=cordierite(data,headers,wantstrctfrm); %outputs the structural formula with cation assignment
        disp('Mineral formula and endmembers was output.')
    else
        [~,APFU]=cordierite(data,headers,wantstrctfrm); %saves only APFU
        disp('The calculated APFU was output without cation assignment.')
    end
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        if strcmp(wantstrctfrm, 'y') %if the structural formula option was chosen, a text file for the structural formula is saved
            writetable(StrctFrm,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
        else %if the APFU only was selected, then the APFU output is chosen
            writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_APFU.txt'),'Delimiter','\t');
        end
    end
        
end

%sulfide
if strcmp(wantformula,'sulfide')
    APFU=sulfide(data,headers); %calls sulfide calculation
    
    %prompts you to save the output
    prompt5='Do you wish to save the recalculated data? (y|n):';
    save=input(prompt5,'s');
    if strcmp(save,'y')
        writetable(APFU,strcat(PathName,FileName(1:size(FileName,2)-4),'_structuralformula.txt'),'Delimiter','\t');
    end
end
