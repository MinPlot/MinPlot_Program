%epidote structural formula 
function [StrctFrm, APFU]=epidote(data,headers,wantstrctfrm)

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

%Makes Fe2O3 optional
if strcmp(headers,'Fe2O3')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'Fe2O3'));
end

%Makes FeO optional
if strcmp(headers,'FeO')==zeros(1,length(headers))
    I(1,6)=0;
else
    I(1,6)=find(strcmp(headers,'FeO'));
end

I(1,7)=find(strcmp(headers,'MnO'));
I(1,8)=find(strcmp(headers,'MgO'));
I(1,9)=find(strcmp(headers,'CaO'));

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

%% organize the input wt %'s
%create columns of zeros if optional data are not included 

%calculates for FeO if it is included in the analysis 
D(:,1)=data(:,I(1,1)); % for SiO2

%For TiO2
if I(1,2)==0
    D(:,2)=zeros(m,1); %makes column of zeros if Ti isn't included
else
    D(:,2)=data(:,I(1,2));
end 

D(:,3)=data(:,I(1,3)); % for Al2O3

%For Cr2O3
if I(1,4)==0
    D(:,4)=zeros(m,1); %makes column of zeros if Cr isn't included
else
    D(:,4)=data(:,I(1,4));
end 

%For Fe2O3
if I(1,5)>0
    D(:,5)=data(:,I(1,5));
else 
     D(:,5)=data(:,I(1,6))*(159.6874./(2*71.8442));
end 

D(:,6)=data(:,I(1,7))*(157.873./(2*70.937)); % converts MnO to Mn2O3
D(:,7)=data(:,I(1,8)); % for MgO
D(:,8)=data(:,I(1,9)); % for CaO

%For Na2O
if I(1,10)==0
    D(:,9)=zeros(m,1); %makes column of zeros if Na2O isn't included
else
    D(:,9)=data(:,I(1,10));
end 

%For K2O
if I(1,11)==0
    D(:,10)=zeros(m,1); %makes column of zeros if K2O isn't included
else
    D(:,10)=data(:,I(1,11));
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
Mn2O3=157.873;
MgO=40.304;
CaO=56.0774;
Na2O=61.979;
K2O=94.195;
BaO=153.329;
F=18.998;
Cl=35.45;

W=[SiO2,TiO2,Al2O3,Cr2O3,Fe2O3,Mn2O3,MgO,CaO,Na2O,K2O];


%% Calculate cations units

MC=D./W; 

%% Calculate Oxygen Units

O2(:,1)=MC(:,1).*2; %for SiO2
O2(:,2)=MC(:,2).*2; %for TiO2
O2(:,3)=MC(:,3).*3; %for Al2O3
O2(:,4)=MC(:,4).*3; %for Cr2O3
O2(:,5)=MC(:,5).*3; %for Fe2O3
O2(:,6)=MC(:,6).*3; %for Mn2O3
O2(:,7)=MC(:,7); %for MgO
O2(:,8)=MC(:,8); %for CaO
O2(:,9)=MC(:,9); %for Na2O
O2(:,10)=MC(:,10); %for K2O

O2_Sum=sum(O2,2); %O2 totals
O2_N=(12.5)./O2_Sum; %normalization factor

%% Normalized Oxygen Units

N_Ox=O2.*O2_N;


%% Atoms pfu

APFU(:,1)=N_Ox(:,1)./2; %Si
APFU(:,2)=N_Ox(:,2)./2; %Ti
APFU(:,3)=N_Ox(:,3).*(2/3); %Al
APFU(:,4)=N_Ox(:,4).*(2/3); %Cr
APFU(:,5)=N_Ox(:,5).*(2/3); %Fe3+
APFU(:,6)=N_Ox(:,6).*(2/3); %Mn3+
APFU(:,7)=N_Ox(:,7); %Mg
APFU(:,8)=N_Ox(:,8); %Ca
APFU(:,9)=N_Ox(:,9).*2; %Na
APFU(:,10)=N_Ox(:,10).*2; %K
APFU(:,11)=sum(APFU,2);

%% Structural Formula
%T site - Si, Al (sums to 3)
%M sites - Ti, Al, Cr3+, Fe3+, Mn3+
%A site - Ca, Na, K, Mg, Fe2+ (Sums to 2)

%T SITE
%Si 
for c=1:m
    if APFU(c,1)<3.000 
        StrctFrm(c,1)=APFU(c,1); %If Si < 3, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=3; %If Si is in excess, then Si(T) = 3
    end
end

%Al(T)
for c=1:m
    if 3-StrctFrm(c,1)>0 %Is 3-Si > 0? If y, then some Al goes into T
        StrctFrm(c,2)=3-StrctFrm(c,1);
    else
        StrctFrm(c,2)=0; %if Si=3, then no Al goes into T
    end
end


StrctFrm(:,3)=StrctFrm(:,2)+StrctFrm(:,1); %T site sum


StrctFrm(:,4)=APFU(:,3)-StrctFrm(:,2); %Al (M)
StrctFrm(:,5)=APFU(:,2);%Ti (M)
StrctFrm(:,6)=APFU(:,4); %Cr (M)
StrctFrm(:,7)=APFU(:,5); %Fe3+ (M)
StrctFrm(:,8)=APFU(:,6); %Mn3+ (M)
StrctFrm(:,9)=sum(StrctFrm(:,4:8),2); %M site sum

StrctFrm(:,10)=APFU(:,7); %Mg (A)
StrctFrm(:,11)=APFU(:,8); %Ca (A)
StrctFrm(:,12)=APFU(:,9); %Na (A)
StrctFrm(:,13)=APFU(:,10); %K2O (A)
StrctFrm(:,14)=sum(StrctFrm(:,11:13),2); %A site sum

%endmembers
Endmembers(:,1)=(StrctFrm(:,4)-2)./(APFU(:,5)+StrctFrm(:,4)+APFU(:,4)+APFU(:,6)-2); %XCzo = (Al-2)/(Fe3 + Al + Cr +Mn -2)
Endmembers(:,2)=(APFU(:,5))./(APFU(:,5)+StrctFrm(:,4)+APFU(:,4)+APFU(:,6)-2); %XEp = (Fe3+)/(Fe3 + Al + Cr + Mn -2)
Endmembers(:,3)=(APFU(:,6))./(APFU(:,5)+StrctFrm(:,4)+APFU(:,4)+APFU(:,6)-2); %XPie = (Mn3+)/(Fe3 + Al + Cr + Mn -2)
Endmembers(:,4)=(APFU(:,4))./(APFU(:,5)+StrctFrm(:,4)+APFU(:,4)+APFU(:,6)-2); %XTaw = (Cr3+)/(Fe3 + Al + Cr + Mn-2)

all=[StrctFrm Endmembers];
StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Sum_T','Al_M','Ti_M','Cr_M','Fe3_M','Mn3_M','M_Sum','Mg_A','Ca_A','Na_A','K_A','A_Sum','XCzo','XEp','XPmt','XTaw'});
APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Cr3','Fe3','Mn3','Mg','Ca','Na','K','Cation_sum'});

if strcmp(wantstrctfrm, 'y')
    %prompts the user if they wish to plot
    prompt1='Do you wish to plot the epidote-clinozoisite binary? (y|n): ';
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
        
        prompt6='Lower XEp limit (0-1):';
        EpLow=input(prompt6);
        
        prompt7='Upper XEp limit (0-1):';
        EpHi=input(prompt7);
        
        figure('Name','Epidote-Clinozoisite Binary')
        y(1:length(Endmembers(:,2)),1)=0;
        scatter(Endmembers(:,2),y,symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
        ylim([-0.5 0.5]);
        stp=(EpHi-EpLow).*0.05; %adjusts the ends of the X limits so that the
        %numerical values are all shown (e.g., for X= 0 to 100, 0 and 100 would be cut off otherwise)
        xlim([EpLow-stp EpHi+stp]); %changes X relative to the Fo limits chosen
        
        %changes the shape of the Y axis to a smaller box and puts
        %the plot in the middle of the figure
        ax=gca;
        sz=ax.OuterPosition;
        ax.OuterPosition=[0 0.5 1 0.15];
        %changes the location of the X axis to the origin
        xloc=ax.XAxisLocation;
        ax.XAxisLocation='origin';
        %remove Y axis
        set(gca,'ytick',[])
        ax.YAxis.Visible = 'off';
        %make the x axis ticks stick out in both directions
        ax.TickDir='both';
        %change the position of the X axis
        Xlb=((EpHi-EpLow)/2)+EpLow;
        xlabel('XEp','Position',[Xlb -0.6], 'VerticalAlignment','Top','HorizontalAlignment','center')
                
    end  

end
    
    

