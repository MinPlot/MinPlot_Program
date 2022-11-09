%%olivine Structural Formula with Fe3+ estimationn

function [StrctFrm, APFU]=olivine_fe3(data,headers,wantstrctfrm)

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

%Makes Al2O3 optional
if strcmp(headers,'Al2O3')==zeros(1,length(headers))
    I(1,3)=0;
else
    I(1,3)=find(strcmp(headers,'Al2O3'));
end

%Makes Cr2O3 optional
if strcmp(headers,'Cr2O3')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'Cr2O3'));
end

%Makes NiO optional
if strcmp(headers,'NiO')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'NiO'));
end

I(1,6)=find(strcmp(headers,'FeO'));
I(1,7)=find(strcmp(headers,'MnO'));
I(1,8)=find(strcmp(headers,'MgO'));

%Makes CaO optional
if strcmp(headers,'CaO')==zeros(1,length(headers))
    I(1,9)=0;
else
    I(1,9)=find(strcmp(headers,'CaO'));
end

cat=3.0; %cations per formula unit
Opfu=4.0; %oxygens per formula unit

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

W=[SiO2,TiO2,Al2O3,Cr2O3,NiO,FeO,MnO,MgO,CaO];

%% Calculate cations units

MC(:,1)=data(:,I(1,1))./W(:,1); %for SiO2

%calculates for TiO2 if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1); 
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for TiO2
end 

%calculates for Al2O3 if it is included in the analysis 
if I(1,3)==0
    MC(:,3)=zeros(m,1); 
else
    MC(:,3)=(data(:,I(1,3))./W(:,3)).*2; %for Al2O3
end 

%calculates for Cr2O3 if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=(data(:,I(1,4))./W(:,4)).*2; %for Cr2O3
end

%calculates for NiO if it is included in the analysis 
if I(1,5)==0
    MC(:,5)=zeros(m,1);
else
    MC(:,5)=data(:,I(1,5))./W(:,5); %for NiO
end

MC(:,6)=data(:,I(1,6))./W(:,6); %for FeO
MC(:,7)=data(:,I(1,7))./W(:,7); %for MnO
MC(:,8)=data(:,I(1,8))./W(:,8); %for MgO

%calculates for CaO if it is included in the analysis 
if I(1,9)==0
    MC(:,9)=zeros(m,1);
else
    MC(:,9)=data(:,I(1,9))./W(:,9); %for CaO
end

MCnormfact=cat./sum(MC,2); %normalization factor

%% Calculate normalized cations units

MCnorm=MCnormfact.*MC; %creates a matrix of normalized cations

%% Calculate Oxygen Units

O2(:,1)=MCnorm(:,1).*2; %for SiO2
O2(:,2)=MCnorm(:,2).*2; %for TiO2
O2(:,3)=MCnorm(:,3).*(3/2); %for Al2O3
O2(:,4)=MCnorm(:,4).*(3/2); %for Cr2O3
O2(:,5)=MCnorm(:,5); %for NiO
O2(:,6)=MCnorm(:,6); %for FeO
O2(:,7)=MCnorm(:,7); %for MnO
O2(:,8)=MCnorm(:,8); %for MgO
O2(:,9)=MCnorm(:,9); %for CaO

O2total=sum(O2,2); %O2 totals

%% Atoms pfu

APFU(:,1)=MCnorm(:,1); %for Si
APFU(:,2)=MCnorm(:,2); %for Ti
APFU(:,3)=MCnorm(:,3); %for Al
APFU(:,5)=MCnorm(:,4); %for Cr
APFU(:,6)=MCnorm(:,5); %for Ni
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
            APFU(c,4)=2*(Opfu-O2total(c,1)); 
        else
            APFU(c,4)=MCnorm(c,6);
        end
    else
        APFU(c,4)=0;
    end
end

APFU(:,7)=MCnorm(:,6)-APFU(:,4); %the APFU of Fe2+ equals totalFe-Fe3+

APFU(:,11)=sum(APFU,2); %calculations the total, which should be 4

% Oxygen deficiency 
APFU(:,12)=Opfu-O2total; %must be greater than zero

%% structural formula calculation

%T SITE
%Si 
for c=1:m
    if APFU(c,1)<1.000 
        StrctFrm(c,1)=APFU(c,1); %If Si < 1, then Si(T) = the measured Si content
    else
        StrctFrm(c,1)=1; %If Si is in excess, then Si(T) = 1
    end
end

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

%Fe3+(T)
for c=1:m
    if 1-StrctFrm(c,1)-StrctFrm(c,2)>0 %Is 1-(Si+Al) > 0? If y, then some Fe3+ goes into T
        if 1-StrctFrm(c,1)-StrctFrm(c,2)>APFU(c,4) %For low Fe3+ Ol, 1-(Si+Al) may be > Fe3+
            StrctFrm(c,3)=APFU(c,4); %All Fe3+ goes into T
        else
            StrctFrm(c,3)=1-StrctFrm(c,1)-StrctFrm(c,2); %if there isn't enough space in T for all Fe3+, the rest will go to M
        end
    else
        StrctFrm(c,3)=0; %if Si+Al=1, then no Fe3+ goes into T
    end
end

%Sum of T site
StrctFrm(:,4)=StrctFrm(:,1)+StrctFrm(:,2)+StrctFrm(:,3);

%M SITE
%Al
StrctFrm(:,5)=APFU(:,3)-StrctFrm(:,2); %Al(M1) = Total Al - Al(T)

%Ti
StrctFrm(:,6)=APFU(:,2);

%Fe3+
StrctFrm(:,7)=APFU(:,4)-StrctFrm(:,3); %Fe3+(M1) = Total Fe3+ - Fe3+(T)

%Cr
StrctFrm(:,8)=APFU(:,5);

%Ni
StrctFrm(:,9)=APFU(:,6);

%Fe
StrctFrm(:,10)=APFU(:,7);

%Mn
StrctFrm(:,11)=APFU(:,8);

%Mg
StrctFrm(:,12)=APFU(:,9);

%Ca
StrctFrm(:,13)=APFU(:,10);

%M total
StrctFrm(:,14)=sum(StrctFrm(:,5:13),2);

%% Endmember calculation

if strcmp(wantstrctfrm, 'y')
    Endmembers(:,1)=APFU(:,9)./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XFo
    Endmembers(:,2)=(APFU(:,4)+APFU(:,7))./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XFa
    Endmembers(:,3)=(APFU(:,8))./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XTe
    Endmembers(:,4)=(APFU(:,10))./(APFU(:,4)+APFU(:,7)+APFU(:,8)+APFU(:,9)+APFU(:,10)); % XCa-Ol
    
    %prompts the user if they wish to plot the CPX data
    prompt1='Do you wish to plot olivine compositions? (y|n): ';
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
            
            prompt5='Do you wish to make a binary plot? (y|n): ';
            wantbinary=input(prompt5, 's');
            
            if strcmp(wantbinary, 'y')
                prompt6='Lower Fo limit (0-100):';
                FoLow=input(prompt6);
                
                prompt7='Upper Fo limit (0-100):';
                FoHi=input(prompt7);
                
                figure('Name','Olivine Binary')
                y(1:length(Endmembers(:,1)),1)=0;
                scatter(Endmembers(:,1).*100,y,symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                ylim([-0.5 0.5]);
                stp=(FoHi-FoLow).*0.05; %adjusts the ends of the X limits so that the 
                %numerical values are all shown (e.g., for X= 0 to 100, 0 and 100 would be cut off otherwise)
                xlim([FoLow-stp FoHi+stp]); %changes X relative to the Fo limits chosen
                
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
                Xlb=((FoHi-FoLow)/2)+FoLow; 
                xlabel('fo (mol %)','Position',[Xlb -0.6], 'VerticalAlignment','Top','HorizontalAlignment','center')
                
            end
            
            prompt6='Do you wish to a ternary plot? (y|n): ';
            wantternary=input(prompt6, 's');
            
            if strcmp(wantternary, 'y')
                figure('Name','Olivine Ternary')
                pgon=polyshape([0 0.5 1],[0 sqrt(3)/2 0]);
                plot(pgon,'FaceColor','w')
                hold on
                
                %plot grid (spacing of 10 %)
                plot([0.45 0.55], [0.779422863 0.779422863],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.9
                hold on
                plot([0.4 0.6], [0.692820323 0.692820323],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.8
                hold on
                plot([0.35 0.65], [0.606217783 0.606217783],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.7
                hold on
                plot([0.3 0.7], [0.519615242 0.519615242],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.6
                hold on
                plot([0.25 0.75], [0.433012702 0.433012702],'k','linewidth',1.5) %XWo=0.5
                hold on
                plot([0.2 0.8], [0.346410162 0.346410162],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.4
                hold on
                plot([0.15 0.85], [0.259807621 0.259807621],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.3
                hold on
                plot([0.1 0.9], [0.173205081 0.173205081],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XWo=0.2
                hold on
                plot([0.05 0.95], [0.08660254 0.08660254],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XWo=0.1
                hold on
                plot([0.05 0.1], [0.08660254 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.9
                hold on
                plot([0.1 0.2], [0.173205081 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.8
                hold on
                plot([0.15 0.3], [0.259807621 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.7
                hold on
                plot([0.2 0.4], [0.346410162 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.6
                hold on
                plot([0.25 0.5], [0.433012702 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.5
                hold on
                plot([0.3 0.6], [0.519615242 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.4
                hold on
                plot([0.35 0.7], [0.606217783 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.3
                hold on
                plot([0.4 0.8], [0.692820323 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XEn=0.2
                hold on
                plot([0.45 0.9], [0.779422863 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XEn=0.1
                hold on
                plot([0.95 0.9], [0.08660254 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.9
                hold on
                plot([0.9 0.8], [0.173205081 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.8
                hold on
                plot([0.85 0.7], [0.259807621 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.7
                hold on
                plot([0.8 0.6], [0.346410162 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.6
                hold on
                plot([0.75 0.5], [0.433012702 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.5
                hold on
                plot([0.7 0.4], [0.519615242 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.4
                hold on
                plot([0.65 0.3], [0.606217783 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.3
                hold on
                plot([0.6 0.2], [0.692820323 0],'color',[0.3 0.3 0.3],'LineStyle','--','linewidth',0.5) %XFs=0.2
                hold on
                plot([0.55 0.1], [0.779422863 0],'color',[0.3 0.3 0.3],'LineStyle',':','linewidth',0.5) %XFs=0.1
                hold on
                
                %plot boundaries of the triangle
                plot([0 1],[0 0],'k','linewidth',1.5)
                hold on
                plot([0 0.5],[0 sqrt(3)/2],'k','linewidth',1.5)
                hold on
                plot([0.5 1],[sqrt(3)/2 0],'k','linewidth',1.5)
                hold on
                
                %labels
                text(-0.10,0.0,'fo','FontSize',14)
                text(0.51,0.92,'Ca-ol','FontSize',14)
                text(1.02,-0.04,'fa+tep','FontSize',14)
                text(0.16,0.44,'mtc','FontSize',14)
                text(0.77,0.44,'kir','FontSize',14)
                
                %tick labels
                text(-0.02,-0.039,'0.0','FontSize',12,'HorizontalAlignment','center','Rotation',60)
                text(0.18,-0.039,'0.2','FontSize',12,'HorizontalAlignment','center','Rotation',60)
                text(0.38,-0.039,'0.4','FontSize',12,'HorizontalAlignment','center','Rotation',60)
                text(0.58,-0.039,'0.6','FontSize',12,'HorizontalAlignment','center','Rotation',60)
                text(0.78,-0.039,'0.8','FontSize',12,'HorizontalAlignment','center','Rotation',60)
                text(0.98,-0.039,'1.0','FontSize',12,'HorizontalAlignment','center','Rotation',60)
                text(1.02,0.01,'0.0','FontSize',12)
                text(0.92,0.18,'0.2','FontSize',12)
                text(0.82,0.35,'0.4','FontSize',12)
                text(0.72,0.52,'0.6','FontSize',12)
                text(0.62,0.70,'0.8','FontSize',12)
                text(0.52,0.87,'1.0','FontSize',12)
                text(-0.025,0.040,'1.0','FontSize',12,'HorizontalAlignment','center','Rotation',300)
                text(0.075,0.215,'0.8','FontSize',12,'HorizontalAlignment','center','Rotation',300)
                text(0.178,0.386,'0.6','FontSize',12,'HorizontalAlignment','center','Rotation',300)
                text(0.278,0.558,'0.4','FontSize',12,'HorizontalAlignment','center','Rotation',300)
                text(0.382,0.735,'0.2','FontSize',12,'HorizontalAlignment','center','Rotation',300)
                text(0.482,0.901,'0.0','FontSize',12,'HorizontalAlignment','center','Rotation',300)

                X1=0.5.*(Endmembers(:,4))+(Endmembers(:,2)+Endmembers(:,3));
                Y1=(Endmembers(:,4))*(cos(30*pi()/180));
                
                scatter(X1(:),Y1(:),symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                hold off
                
                axis image
                axis off
            end 
            
            prompt7='Do you wish to plot Fo (mol%) vs NiO (wt. %)? (y|n): ';
            wantNi=input(prompt7, 's');
            
            if strcmp(wantNi, 'y')
                
                %set axes limits
                prompt8='Lower Fo limit (0-100):';
                FoLow=input(prompt8);
                
                prompt9='Upper Fo limit (0-100):';
                FoHi=input(prompt9);
                
                prompt10='Lower NiO limit (wt %):';
                NiLow=input(prompt10);
                
                prompt11='Upper NiO limit (wt %):';
                NiHi=input(prompt11);

                figure('Name','fo (mol%) vs NiO (wt. %)')
                Ni_wt=data(:,I(1,5)); %vector with NiO contents in wt %
                scatter(Endmembers(:,1).*100,Ni_wt,symbsize,symb,'filled','MarkerFaceAlpha',3/8,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',fil)
                ylim([NiLow NiHi]);
                xlim([FoLow FoHi]);
                xlabel('fo (mol %)')
                ylabel('NiO (wt %)')

                ax=gca;
                ax.Box = 'on';
            end

    end
            
    all=[StrctFrm Endmembers APFU(:,12)];
    StrctFrm=array2table(all,'VariableNames',{'Si_T','Al_T','Fe3_T','Sum_T','Al_M','Ti_M','Fe3_M','Cr_M','Ni_M','Fe2_M','Mn_M','Mg_M','Ca_M','Sum_M','XFo','XFa','XTep','XLrn','O2_deficiency'});  
    
else
    %outputs APFU table only
    APFU=array2table(APFU,'VariableNames',{'Si','Ti','Al','Fe3','Cr','Ni','Fe2','Mn','Mg','Ca','Sum','O2_deficiency'});
end

end
