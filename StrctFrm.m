%amphibole structural formula 

function Frm=StrctFrm(APFU)

[m,n]=size(APFU); %finds the x and y size of the input data matrix

%T SITE

%Si 
for c=1:m
    if APFU(c,1)<8.000
        Frm(c,1)=APFU(c,1);
    else
        Frm(c,1)=8;
    end
end

%Al(T)
for c=1:m
    if (8-Frm(c,1))>APFU(c,1)
        Frm(c,2)=APFU(c,3);
    else
        Frm(c,2)=8-Frm(c,1);
    end
end


%Sum of T site
Frm(:,3)=Frm(:,1)+Frm(:,2);

%C site 

%Al(C)
Frm(:,4)=APFU(:,3)-Frm(:,2); %Altotal-Al(T)=Al(C)
Frm(:,5)=APFU(:,2); %Ti (C)
Frm(:,6)=APFU(:,4); %Cr (C)
Frm(:,7)=APFU(:,5); %Fe3+ (C)

%Mg (C)
for c=1:m
    if 5-sum(Frm(c,4:1:7),2)>0 %Is 5-(Al(C) + Ti + Cr + Fe3) > 0? If y, then some Mg goes into C
        if 5-sum(Frm(c,4:1:7),2)>APFU(c,8) %is 5-(Al(C) + Ti + Cr + Fe3) > Mg?
            Frm(c,8)=APFU(c,8); %if y, all Mg goes into C
        else
            Frm(c,8)=5-sum(Frm(c,4:1:7),2); %if n, only some Mg goes into C, the rest goes into B
        end
    else
        Frm(c,8)=0; %if C is already filled, all Mg goes into the B site
    end
end

%Fe2+ (C)
for c=1:m
    if 5-sum(Frm(c,4:1:8),2)>0 %Is 5-(Al(C) + Ti + Cr + Fe3 + Mg) > 0? If y, then some Fe2+ goes into C
        if 5-sum(Frm(c,4:1:8),2)>APFU(c,6) %is 5-(Al(C) + Ti + Cr + Fe3 + Mg) > Fe2+?
            Frm(c,9)=APFU(c,6); %if y, all Fe2+ goes into C
        else
            Frm(c,9)=5-sum(Frm(c,4:1:8),2); %if n, only some Fe2+ goes into C, the rest goes into B
        end
    else
        Frm(c,9)=0; %if C is already filled, all Fe2+ goes into the B site
    end
end


%Mn (C)
for c=1:m
    if 5-sum(Frm(c,4:1:9),2)>0 %Is 5-(Al(C) + Ti + Cr + Fe3 + Mg + Fe2) > 0? If y, then some Mn goes into C
        if 5-sum(Frm(c,4:1:9),2)>APFU(c,7) %is 5-(Al(C) + Ti + Cr + Fe3 + Mg + Fe) > Mn?
            Frm(c,10)=APFU(c,7); %if y, All Mn goes into C
        else
            Frm(c,10)=5-sum(Frm(c,4:1:9),2); %if n, only some Mn goes into C, the rest goes into B
        end
    else
        Frm(c,10)=0; %if C is already filled, all Mn goes into the B site
    end
end

%Sum of C site
for c=1:m
    Frm(c,11)=sum(Frm(c,4:1:10));
end

%B Site 

Frm(:,12)=APFU(:,8)-Frm(:,8); %Mg (B)
Frm(:,13)=APFU(:,6)-Frm(:,9); %Fe2+ (B)
Frm(:,14)=APFU(:,7)-Frm(:,10); %Mn(B)

%Ca (B)
for c=1:m
    if 2-sum(Frm(c,12:1:14),2)>0 %Is 2-(Mg(B) + Fe2(B) + Mn(B)) > 0? If y, then some Ca goes into A
        if 2-sum(Frm(c,12:1:14),2)>APFU(c,9) %Is 2-(Mg(B) + Fe2(B) + Mn(B)) > Ca?
            Frm(c,15)=APFU(c,9); %if y, All Ca goes into B
        else
            Frm(c,15)=2-sum(Frm(c,12:1:14),2); %if n, only some Ca goes into B, the rest goes into A
        end
    else
        Frm(c,15)=0; %if B is already filled, all Ca goes into the A site
    end
end

%Na (B)
for c=1:m
    if 2-sum(Frm(c,12:1:15),2)>0 %Is 2-(Mg(B) + Fe2(B) + Mn(B) + Ca(B)) > 0? If y, then some Na goes into B
        if 2-sum(Frm(c,12:1:15),2)>APFU(c,10) %Is 2-(Mg(B) + Fe2(B) + Mn(B) + Ca(B)) > Na?
            Frm(c,16)=APFU(c,10); %if y, All Na goes into B
        else
            Frm(c,16)=2-sum(Frm(c,12:1:15),2); %if n, only some Na goes into B, the rest goes into A
        end
    else
        Frm(c,16)=0; %if B is already filled, all Na goes into the A site
    end
end

%B site sum
Frm(:,17)=Frm(:,12)+Frm(:,13)+Frm(:,14)+Frm(:,15)+Frm(:,16);

%A site 
Frm(:,18)=APFU(:,9)-Frm(:,15); %Ca (A)
Frm(:,19)=APFU(:,10)-Frm(:,16); %Na (A)
Frm(:,20)=APFU(:,11); %K (B)

%Cation total
Frm(:,21)=Frm(:,20)+Frm(:,19)+Frm(:,18)+Frm(:,17)+Frm(:,11)+Frm(:,3);


end