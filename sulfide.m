%generic sulfide calculator

function [APFU]=sulfide(data,headers)

[m,n]=size(data); %finds the x and y size of the input data matrix

I(1,1)=find(strcmp(headers,'S'));

%makes As optional
if strcmp(headers,'As')==zeros(1,length(headers))
    I(1,2)=0;
else
    I(1,2)=find(strcmp(headers,'As'));
end

%makes Fe optional
if strcmp(headers,'Fe')==zeros(1,length(headers))
    I(1,3)=0;
else
    I(1,3)=find(strcmp(headers,'Fe'));
end

%makes Ni optional
if strcmp(headers,'Ni')==zeros(1,length(headers))
    I(1,4)=0;
else
    I(1,4)=find(strcmp(headers,'Ni'));
end

%makes Cu optional
if strcmp(headers,'Cu')==zeros(1,length(headers))
    I(1,5)=0;
else
    I(1,5)=find(strcmp(headers,'Cu'));
end

%makes Co optional
if strcmp(headers,'Co')==zeros(1,length(headers))
    I(1,6)=0;
else
    I(1,6)=find(strcmp(headers,'Co'));
end

%makes Pb optional
if strcmp(headers,'Pb')==zeros(1,length(headers))
    I(1,7)=0;
else
    I(1,7)=find(strcmp(headers,'Pb'));
end

%makes Zn optional
if strcmp(headers,'Zn')==zeros(1,length(headers))
    I(1,8)=0;
else
    I(1,8)=find(strcmp(headers,'Zn'));
end

%% Molecular weights
S=32.06;
Fe=55.8452;
Ni=58.693;
Cu=63.546;
Co=58.933;
Pb=207.2;
Zn=65.382;
As=74.922; 

W=[S,As,Fe,Ni,Cu,Co,Pb,Zn];

%% Moles of elements

MC(:,1)=data(:,I(1,1))./W(:,1); %for S

%calculates for As if it is included in the analysis 
if I(1,2)==0
    MC(:,2)=zeros(m,1);
else
    MC(:,2)=data(:,I(1,2))./W(:,2); %for As
end 

%calculates for Fe if it is included in the analysis 
if I(1,3)==0
    MC(:,3)=zeros(m,1);
else
    MC(:,3)=data(:,I(1,3))./W(:,3); %for Fe
end 

%calculates for Ni if it is included in the analysis 
if I(1,4)==0
    MC(:,4)=zeros(m,1);
else
    MC(:,4)=data(:,I(1,4))./W(:,4); %for Ni
end 

%calculates for Cu if it is included in the analysis 
if I(1,5)==0
    MC(:,5)=zeros(m,1);
else
    MC(:,5)=data(:,I(1,5))./W(:,5); %for Cu
end 

%calculates for Co if it is included in the analysis 
if I(1,6)==0
    MC(:,6)=zeros(m,1);
else
    MC(:,6)=data(:,I(1,6))./W(:,6); %for Co
end 

%calculates for Pb if it is included in the analysis 
if I(1,7)==0
    MC(:,7)=zeros(m,1);
else
    MC(:,7)=data(:,I(1,7))./W(:,7); %for Pb
end 

%calculates for Zn if it is included in the analysis 
if I(1,8)==0
    MC(:,8)=zeros(m,1);
else
    MC(:,8)=data(:,I(1,8))./W(:,8); %for Zn
end 

%% Normalization

%prompts the user if they wish calculate on a cation or anion basis
prompt1='Do you wish to calculate on a cation basis? (y|n): ';
wantcation=input(prompt1, 's');

if strcmp(wantcation, 'y')
    disp('You are calculating on a cation basis.')
    prompt2='How many cations do you wish to normalize to?: ';
    CM=input(prompt2);

    prompt3='Do you wish to treat As as a cation? (y|n): ';
    wantAs=input(prompt3,'s');

    if strcmp(wantAs, 'y')
        disp('As will be treated as a cation.')
        NF=CM./sum(MC(:,2:8),2); %normalization factor w/ As
    else
        disp('As will be treated as an anion.')
        NF=CM./sum(MC(:,3:8),2); %normalization factor w/o As
    end
else
    disp('You are calculating on an anion basis.')
    prompt2='How many anions do you wish to normalize to?: ';
    AM=input(prompt2);

    prompt3='Do you wish to treat As as a cation? (y|n): ';
    wantAs=input(prompt3,'s');

    if strcmp(wantAs, 'y')
        disp('As will be treated as a cation.')
        NF=AM./MC(:,1); %normalization factor w/o As
    else
        disp('As will be treated as an anion.')
        NF=AM./(MC(:,1)+MC(:,2)); %normalization factor w/ As
    end
end

APFU=MC.*NF; %Normalized moles = moles of elements * normalization factor

APFU=array2table(APFU,'VariableNames',{'S','As','Fe','Ni','Cu','Co','Pb','Zn'});

end






