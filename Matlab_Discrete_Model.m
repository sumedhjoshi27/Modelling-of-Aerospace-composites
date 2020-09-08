clc
clear 

%INPUT DATA
%*************************************************************************%
%*************************************************************************%
lengthX=100;
lengthY=200;
thickness = 10;

N_col = 10; %X direction
N_row = 20; %Y direction

%*************************************************************************%
%*************************************************************************%

D = lengthX/N_col;

Inp=fopen('Output.inp','wt+');
fprintf(Inp,'%s\n','*Node');
StrNodeFormat='%d,  %d,  %d,  %d';
StrShellFormat='%d,  %d,  %d,  %d,  %d';
StrTriFormat='%d, %d, %d, %d';
StrBeamFormat='%d,  %d,  %d';
StrEsetFormat = '%d,  %d,   1';

%***************NODE GENERATION*****************%

TotalNode = (N_col+1)*(N_row +2) + (N_col+2)*(N_row + 1);

for i=1:(2*N_row+3)
    if mod(i,2)~=0
       for j=1:(N_col+1)
           NodeNumber = j+  (2*N_col+3)*(i-1)/2;
           index{NodeNumber} = NodeNumber;
           x{NodeNumber} =  (j-1)*D;
           y{NodeNumber} = (i-2)*(D/2);
           z{NodeNumber} = 0;
       end
    else
        for j=1:(N_col+2)
           NodeNumber = j+ (N_col+1) + (2*N_col+3)*(i-2)/2;
           index{NodeNumber} = NodeNumber;
           x{NodeNumber} =  (2*j-3)*(D/2);
           y{NodeNumber} = (i-2)*(D/2);
           z{NodeNumber} = 0;
       end  
    end
end

for i=1: TotalNode
    Node_M{i}=sprintf(StrNodeFormat,i,x{i},y{i},z{i});
end
for i=1 : (TotalNode)
    fprintf(Inp,'%s\n',Node_M{i});
end    

%***************SHELL ELEMENT GENERATION*****************%


fprintf(Inp,'%s\n','*Element, type=S4R');

E_index = 0;
for i=1:(2*N_row+1)
    if mod(i,2)~=0
       for j=(1+(2*N_col+3)*(i-1)/2) : (N_col+1  +(2*N_col+3)*(i-1)/2)
           E_index = E_index +1;
           Shell{E_index} = sprintf(StrShellFormat,E_index,j,(j+N_col+2),(j+2*N_col+3),(j+N_col+1));
           fprintf(Inp,'%s\n',Shell{E_index});
       end
    else
        for j=(1+ (N_col+2) + (2*N_col+3)*(i-2)/2) : (N_col + (N_col+2) + (2*N_col+3)*(i-2)/2)
           E_index = E_index +1;
           Shell{E_index} = sprintf(StrShellFormat,E_index,j,(j+N_col+2),(j+2*N_col+3),(j+N_col+1));
           fprintf(Inp,'%s\n',Shell{E_index});
       end  
    end
end
TotalElement = E_index;


fprintf(Inp,'%s\n','*Elset, elset=QuadshellElem, generate');
QuadSet = sprintf(StrEsetFormat,1,TotalElement);
fprintf(Inp,'%s\n',QuadSet);



%---------------Generation of Triangular shell elements--------------%
%--------------------------------------------------------------------%


fprintf(Inp,'%s\n','*Element, type=S3R');

Tri_index=0;
for i=[1 (2*N_row+3)]
        if (i==1)
            for j=1:N_col
                Tri_index=Tri_index+1;
                Triele{Tri_index}=sprintf(StrTriFormat,Tri_index+TotalElement,j,j+N_col+2,j+1);
            end
        end
              
        if (i==(2*N_row+3))
            for j=(1+(2*N_col+3)*(i-1)/2) : (N_col+(2*N_col+3)*(i-1)/2)
                   Tri_index=Tri_index+1;
                   Triele{Tri_index}=sprintf(StrTriFormat,Tri_index+TotalElement,j,j-N_col-1,j+1);
            end
        end
end



for j= [1 2]
    if (j==1)
        for i=1:N_row
            Tri_index=Tri_index+1;
            k=((2*N_col+3)*(i-1));
            Triele{Tri_index}=sprintf(StrTriFormat,Tri_index+TotalElement,N_col+2+k,(2*N_col)+4+k,(3*N_col)+5+k);
        end
    end
    
    if (j==2)
        for i=1:N_row
            Tri_index=Tri_index+1;
            k=((2*N_col+3)*(i-1));
            Triele{Tri_index}=sprintf(StrTriFormat,Tri_index+TotalElement,(2*N_col)+3+k,(3*N_col)+4+k,(4*N_col)+6+k);
        end
    end
end


for i=1:Tri_index
    fprintf(Inp,'%s\n', Triele{i});
end


fprintf(Inp,'%s\n','*Elset, elset=TriShellElement, generate');
TriSet = sprintf(StrEsetFormat,1+TotalElement,Tri_index+TotalElement);
fprintf(Inp,'%s\n',TriSet);



%----------------Beam Element Generation------------------%


fprintf(Inp,'%s\n','*Element, type=B31');

 B_index=-1;
 for i= 1:(2*N_row+2)
     if mod(i,2)~=0
         for j=(1+(2*N_col+3)*(i-1)/2) : (N_col+1  +(2*N_col+3)*(i-1)/2)
             B_index=B_index+2;
             Beam{B_index}=sprintf(StrBeamFormat,B_index+Tri_index+TotalElement,j,j+N_col+1);
             Beam{B_index+1}=sprintf(StrBeamFormat, B_index+1+Tri_index+TotalElement,j,j+N_col+2);
         end
     else
      
         for j=(N_col+2 + (2*N_col+3)*(i-2)/2) : ((2*N_col+2) + (2*N_col+3)*(i-2)/2)
              B_index=B_index+2;
              Beam{B_index}=sprintf(StrBeamFormat,B_index+Tri_index+TotalElement,j,j+N_col+2);
              Beam{B_index+1}=sprintf(StrBeamFormat,B_index+1+Tri_index+TotalElement,j+N_col+2,j+1);
         end
     end
 end
B_index=B_index+1;
for i=[1 (2*N_row+3)]
        if (i==1)
            for j=1:N_col
                B_index=B_index+1;
                Beam{B_index}=sprintf(StrBeamFormat,B_index+Tri_index+TotalElement,j,j+1);
            end
        end
        
        if (i==(2*N_row+3))
            for j=(1+(2*N_col+3)*(i-1)/2) : (N_col+(2*N_col+3)*(i-1)/2)
                B_index=B_index+1;
                Beam{B_index}=sprintf(StrBeamFormat,B_index+Tri_index+TotalElement,j,j+1);
            end
        end
end

for j= [1 2]
    if (j==1)
        for i=1:N_row  
            B_index=B_index+1;
            k=((2*N_col+3)*(i-1));
            Beam{B_index}=sprintf(StrBeamFormat,B_index+Tri_index+TotalElement,N_col+2+k,(3*N_col)+5+k);
        end
    end
    
    if (j==2)
        for i=1:N_row  
            B_index=B_index+1;
            k=((2*N_col+3)*(i-1));
            Beam{B_index}=sprintf(StrBeamFormat,B_index+Tri_index+TotalElement,(2*N_col)+3+k,(4*N_col)+6+k);
        end
    end
end
 
            
 for i=1:B_index
     fprintf(Inp,'%s\n', Beam{i});
 end
 
fprintf(Inp,'%s\n','*Elset, elset=BeamLayer, generate');
Beam1 = sprintf(StrEsetFormat,1+Tri_index+TotalElement,B_index+Tri_index+TotalElement);
fprintf(Inp,'%s\n',Beam1);


fclose(Inp);

        