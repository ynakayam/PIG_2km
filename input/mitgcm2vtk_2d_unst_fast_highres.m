%MITgcm2vtk written by Yoshihiro Nakayama 04/05/2016

clear;
%Set output time
 addpath /u/ynakaya2/slr/components/mitgcm/install/utils/matlab/cs_grid/read_cs/

runid='Latlon_run13_2d_fast_ver3'
time_start=30*1;
time_end  =30*24*100;
time_int  =30;

hFacC=rdmds('hFacC');
DD=size(50,60);

%Read Grids
XC=rdmds('XC');
YC=rdmds('YC');
ZC=rdmds('RC');
DEPTH=readbin('topo_2000m.bin', [50 60]);
ZCC=reshape(ZC(1,1,:),1,max(size(ZC)));

%Changing ratio for Z
XCm=XC+360;
YCm=YC.*3;%5;
ZCCm=ZCC./400.;
cells2=0; 


%Prepare data for vtk
lo=XCm(:);la=YCm(:);

np=length(la)*1;  % number of points
ncell=(size(XC,1)-1)*(size(XC,2)-1); %*(max(size(ZCCm(1,:)))-1);
celltype(1:ncell,1)=8;

sl=size(XC,1);
ce=[];
sss=max(size(XCm(:,1)))*max(size(YCm(1,:)))
%for k=1:(max(size(ZCCm(1,:)))-1)
k=1;
for i=1:size(XC,2)-1;
cc=[[1+(sl*(i-1))+(k-1)*sss:sl-1+(sl*(i-1))+(k-1)*sss]',...
    [2+(sl*(i-1))+(k-1)*sss:sl+(sl*(i-1))+(k-1)*sss]',...
    [[2+(k-1)*sss:sl+(k-1)*sss]+sl+(sl*(i-1))]',...
    [[1+(k-1)*sss:sl-1+(k-1)*sss]+sl+(sl*(i-1))]'  ];
ce=[ce;cc];
end
%end
cells(1:size(ce,1),1)=4;

cells2=cells2+max(size(XCm(:,1)))*max(size(YCm(1,:)));
ce=ce-1;
%ce2=ce+max(size(XCm(:,1)))*max(size(YCm(1,:)));
cells=[cells,ce];
ncell2=size(cells(:),1);

for k=1:1; %max(size(ZCCm(1,:)))
for kk=1:np
ZCCCm(kk+(k-1)*np)=ZCCm(1,k);    
lok(kk+(k-1)*np)=lo(kk);
lak(kk+(k-1)*np)=la(kk);
end;
end;


for t=time_start:time_int:time_end;

time=t;
  T=rdmds('T',time);
  S=rdmds('S',time);
  U=rdmds('U',time);
  V=rdmds('V',time);
  W=rdmds('W',time);

Mflux=rdmds('SHICE_fwFlux',time);
icetopo=readbin('itopo_2000m.bin', [50 60]);

   Eta    =rdmds('ETAtave',time);
   siAREA =rdmds('AREAtave',time);
   siHEFF =rdmds('HEFFtave',time);

T105   =T(:,:,10);
S105   =S(:,:,10);
U105   =U(:,:,10);
V105   =V(:,:,10);
W105   =W(:,:,10);
T222   =T(:,:,20);
S222   =S(:,:,20);
U222   =U(:,:,20);
V222   =V(:,:,20);
W222   =W(:,:,20);
T299   =T(:,:,30);
S299   =S(:,:,30);
U299   =U(:,:,30);
V299   =V(:,:,30);
W299   =W(:,:,30);
T409   =T(:,:,40);
S409   =S(:,:,40);
U409   =U(:,:,40);
V409   =V(:,:,40);
W409   =W(:,:,40);
T477   =T(:,:,50);
S477   =S(:,:,50);
U477   =U(:,:,50);
V477   =V(:,:,50);
W477   =W(:,:,50);
T552   =T(:,:,60);
S552   =S(:,:,60);
U552   =U(:,:,60);
V552   =V(:,:,60);
W552   =W(:,:,60);

TBot=T105.*0.0;
SBot=S105.*0.0;
ccount=0;

% search good value

for i=1:50;
for j=1:60;
set=0;
DD(i,j)= 0;
for k=130:-1:1
if(S(i,j,k)>1.0 & set==0);
TBot(i,j)   =T(i,j,k);
SBot(i,j)   =S(i,j,k);
set=1;
DD(i,j)= k;
end
end
end
end
save DD.mat DD

load DD;
%DD(1977:2016,1:768)=0; %just because i forgot to preset DD (initially)

for i=1:50;
for j=1:60;
    if(DD(i,j)>0);
TBot(i,j)=T(i,j,DD(i,j));
SBot(i,j)=S(i,j,DD(i,j));
    end
end
end

for k=1:1;
for j=1:max(size(YCm(1,:)));
for i=1:max(size(XCm(:,1)));

    tot=i+(j-1)*max(size(XCm(:,1)))+(k-1)*max(size(XCm(:,1)))*max(size(YCm(1,:)));
    Et(tot)=Eta(i,j);
    deP(tot)=DEPTH(i,j);
    siA(tot)=siAREA(i,j);
    siHE(tot)=siHEFF(i,j);

siT105(tot)=T105(i,j);
siS105(tot)=S105(i,j);
siU105(tot)=U105(i,j);
siV105(tot)=V105(i,j);
siW105(tot)=W105(i,j);

siT222(tot)=T222(i,j);
siS222(tot)=S222(i,j);
siU222(tot)=U222(i,j);
siV222(tot)=V222(i,j);
siW222(tot)=W222(i,j);

siT299(tot)=T299(i,j);
siS299(tot)=S299(i,j);
siU299(tot)=U299(i,j);
siV299(tot)=V299(i,j);
siW299(tot)=W299(i,j);

siT409(tot)=T409(i,j);
siS409(tot)=S409(i,j);
siU409(tot)=U409(i,j);
siV409(tot)=V409(i,j);
siW409(tot)=W409(i,j);

siT477(tot)=T477(i,j);
siS477(tot)=S477(i,j);
siU477(tot)=U477(i,j);
siV477(tot)=V477(i,j);
siW477(tot)=W477(i,j);

siT552(tot)=T552(i,j);
siS552(tot)=S552(i,j);
siU552(tot)=U552(i,j);
siV552(tot)=V552(i,j);
siW552(tot)=W552(i,j);

siTBot(tot)=TBot(i,j);
siSBot(tot)=SBot(i,j);
siMflux(tot)=Mflux(i,j);
    
if(icetopo(i,j)==0 & S(i,j,1)>10.)
  icetopo(i,j)=1.;
  icet(tot)=icetopo(i,j);
elseif(icetopo(i,j)==0 & S(i,j,1)<10.)
  icetopo(i,j)=0.;
  icet(tot)=icetopo(i,j);
elseif(icetopo(i,j)~=0 )
  icet(tot)=-100;
end

end;
end;
end;

%%%convert VTK
%  writeVTK_scalar_test(filename,XCm(:,1),YCm(1,:),ZCCm(1,:),'hFacC',hFacC,'T',T,'S',S,'RHO',R,'U',U,'V',V,'W',W)

%%%VTK converter unstractured mesh %%%
filename=[runid,'_',num2str(time),'.vtk']
fid=fopen(filename,'w');
%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from RTOPO\n');
fprintf(fid, 'BINARY\n');
%fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

fprintf(fid, ['POINTS ', num2str(np) ,' float\n']);
    
fwrite(fid, [reshape(lok,1,np*1); reshape(lak,1,np*1); reshape(ZCCCm,1,np*1)],'float','b');


%for k=1:max(size(ZCCm(1,:)))
%for i=1:np
%	fprintf(fid,'%.2f %.2f %.2f\n' ,[lo(i); la(i); ZCCm(1,k)]);
%end
%end
fprintf(fid, '\n\n');

index1=cells(:,1);index2=cells(:,2);index3=cells(:,3);index4=cells(:,4);
index5=cells(:,5);%index6=cells(:,6);index7=cells(:,7);index8=cells(:,8);index9=cells(:,9);
numcell=max(size(cells));

fprintf(fid, ['CELLS ',num2str(ncell),' ',num2str(ncell2) '\n']);
fwrite(fid, [reshape(index1,1,numcell); reshape(index2,1,numcell); reshape(index3,1,numcell); reshape(index5,1,numcell);...
             reshape(index4,1,numcell)],'int','b');
%fprintf(fid,'%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n' ,cells');

fprintf(fid, '\n\n');

fprintf(fid, ['CELL_TYPES ',num2str(ncell),'\n']);
fwrite(fid, reshape(celltype,1,max(size(celltype))),'int','b');
%fprintf(fid,'%.0f\n ' ,celltype);
fprintf(fid, '\n\n');


fprintf(fid, ['POINT_DATA ' num2str(np*1) ' \n']);
fprintf(fid, ['SCALARS Eta FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(Et,1,np*1),'float','b'); %binary data

fprintf(fid, ['SCALARS siAREA FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(siA,1,np*1),'float','b'); %binary data

fprintf(fid, ['SCALARS siHEFF FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(siHE,1,np*1),'float','b'); %binary data

fprintf(fid, ['SCALARS DEPTH FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(deP,1,np*1),'float','b')

fprintf(fid, ['SCALARS icetopo FLOAT 1\n']); 
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(icet,1,np*1),'float','b') 

fprintf(fid, ['SCALARS DEPTH FLOAT 1\n']); 
fprintf(fid, ['LOOKUP_TABLE default\n']);
fwrite (fid, reshape(deP,1,np*1),'float','b') 

fprintf(fid, ['SCALARS siMflux FLOAT 1\n']); 
fprintf(fid, ['LOOKUP_TABLE default\n']);
fwrite (fid, reshape(siMflux,1,np*1),'float','b')         

fprintf(fid, ['SCALARS T100 FLOAT 1\n']);            
fprintf(fid, ['LOOKUP_TABLE default\n']);                  
fwrite (fid, reshape(siT105,1,np*1),'float','b') 

fprintf(fid, ['SCALARS S100 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siS105,1,np*1),'float','b')  

fprintf(fid, ['SCALARS T200 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siT222,1,np*1),'float','b')  

fprintf(fid, ['SCALARS S200 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siS222,1,np*1),'float','b')  

fprintf(fid, ['SCALARS T300 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siT299,1,np*1),'float','b')  

fprintf(fid, ['SCALARS S300 FLOAT 1\n']);          
fprintf(fid, ['LOOKUP_TABLE default\n']); 
fwrite (fid, reshape(siS299,1,np*1),'float','b')       

fprintf(fid, ['SCALARS T400 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siT409,1,np*1),'float','b')  

fprintf(fid, ['SCALARS S400 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siS409,1,np*1),'float','b')  

fprintf(fid, ['SCALARS T500 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siT477,1,np*1),'float','b')  

fprintf(fid, ['SCALARS S500 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siS477,1,np*1),'float','b')  

fprintf(fid, ['SCALARS T600 FLOAT 1\n']);                                
fprintf(fid, ['LOOKUP_TABLE default\n']);                                     
fwrite (fid, reshape(siT552,1,np*1),'float','b')  

fprintf(fid, ['SCALARS S600 FLOAT 1\n']);                                       
fprintf(fid, ['LOOKUP_TABLE default\n']);                                       
fwrite (fid, reshape(siS552,1,np*1),'float','b')   

fprintf(fid, ['SCALARS TBot FLOAT 1\n']);                                       
fprintf(fid, ['LOOKUP_TABLE default\n']);                                       
fwrite (fid, reshape(siTBot,1,np*1),'float','b')   

fprintf(fid, ['SCALARS SBot FLOAT 1\n']);  
fprintf(fid, ['LOOKUP_TABLE default\n']);
fwrite (fid, reshape(siSBot,1,np*1),'float','b')   


for k=1:np*1                                                       
     vdata(3*k-2)= siU105(k);                             
     vdata(3*k-1)= siV105(k);                             
     vdata(3*k  )= siW105(k);                                      
end      

fprintf(fid, ['VECTORS v_100m float\n']);                    
fwrite (fid, reshape(vdata,1,np*1*3),'float','b'); %binary data  

for k=1:np*1                                                       
     vdata(3*k-2)= siU222(k);                             
     vdata(3*k-1)= siV222(k);                             
     vdata(3*k  )= siW222(k);                                      
end                                                                
                                                                   
fprintf(fid, ['VECTORS v_200m float\n']);                    
fwrite (fid, reshape(vdata,1,np*1*3),'float','b'); %binary data                              
                                      
for k=1:np*1                                                       
     vdata(3*k-2)= siU299(k);                             
     vdata(3*k-1)= siV299(k);                             
     vdata(3*k  )= siW299(k);                                      
end                                                                
                                                                   
fprintf(fid, ['VECTORS v_300m float\n']);                    
fwrite (fid, reshape(vdata,1,np*1*3),'float','b'); %binary data  

for k=1:np*1                                                       
     vdata(3*k-2)= siU409(k);                             
     vdata(3*k-1)= siV409(k);                             
     vdata(3*k  )= siW409(k);                                      
end                                                                
                                                                   
fprintf(fid, ['VECTORS v_400m float\n']);                    
fwrite (fid, reshape(vdata,1,np*1*3),'float','b'); %binary data  


for k=1:np*1                                                       
     vdata(3*k-2)= siU477(k);                             
     vdata(3*k-1)= siV477(k);                             
     vdata(3*k  )= siW477(k);                                      
end                                                                
                                                                   
fprintf(fid, ['VECTORS v_500m float\n']);                    
fwrite (fid, reshape(vdata,1,np*1*3),'float','b'); %binary data  

for k=1:np*1     
     vdata(3*k-2)= siU552(k); 
     vdata(3*k-1)= siV552(k);
     vdata(3*k  )= siW552(k); 
end   
fprintf(fid, ['VECTORS v_600m float\n']);
fwrite (fid, reshape(vdata,1,np*1*3),'float','b'); %binary data       

fclose(fid);

end


