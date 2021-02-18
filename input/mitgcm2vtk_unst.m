% MITgcm2vtk written by Yoshihiro Nakayama 04/05/2016
% 
clear all;

%Set output time
runid='ABS_latlon_3D_run12'
time_start=360;
time_end  =360*24*30;
time_int  =360*12;
addpath /u/ynakaya2/sw_lib;

hFacC=rdmds('hFacC');
hFacW=rdmds('hFacC');
hFacS=rdmds('hFacC');

%Read Grids
XC=rdmds('XC');
YC=rdmds('YC');
ZC=rdmds('RC');
ZCC=reshape(ZC(1,1,:),1,max(size(ZC)));

%Changing ratio for Z
XCm=XC+360;
YCm=YC.*3;%5;
ZCCm=ZCC./400.;

%Prepare data for vtk
%np=max(size(XCm(:,1)))*max(size(YCm(1,:)))*max(size(ZCCm(1,:)))
lo=XCm(:);la=YCm(:);

np=length(la)*1;  % number of points
ncell=(size(XC,1)-1)*(size(XC,2)-1)*(max(size(ZCCm(1,:)))-1);
celltype(1:ncell,1)=12;

sl=size(XC,1);
ce=[];
sss=max(size(XCm(:,1)))*max(size(YCm(1,:)))
for k=1:(max(size(ZCCm(1,:)))-1)
for i=1:size(XC,2)-1;
cc=[[1+(sl*(i-1))+(k-1)*sss:sl-1+(sl*(i-1))+(k-1)*sss]',...
    [2+(sl*(i-1))+(k-1)*sss:sl+(sl*(i-1))+(k-1)*sss]',...
    [[2+(k-1)*sss:sl+(k-1)*sss]+sl+(sl*(i-1))]',...
    [[1+(k-1)*sss:sl-1+(k-1)*sss]+sl+(sl*(i-1))]'  ];
ce=[ce;cc];
end
end
cells(1:size(ce,1),1)=8;
cells2=cells;
cells2=cells2+max(size(XCm(:,1)))*max(size(YCm(1,:)));
ce=ce-1;
ce2=ce+max(size(XCm(:,1)))*max(size(YCm(1,:)));
cells=[cells,ce,ce2];
ncell2=size(cells(:),1);

for k=1:max(size(ZCCm(1,:)))
for kk=1:np
ZCCCm(kk+(k-1)*np)=ZCCm(1,k);    
lok(kk+(k-1)*np)=lo(kk);
lak(kk+(k-1)*np)=la(kk);
end;
end;

for t=time_start:time_int:time_end;

time=t
U=rdmds('U',t);
V=rdmds('V',t);
W=rdmds('W',t);
T=rdmds('T',t);
S=rdmds('S',t);

for i=1:max(size(XCm(:,1)));
  for j=1:max(size(YCm(1,:)));
    for k=1:max(size(ZCCm(1,:)));
     if(S(i,j,k)>10.)
      R(i,j,k)=sw_dens(S(i,j,k),T(i,j,k),0);
     else
      R(i,j,k)=0;
     end
    end
  end
end
%A=fld8(:,:,:,1);
%B=fld8(:,:,:,2);
%C=fld8(:,:,:,3);
%D=fld8(:,:,:,4);
%E=fld8(:,:,:,5);
%F=fld8(:,:,:,6);
%G=fld8(:,:,:,7);
%H=fld8(:,:,:,8);
%I=fld8(:,:,:,9);
%J=fld8(:,:,:,10);

filename=[runid,'_',num2str(time),'.vtk']

for k=1:max(size(ZCCm(1,:)));
for j=1:max(size(YCm(1,:)));
for i=1:max(size(XCm(:,1)));
    tot=i+(j-1)*max(size(XCm(:,1)))+(k-1)*max(size(XCm(:,1)))*max(size(YCm(1,:)));
    hF(tot)=hFacC(i,j,k);
    hFS(tot)=hFacS(i,j,k);
    hFW(tot)=hFacW(i,j,k);
    te(tot)=T(i,j,k);
    sa(tot)=S(i,j,k);
    rh(tot)=R(i,j,k);
    uv(tot)=U(i,j,k)*100; % m/s -> cm/s
    vv(tot)=V(i,j,k)*100;
    wv(tot)=W(i,j,k)*100;

%    aa(tot)=A(i,j,k);
%    bb(tot)=B(i,j,k);
%    ccc(tot)=C(i,j,k);
%    dd(tot)=D(i,j,k);
%    ee(tot)=E(i,j,k);
%    ff(tot)=F(i,j,k);
%    gg(tot)=G(i,j,k);
%    hh(tot)=H(i,j,k);
%    ii(tot)=I(i,j,k);
%    jj(tot)=J(i,j,k);
end;
end;
end;

%%%convert VTK
%  writeVTK_scalar_test(filename,XCm(:,1),YCm(1,:),ZCCm(1,:),'hFacC',hFacC,'T',T,'S',S,'RHO',R,'U',U,'V',V,'W',W)
%%%VTK converter unstractured mesh %%%

fid=fopen(filename,'w');
%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from RTOPO\n');
fprintf(fid, 'BINARY\n');
%fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

fprintf(fid, ['POINTS ', num2str(np*max(size(ZCCm(1,:)))) ,' float\n']);
d    
fwrite(fid, [reshape(lok,1,np*max(size(ZCCm(1,:)))); reshape(lak,1,np*max(size(ZCCm(1,:)))); reshape(ZCCCm,1,np*max(size(ZCCm)))],'float','b');


%for k=1:max(size(ZCCm(1,:)))
%for i=1:np
%	fprintf(fid,'%.2f %.2f %.2f\n' ,[lo(i); la(i); ZCCm(1,k)]);
%end
%end
fprintf(fid, '\n\n');

index1=cells(:,1);index2=cells(:,2);index3=cells(:,3);index4=cells(:,4);
index5=cells(:,5);index6=cells(:,6);index7=cells(:,7);index8=cells(:,8);index9=cells(:,9);
numcell=max(size(cells));

fprintf(fid, ['CELLS ',num2str(ncell),' ',num2str(ncell2) '\n']);
fwrite(fid, [reshape(index1,1,numcell); reshape(index2,1,numcell); reshape(index3,1,numcell); reshape(index4,1,numcell);...
             reshape(index5,1,numcell); reshape(index6,1,numcell); reshape(index7,1,numcell); reshape(index8,1,numcell);...
             reshape(index9,1,numcell)],'int','b');
%fprintf(fid,'%.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f %.0f\n' ,cells');

fprintf(fid, '\n\n');

fprintf(fid, ['CELL_TYPES ',num2str(ncell),'\n']);
fwrite(fid, reshape(celltype,1,max(size(celltype))),'int','b');
%fprintf(fid,'%.0f\n ' ,celltype);
fprintf(fid, '\n\n');

fprintf(fid, ['POINT_DATA ' num2str(np*max(size(ZCCm(1,:)))) ' \n']);
fprintf(fid, ['SCALARS hFacS FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n']);
fwrite (fid, reshape(hFS,1,np*max(size(ZCCm(1,:)))),'float','b'); %binary data  

fprintf(fid, ['SCALARS Salt FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(sa,1,np*max(size(ZCCm(1,:)))),'float','b'); %binary data

fprintf(fid, ['SCALARS Temp FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(te,1,np*max(size(ZCCm(1,:)))),'float','b'); %binary data

fprintf(fid, ['SCALARS RHO FLOAT 1\n']);                                              
fprintf(fid, ['LOOKUP_TABLE default\n'])                                                
fwrite (fid, reshape(rh,1,np*max(size(ZCCm(1,:)))),'float','b'); %binary data   

fprintf(fid, ['SCALARS hFacC FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(hF,1,np*max(size(ZCCm(1,:)))),'float','b'); %binary data

%fprintf(fid, ['SCALARS hFacW FLOAT 1\n']);
%fprintf(fid, ['LOOKUP_TABLE default\n']);
%fwrite (fid, reshape(hFW,1,np*max(size(ZCCm(1,:)))),'float','b'); %binary data     

%{
fprintf(fid, ['SCALARS WTHMASS FLOAT 1\n']);
fprintf(fid, ['LOOKUP_TABLE default\n'])
fwrite (fid, reshape(aa,1,np*max(size(ZCCm(1,:)))),'float','b'); 

fprintf(fid, ['SCALARS WSLTMASS FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(bb,1,np*max(size(ZCCm(1,:)))),'float','b'); 

fprintf(fid, ['SCALARS ADVr_TH FLOAT 1\n']);            
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(ccc,1,np*max(size(ZCCm(1,:)))),'float','b'); 

fprintf(fid, ['SCALARS DFrE_TH FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(dd,1,np*max(size(ZCCm(1,:)))),'float','b');
 
fprintf(fid, ['SCALARS DFrI_TH FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(ee,1,np*max(size(ZCCm(1,:)))),'float','b');
 
fprintf(fid, ['SCALARS ADVr_SLT FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(ff,1,np*max(size(ZCCm(1,:)))),'float','b'); 

fprintf(fid, ['SCALARS DFrE_SLT FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(gg,1,np*max(size(ZCCm(1,:)))),'float','b'); 

fprintf(fid, ['SCALARS DFrI_SLT FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(hh,1,np*max(size(ZCCm(1,:)))),'float','b'); 

fprintf(fid, ['SCALARS GM_Kwz FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(ii,1,np*max(size(ZCCm(1,:)))),'float','b'); 

fprintf(fid, ['SCALARS KPPdiffT FLOAT 1\n']);                           
fprintf(fid, ['LOOKUP_TABLE default\n'])                           
fwrite (fid, reshape(jj,1,np*max(size(ZCCm(1,:)))),'float','b'); 
%}
%{
fprintf(fid, ['VECTORS vector FLOAT\n']);

for i=1:np*max(size(ZCCm(1,:)))
fprintf(fid,'%.4f %.4f %.4f\n', [uv(i); vv(i); wv(i)]); 
end
fprintf(fid, '\n\n');
%}
for k=1:np*max(size(ZCCm(1,:)))
     vdata(3*k-2)=uv(k);
     vdata(3*k-1)=vv(k);
     vdata(3*k  )=wv(k);
end

fprintf(fid, ['VECTORS velocity float\n']);
fwrite (fid, reshape(vdata,1,np*max(size(ZCCm(1,:)))*3),'float','b'); %binary data
%clear cells;

fclose(fid);

end


 

