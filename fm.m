clear;
fname='D:\MinGW\GPSS\R_12.TXT';
frmt='%d  %d %f %f %d %f %f %d %f %f %d %d %f %f %d %d';
fid=fopen(fname);
ncx=2;
ncy=7;
%ncz=9;
%ncz=13;
%ncz=14;
ix=0;
iy=0;
while (feof(fid)==0)
     s=fgets(fid);
    [A,count] = sscanf(s,frmt);
    if (ix==0)
       X=A(ncx);
    else
       m=numel(X);
       finded=false;
       for ixx=1:m
           if (X(ixx)==A(ncx))
              finded=true;
              break; 
           end
       end 
       if (finded==false)
          X=[X,A(ncx)];
       end
       
    end    
    if (iy==0)
       Y=A(ncy);
    else
       m=numel(Y);
       finded=false;
       for iyy=1:m
           if (abs(abs(Y(iyy))-abs(A(ncy)))<0.005)
              finded=true;
              break; 
           end
       end 
       if (finded==false)
          Y=[Y,A(ncy)];
       end
    end    
    ix=ix+1;
    iy=iy+1;
end
f=fclose(fid);
fid=fopen(fname);
iz=0;
iy=0;
mr=numel(X);
mc=numel(Y);
Z=zeros(mr,mc);
XX=sort(X);
YY=sort(Y);
for i=1:mr
    for j=1:mc
        Z(i,j)=-1;
    end    
end
while (feof(fid)==0)
    s=fgets(fid);
    [A,count] = sscanf(s,frmt);
    m=mr;
    ir=0;
    for ixx=1:m
        if (XX(ixx)==A(ncx))
           finded=true;
           ir=ixx;
           break; 
        end
    end 
    if (ir>0)
        n=mc;
        ic=0;
        for iyy=1:n
            if (abs(abs(YY(iyy))-abs(A(ncy)))<0.005)
               finded=true;
               ic=iyy;
               break; 
            end     
        end    
        if ((ir>0) & (ic>0))
           if (Z(ir,ic) <A(ncz))
               Z(ir,ic) =A(ncz);
           end    
        end
    end 
end
for i=1:mr
    for j=1:mc
        if (Z(i,j)==-1)
            Z(i,j)=0;
        end
    end    
end
f=fclose(fid);
[XXX,YYY]=meshgrid(XX,YY);
ZZ=Z';
colormap(white);
surfc(XXX,YYY,ZZ);
grid ON;

