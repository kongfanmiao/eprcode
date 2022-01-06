% operators
clear;

sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];
one = eye(2);

sx = 1/2*sigma_x;
sy = 1/2*sigma_y;
sz = 1/2*sigma_z;
ix = 1/2*sigma_x;
iy = 1/2*sigma_y;
iz = 1/2*sigma_z;

I = kron(one,one);
Sx = kron(sx,one); 
Sy = kron(sy,one); 
Sz = kron(sz,one);
Ix = kron(one,ix); 
Iy = kron(one,iy); 
Iz = kron(one,iz);
SzIz = kron(sz,iz);
SzIx = kron(sz,ix);
SzIy = kron(sz,iy);
SxIz = kron(sx,iz);
SyIz = kron(sy,iz);
SxIx = kron(sx,ix);
SyIy = kron(sy,iy);
SxIy = kron(sx,iy);
SyIx = kron(sy,ix);

basis = {I,Sx,Sy,Sz,Ix,Iy,Iz,...
    SzIz,SxIz,SyIz,SzIx,SzIy,SxIx,SyIy,SxIy,SyIx};

