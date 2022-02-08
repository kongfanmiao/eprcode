% operators
clear;

sx = sop(1/2,'x');
sy = sop(1/2,'y');
sz = sop(1/2,'z');
one = eye(2);

ix = sx;
iy = sy;
iz = sz;

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

