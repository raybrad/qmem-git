function displaySlice(variable,sliceNr,sliceDim)

%%% display variables on a given slide of the structure %%%%%%%%
%%% variable: Vs, ns, etc
%%% sliceNr: the number of the slice
%%% sliceDim: direction of the slice ('x','y' or 'z')

global kx ky kz;
global nodes;
switch sliceDim
    case 'x'
        sliceNodes = defBrickNodes([sliceNr,1,1],[sliceNr,ky+1,kz+1]);
        x = nodes(sliceNodes,2);
        y = nodes(sliceNodes,3);
        nx = ky+1; ny = kz+1;
    case 'y'
        sliceNodes = defBrickNodes([1,sliceNr,1],[kx+1,sliceNr,kz+1]);
        x = nodes(sliceNodes,1);
        y = nodes(sliceNodes,3);
        nx = kx+1; ny = kz+1;
    case 'z'
        sliceNodes = defBrickNodes([1,1,sliceNr],[kx+1,ky+1,sliceNr]);
        x = nodes(sliceNodes,1);
        y = nodes(sliceNodes,2);
        nx = kx+1; ny = ky+1;
end
x = reshape(x,nx,ny);
y = reshape(y,nx,ny);
sliceNodes_ = reshape(sliceNodes,nx,ny);
values = variable(sliceNodes_);
figure;
pcolor(x,y,values);
shading interp;

