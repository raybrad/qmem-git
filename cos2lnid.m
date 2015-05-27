function lnid = cos2lnid(nco1,nco2)

global kx ky kz;

global links;
global Nlink;

n1 = co2id(nco1);
n2 = co2id(nco2); 

for i = 1:Nlink
    if n1 == links(i,1) && n2 == links(i,2)
       lnid =i;
    end
end
