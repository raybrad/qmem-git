function  printVoltage(V,currtime,dt)

global qmx1 qmx2 qmy1 qmy2 qmz1 qmz2;
%
leftnodes=defBrickNodes([qmx1,qmy1,qmz1],[qmx1,qmy2,qmz2]);
rghtnodes=defBrickNodes([qmx2,qmy1,qmz1],[qmx2,qmy2,qmz2]);
topnodes =defBrickNodes([qmx1,qmy1,qmz2],[qmx2,qmy2,qmz2]);

ltv=mean(V(leftnodes));

rtv=mean(V(rghtnodes));

ttv=mean(V(topnodes));

fp=fopen('voltage.dat','a+');
fprintf(fp,'%f %13.5e %13.5e %13.5e \n',currtime,ltv,rtv,ttv);
fclose(fp);

