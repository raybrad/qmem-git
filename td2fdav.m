function [frev,frea] = td2fdav(nsteps,nstat,ffttn,ele)

global Nnode Nlink;

[VT,nT,pT,AT,HT]=tdvara(nsteps);

vtt=VT';
att=AT';

vabs=[];
aabs=[];
vang=[];
aang=[];

ffta=[];
fftv=[];

stata=[];
statv=[];

temp1=0;
temp2=0;

for i=1:Nnode
  
    fftv=fft(vtt(nstat:nsteps,i),ffttn);
    
    temp1=abs(fftv(ele))/ffttn;
    temp2=abs(fftv(1))/ffttn;

    vabs=[vabs;4*temp1];
    vang=[vang;180*angle(fftv(ele))/pi];

    statv=[statv;(temp2-2*temp1)];
    
end

for i=1:Nlink

    ffta=fft(att(nstat:nsteps,i),ffttn);

    temp1=abs(ffta(ele))/ffttn;
    temp2=abs(ffta(1))/ffttn;

    aabs=[aabs;4*temp1];
    aang=[aang;180*angle(ffta(ele))/pi];

    stata=[stata;(temp2-2*temp1)];

end

frev=[statv,vabs,vang];
frea=[stata,aabs,aang];

