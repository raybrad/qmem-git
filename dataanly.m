function nnnttt = dataanly(nsteps)
global nnnttt;
nnnttt=[];
for tttt=1:nsteps
   if mod(tttt,10)==1
    nnnttt=[nnnttt;tttt];
   end
end
