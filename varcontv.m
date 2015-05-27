function varcontv = varcontv(dc1,ac1,dc2,ac2)

global contacts;
global dcVolDirNodes acVolDirNodes;

contacts{1}.dcVol=ones(length(contacts{1}.nodes),1)*dc1;
contacts{2}.dcVol=ones(length(contacts{2}.nodes),1)*dc2;
contacts{1}.acVol=ones(length(contacts{1}.nodes),1)*ac1;
contacts{2}.acVol=ones(length(contacts{2}.nodes),1)*ac2;

dcVolDirNodes = [contacts{1}.dcVol;contacts{2}.dcVol];
acVolDirNodes = [contacts{1}.acVol;contacts{2}.acVol];
