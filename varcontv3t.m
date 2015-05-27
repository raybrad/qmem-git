function varcontv = varcontv(dc1,ac1,dc2,ac2,dc3,ac3)

global contacts;
global dcVolDirNodes acVolDirNodes;

contacts{1}.dcVol=ones(length(contacts{1}.nodes),1)*dc1;
contacts{2}.dcVol=ones(length(contacts{2}.nodes),1)*dc2;
contacts{1}.acVol=ones(length(contacts{1}.nodes),1)*ac1;
contacts{2}.acVol=ones(length(contacts{2}.nodes),1)*ac2;

contacts{3}.dcVol=ones(length(contacts{3}.nodes),1)*dc3;
contacts{3}.acVol=ones(length(contacts{3}.nodes),1)*ac3;

dcVolDirNodes = [contacts{1}.dcVol;contacts{2}.dcVol;contacts{3}.dcVol];
acVolDirNodes = [contacts{1}.acVol;contacts{2}.acVol;contacts{3}.acVol];
