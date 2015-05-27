function contact = defContacts(contnodes,dcVol,acVol)

contact = struct('nodes',[],'dcVol',[],'acVol',[]);
contact.nodes = contnodes;
contact.dcVol = ones(length(contact.nodes),1)*dcVol;
contact.acVol = ones(length(contact.nodes),1)*acVol;
