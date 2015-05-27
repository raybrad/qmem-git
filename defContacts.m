function contact = defContacts(p1,p2,dcVol,acVol)

contact = struct('nodes',[],'dcVol',[],'acVol',[]);
contact.nodes = defBrickNodes(p1,p2); % return nodes index
contact.dcVol = ones(length(contact.nodes),1)*dcVol;
contact.acVol = ones(length(contact.nodes),1)*acVol;
