function contact = defMetalContacts(tmpContact,p1,p2,dcVol,acVol)

contact = struct('nodes',[],'dcVol',[],'acVol',[]);
contact.nodes = tmpContact.nodes;
contact.dcVol = tmpContact.dcVol;
contact.acVol = tmpContact.acVol;
tmpnodes = defBrickNodes(p1,p2); % return nodes index
[~, ind] = ismember(tmpnodes,contact.nodes);

contact.dcVol(ind)=dcVol;
contact.acVol(ind)=acVol;
