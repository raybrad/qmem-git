function label(locType,var,value)

%%%% label nodes or links with numbers %%%%%%%%%%%%%

global nodes linkCenter surfCenter;


if length(var) ~= length(value), error('numbers of locations and values must be equal'); end
switch locType
    case 'node'
        loc = nodes(var,:);
    case 'link'
        loc = linkCenter(var,:);
    case 'surf'
        loc = surfCenter(var,:);
    otherwise
        error('No such type of variable');
end
for i = 1:size(loc,1)
    lc = loc(i,:);
    text(lc(1),lc(2),lc(3),num2str(value(i)),'Fontsize',4,'Color','r');
end
    