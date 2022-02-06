function rgb_demo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
clf
cols = get_cols;
cols = {cols{:,1}}'; %#ok<CCAT1>
cols = { cols{:}, ... %#ok<CCAT1>
	'k', ...
	'r', ...
	'g', ...
	'b', ...
	'y', ...
	'm', ...
	'c', ...
	'w', ...
	'', ...
	'extremely dark green', ...
	'very dark green', ...
	'dark green', ...
	'slightly dark green', ...
	'green', ...
	'slightly pale green', ...
	'pale green', ...
	'very pale green', ...
	'extremely pale green', ...
}; %#ok<CCAT>

height = 9;
x = 0;
y = 0;
for n = 1:length(cols)
	rect(x,y,cols{n})
	y = y+1;
	if y == height
		x = x+2;
		y = 0;
	end
end
if y == 0
    x = x-2;
end
axis([0 (x+2) 0 height])
title('names on different rows are alternates')

end % subfunction {rgb_demo}


function rect(x,y,col)
if isempty(col)
    return
end
r = rectangle('position',[x+0.1 y+0.1 1.8 0.8]);
col_ = col;
if iscell(col)
    col = col{1};
end
colrgb = rgb(col);
if strcmp(col(1),'u') && length(col)==2
	t = text(x+1,y+0.5,{'unnamed',['colour (' col(2) ')']});
	set(r,'facecolor',colrgb);
else
	t = text(x+1,y+0.5,col_);
	set(r,'facecolor',colrgb);
	if sum(colrgb)<1.5
        set(t,'color',[1 1 1]);
	end
end
set(t,'horizontalalignment','center')
set(t,'fontsize',10)

end % subfunction {rect}


function cols = get_cols

cols = {
	'black', [0 0 0]; ...
	'navy', [0 0 0.5]; ...
	'blue', [0 0 1]; ...
	'u1', [0 0.5 0]; ...
	{'teal','turquoise'}, [0 0.5 0.5]; ...
	'slateblue', [0 0.5 1]; ...
	{'green','lime'}, [0 1 0]; ...
	'springgreen', [0 1 0.5]; ...
	{'cyan','aqua'}, [0 1 1]; ...
	'maroon', [0.5 0 0]; ...
	'purple', [0.5 0 0.5]; ...
	'u2', [0.5 0 1]; ...
	'olive', [0.5 0.5 0]; ...
	{'gray','grey'}, [0.5 0.5 0.5]; ...
	'u3', [0.5 0.5 1]; ...
	{'mediumspringgreen','chartreuse'}, [0.5 1 0]; ...
	'u4', [0.5 1 0.5]; ...
	'sky', [0.5 1 1]; ...
	'red', [1 0 0]; ...
	'u5', [1 0 0.5]; ...
	{'magenta','fuchsia'}, [1 0 1]; ...
	'orange', [1 0.5 0]; ...
	'u6', [1 0.5 0.5]; ...
	'u7', [1 0.5 1]; ...
	'yellow', [1 1 0]; ...
	'u8', [1 1 0.5]; ...
	'white', [1 1 1]; ...
	};

for n = 1:size(cols,1)
	if ~iscell(cols{n,1})
        cols{n,1} = {cols{n,1}}; %#ok<CCAT1>
	end
end

end