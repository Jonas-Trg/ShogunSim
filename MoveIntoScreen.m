function pos = MoveIntoScreen(pos)
% Get monitor positions: each row = [left bottom width height]
monitors = get(0, 'MonitorPositions');
% add title bar height
titlebarheight = 40;
pos(3) = pos(3) + titlebarheight;
% Choose the monitor where most of the figure currently lies
overlapArea = zeros(size(monitors,1),1);
for i = 1:size(monitors, 1)
    mon = monitors(i, :);
    % Compute overlap rectangle
    overlapWidth  = max(0, min(pos(1) + pos(3), mon(1) + mon(3)) - max(pos(1), mon(1)));
    overlapHeight = max(0, min(pos(2) + pos(4), mon(2) + mon(4)) - max(pos(2), mon(2)));
    overlapArea(i) = overlapWidth * overlapHeight;
end
[~, bestMonIdx] = max(overlapArea);
mon = monitors(bestMonIdx,:);

% Clamp figure coordinates to stay inside chosen monitor
newLeft   = max(mon(1), min(pos(1), mon(1)+mon(3)-pos(3)));
newBottom = max(mon(2), min(pos(2), mon(2)+mon(4)-pos(4)));

% Apply new position
pos = [newLeft newBottom pos(3) pos(4)];
pos(3) = pos(3) - titlebarheight;
end % MoveIntoScreen
