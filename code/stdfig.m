function r=stdfig(type, figno, name, sb, offset,screen)
% stdfig
% This function ease the creation of figures
%   r=stdfig(type, figno, name, sb)
%
%   Input
%       type - this defines the dimension of the figure itself. it can
%         either be the size of the plotting area in pixels e.g. [350 200]
%         for 350 width and 200 height, or in part of the screen [1/2 1/3]
%         for a figure that takes half height and a third of the screen 
%         wideness.
%         
%       figno - gives the number of the figure. If the figure already
%         exists, it will reshape this figure according to the inputs.
%       name - gives the title to be displayed in the window bar.
%       sb - defines the axes of the subplots. If "1" only one axes, 
%         otherwise the value shall be similar to subplot command: 
%           sb=[row, cols, selection]
%       offset - Defines the offset of the figure. E.g. 2 x 1 containing 
%         the offset in number of figures (e.g. [1 0] is plotted to the 
%         right of [0 0] side by side.
%       screen - contains the coordinates for the screen you want to use 
%         for the plots. If omitted the primary screen is assumed.
%
%   Output
%       r - if the figure is created, r will return the figure number.
%       Otherwise r will return the created axes handle.
%
scrsz = get(0,'ScreenSize');
if nargin==0
    type=1;
end
if nargin<4
    sb=0;
end
if nargin<5
    offset=[0 0];
end
if nargin>=6
    scrsz=screen;
end
if size(type,2)==2
    if any(type<=1)
        w=(scrsz(3)-5)*type(2)-5;
        h=(scrsz(4)-45)*type(1)-85;
    else
        w=80*type(2);
        h=80*type(1);
    end
else
    w=560;
    h=420;
end
if nargin<2
    r=figure('Position',[5+offset(2)*(w+5)+scrsz(1) scrsz(2)+scrsz(4)-h-85-offset(1)*(h+85) w h]);
elseif figno
    if ishandle(figno)
        r=figno;
        set(r,'Position',[5+offset(2)*(w+5)+scrsz(1) scrsz(2)+scrsz(4)-h-85-offset(1)*(h+85) w h]);
        if size(sb,2)==1
            b=findobj(r,'Type','Axes');
            for ii=1:size(b,1)
                x=get(b(ii),'Position');
                if x(1)<0.5
                    x(1)=64/w;
                else
                    x(1)=64/w+0.5;
                end
                if x(2)<0.5
                    x(2)=48/h;
                else
                    x(2)=48/h+0.5;
                end
                if x(3)>0.5
                    x(3)=1-96/w;
                else
                    x(3)=0.5-96/w;
                end
                if x(4)>0.5
                    x(4)=1-80/h;
                else
                    x(4)=0.5-80/h;
                end
                set(b(ii),'Position',x);
            end
            %         if size(b,1)==1
            %             set(b,'Position',[64/w 48/h 1-96/w 1-80/h])
            %         end
        end
        figure(r)
    elseif figno>0
        r=figure(figno);
        set(r,'Position',[5+offset(2)*(w+5)+scrsz(1) scrsz(2)+scrsz(4)-h-85-offset(1)*(h+85) w h]);
    else
        r=figure;
        set(r,'Position',[5+offset(2)*(w+5)+scrsz(1) scrsz(2)+scrsz(4)-h-85-offset(1)*(h+85) w h]);
    end
else
    r=figure;
    set(r,'Position',[5+offset(2)*(w+5)+scrsz(1) scrsz(2)+scrsz(4)-h-85-offset(1)*(h+85) w h]);
end
if sb
    if size(sb,2)==3 %rows,cols,sel
        sel=sb(3);
        nrows=sb(1);
        ncols=sb(2);
        col=mod(sel,ncols);
        row=nrows-floor((sel-1)/ncols)-1;
        r=subplot('Position',[64/w+col/ncols 48/h+row/nrows 1/ncols-96/w 1/nrows-80/h]);
%         r=subplot('Position',[64/w+col/ncols 48/h+row/nrows (col+1)/ncols-96/w (row+1)/nrows-80/h]);
    else
        subplot('Position',[64/w 48/h 1-96/w 1-80/h])
    end
end
if exist('name','var') && ~isempty(name)
    set(gcf,'Name',name)
end
return
