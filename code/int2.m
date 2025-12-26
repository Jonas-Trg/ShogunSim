function z_out = int2(X,Y,Z,xIn,yIn)
% int2
%
if size(xIn)>1
    z_out = interp2(X,Y,Z,xIn,yIn);
elseif isempty(xIn)
    z_out = [];
else
    xb=find(abs(xIn)>=abs(X(1,1:end-1)),1,'last');
    yb=find(yIn>=Y(1:end-1,1),1,'last');
    if ~isempty(xb) || ~isempty(yb)
        xp=(xIn-X(1,xb))/(X(1,xb+1)-X(1,xb));
        yp=(yIn-Y(yb,1))/(Y(yb+1,1)-Y(yb,1));
        z1=Z(yb,  xb)+xp*(Z(yb,  xb+1)-Z(yb,  xb));
        z2=Z(yb+1,xb)+xp*(Z(yb+1,xb+1)-Z(yb+1,xb));
        z_out=z1+yp*(z2-z1);
    else
        z_out=[];
    end
end
return
