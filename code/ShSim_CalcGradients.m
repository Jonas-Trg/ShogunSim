function slope = ShSim_CalcGradients(vehicle, slopedata, do_echo)
% ShSim_calc_slope
% 0.1.1.0
%
if do_echo
    fprintf('Calculating slope matrix\n');
end
vLen = round(vehicle.unitLength * vehicle.nUnits);
newslope = round([slopedata(:,1)-vLen;slopedata(:,1)]);
newslope(end,3)=0;
newslope=sortrows(newslope, 1);
for ii=1:size(newslope,1)
    rear=newslope(ii,1);
    front=rear+vLen;
    pf=find(front>=slopedata(:,1));
    pf=pf(end);
    pr=find(rear<=slopedata(:,1));
    if isempty(pr)
        pr=1;
    else
        pr=pr(1);
    end
    ll=0;
    for jj=pr:pf-1 %alla hel sectioner
        ll=ll+(slopedata(jj+1,1)-slopedata(jj,1))*slopedata(jj,2);
    end
    if rear>0 & pr>1 %#ok<AND2> %slope on track below 0 is always zero!
        ll=ll+(slopedata(pr,1)-rear)*slopedata(pr-1,2);
    end
    ll=ll+(front-slopedata(pf,1))*slopedata(pf,2);
    slope(ii,:)=[front ll/vLen 0]; %#ok<AGROW>
end
for ii=1:size(slope,1)-1
    slope(ii,3)=(slope(ii+1,2)-slope(ii,2))/(slope(ii+1,1)-slope(ii,1));
end
%ShSim_calc_slope
