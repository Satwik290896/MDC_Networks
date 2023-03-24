%function [out] = I1(Tt,r,alp_0,d,alp_1)
%    out = hyp_geo((alp_0/d),1/(Tt*(r^(alp_0/d)))) + hyp_geo((-alp_1/d),(Tt*(r^(alp_0/d)))) - (r*hyp_geo((alp_0/d),1/Tt)) + r - 1;
%end

function [out] = I1(Tt,r,alp_0,d,alp_1)
    
    % Tt Threshold > vector 
    % r : integrand 
    % r is the integrand, everything related to it should have dot operator
    % alp_0, d, alp_1 : scalars

    %out = hyp_geo((alp_0./d),1./(Tt.*(r.^(alp_0/d)))) + hyp_geo((-alp_1/d),(Tt.*(r.^(alp_0/d)))) - (r.*hyp_geo((alp_0/d),1./Tt)) + r - 1;
    
    %out = hyp_geo((alp_0/d),1/(Tt*(r^(alp_0/d)))) + hyp_geo((-alp_1/d),(Tt*(r^(alp_0/d)))) - (r*hyp_geo((alp_0/d),1/Tt)) + r - 1;

    out = hyp_geo((alp_0/d),1./(Tt*(r.^(alp_0/d)))) + hyp_geo((-alp_1/d),(Tt*(r.^(alp_0/d)))) - (r.*hyp_geo((alp_0/d),1/Tt)) + r - 1;


end