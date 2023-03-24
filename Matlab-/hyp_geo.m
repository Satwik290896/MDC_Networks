function [C] = hyp_geo(b,z)
    C = hypergeom([1,1/b],1+1/b,-z);
end