function [c] = GetInitialVelocity(c_true, value, analysis_to_make)
    nz = length(c_true);
    if(analysis_to_make==0)
        c = value * ones(nz, 1);
    elseif(analysis_to_make==3)
        c = smooth(c_true, value);
    else
        c = 2000*ones(nz,1);
    end
end

