function [c_true,cmin,cmax] = GetTrueVelocity(nz, value, analysis_to_make)
    if(analysis_to_make ~= 2)
        value = 0;
    end
    c_true = 2000*ones(nz, 1);
    c_true(51 + value: 80 + value) = 1600;
    c_true(121:150) = 2500;
    cmin=min(c_true);
    cmax=max(c_true);
end

