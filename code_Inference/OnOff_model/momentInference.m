function theta = momentInference(vals)

theta = [];
num_samples = length(vals);
m1 = mean(vals);
m2 = sum(vals.*(vals - 1))/num_samples;
m3 = sum(vals.*(vals - 1).*(vals - 2))/num_samples;

if sum(vals) == 0 || m1 == 0 || m2 == 0
    return
else
    r1 = m1;
    r2 = m2/m1;
    r3 = m3/m2;
    if (r1*r2 - 2*r1*r3 + r2*r3) == 0 || ((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3)) == 0 || (r1 - 2*r2 + r3) == 0
        return
    else
        kon_est = (2*r1*(r3-r2))/(r1*r2-2*r1*r3 + r2*r3);
        koff_est = (2*(r3-r2)*(r1-r3)*(r2-r1))/((r1*r2 - 2*r1*r3 + r2*r3)*(r1-2*r2+r3));
        ksyn_est = (2*r1*r3 - r1*r2 - r2*r3)/(r1 - 2*r2 + r3);
        theta = [kon_est, koff_est, ksyn_est];
    end
end
end
