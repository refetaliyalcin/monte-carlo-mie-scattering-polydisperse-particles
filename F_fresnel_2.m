function theta_2=F_fresnel_2(n2,k2,theta_1)
    temp=n2^2-k2^2-sin(theta_1)^2;
    p_sq=0.5*(sqrt(temp^2+4*n2^2*k2^2)+temp);
    theta_2=atan(sin(theta_1)/sqrt(p_sq));
end
 

