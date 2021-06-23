function result = Mie_phasefn(m, x, nsteps)


    teta=linspace(0,pi,nsteps)';
    P=zeros(length(teta),1);
    for j = 1:nsteps

        u=cos(teta(j));

        a=Mie_S12(m,x,u);

        S1=norm(a(1,:));

        S2=norm(a(2,:));
        
        P(j)=S1^2+S2^2;

    end

    
%     figure
%     
%     polarplot(teta,P)


result=P; 

