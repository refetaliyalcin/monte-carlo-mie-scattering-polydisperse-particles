function [r_tot,t_tot,a_tot]=monte_carlo(photon_number,s_ref,cos_gecen,h,scat_prob,mu_tot,n_medium,k_medium,n_subs,k_subs,g)
r_tot=0;
t_tot=0;
a_tot=0;
for i=1:photon_number
    r_no=0;
    t_no=0;
    a_no=0;
    Rastgele=rand();
    if Rastgele>s_ref
        %ray penetrates to the medium
        alive=1;
    else
        %specular reflection
        alive=0;
        r_no=1;
    end
    x=0;%x component of position vector
    y=0;%y component of position vector
    z=0;%z component of position vector
    s_x=0;%x component of direction vector
    s_y=sqrt(1-cos_gecen*cos_gecen);%y component of direction vector. could be defined better along with s_x but it is a 1D code anyway, nothing will change
    s_z=cos_gecen; %z component of direction vector 
    l_beta=-log(rand())/mu_tot; %excitation length
    while alive   
        if (s_z>0) %going down to substrate
            l_w = (h - z)/s_z; %distance to lower boundary
        else %going up to air
            l_w = -z/s_z; %distance to upper boundary
        end
        if l_w<l_beta %check if ray reaches to boundary or extinct?
            min_index=1; %reach boundary
            min_l=l_w;
        else
            min_index=2; %extinct (will absorbed or scatter? we will check below)
            min_l=l_beta;
        end
        x=x+min_l*s_x;% move
        y=y+min_l*s_y;% move
        z=z+min_l*s_z;% move
        if (min_index==1)
    %            disp('hit boundary');
            alive=snell(s_z,n_medium,k_medium,n_subs,k_subs);% check if the ray can leave the medium or not
            if (alive==0)
                if s_z>0
                    t_no=1;%it left from bottom so transmitted
                else
                    r_no=1;%it left from top so reflected
                end
            else
                l_beta=l_beta-l_w;%ray did not extinct, so don't behave like it starts as a new ray and reduce l_beta by l_w
                s_z=-s_z; %specularly reflected and direction changed
            end
        else
            random_no=rand();
            if random_no<scat_prob
    %               disp('scattering');
                    [s_x,s_y,s_z]=scatter_hg(g,s_x,s_y,s_z);%find new trajectory by henyey greenstein phase function

                l_beta=-log(rand())/mu_tot;%it extincted so don't keep the old l_beta create new for new event
            else
    %               disp('absorption');
                alive=0;%it is aborbed, game over for this ray. exit this loop
                a_no=1;
            end
        end
    end
    r_tot=r_tot+r_no;
    t_tot=t_tot+t_no;
    a_tot=a_tot+a_no;
end
r_tot=r_tot/photon_number;
t_tot=t_tot/photon_number;
a_tot=a_tot/photon_number;