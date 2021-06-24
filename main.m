clc
clear all
close all

%this code considers the size distribution of spherical particles

%solution of radiative transfer equation in one layer pigmented plane parallel medium
%ray is incident from air to coating. coating is coated on a substrate
%substrate could be air or other material such as silver, glass etc.
%the code estimates spectral hemispherical reflectance, transmittance and absorptance
%can handle; independent scattering, boundary reflections, absorption in
%medium. can't handle; coherent backscattering, dependent scattering and polarized ray tracing.
%while calculating the scattering direction the code uses cumulative inverse
%relation of exact single scattering pahse function or henyey greenstein 
%function approximation depending on choice 

lamda=(300:10:700)*10^-9; %freespace wavelength of incident ray in meter 
thickness=10*10^-6;  %thickness of coating in meter 

f_v=0.02; %volume fraction. 0.01 corresponds to 1% 
polar_angle=0;%linspace(0,89.99999,9); %incident angles. 0 = perpendicular to slab face. 90 parallel and should be avoided.
use_HG=0; %if 0 use exact scattering phase function, if 1 uses henyey greenstein phase function approximation


r_vector_nm=(1:200); %size range to be considered
weight_vector = lognpdf(r_vector_nm,3.8,0.5); %the distribution
r_vector=r_vector_nm*10^-9;
weight_vector=weight_vector/trapz(r_vector,weight_vector); 

% r_vector=49:1:51; %size range to be considered
% weight_vector = [0,1,0]; %the distribution
% r_vector=(49:1:51)*10^-9;
% weight_vector=weight_vector/trapz(r_vector,weight_vector); 

% figure %show the size distribution
% plot(10^9*r_vector,weight_vector,'LineWidth',2)
% ylabel('Frequency')
% xlabel('Radius [nm]')



n_pigment=sio2_n(lamda);  %real refractive index of pigment
k_pigment=sio2_k(lamda);  %imaginary refractive index of pigment
n_medium=ones(length(lamda),1); %real refractive index of medium
k_medium=zeros(length(lamda),1); %imaginary refractive index of medium
n_substrat=ones(length(lamda),1); %real refractive index of substrate
k_substrat=zeros(length(lamda),1); %imaginary refractive index of substrate

photon_number=10^5; %number of rays that will be traced, higher the number more accurate the result
n_cdf_random=1000; %how many pieces will be between (0,1) in random number relation with inverse cumulative function, higher the number accurate the phase function. useless if HG is used
nang_gid=1000; %how many pieces will be between (0,pi) in cumulative distribution function, higher the number more accurate the result. useless if HG is used


lamda_nm=lamda*10^9; %for plot
polar_angle_rad=polar_angle*pi/180;

% calculate surface reflection from air to medium. 
% medium to air and medium to substrate is calculated within snell.m.
% air to medium is calculated here seperately since we need refraction angle
teta_prime=zeros(length(lamda),length(polar_angle));
sur_reflection=zeros(length(lamda),length(polar_angle));
for i=1:length(lamda)
    for j=1:length(polar_angle_rad)
        teta_prime(i,j)=F_fresnel_2(n_medium(i),k_medium(i),polar_angle_rad(j))*180/pi;
        cos_teta=cosd(polar_angle(j));
        sin_teta=sqrt(1-cos_teta*cos_teta);
        carpan2=1/(n_medium(i)-1i*k_medium(i));
        sin_x2=sin_teta*carpan2;
        cos_x2=sqrt(1-sin_x2*sin_x2);
        carpan1=cos_teta/cos_x2;
        carpan3=cos_x2/cos_teta;
        E_parallel=(carpan1-carpan2)/(carpan1+carpan2);
        R_parallel=E_parallel*conj(E_parallel);
        E_orth=(carpan3-carpan2)/(carpan3+carpan2);
        R_orth=E_orth*conj(E_orth);
        reflectance=real(R_parallel+R_orth)*0.5;
        sur_reflection(i,j)=reflectance;
    end
end

%initialize variables
ref_lamda=zeros(length(lamda),length(polar_angle));
tra_lamda=zeros(length(lamda),length(polar_angle));
abs_lamda=zeros(length(lamda),length(polar_angle));


mu_tot_arr_sigma=zeros(length(lamda),length(r_vector));
beta_sigma=zeros(length(lamda),length(r_vector));
alfa_sigma=zeros(length(lamda),length(r_vector));
C_sca_sigma=zeros(length(lamda),length(r_vector));
C_abs_sigma=zeros(length(lamda),length(r_vector));
g_arr_sigma=zeros(length(lamda),length(r_vector));
for z=1:length(r_vector)
    Area=pi*r_vector(z)^2;
    V=(4/3)*pi*r_vector(z)^3;
    for i=1:length(lamda)
        x=2*pi*r_vector(z)*n_medium(i)/lamda(i);
        m=(n_pigment(i)+1i*k_pigment(i))/n_medium(i);
        fonksiyon=Mie(m,x);
        Qabs=fonksiyon(3);
        Qsca=fonksiyon(2);
        C_sca_sigma(i,z)=Qsca*Area;
        C_abs_sigma(i,z)=Qabs*Area;
        g_arr_sigma(i,z)=fonksiyon(5);
    end
end

r_3=trapz(r_vector,(weight_vector.*r_vector.^3)');
V_avg=(4/3)*pi*r_3;
C_sca=trapz(r_vector,(weight_vector.*C_sca_sigma)');
C_abs=trapz(r_vector,(weight_vector.*C_abs_sigma)');
ext_tot=(C_sca+C_abs)*f_v/V_avg;
g=trapz(r_vector,(weight_vector.*g_arr_sigma.*C_sca_sigma)')./trapz(r_vector,(weight_vector.*C_sca_sigma)');
scat_prob=trapz(r_vector,(weight_vector.*C_sca_sigma)')./trapz(r_vector,(weight_vector.*(C_sca_sigma+C_abs_sigma))');


tic
%loop the monte carlo code for lamda and polar_angle
for j=1:length(polar_angle)
    for i=1:length(lamda)
        [ref_lamda(i,j),tra_lamda(i,j),abs_lamda(i,j)]=monte_carlo(photon_number,sur_reflection(i,j),cosd(teta_prime(i,j)),thickness,scat_prob(i),ext_tot(i),n_medium(i),k_medium(i),n_substrat(i),k_substrat(i),g(i));
    end
end
toc
figure %draw normal to diffuse R, T and A for normal incidence (first index in my case)
plot(lamda_nm,ref_lamda(:,1),lamda_nm,tra_lamda(:,1),lamda_nm,abs_lamda(:,1),'LineWidth',2)
ylim([0 1])
xlim([min(lamda_nm) max(lamda_nm)])
legend('Reflectance','Transmittance','Absorptance','Location', 'Best')
xlabel('Wavelength [nm]')
ylabel({'Normal to Hemispherical';'Reflectance, Transmittance, Absorptance'})
