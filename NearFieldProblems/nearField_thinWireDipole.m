clc;clear;close all;
%% Dipole properties
lambda = 4;
wireLength = lambda/2;
a = 0.005*lambda;   %wire radius
%% Constants definitions
eps0 = 8.854e-12;
mu0 = pi*4e-7;
c0 = 1/sqrt(eps0*mu0);
freq = c0/lambda;
omega = 2*pi*freq;
k = 2*pi/lambda;
%% Segmenting the dipole
noOfSeg = 52;
delta = wireLength/noOfSeg;
basisPoint_z = -wireLength/2:delta:wireLength/2;
%% Voltage source at the center of the dipole
V = zeros(noOfSeg-1,1);
V(ceil(noOfSeg/2),1) = 1;
%% Finding Impedance matrix and current
s = linspace(0,delta,5);delta_s = s(2)-s(1);
factor1 = 1j*omega*mu0*delta/(4*pi);
factor2 = 1/(1j*omega*eps0*4*pi);
basisPoint_n_minus = zeros;basisPoint_m_plus = zeros;basisPoint_m_minus = zeros;Z = zeros;

for m = 2:length(basisPoint_z)-1
    for n = 2:length(basisPoint_z)-1

        basisPoint_n_minus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
        F1 = exp(-1j*k*sqrt((basisPoint_z(m)-basisPoint_n_minus(n)-s).^2+ a^2))./...
            sqrt((basisPoint_z(m)-basisPoint_n_minus(n)-s).^2+ a^2);
        psi_m_n = (1/delta)*simp(F1.',delta_s);
        
        basisPoint_m_plus(m) = (basisPoint_z(m+1)+basisPoint_z(m))/2;
        F2 = exp(-1j*k*sqrt((basisPoint_m_plus(m)-basisPoint_z(n)-s).^2+ a^2))./...
            sqrt((basisPoint_m_plus(m)-basisPoint_z(n)-s).^2+ a^2);
        psi_mPlus_nPlus = (1/delta)*simp(F2.',delta_s);
        
        F3 = exp(-1j*k*sqrt((basisPoint_m_plus(m)-basisPoint_z(n-1)-s).^2+ a^2))./...
            sqrt((basisPoint_m_plus(m)-basisPoint_z(n-1)-s).^2+ a^2);
        psi_mPlus_nMinus = (1/delta)*simp(F3.',delta_s);
        
        basisPoint_m_minus(m) = (basisPoint_z(m)+basisPoint_z(m-1))/2;
        F4 = exp(-1j*k*sqrt((basisPoint_m_minus(m)-basisPoint_z(n)-s).^2+ a^2))./...
            sqrt((basisPoint_m_minus(m)-basisPoint_z(n)-s).^2+ a^2);
        psi_mMinus_nPlus = (1/delta)*simp(F4.',delta_s);
        
        F5 = exp(-1j*k*sqrt((basisPoint_m_minus(m)-basisPoint_z(n-1)-s).^2+ a^2))./...
            sqrt((basisPoint_m_minus(m)-basisPoint_z(n-1)-s).^2+ a^2);
        psi_mMinus_nMinus = (1/delta)*simp(F5.',delta_s);
        
        Z(m,n) = factor1*delta*psi_m_n + factor2*(psi_mPlus_nPlus + psi_mMinus_nMinus - psi_mMinus_nPlus - psi_mPlus_nMinus);
    end
end
Z(1,:) = [];Z(:,1) = [];
I = Z\V;

I1 = I.'*1e3;
%% Plotting the current along wire length
fig1 = figure('Color','w');
yyaxis left
plot(basisPoint_z,[0 abs(I1) 0],'LineWidth',1.4);
ylabel('Magnitude of Current (mA)')
ylim([0 12])
hold on;
yyaxis right
plot(basisPoint_z,[rad2deg(angle(I1(1))) rad2deg(angle(I1)) rad2deg(angle(I1(noOfSeg-1)))],'--','LineWidth',1.4);
ylim([-180 180])
ylabel('Phase of Current (Degree)')
hold on;
xlim([-wireLength/2 wireLength/2])
xlabel('Wire Length (m)')
legend('Magnitude','Phase');
hold off;
%% Near-field along z and rho directions
z_vector = 0:0.001*lambda:0.3*lambda;
rho_vector = [0.01 0.02 0.03 0.05 0.1]*lambda;
s1 = linspace(-delta,delta,5);delta_s1 = s1(2)-s1(1);
triFun = 1-abs(s1)/delta;
vectorPot = zeros;scalarPot_z = zeros;E_z = zeros;E_y = zeros;scalarPot_y=zeros;
for o = 1:length(rho_vector)
    for m = 1:length(z_vector)
        for n = 2:length(basisPoint_z)-1
            distance_z = sqrt(a^2 + rho_vector(o)^2 + (basisPoint_z(n)-z_vector(m)+s1).^2);
            distance1_z = sqrt(a^2 + rho_vector(o)^2 + (basisPoint_z(n-1)-z_vector(m)+s).^2);
            distance2_z = sqrt(a^2 + rho_vector(o)^2 + (basisPoint_z(n)-z_vector(m)+s).^2);
            distance1_rho = sqrt(a^2 + rho_vector(o)^2 + (basisPoint_z(n-1)-z_vector(m)+s).^2);
            distance2_rho = sqrt(a^2 + rho_vector(o)^2 + (basisPoint_z(n)-z_vector(m)+s).^2);
            y_term1 = (rho_vector(o))./distance1_rho;
            y_term2 = (rho_vector(o))./distance2_rho;
            A2 = (1/(4*pi))*exp(-1j*k*distance_z)./distance_z;
            A2_TnIn = A2.*triFun*I(n-1,1);
            vectorPot(n) = (-1j*omega*mu0)*simp(A2_TnIn.',delta_s1);
            z_term1 = (basisPoint_z(n-1)-z_vector(m)+s)./distance1_z;
            z_term2 = (basisPoint_z(n)-z_vector(m)+s)./distance2_z;
            Phi2_1_z = (1/delta)*(exp(-1j*k*distance1_z)./distance1_z).*(1j*k + 1./distance1_z).*z_term1;
            Phi2_2_z = -(1/delta)*(exp(-1j*k*distance2_z)./distance2_z).*(1j*k + 1./distance2_z).*z_term2;
            Phi2_z = (1/(1j*omega*eps0*4*pi))*(Phi2_1_z + Phi2_2_z);
            Phi2_In_z = I(n-1,1)*Phi2_z;
            scalarPot_z(n) = simp(Phi2_In_z.',delta_s);
            Phi2_1_y = (1/delta)*(exp(-1j*k*distance1_rho)./distance1_rho).*(1j*k + 1./distance1_rho).*y_term1;
            Phi2_2_y = -(1/delta)*(exp(-1j*k*distance2_rho)./distance2_rho).*(1j*k + 1./distance2_rho).*y_term2;
            Phi2_y = (1/(1j*omega*eps0*4*pi))*(Phi2_1_y + Phi2_2_y);
            Phi2_In_y = I(n-1,1)*Phi2_y;
            scalarPot_y(n) = simp(Phi2_In_y.',delta_s);
        end
        E_z(o,m) = (lambda/2)*(sum(vectorPot) + sum(scalarPot_z))/(1.414);
        E_y(o,m) = (lambda/2)*sum(scalarPot_y)/(1.414);
        E_rho(o,m) = E_y(o,m)*1.414;
    end
end

fig2 = figure('Color','w');
for p = 2:length(rho_vector)
    plot(z_vector/lambda,abs(E_z(p,:)),'LineWidth',1.4);hold on;
end
fig3 = figure('Color','w');
for p = 1:length(rho_vector)
    plot(z_vector/lambda,abs(E_rho(p,:)),'LineWidth',1.4);hold on;
end
save nearField_thinWireDipole   I ...
                                E_z ...
                                E_rho
%% Simpsons function of numerical integration
function q=simp(y,dx)
N=length(y);

mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;

q=(dx/3)*sum(y.'.*mul);
end
