clc;clear;close all;
tic
format short;
noOfSeg = 162;        %Must be an even number
noOfBasisPoints = noOfSeg+1;
lambda = 1;
eps0 = 8.854e-12;mu0 = pi*4e-7;c0 = 1/sqrt(eps0*mu0);
freq = c0/lambda;
omega = 2*pi*freq;
wireLength = lambda/2;
a = 0.005*lambda;       %wire radius
delta = wireLength/noOfSeg;
k = 2*pi/lambda;

V = zeros(noOfSeg-1,1);
V(ceil(noOfSeg/2),1) = 1;

basisPoint_z = -wireLength/2:delta:wireLength/2;
basisPoint_y = zeros(length(basisPoint_z));
basisPoint_x = zeros(length(basisPoint_z));

s = linspace(0,delta,5);delta_s = s(2)-s(1);
factor = -1/(1j*4*pi*omega*eps0);
basisPoint_n_minus = zeros;basisPoint_m_plus = zeros;basisPoint_m_minus = zeros;Z = zeros;
for m = 2:length(basisPoint_z)-1
    for n = 2:length(basisPoint_z)-1
        basisPoint_n_minus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
        F1 = exp(-1j*k*sqrt((basisPoint_z(m)-basisPoint_n_minus(n)-s).^2+ a^2))./...
            sqrt((basisPoint_z(m)-basisPoint_n_minus(n)-s).^2+ a^2);
        psi_m_n = simp(F1.',delta_s);
        
        basisPoint_m_plus(m) = (basisPoint_z(m+1)+basisPoint_z(m))/2;
        F2 = exp(-1j*k*sqrt((basisPoint_m_plus(m)-basisPoint_z(n)-s).^2+ a^2))./...
            sqrt((basisPoint_m_plus(m)-basisPoint_z(n)-s).^2+ a^2);
        psi_mPlus_nPlus = simp(F2.',delta_s);
        
        F3 = exp(-1j*k*sqrt((basisPoint_m_plus(m)-basisPoint_z(n-1)-s).^2+ a^2))./...
            sqrt((basisPoint_m_plus(m)-basisPoint_z(n-1)-s).^2+ a^2);
        psi_mPlus_nMinus = simp(F3.',delta_s);
        
        basisPoint_m_minus(m) = (basisPoint_z(m)+basisPoint_z(m-1))/2;
        F4 = exp(-1j*k*sqrt((basisPoint_m_minus(m)-basisPoint_z(n)-s).^2+ a^2))./...
            sqrt((basisPoint_m_minus(m)-basisPoint_z(n)-s).^2+ a^2);
        psi_mMinus_nPlus = simp(F4.',delta_s);
        
        F5 = exp(-1j*k*sqrt((basisPoint_m_minus(m)-basisPoint_z(n-1)-s).^2+ a^2))./...
            sqrt((basisPoint_m_minus(m)-basisPoint_z(n-1)-s).^2+ a^2);
        psi_mMinus_nMinus = simp(F5.',delta_s);
        
        vectorPot = k^2*delta*psi_m_n;
        
        Z(m,n) = factor*(vectorPot-(psi_mPlus_nPlus/delta)-(psi_mMinus_nMinus/delta)+...
            (psi_mMinus_nPlus/delta)+(psi_mPlus_nMinus/delta));
    end
end
Z(1,:) = [];Z(:,1) = [];

I = Z\V;
I1 = I.'*1e3;

fig1 = figure('Color','w');
yyaxis left
plot(basisPoint_z,[0 abs(I1) 0],'LineWidth',1.4);
ylabel('Magnitude of Current (mA)')
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

z_vector = 0.001*lambda:0.001:0.3*lambda;
rho_vector = [0.01 0.02 0.05]*lambda;
tempArray = linspace(-delta/2,delta/2,5);delta_ta = tempArray(2)-tempArray(1);

I = I.';
factor2 = 1/(1j*4*pi*omega*eps0);
for o = 1:length(rho_vector)
    for m = 1:length(z_vector)
        for n = 2:length(basisPoint_z)-1
            distance = sqrt(a^2 + rho_vector(o)^2 + (z_vector(m)-basisPoint_z(n)-tempArray).^2);
            distance1 = sqrt(a^2 + rho_vector(o)^2 + (z_vector(m)-basisPoint_z(n)).^2);
%             distance1 = sqrt(a^2 + rho_vector(o)^2 + (z_vector(m)-basisPoint_z(n)-(linspace(-delta,0,5))).^2);
%             distance2 = sqrt(a^2 + rho_vector(o)^2 + (z_vector(m)-basisPoint_z(n)-(linspace(0,delta,5))).^2);
            firstTerm(n) = k*k*I(n-1)*simp((exp(-1j*k*distance)./distance).',delta_ta);
            secondTerm(n) = I(n-1)*(z_vector(m)-basisPoint_z(n)).*((1./(distance1.^5)) + ((1j*k)./(distance1.^4))).*exp(-1j*k*distance1);
%             secondTerm(n) = secondTerm_part1(n)+secondTerm_part2(n);
        end
        firstTerm1 = sum(firstTerm);secondTerm1 = sum(secondTerm);
        E_z(o,m) = factor2*(firstTerm1+secondTerm1)/1.414;
    end
end
fig2 = figure('Color','w');
for o = 1:length(rho_vector)
plot(z_vector, abs(E_z(o,:)),'LineWidth',1.4);hold on;
end



function q=simp(y,dx)
N=length(y);

mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;

q=(dx/3)*sum(y.'.*mul);
end
