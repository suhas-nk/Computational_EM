clc;clear;close all;
format long;
noOfSeg = 24;        %Must be an even number
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

basisPoint = -wireLength/2:delta:wireLength/2;

factor = -1/(1j*4*pi*omega*eps0);
basisPoint_n_minus = zeros;basisPoint_m_plus = zeros;basisPoint_m_minus = zeros;Z = zeros;
for m = 2:length(basisPoint)-1
    for n = 2:length(basisPoint)-1
        basisPoint_n_minus(n) = (basisPoint(n)+basisPoint(n-1))/2;
        F1 = @(s) exp(-1j*k*sqrt((basisPoint(m)-basisPoint_n_minus(n)-s).^2+ a^2))./...
            sqrt((basisPoint(m)-basisPoint_n_minus(n)-s).^2+ a^2);
        psi_m_n = integral(F1,0,delta);
        
        basisPoint_m_plus(m) = (basisPoint(m+1)+basisPoint(m))/2;
        F2 = @(s) exp(-1j*k*sqrt((basisPoint_m_plus(m)-basisPoint(n)-s).^2+ a^2))./...
            sqrt((basisPoint_m_plus(m)-basisPoint(n)-s).^2+ a^2);
        psi_mPlus_nPlus = integral(F2,0,delta);
        
        F3 = @(s) exp(-1j*k*sqrt((basisPoint_m_plus(m)-basisPoint(n-1)-s).^2+ a^2))./...
            sqrt((basisPoint_m_plus(m)-basisPoint(n-1)-s).^2+ a^2);
        psi_mPlus_nMinus = integral(F3,0,delta);
        
        basisPoint_m_minus(m) = (basisPoint(m)+basisPoint(m-1))/2;
        F4 = @(s) exp(-1j*k*sqrt((basisPoint_m_minus(m)-basisPoint(n)-s).^2+ a^2))./...
            sqrt((basisPoint_m_minus(m)-basisPoint(n)-s).^2+ a^2);
        psi_mMinus_nPlus = integral(F4,0,delta);
        
        F5 = @(s) exp(-1j*k*sqrt((basisPoint_m_minus(m)-basisPoint(n-1)-s).^2+ a^2))./...
            sqrt((basisPoint_m_minus(m)-basisPoint(n-1)-s).^2+ a^2);
        psi_mMinus_nMinus = integral(F5,0,delta);
        
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
plot(basisPoint,[0 abs(I1) 0],'LineWidth',1.4);
ylabel('Magnitude of Current (mA)')
hold on;
yyaxis right
plot(basisPoint,[rad2deg(angle(I1(1))) rad2deg(angle(I1)) rad2deg(angle(I1(noOfSeg-1)))],'--','LineWidth',1.4);
ylim([-180 180])
ylabel('Phase of Current (Degree)')
hold on;
xlim([-wireLength/2 wireLength/2])
xlabel('Wire Length (m)')
legend('Magnitude','Phase');
hold off;

z_vector = 0.001*lambda:0.001:0.3*lambda;
testingDipole_halfLength = 0.001*lambda/2;
rho_vector = [0.02 0.03 0.05 0.1]*lambda;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector)
        
        basisPoint(noOfBasisPoints+1) = z_vector(p)-testingDipole_halfLength;
        basisPoint(noOfBasisPoints+2) = z_vector(p);
        basisPoint(noOfBasisPoints+3) = z_vector(p)+testingDipole_halfLength;
        testingDipole_lowerMidPoint = (basisPoint(noOfBasisPoints+1)+basisPoint(noOfBasisPoints+2))/2;
        testingDipole_upperMidPoint = (basisPoint(noOfBasisPoints+2)+basisPoint(noOfBasisPoints+3))/2;
        diffz = testingDipole_upperMidPoint - testingDipole_lowerMidPoint;
        
        esum = 0;
        for n = 2:noOfBasisPoints-1
            basisPoint_n_minus(n) = (basisPoint(n)+basisPoint(n-1))/2;
            F1 = @(s) exp(-1j*k*sqrt((basisPoint(noOfBasisPoints+2)-basisPoint_n_minus(n)-s).^2 + rho_vector(q)^2 + a^2))./...
                sqrt((basisPoint(noOfBasisPoints+2)-basisPoint_n_minus(n)-s).^2 + rho_vector(q)^2 + a^2);
            psi_m_n = integral(F1,0,delta);
            
            F2 = @(s) exp(-1j*k*sqrt((testingDipole_upperMidPoint-basisPoint(n)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_upperMidPoint-basisPoint(n)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mPlus_nPlus = integral(F2,0,delta);
            
            F3 = @(s) exp(-1j*k*sqrt((testingDipole_upperMidPoint-basisPoint(n-1)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_upperMidPoint-basisPoint(n-1)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mPlus_nMinus = integral(F3,0,delta);
            
            F4 = @(s) exp(-1j*k*sqrt((testingDipole_lowerMidPoint-basisPoint(n)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_lowerMidPoint-basisPoint(n)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mMinus_nPlus = integral(F4,0,delta);
            
            F5 = @(s) exp(-1j*k*sqrt((testingDipole_lowerMidPoint-basisPoint(n-1)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_lowerMidPoint-basisPoint(n-1)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mMinus_nMinus = integral(F5,0,delta);
            
            vectorPot = k^2*diffz*psi_m_n;
            matrix = factor*(vectorPot-(psi_mPlus_nPlus/delta)-(psi_mMinus_nMinus/delta)+...
                (psi_mMinus_nPlus/delta)+(psi_mPlus_nMinus/delta));
            esum = esum+matrix*I(n-1);
        end
        efield(q,p) = abs(-esum/0.001/1.414);
    end
end

fig2 = figure('Color','w');
for p = 1:length(rho_vector)
    plot(z_vector,efield(p,:),'LineWidth',1.4);hold on;
end




























