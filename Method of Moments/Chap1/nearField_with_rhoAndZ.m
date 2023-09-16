clc;clear;close all;
tic
format short;
noOfSeg = 32;        %Must be an even number
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
testingDipole_halfLength = 0.001*lambda/2;
rho_vector = [0.02 0.03 0.05 0.1]*lambda;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector)
        
        basisPoint_z(noOfBasisPoints+1) = z_vector(p)-testingDipole_halfLength;
        basisPoint_z(noOfBasisPoints+2) = z_vector(p);
        basisPoint_z(noOfBasisPoints+3) = z_vector(p)+testingDipole_halfLength;
        testingDipole_lowerMidPoint = (basisPoint_z(noOfBasisPoints+1)+basisPoint_z(noOfBasisPoints+2))/2;
        testingDipole_upperMidPoint = (basisPoint_z(noOfBasisPoints+2)+basisPoint_z(noOfBasisPoints+3))/2;
        diffz = testingDipole_upperMidPoint - testingDipole_lowerMidPoint;
        
        esum = 0;
        for n = 2:noOfBasisPoints-1
            basisPoint_n_minus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
            F1 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_n_minus(n)-s).^2 + rho_vector(q)^2 + a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_n_minus(n)-s).^2 + rho_vector(q)^2 + a^2);
            psi_m_n = simp(F1.',delta_s);
            
            F2 = exp(-1j*k*sqrt((testingDipole_upperMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_upperMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mPlus_nPlus = simp(F2.',delta_s);
            
            F3 = exp(-1j*k*sqrt((testingDipole_upperMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_upperMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mPlus_nMinus = simp(F3.',delta_s);
            
            F4 = exp(-1j*k*sqrt((testingDipole_lowerMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_lowerMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mMinus_nPlus = simp(F4.',delta_s);
            
            F5 = exp(-1j*k*sqrt((testingDipole_lowerMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2+ a^2))./...
                sqrt((testingDipole_lowerMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2+ a^2);
            psi_mMinus_nMinus = simp(F5.',delta_s);
            
            vectorPot = k^2*diffz*psi_m_n;
            matrix = factor*(vectorPot-(psi_mPlus_nPlus/delta)-(psi_mMinus_nMinus/delta)+...
                (psi_mMinus_nPlus/delta)+(psi_mPlus_nMinus/delta));
            esum = esum+matrix*I(n-1);
        end
        efield_Z(q,p) = abs(-esum/0.001/1.414);
    end
end

fig2 = figure('Color','w');
for p = 1:length(rho_vector)
    plot(z_vector,efield_Z(p,:),'LineWidth',1.4);hold on;
end

rho_vector = [0.01 0.02 0.05]*lambda;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector)
        
        basisPoint_z(noOfBasisPoints+1) = z_vector(p);
        basisPoint_z(noOfBasisPoints+2) = z_vector(p);
        basisPoint_z(noOfBasisPoints+3) = z_vector(p);
        testingDipole_lowerPoint_y = rho_vector(q)-testingDipole_halfLength;
        testingDipole_upperPoint_y = rho_vector(q)+testingDipole_halfLength;
        testingDipole_mainPoint_y = rho_vector(q);
        testingDipole_lowerMidPoint = (testingDipole_mainPoint_y+testingDipole_lowerPoint_y)/2;
        testingDipole_upperMidPoint = (testingDipole_mainPoint_y+testingDipole_upperPoint_y)/2;
        diffy = testingDipole_upperMidPoint - testingDipole_lowerMidPoint;
        
        esum = 0;
        for n = 2:noOfBasisPoints-1
            basisPoint_n_minus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
            
            F2 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_upperMidPoint^2+ a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_upperMidPoint^2+ a^2);
            psi_mPlus_nPlus = simp(F2.',delta_s);
            
            F3 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_upperMidPoint^2+ a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_upperMidPoint^2+ a^2);
            psi_mPlus_nMinus = simp(F3.',delta_s);
            
            F4 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_lowerMidPoint^2+ a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_lowerMidPoint^2+ a^2);
            psi_mMinus_nPlus = simp(F4.',delta_s);
            
            F5 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_lowerMidPoint^2+ a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_lowerMidPoint^2+ a^2);
            psi_mMinus_nMinus = simp(F5.',delta_s);
            
            matrix = factor*(-(psi_mPlus_nPlus/delta)-(psi_mMinus_nMinus/delta)+...
                (psi_mMinus_nPlus/delta)+(psi_mPlus_nMinus/delta));
            esum = esum+matrix*I(n-1);
        end
        efield_Y(q,p) = abs(-esum/0.001/1.414);
    end
end

for q = 1:length(rho_vector)
    for p = 1:length(z_vector)
        
        basisPoint_z(noOfBasisPoints+1) = z_vector(p);
        basisPoint_z(noOfBasisPoints+2) = z_vector(p);
        basisPoint_z(noOfBasisPoints+3) = z_vector(p);
        testingDipole_lowerPoint_x = -testingDipole_halfLength;
        testingDipole_upperPoint_x = testingDipole_halfLength;
        testingDipole_mainPoint_x = 0;
        testingDipole_lowerMidPoint = (testingDipole_mainPoint_x+testingDipole_lowerPoint_x)/2;
        testingDipole_upperMidPoint = (testingDipole_mainPoint_x+testingDipole_upperPoint_x)/2;
        diffx = testingDipole_upperMidPoint - testingDipole_lowerMidPoint;
        
        esum = 0;
        for n = 2:noOfBasisPoints-1
            basisPoint_n_minus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
            
            F2 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_upperMidPoint^2+ rho_vector(q)^2 + a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_upperMidPoint^2+ rho_vector(q)^2 + a^2);
            psi_mPlus_nPlus = simp(F2.',delta_s);
            
            F3 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_upperMidPoint^2+ rho_vector(q)^2 + a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_upperMidPoint^2+ rho_vector(q)^2 + a^2);
            psi_mPlus_nMinus = simp(F3.',delta_s);
            
            F4 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_lowerMidPoint^2+ rho_vector(q)^2 + a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n)-s).^2 + testingDipole_lowerMidPoint^2+ rho_vector(q)^2 + a^2);
            psi_mMinus_nPlus = simp(F4.',delta_s);
            
            F5 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_lowerMidPoint^2+ rho_vector(q)^2 + a^2))./...
                sqrt((basisPoint_z(noOfBasisPoints+2)-basisPoint_z(n-1)-s).^2 + testingDipole_lowerMidPoint^2+ rho_vector(q)^2 + a^2);
            psi_mMinus_nMinus = simp(F5.',delta_s);
           
            matrix = factor*(-(psi_mPlus_nPlus/delta)-(psi_mMinus_nMinus/delta)+...
                (psi_mMinus_nPlus/delta)+(psi_mPlus_nMinus/delta));
            esum = esum+matrix*I(n-1);
        end
        efield_X(q,p) = abs(-esum/0.001/1.414);
    end
end

for q = 1:length(rho_vector)
    efield_rho(q,:) = sqrt(efield_X(q,:).^2 + efield_Y(q,:).^2);
end

fig3 = figure('Color','w');
for p = 1:length(rho_vector)
    plot(z_vector,efield_rho(p,:),'LineWidth',1.4);hold on;
end
toc
function q=simp(y,dx)
N=length(y);

mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;

q=(dx/3)*sum(y.'.*mul);
end