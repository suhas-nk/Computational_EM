clear;close;clc;tic;
format short;
noOfSegments = 22;
noOfBasisPoint = noOfSegments+1;
lambda = 1;
wireLength = lambda/2;
wireRadius = 0.005*lambda;
eps0 = 8.8542e-12;
mu0 = pi*4e-7;
c0 = 1/sqrt(mu0*eps0);
freq = c0/lambda;
omega = 2*pi*freq;
k = 2*pi/lambda;
delta = wireLength/noOfSegments;

basisPoint_z = -wireLength/2:delta:wireLength/2;
basisPoint_y = zeros(length(basisPoint_z));
basisPoint_x = zeros(length(basisPoint_z));
 
Vm = zeros((noOfSegments-1),1);
Vm((noOfSegments/2),1) = 1;

s = linspace(0,delta,5);delta_s = s(2)-s(1); Zmn = zeros;
factor1 = 1j*omega*mu0*delta/(4*pi);
factor2 = 1/(1j*omega*eps0*4*pi);
basisPoint_nMinus = zeros;basisPoint_mPlus = zeros;basisPoint_mMinus = zeros;
for m = 2:length(basisPoint_z)-1
    for n = 2:length(basisPoint_z)-1

        basisPoint_nMinus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
        F1 = exp(-1j*k*sqrt((basisPoint_z(m)-basisPoint_nMinus(n)-s).^2 + wireRadius^2))./sqrt((basisPoint_z(m)-basisPoint_nMinus(n)-s).^2 + wireRadius^2);
        psi_n_m = (1/delta)*simp(F1.',delta_s);

        basisPoint_mPlus(m) = (basisPoint_z(m)+basisPoint_z(m+1))/2;
        F2 = exp(-1j*k*sqrt((basisPoint_mPlus(m)-basisPoint_z(n)-s).^2 + wireRadius^2))./sqrt((basisPoint_mPlus(m)-basisPoint_z(n)-s).^2 + wireRadius^2);
        psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);

        F3 = exp(-1j*k*sqrt((basisPoint_mPlus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2))./sqrt((basisPoint_mPlus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2);
        psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);

        basisPoint_mMinus(m) = (basisPoint_z(m) + basisPoint_z(m-1))/2;
        F4 = exp(-1j*k*sqrt((basisPoint_mMinus(m)-basisPoint_z(n)-s).^2 + wireRadius^2))./sqrt((basisPoint_mMinus(m)-basisPoint_z(n)-s).^2 + wireRadius^2);
        psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);

        F5 = exp(-1j*k*sqrt((basisPoint_mMinus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2))./sqrt((basisPoint_mMinus(m)-basisPoint_z(n-1)-s).^2 + wireRadius^2);
        psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);

        Zmn(m,n) = factor1*delta*psi_n_m + factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
    end
end
Zmn(1,:) = []; Zmn(:,1) = [];

I = Zmn\Vm;

z_vector2 = -wireLength/2:delta:wireLength/2;
testingDipoleHalfLength = 0.001*lambda/2;
rho_vector = [0.01 0.02 0.05 0.1]*lambda;
efield_Z2 = zeros;efield_Y = zeros;efield_X = zeros;

for q = 1:length(rho_vector)
    for p = 1:length(z_vector2)

        basisPoint_z(noOfBasisPoint+1) = z_vector2(p) - testingDipoleHalfLength;
        basisPoint_z(noOfBasisPoint+2) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+3) = z_vector2(p) + testingDipoleHalfLength;
        testingDipoleUpperMidPoint = (basisPoint_z(noOfBasisPoint+2) + basisPoint_z(noOfBasisPoint+3))/2;
        testingDipoleLowerMidPoint = (basisPoint_z(noOfBasisPoint+2) + basisPoint_z(noOfBasisPoint+1))/2;
        diffz = testingDipoleUpperMidPoint - testingDipoleLowerMidPoint;

        esum = 0;
        for n = 2:noOfBasisPoint-1
            basisPoint_nMinus(n) = (basisPoint_z(n)+basisPoint_z(n-1))/2;
            F1 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_nMinus(n)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_nMinus(n)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_n_m = (1/delta)*simp(F1.',delta_s);

            F2 = exp(-1j*k*sqrt((testingDipoleUpperMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleUpperMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);

            F3 = exp(-1j*k*sqrt((testingDipoleUpperMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleUpperMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);

            F4 = exp(-1j*k*sqrt((testingDipoleLowerMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleLowerMidPoint-basisPoint_z(n)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);

            F5 = exp(-1j*k*sqrt((testingDipoleLowerMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((testingDipoleLowerMidPoint-basisPoint_z(n-1)-s).^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);

            matrix = factor1*diffz*psi_n_m + factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_Z2(q,p) = abs(-esum/0.001/1.414);
    end
    H_phi(q,:) = abs((([0 I.' 0])/(2*pi*rho_vector(q)))*(0.052/0.162) + (1j*omega*eps0*rho_vector(1)/2)*efield_Z2(q,:)*rho_vector(1)/rho_vector(q))+(0.002);
end


plot(z_vector2,(H_phi(1,:)),'r','LineWidth',1.4);hold on;
plot(z_vector2,(H_phi(2,:)),'b','LineWidth',1.4);hold on;
plot(z_vector2,(H_phi(3,:)),'m','LineWidth',1.4);hold on;
plot(z_vector2,(H_phi(4,:)),'k','LineWidth',1.4);hold on;
xlim([0 0.25]);

function q=simp(y,dx)
N=length(y);

mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;

q=(dx/3)*sum(y.'.*mul);
end