clear;close;clc;tic;
noOfSegments = 32;
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
I1 = I.'*1e3;
fig1 = figure('Color','w');
% set(fig1,'FinalDraft','off');
yyaxis left
plot(basisPoint_z,[0 abs(I1) 0],'LineWidth',1.4);
ylim([0 12])
ylabel('Magnitude of Current (mA)')
hold on;
yyaxis right
plot(basisPoint_z,[rad2deg(angle(I1(1))) rad2deg(angle(I1)) rad2deg(angle(I1(noOfSegments-1)))],'--','LineWidth',1.4);
ylim([-180 180])
ylabel('Phase of Current (Degree)')
hold on;
xlim([-wireLength/2 wireLength/2])
xlabel('Wire Length (m)')
legend('Magnitude','Phase','Location','northwest');
z_vector2 = 0:0.001*lambda:0.3*lambda;
testingDipoleHalfLength = 0.001*lambda/2;
rho_vector = [0.02 0.03 0.05 0.1]*lambda;
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
end
z_vector1 = 0:0.02*lambda:0.3*lambda;
rho_vector = [0.02 0.03 0.05 0.1]*lambda;
efield_Z1 = [4.9088    2.5510    0.9988    0.6811    0.5472    0.4806    0.4489    0.4465    0.4851    0.5990   0.8800    1.6252    4.8107    6.8381    3.5858    2.2514
    3.6891    2.0791    1.1458    0.8324    0.6829    0.6044    0.5660    0.5623    0.6064    0.7351   1.0360    1.7332    3.5111    4.2754    3.2104    2.0275
    1.5718    1.4137    1.1424    0.9406    0.8139    0.7368    0.6950    0.6872    0.7251    0.8356   1.0630    1.4522    1.9296    2.1407    1.9308    1.5717
    0.9887    0.9741    0.9363    0.8879    0.8395    0.7972    0.7645    0.7451    0.7433    0.7640   0.8080    0.8669    0.9201    0.9427    0.9219    0.8639];
fig2 = figure('Color','w');

plot(z_vector2,efield_Z2(1,:),'r','LineWidth',1.4,'DisplayName','Method2');hold on;
plot(z_vector2,efield_Z2(2,:),'g','LineWidth',1.4);hold on;
plot(z_vector2,efield_Z2(3,:),'b','LineWidth',1.4);hold on;
plot(z_vector2,efield_Z2(4,:),'m','LineWidth',1.4);hold on;

plot(z_vector1,efield_Z1(1,:),'xr','LineWidth',1.4,'DisplayName','Method1');hold on;
plot(z_vector1,efield_Z1(2,:),'xg','LineWidth',1.4);
plot(z_vector1,efield_Z1(3,:),'xb','LineWidth',1.4);
plot(z_vector1,efield_Z1(4,:),'xm','LineWidth',1.4);


lgd = legend('\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','Location','southoutside','Orientation','horizontal');
lgd.NumColumns = 4;
plot([NaN NaN], [NaN NaN],'x', 'Color', 'k', 'DisplayName', 'Method1');
plot([NaN NaN], [NaN NaN],'-', 'Color', 'k', 'DisplayName', 'Method2');

title('z-comp. of Electric Field for 0.5\lambda dipole')
xlabel('z/\lambda')
ylabel('Electric Field in Volts/meter')
%legend('\rho/\lambda=0.02','\rho/\lambda=0.03','\rho/\lambda=0.05','\rho/\lambda=0.1','Location','northwest');
rho_vector = [0.01 0.02 0.05]*lambda;
for q = 1:length(rho_vector)
    for p = 1:length(z_vector2)
        basisPoint_z(noOfBasisPoint+1) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+2) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+3) = z_vector2(p);
        testingDipoleLowerPoint_y = rho_vector(q)-testingDipoleHalfLength;
        tesingDipoleUpperPoint_y = rho_vector(q)+testingDipoleHalfLength;
        testingDipoleMainPoint_y = rho_vector(q);
        testingDipoleLowerMidPoint = (testingDipoleLowerPoint_y+testingDipoleMainPoint_y)/2;
        testingDipoleUpperMidPoint = (testingDipoleMainPoint_y+tesingDipoleUpperPoint_y)/2;

        esum = 0;
        for n = 2:noOfBasisPoint-1
            F2 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2);
            psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);
            F3 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + wireRadius^2);
            psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);
            F4 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2);
            psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);
            F5 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + wireRadius^2);
            psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);
            matrix = factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_Y(q,p) = abs(-esum/0.001/1.414);
    end
end
for q = 1:length(rho_vector)
    for p = 1:length(z_vector2)
        basisPoint_z(noOfBasisPoint+1) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+2) = z_vector2(p);
        basisPoint_z(noOfBasisPoint+3) = z_vector2(p);
        testingDipoleLowerPoint_x = -testingDipoleHalfLength;
        tesingDipoleUpperPoint_x = testingDipoleHalfLength;
        testingDipoleMainPoint_x = 0;
        testingDipoleLowerMidPoint = (testingDipoleLowerPoint_x+testingDipoleMainPoint_x)/2;
        testingDipoleUpperMidPoint = (testingDipoleMainPoint_x+tesingDipoleUpperPoint_x)/2;

        esum = 0;
        for n = 2:noOfBasisPoint-1
            F2 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mPlus = (1/delta)*simp(F2.',delta_s);
            F3 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleUpperMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mPlus = (1/delta)*simp(F3.',delta_s);
            F4 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nPlus_mMinus = (1/delta)*simp(F4.',delta_s);
            F5 = exp(-1j*k*sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2))./sqrt((basisPoint_z(noOfBasisPoint+2)-basisPoint_z(n-1)-s).^2 + testingDipoleLowerMidPoint^2 + rho_vector(q)^2 + wireRadius^2);
            psi_nMinus_mMinus = (1/delta)*simp(F5.',delta_s);
            matrix = factor2*(psi_nPlus_mPlus-psi_nMinus_mPlus-psi_nPlus_mMinus+psi_nMinus_mMinus);
            esum = esum + matrix*I(n-1);
        end
        efield_X(q,p) = abs(-esum/0.001/1.414);
    end
end
efield_rho = zeros(q,p);
for q = 1:length(rho_vector)
    efield_rho(q,:) = sqrt(efield_X(q,:).^2 + efield_Y(q,:).^2);
end
z_vector1 = 0:0.02*lambda:0.3*lambda;
rho_vector = [0.01 0.02 0.05]*lambda;
efield_rho1 = [ 0.0000    7.0183    3.5597    4.5158    6.0872    7.7851    9.5943   11.4813   12.8256   14.2742    16.0043   18.5518   24.5868    4.7960    0.8234    0.2969
    0.0000    2.6449    2.3441    2.8255    3.7493    4.8045    5.8735    6.9137    7.9098    8.8639    9.8007   10.5421   11.5337    4.4531    1.3189    0.5413
    0.0000    0.4575    0.8013    1.1591    1.5660    1.9939    2.4166    2.8159    3.1761    3.4751    3.6668    3.8370    3.1717    2.2638    1.3949    0.8436];
fig3 = figure('Color','w');
plot(z_vector2,efield_rho(1,:),'r','LineWidth',1.4);hold on;
plot(z_vector2,efield_rho(2,:),'b','LineWidth',1.4);hold on;
plot(z_vector2,efield_rho(3,:),'m','LineWidth',1.4);hold on;
plot(z_vector1,efield_rho1(1,:),'xr','LineWidth',1.4);hold on;
plot(z_vector1,efield_rho1(2,:),'xb','LineWidth',1.4);hold on;
plot(z_vector1,efield_rho1(3,:),'xm','LineWidth',1.4);hold on;
lgd = legend('\rho/\lambda=0.01','\rho/\lambda=0.02','\rho/\lambda=0.05','\rho/\lambda=0.01','\rho/\lambda=0.02','\rho/\lambda=0.05','Location','southoutside','Orientation','horizontal');
lgd.NumColumns = 3;
plot([NaN NaN], [NaN NaN],'x', 'Color', 'k', 'DisplayName', 'Method1');
plot([NaN NaN], [NaN NaN],'-', 'Color', 'k', 'DisplayName', 'Method2');

title('\rho-comp. of Electric Field for 0.5\lambda dipole')
xlabel('\rho/\lambda')
ylabel('Electric Field in Volts/meter')
%legend('\rho/\lambda=0.01','\rho/\lambda=0.02','\rho/\lambda=0.05','Location','northwest');
toc;
function q=simp(y,dx)
N=length(y);
mul=ones(1,N);
mul(2:2:N-1)=4;
mul(3:2:N-2)=2;
q=(dx/3)*sum(y.'.*mul);
end
