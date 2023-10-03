close all
clear all

t = load('data/mesh_points_inc.mat');
p = t.u;

I = 0;
fmin = 33;
fmax = 81;
Nf = 25;
freqs = linspace(fmin, fmax, Nf);

c0 = 1500;
f0 = 75;
lambda_0 = c0/f0;
Dm = 20*lambda_0;


Nr = 41;
h = Dm/(Nr-1);
y_a = linspace(h, Dm-h, Nr);

for  i = 1 : Nf
    i
    % Load incident field.
    ui = load(['data/inc_f',num2str(freqs(i)),'.0.mat']);
    uinc = ui.u;

    % Load total field.
    ut = load(['data/tot_f',num2str(freqs(i)),'.0.mat']);
    utot = ut.u;

    % Create scattered field.
    usc = utot - uinc;

    Gf = load(['data/green_f',num2str(freqs(i)),'.0.mat']);
    G = Gf.u;

% scind = 6843;
% usc = G(:,scind)*G(:,scind).';
    omega = 2.*pi*freqs(i);
    NPM = floor((2 * Dm * freqs(i))/c0);
    kn = omega/c0;
    m = 1:NPM;
    VV = sqrt(2/Dm) *sin((pi*y_a'*m)/Dm);
    lambdam = (m*pi / Dm).^2;
    betam = sqrt(kn*kn - lambdam);
    Dbinv = diag(betam);
    B = h*(VV'*VV);

    [Ua,Sa,Va] = svd(B);
    Dvinv = diag(1./diag(Sa));
    SJ = Dvinv*Ua'*VV';
    Shat  = SJ*usc*SJ';
    pproj = Dbinv*Ua*Shat*Ua'*Dbinv;
%     pproj = h^2*Dbinv*VV'*usc*VV*Dbinv;
    Gp = h*VV'*G;

for m = 1 : NPM
    for n = 1 : NPM
        I = I + conj(pproj(m,n)).*Gp(m,:).*Gp(n,:);
    end
end
end
figure(2)
scatter(p(:,1),p(:,2),35,abs(I),'filled')
hold on
% plot(p(scind,1),p(scind,2),'w*')
circle([390,100],20,32,'w');

axis equal
axis image
% axis([300 500 0 400])
colormap jet
colorbar
set(gca,'Ydir','reverse')
drawnow()
shg
pause(0.5)
shg

