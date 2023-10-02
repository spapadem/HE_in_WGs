t = load('data/mesh_points.mat');
p = t.u;

I = 0;
fmin = 33;
fmax = 81;
Nf = 25;
freqs = linspace(fmin, fmax, Nf);

c0 = 1500;
f0 = 75;
lambda_0 = c0/f0;
Dm = 10*lambda_0;


Nr = 41;
h = Dm/(Nr-1);
y_a = linspace(h, Dm-h, Nr);
freqs = 73;
for  i = 1 : 1
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

for m = 1 : size(usc,1)
    for n = 1 : size(usc,2)
        usc(m,n) = G(m,3000)*G(n,3000);
    end
end


    omega = 2.*pi*freqs(i);
    NPM = floor((2 * Dm * freqs(i))/c0);
    kn = omega/c0;
    m = 1:NPM;
    VV = sqrt(2/Dm) *sin((pi*y_a'*m)/Dm);
    lambdam = (m*pi / Dm).^2;
    betam = sqrt(kn*kn - lambdam);
    Dbinv = diag(betam);
%     B = h*(VV'*VV);
% 
%     [Ua,Sa,Va] = svd(B);
%     Dvinv = diag(1./diag(Sa));
%     SJ = Dvinv*Ua'*VV';
%     Shat  = SJ*usc*SJ';
%     pproj = Ua*Shat*Ua';
    pproj = h^2*Dbinv*VV'*usc*VV*Dbinv;
    Gp = h*VV'*G;

for m = 1 : NPM
    for n = 1 : NPM
        I = I + conj(pproj(m,n)).*Gp(m,:).*Gp(n,:);
    end
end
end
figure(2)
scatter(p(:,1),p(:,2),45,abs(I),'filled')
hold on
circle([380,140],40,32,'w');
axis equal
axis image
colormap jet
colorbar
drawnow()
shg

