close all
clear


filename = ['mesh_inc_var_depth'];

gmsh_to_fem(filename)
nodes = load([filename,'_nodes.txt']);
triangles = load([filename,'_elements.txt']);




I = 0;
fmin = 41;
fmax = 89;
h = 2;
Nf = (fmax-fmin)/h + 1;
freqs = linspace(fmin, fmax, Nf);

c0 = 1500;
f0 = 75;
lambda_0 = c0/f0;
Dm = 20*lambda_0;

Nx = 201;
Ny = 201;

x_sc = 19.5*lambda_0;
y_sc = 5*lambda_0;
x = linspace(x_sc-4*lambda_0,x_sc+4*lambda_0,Nx);
y = linspace(y_sc-4*lambda_0,y_sc+4*lambda_0,Ny);

[X,Y] = meshgrid(x,y);

Nr = 41;
h = Dm/(Nr-1);
y_a = linspace(h, Dm-h, Nr);

for  i = 1 : Nf
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
    pproj = Ua*Shat*Ua';
%     pproj = h^2*Dbinv*VV'*usc*VV*Dbinv;
    Gp = h*VV'*G;

for m = 1 : NPM
    for n = 1 : NPM
        I = I + conj(pproj(m,n)).*Gp(m,:).*Gp(n,:);
    end
end


figure(2)
trisurf(triangles,nodes(:,1),nodes(:,2),abs(I),linestyle="none");
shading interp
% title('TPs','interpreter','latex','FontSize',16)
colorbar
axis equal; axis image;
colormap jet
view([0 0 90])
drawnow;
shg
% imagesc(x,y,abs(reshape(I,Nx,Ny)))
hold on
circle([390,100],20,32,'w');
axis equal
axis image
colormap jet
colorbar
set(gca,'Ydir','reverse')
drawnow()
shg
pause(0.1)
shg
end



