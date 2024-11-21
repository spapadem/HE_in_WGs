close all
clear


filename = ['mesh_inc_var_depth'];

gmsh_to_fem(filename)
nodes = load([filename,'_nodes.txt']);
triangles = load([filename,'_elements.txt']);




I = 0;
fmin = 70.5;
fmax = 79.5;
h = 1;
Nf = (fmax-fmin)/h + 1;
freqs = linspace(fmin, fmax, Nf);

c0 = 1500;
f0 = 75;
lambda_0 = c0/f0;
Dm = 20*lambda_0;

Nx = 201;
Ny = 201;

x_sc = 18*lambda_0;
y_sc = 3.5*lambda_0;
x = linspace(x_sc-4*lambda_0,x_sc+4*lambda_0,Nx);
y = linspace(y_sc-4*lambda_0,y_sc+4*lambda_0,Ny);

[X,Y] = meshgrid(x,y);

Nr = 81;
h = Dm/(Nr-1);
y_a = linspace(h, Dm-h, Nr)';
b_a = y_a(end) - y_a(1);
x_a = 2*lambda_0*ones(Nr,1);
for  i = 1:1
    % Load incident field.
    ui = load(['data/inc_f',num2str(freqs(i)),'.mat']);
    uinc = ui.u;

    % Load total field.
    ut = load(['data/tot_f',num2str(freqs(i)),'.mat']);
    utot = ut.u;

    % Create scattered field.
    usc = utot - uinc;

    Gf = load(['data/green_f',num2str(freqs(i)),'.mat']);
    G = Gf.u;

    omega = 2.*pi*freqs(i);
    NPM = floor((2 * Dm * freqs(i))/c0);
    kn = omega/c0;
    m = 1:NPM;
    VV = sqrt(2/Dm) *sin((pi*y_a*m)/Dm);
    lambdam = (m*pi / Dm).^2;
    betam = sqrt(kn*kn - lambdam);
    if betam(end) <=1e-8
        fprintf('Warning: betam is too small for i = %d\n',i);
    end
    Dbinv = diag(betam);

    A = h*VV'*(VV);

    for m = 1 : length(y_a)
%     figure(2)
% trisurf(triangles,nodes(:,1),nodes(:,2),abs(G(m,:)),linestyle="none");
% shading interp
% % title('TPs','interpreter','latex','FontSize',16)
% colorbar
% axis equal; axis image;
% colormap jet
% view([0 0 90])
% % drawnow;
% % shg
% pause(0.2)
    end


    [Ua,Sa,Va] = svd(A);
    Dvinv = diag(1./diag(Sa));
    SJ = zeros(Nr,NPM);

    for j = 1 : NPM
        for ii = 1 : NPM
            SJ(:,j) = SJ(:,j) + Ua(ii,j)*VV(:,ii);
        end
        % SJ(:,j) = SJ(:,j)/Sa(j,j);
    end
    Shat  = SJ.'*usc*SJ;
    pproj = Dbinv*Ua*Shat*Ua'*Dbinv;
    % pproj = h^2*Dbinv*VV'*usc*VV*Dbinv;
    Gp = h*Dbinv*VV'*G;

    for m = 1 : NPM
        for n = 1 : NPM
            I= I + conj(pproj(m,n))*Gp(m,:).*Gp(n,:);
        end  
    end
% end
% return
figure(2)
trisurf(triangles,nodes(:,1),nodes(:,2),abs(I),linestyle="none");
shading interp
% title('TPs','interpreter','latex','FontSize',16)
colorbar
axis equal; axis image;
colormap jet
view([0 0 90])
% drawnow;
% shg
% imagesc(x,y,abs(reshape(I,Nx,Ny)))
hold on
circle([x_sc,y_sc],20,32,'w');
axis equal
axis image
% colormap jet
colorbar
set(gca,'Ydir','reverse')
drawnow()
shg
pause(0.1)
shg
end



