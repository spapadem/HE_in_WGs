t = load('data/mesh_points_inc.mat');
p = t.u;
t = load('data/mesh_points_tot.mat');
pt = t.u;

I = 0;
fmin = 33;
fmax = 81;
Nf = 25;
freqs = linspace(fmin, fmax, Nf);
for  i = 1 : 10
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

    Gt = load(['data/tot_f_mesh',num2str(freqs(i)),'.0.mat']);
    GT = Gt.u;
% for m = 1 : size(usc,1)
%     for n = 1 : size(usc,2)
%         usc(m,n) = G(m,3000)*G(n,3000);
%     end
% end

% usc = G(:,3000)*G(:,3000).';

for m = 1 : size(usc,1)
    for n = 1 : size(usc,2)
        I = I + (usc(m,n)).*conj(G(m,:)).*conj(G(n,:));
    end
end
end
figure(3)
scatter(p(:,1),p(:,2),45,abs(I),'filled')
hold on
circle([380,140],20,32,'w');
axis equal
axis image
colormap jet
colorbar
set(gca,'Ydir','reverse')
drawnow()
shg

