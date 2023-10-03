clear all
t = load('data/mesh_points_inc.mat');
p = t.u;
t = load('data/mesh_points_tot.mat');
pt = t.u;

I = 0;
fmin = 33;
fmax = 81;
Nf = 25;
freqs = linspace(fmin, fmax, Nf);
for  i = 1 : 19
    i
    % Load incident field.
    ui = load(['data/inc_f',num2str(freqs(i)),'.0.mat']);
    uinc = ui.u;

%     Load total field.
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
% scind = 6843;
% usc = G(:,scind)*G(:,scind).';

for m = 1 : size(usc,1)
    for n = 1 : size(usc,2)
        I = I + conj(usc(m,n)).*G(m,:).*G(n,:);
    end
end

end

figure(3)
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
