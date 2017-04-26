function plot_field_2D(fname)

field2D = load(fname);
ncnt = length(field2D);

% Assume integer coordinates
min_i = min(field2D(:,1));
min_j = min(field2D(:,2));
max_i = max(field2D(:,1));
max_j = max(field2D(:,2));

ni = max_i - min_i + 1
nj = max_j - min_j + 1

% Allocate 2D arrays for the flow field
is = zeros(nj,ni);
js = zeros(nj,ni);
dns = zeros(nj,ni);
uxs = zeros(nj,ni);
uys = zeros(nj,ni);
us = zeros(nj,ni);
vorts = zeros(nj,ni);
ekins = zeros(nj,ni);
%is = zeros(ni,nj);
%js = zeros(ni,nj);
%dns = zeros(ni,nj);
%uxs = zeros(ni,nj);
%uys = zeros(ni,nj);
%us = zeros(ni,nj);
%vorts = zeros(ni,nj);
%ekins = zeros(ni,nj);

% Read and store flow field into 2D arrays
for n=1:ncnt
    i = field2D(n,1) - min_i + 1;
    j = field2D(n,2) - min_j + 1;
    dn = field2D(n,3);
    ux = field2D(n,4);
    uy = field2D(n,5);
    u2 = ux^2 + uy^2;
    u = sqrt(u2);   
    ekin = 0.5*dn*u2;
    
    is(j,i) = i;
    js(j,i) = j;
    dns(j,i) = dn;
    uxs(j,i) = ux;
    uys(j,i) = uy;
    us(j,i)  = u;
    ekins(j,i) = ekin;
end

% Compute vorticity
for i=2:(ni-1)
    for j=2:(nj-1)
        dx_uy = 0.5*(uys(j,i+1) - uys(j,i-1));
        dy_ux = 0.5*(uxs(j+1,i) - uxs(j-1,i));
        vorts(j,i) = dx_uy - dy_ux;
     end
 end

figure;
surface(dns);
%contourf(dns);
title('Simulated den');

figure;
surface(uxs);
title('Simulated ux');

figure;
contourf(uxs);
title('Simulated ux');

figure;
surface(uys);
title('Simulated uy');

figure;
contourf(uys);
title('Simulated uy');

figure;
%surface(us);
contourf(us);
title('Simulated flow speed');

%figure;
%surface(ekins);
%contourf(ekins);
%title('Simulated kinetic energy');

%figure;
%surface(vorts);
%contourf(vorts);
%title('Simulated vorticity');

figure;
title('Simulated velocity');
%step = 2
%[X,Y] = meshgrid(1:step:nj, 1:step:ni);
%quiver(X,Y,uxs(1:step:nj,1:step:ni),uys(1:step:nj,1:step:ni),6);
quiver(uxs,uys,6);

figure;
title('Simulated velocity');
%[X,Y] = meshgrid(1:(nj/5):nj, 1:(ni/5):ni);
%streamline(uxs,uys,X,Y);
streamslice(uxs,uys,2.5,'cubic');
