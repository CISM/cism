
%% do some interpolation / massaging of Goldberg input .mat file for use in CISM-based
%% setup of Goldberg stream-shelf test case

%% this version used for patching Dan's fields into .nc input that works w/ native CISM dycore

clear all; close all
load Goldberg-input.mat         % Dan's original file, w/ a few regular meshes added

nz = 5;                 % no. of vert levels (needs to be the same as in .config file)

% interpolate thickness, elev., draft onto appropriate mesh

h = interp2( x0g, y0g, hmid, xg, yg, 'spline');    % ice thickness?

u = interp2( xg, yg, U*120, x0g, y0g, 'linear' );      % x comp of vel
v = interp2( xg, yg, V*120, x0g, y0g, 'linear' );      % y comp of vel

% zero vels on edges? Not sure what else to do for now ...
u(1,:) = 0; u(end,:) = 0; u(:,1) = 0; u(:,end) = 0;
v(1,:) = 0; v(end,:) = 0; v(:,1) = 0; v(:,end) = 0;

% place inside of CISM grid, which includes cell buffers at margins for remapping

hh = zeros( size( yyg ) ); 
hh(3:end-2,3:end-2) = h;

uvel = zeros( nz, length(yy0g(:,1)), length(yy0g(1,:)) );
vvel = zeros( nz, length(yy0g(:,1)), length(yy0g(1,:)) );

% extrap vel fields to vertical dim
for j = 1:nz
    uu(j,:,:) = u;
    vv(j,:,:) = v;
end

thck = hh;


% set up topog as well (NOTE: this is now done correctly in the original .py setup script)
R = -( 300 + 600 * sin( pi * ((y-1)*1e3 / 50e3 ) ) );  % across-flow profile
topg = repmat( R, 1, 151+4 );

% pad the n/s sides w/ const. elevation
temp = zeros( size( thck ) );
temp(3:end-2,:) = topg;
temp(1:2,:) = repmat( temp(3,:), 2, 1);
temp(end-1:end,:) = repmat( temp(end-2,:), 2, 1);
topg = temp; clear temp


for j = 1:nz
    uvel(j,3:end-2,3:end-2) = uu(j,:,:);
    vvel(j,3:end-2,3:end-2) = vv(j,:,:);
end

% % crop the last few thickness cells
% thck(:,end-3:end) = 0;

% smooth out upstream-most thickness profiles      
thck(:,3:4) = repmat( thck(:,5), 1, 2 );
% smooth out thickness on north/south edges as well
thck(3,:) = thck(4,:) - (topg(3,:) - topg(4,:) );
thck(end-2,:) = thck(end-3,:) - (topg(end-2,:) - topg(end-3,:) );

% apply a constant initial thickness field
% ind = find( thck > 0 ); thck(ind) = 1200;

% make sure all buffer area thickness is zero
thck(1:2,:) = 0; thck(end-1:end,:) = 0; thck(:,1:2) = 0; thck(:,end-3:end) = 0;
% thck(1:2,:) = 0; thck(end-1:end,:) = 0; thck(:,1:2) = 0; thck(:,end-1:end) = 0;

save stream-shelf.mat thck uvel vvel 

