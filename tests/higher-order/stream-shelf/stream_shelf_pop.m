
%% do some interpolation / massaging of Goldberg input .mat file for use in CISM-based
%% setup of Goldberg stream-shelf test case

%% this version used for patching Dan's fields into .nc input that works w/ POP2x and Bisicles dycore

clear all; close all
load Goldberg-input-pop.mat         % Dan's original file, w/ a few regular meshes added

nz = 5;                 % no. of vert levels (needs to be the same as in .config file)

% interpolate thickness, elev., draft onto appropriate mesh

h = interp2( x0g, y0g, hmid, xg, yg, 'spline');    % ice thickness?

u = interp2( x1g, y1g, U*120, xg, yg, 'linear' );      % x comp of vel
v = interp2( x1g, y1g, V*120, xg, yg, 'linear' );      % y comp of vel

% extrap vel fields to vertical dim
for j = 1:nz
    uu(j,:,:) = u;
    vv(j,:,:) = v;
end

thck = h;

% set up topog as well (NOTE: this is now done correctly in the original .py setup script)
R = -( 300 + 600 * sin( pi * ( y / 50e3 ) ) );  % across-flow profile
topg = repmat( R, 1, 168 );

uvel = uu;
vvel = vv;

% % crop the last few thickness cells
% thck(:,end-3:end) = 0;

% smooth out upstream-most thickness profiles      
thck(:,1:2) = repmat( thck(:,3), 1, 2 );
% smooth out thickness on north/south edges as well
thck(1,:) = thck(2,:) - (topg(1,:) - topg(2,:) );
thck(end,:) = thck(end-1,:) - (topg(end,:) - topg(end-1,:) );


save stream-shelf-pop.mat thck uvel vvel 

