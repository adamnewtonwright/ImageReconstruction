%% Adam Newton Wright
... July 28, 2017
... Laboratory for Laser Energetics

%% This code is a simulation of how one would go about laser etching a crater on a sample that will be analyzed under microscope
... Using the crater information in image space, one can recover the a focused image of the object even if the image taken was defocused
    
clear all, close all, clc

%% We start by creating an arbitrary object (supposed to look like an etched CR39 sample)
% Black background
o(1:128, 1:128) = 0;

% Place white circles on object
centerx = 64;
centery = 64;
radius = 8;
for ix = 1:128
    for iy = 1:128
        for rf = 1:9
            if (ix-centerx)^2+(iy-centery)^2<=rf^2 && (ix-centerx)^2+(iy-centery)^2>(rf-1)^2
                o(ix,iy) = .1*rf - .1;
            end
        end
    end
end

centerx = 20;
centery = 20;
radius = 4;
for ix = 1:128
    for iy = 1:128
        for rf = 1:9
            if (ix-centerx)^2+(iy-centery)^2<=rf^2 && (ix-centerx)^2+(iy-centery)^2>(rf-1)^2
                o(ix,iy) = .1*rf - .1;
            end
        end
    end
end

centerx = 20;
centery = 80;
radius = 6;
for ix = 1:128
    for iy = 1:128
        for rf = 1:9
            if (ix-centerx)^2+(iy-centery)^2<=rf^2 && (ix-centerx)^2+(iy-centery)^2>(rf-1)^2
                o(ix,iy) = .1*rf - .1;
            end
        end
    end
end

% Place white dot on upper right hand corner @ (8,120) = (row, column)
% This would be our laser ethced crater
o(8,120) = 1;

O = ifftshift(fftn(fftshift(o)));

%% We will then create a point spread function to image this object. The PsF will contain blur
% Pupil Function
ex = linspace(-2,2-4/128,128);
ey = linspace(-2,2-4/128,128);

[xx,yy] = meshgrid(ex,ey);
rr = sqrt((xx).^2 +(yy).^2);

a = zeros(128,128);

a(rr<=.5) = 10;

a2 = zeros(128,128);
a2(rr<1) = 10;
%A = ifftshift(fftn(fftshift(a)));

% Create a blur with zernikes
x = linspace(-2,2-4/128,128);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
one = ones(size(r));

c1(1:11) = 0;
c1(1) = 0;
c1(2) = 0;
c1(3) = 0;
c1(4) = .01;
c1(5) = 0;

c2(1:11) = 0;
c2(1) = 0;
c2(2) = 0;
c2(3) = 0;
c2(4) = 53;
c2(5) = 0;

Theta1 = zernikeps(c1,r,theta,one);
Theta2 = zernikeps(c2,r,theta,one);

H1 = a.*exp(1i.*Theta1);
H2 = a.*exp(1i.*Theta2);

h1 = ifftshift(ifftn(fftshift(H1)));
h2 = ifftshift(ifftn(fftshift(H2)));

p1 = abs(h1).^2;
p2 = abs(h2).^2;

P1 = ifftshift(fftn(fftshift(p1)));
P2 = ifftshift(fftn(fftshift(p2)));

%% Image the objects through point spread functions

O1 = O .* P1;
O2 = O .* P2;

o1 = ifftshift(ifftn(fftshift(O1)));
o2 = ifftshift(ifftn(fftshift(O2)));

%% Extract the observed point spread functions from the images
tic
pobs1 = o1(1:16, 113:128);
pobs2 = o2(1:16, 113:128);

% Add black matrix around it so that dimensions agree
long(1:56,1:128) = 0;
short(1:16,1:56) = 0;

pobs1 = [long; short pobs1 short; long];
pobs2 = [long; short pobs2 short; long];

Pobs1 = ifftshift(fftn(fftshift(pobs1)));
Pobs2 = ifftshift(fftn(fftshift(pobs2)));

%% Reconstruct
% We will reconstruct using the blurred images of the dots as out point
% spread function

Or = ( conj(Pobs1) .* O1 + conj(Pobs2) .*  O2 ) ./ ( abs(Pobs1).^2 + abs(Pobs2).^2 + .0000000001E-300);

or = ifftshift(ifftn(fftshift(O)));

figure
subplot(2,2,1)
imshow(o)
title('Object')

subplot(2,2,2)
imshow(o1)
title('Blurred Image')

subplot(2,2,3)
surf(pobs1)
title('Observed Point Spread Function')

subplot(2,2,4)
imshow(or)
title('Reconstructed Image')