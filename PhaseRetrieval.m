%% Adam Newton Wright
... July 6, 2017
... Phase Retrieval and Reconstruction Simulation
... Laboratory for Laser Energetics

%% Code is based on a paper written by Gonsalves, 1982, "Phase Retrieval and diversity in Adaptive Optics"

clear all, close all, clc


%% We first create an artificial, simulated object; it will be a white circle on a black background
% Gray background
o(1:128, 1:128) = 0;

% Place white circle on image
centerx = 64;
centery = 64;
radius = 4;
for ix = 1:128
    for iy = 1:128
        if (ix-centerx)^2+(iy-centery)^2<=radius^2
            o(ix,iy) = 1;
        end
    end
end


%% We then create a pupil function, which will be the lens the "object" is imaged through
x = linspace(-2,2-4/128,128);
y = linspace(-2,2-4/128,128);

[xx,yy] = meshgrid(x,y);
rr = sqrt((xx).^2 +(yy).^2);
a   = zeros(128,128);

a(rr<=.5) = 10;
A = ifftshift(fftn(fftshift(a)));

%% We then create an arbitrary phase that the object will encounter as it passes through the lens
... We do this using a few of the Zernike polynomials
... We are interested in the fourth polynomial, which is defocus
x = linspace(-1,1,128);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
one = ones(size(r));
ctest(1:11) = 0;
ctest(1) = .2;
ctest(2) = .1;
ctest(3) = 0.1;
ctest(4) = 50;
ctest(5) = 0.1;
ctest(6) = -.7;
Thetatest = zernikeps(ctest,r,theta,one);


%% Using the pupil function, and the phase, we can create a point spread function
Htest = a.*exp(1i.*Thetatest);
htest = ifftshift(ifftn(fftshift(Htest)));
ptest = conj(htest).*htest;
Ptest = ifftshift(fftn(fftshift(ptest)));


%% We can also create a point spread function with defocus using the following information
% Complex exponential function
f = 2.24;% 586.5; %pixels
L = 550E-9; %meters %wavelength

dz1 = 1E-6; %meters
g1=exp(-1i*pi/L*dz1*(xx.^2+yy.^2)/f^2);

dz2 = 9E-6;
g2 = exp(-1i*pi/L*dz2*(xx.^2+yy.^2)/f^2);

% Multiply pupil function and G
pg1 = a .* g1;
pg2 = a .* g2;

% Fourier transform pupilG
PG1 = ifftshift( fftn( fftshift( pg1)));
PG2 = ifftshift( fftn( fftshift( pg2)));

% PSF
p1=abs(PG1).^2;
p2 = abs(PG2).^2;


%% We then convolve the point spread function with our object to create our image
... convolution is simply multiplication in the fourier domain
O = ifftshift(fftn(fftshift(o)));
P1 = ifftshift(fftn(fftshift(p1)));
P2 = ifftshift(fftn(fftshift(p2)));

Z1 = O .* P1;
Z2 = O .* P2;
Ztest = O .* Ptest;

z1 = ifftshift(ifftn(fftshift(Z1)));
z2 = ifftshift(ifftn(fftshift(Z2)));
ztest = real(ifftshift(ifftn(fftshift(Ztest))));


%% Now that we have an image, we want to use that image to find the point spread function that was used in the convolution
... without any information except the image itself
... We start out by making a guess of the phase and the point spread function

% Initial C vector
c(1:11) = 0;

% Zernike Functions & initial Wavefront
Theta = zernikeps(c,r,theta,one);

Hhat = a.*exp(1i.*Theta);

hhat = ifftshift(ifftn(fftshift(Hhat)));

phat = conj(hhat).*hhat;

% Error in initial guess
error = sum(sum((ptest - phat).^2));

q1 = 1;
q2 = 1;
q3 = 1;
q4 = 1;
q5 = 1;
q6 = 1;
cg(1:11) = 0;

figure
movegui('center')
set(gcf,'Units','inches',...
 'Position',[0 0 20 12]);

errorlist = [error];
%% What we will now do is change the coefficients of each zernike polynomial by some amount, create a new point spread function from that
... and then calculate a new error.
... If this error is smaller than the previous error, we will continue to move that polynomial's coefficient in that direction
... If this error is larger than the previous error, we will move the coefficient in the other direction
... We continue until the error has converged to a small enough value
while 1,
    %FIRST
    c = cg;
    c(1) = c(1) +q1*.01;

    Thetag = zernikeps(c,r,theta,one);
    
    hhatg = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetag))));
    
    phatg = conj(hhatg).*hhatg;
    
    errorg = sum(sum((ptest - phatg).*(ptest-phatg)));
    if errorg < error
        cg(1) = c(1);
    elseif errorg > error
        q1 = -1;
    else
        q1 = 0;
        cg(1) = c(1);
    end

 

    %SECOND
    c = cg;
    c(2) = c(2) +q2*.01;

    Thetag = zernikeps(c,r,theta,one);
    
    hhatg = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetag))));
    
    phatg = conj(hhatg).*hhatg;
    
    errorg = sum(sum((ptest - phatg).*(ptest-phatg)));
    
    if errorg < error
        cg(2) = c(2);
    elseif errorg > error
        q2 = -1;
    else
        q2 = 0;
        cg(2) = c(2);
    end

    
    %THIRD
    c = cg;
    c(3) = c(3) +q3*.01;

    Thetag = zernikeps(c,r,theta,one);
    
    hhatg = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetag))));
    
    phatg = conj(hhatg).*hhatg;
    
    errorg = sum(sum((ptest - phatg).*(ptest-phatg)));
    
    if errorg < error
        cg(3) = c(3);
    elseif errorg > error
        q3 = -1;
    else
        q3 = 0;
        cg(3) = c(3);
    end

    
    %FOURTH
    c = cg;
    c(4) = c(4) +q4*1;
    
    Thetag = zernikeps(c,r,theta,one);
    
    hhatg = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetag))));
    
    phatg = conj(hhatg).*hhatg;
    
    errorg = sum(sum((ptest - phatg).*(ptest-phatg)));
    
    if errorg < error
        cg(4) = c(4);
    elseif errorg > error
        q4 = -1;
    else
        q4 = 0;
        cg(4) = c(4);
    end

    
    %FIFTH

    c = cg;
    c(5) = c(5) +q5*.01;
    Thetag = zernikeps(c,r,theta,one);
    
    hhatg = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetag))));
    
    phatg = conj(hhatg).*hhatg;
    
    errorg = sum(sum((ptest - phatg).*(ptest-phatg)));
    
    if errorg < error
        cg(5) = c(5);
    elseif errorg > error
        q5 = -1;
    else
        q5 = 0;
        cg(5) = c(5);
    end

    %SIXTH

    c = cg;
    c(6) = c(6) +q6*.01;
    Thetag = zernikeps(c,r,theta,one);
    
    hhatg = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetag))));
    
    phatg = conj(hhatg).*hhatg;
    
    errorg = sum(sum((ptest - phatg).*(ptest-phatg)));
    
    if errorg < error
        cg(6) = c(6);
    elseif errorg > error
        q6 = -1;
    else
        q6 = 0;
        cg(6) = c(6);
    end
    %surf(phatg);
    %pause(.1)
    error = errorg;
    errorlist = [errorlist, error];
    if error <.00001
        break
    end
    subplot(2,2,2)
    plot(errorlist,'.')
    grid on
    title('Error Convergence', 'FontSize',24)
    xlabel('Iterations','FontSize',20)
    ylabel('Error','FontSize',20)
    
    
    subplot(2,2,1)
    surf(ptest)
    title('Simulated PSF', 'FontSize',24)
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    subplot(2,2,3)
    surf(phat)
    title('Estimated PSF', 'FontSize',24)
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    
    
    subplot(2,2,4)
    surf(phatg)
    title('Recovered PSF', 'FontSize',24)
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    pause(.1)
end

%% Now that we have found our best c-vector coefficients, we will calculate our point spread function
... and our reconstructed object
c(1:6)
Thetag = zernikeps(cg,r,theta,one);

hhatg = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetag))));

phatg = conj(hhatg).*hhatg;

Phatg = ifftshift(fftn(fftshift(phatg)));

Zrec = O .* Phatg;

%sharp point spread function
Thetarec = Thetatest - Thetag;

hhatrec = ifftshift(ifftn(fftshift(a.*exp(1i.*Thetarec))));

phatrec = conj(hhatrec).*hhatrec;

% figure
% subplot(2,2,2)
% plot(errorlist,'.')
% grid on
% title('Error Convergence', 'FontSize',24)
% xlabel('Iterations','FontSize',20)
% ylabel('Error','FontSize',20)
% 
% 
% subplot(2,2,1)
% surf(ptest)
% title('Original PSF', 'FontSize',24)
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% 
% subplot(2,2,3)
% surf(phat)
% title('Estimated PSF', 'FontSize',24)
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% 
% 
% subplot(2,2,4)
% surf(phatg)
% title('Recovered PSF', 'FontSize',24)
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])






