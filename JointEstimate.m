%% Adam Newton Wright
... Joint Estimate Algorithm
... July 28, 2017
... Laboratory for Laser Energetics
    
clear all, close all, clc

%% Code is based on a paper written by Gonsalves, 1982, "Phase Retrieval and diversity in Adaptive Optics"
... And is aimed to take two simulated images that are defocused, estimated their phase and point spread functions, and then restore the images

%% Create a simulated object
% Gray background
o(1:128, 1:128) = .5;

% Randomly create number of circles
num = randi([20,40],1,1);

% Place black circles on image
for k = 1:num
    centerx = randi([5,125],1,1); % Random x-center
    centery = randi([5,125],1,1); % Random y-center
    radius = randi([1,4],1,1);     % Random radius
    for ix = 1:128
        for iy = 1:128
            if (ix-centerx)^2+(iy-centery)^2<=radius^2
                o(ix,iy) = 0;
            %elseif (x-centerx)^2+(y-centery)^2>radius^2 && (x-centerx)^2+(y-centery)^2<=radius^2+(floor(radius*.6))^2
             %   I(x,y) = .7;
            end
        end
    end
end


% Fourier Transform of object
O = ifftshift(fftn(fftshift(o)));


%% Create a Pupil Function
ex = linspace(-2,2-4/128,128);
ey = linspace(-2,2-4/128,128);

[xx,yy] = meshgrid(ex,ey);
rr = sqrt((xx).^2 +(yy).^2);

a = zeros(128,128);
a(rr<=1) = 1;

% Create a defocued phase with zernike polynomials
x = linspace(-2,2-4/128,128);
[X,Y] = meshgrid(x,x);
[theta,r] = cart2pol(X,Y);
one = ones(size(r));

% Initial C-vector
c1(1:11) = 0;
c1(1) = 0;
c1(2) = 0;
c1(3) = 0;
c1(4) = .5; %Blur
c1(5) = 0;
c1(6) = 0;

Theta1 = zernikeps(c1,r,theta,one); % returns fourth zernike polynomial: c(4).*sqrt(3).*(-one + 2.*r.^2)

%% Also can create defocus based on some "camera" parameters
lambda = 550E-9;
deltaz = 20E-6;
radius = 32; %pixels
focal = 73; %pixels for na = .4

phase_defocus=-pi*deltaz/lambda*radius^2/focal^2*r.^2;
phi = exp(1i*phase_defocus);

H1 = a.*exp(2*pi*1i.*Theta1);
% H2 = a.*exp(1i.*Theta1)*exp(-1i*pi*(deltaz/lambda)*(radius^2/focal^2)); %intentional blur
H2 = a.*exp(2*pi*1i.*Theta1).*phi; %intentional blur

h1 = ifftshift(ifftn(fftshift(H1)));
h2 = ifftshift(ifftn(fftshift(H2)));

p1 = abs(h1).^2;
p2 = abs(h2).^2;

P1 = ifftshift(fftn(fftshift(p1)));
P2 = ifftshift(fftn(fftshift(p2)));

Z1 = O .* P1;
Z2 = O .* P2;

z1 = real(ifftshift(ifftn(fftshift(Z1))));
z2 = real(ifftshift(ifftn(fftshift(Z2))));


%% Algorithm takes in two defocused images and an initial guess about the
... phase, point spread function, and an guess about the object

%Initial Guess of zernike polynomials  % "hat" is just variable to distinguish estimate
c(1:11) = 0;

Theta = zernikeps(c,r,theta,one);

Hhat1 = a.*exp(2*pi*1i.*Theta);
Hhat2 = a.*exp(2*pi*1i.*(Theta)).*phi;

hhat1 = ifftshift(ifftn(fftshift(Hhat1)));
hhat2 = ifftshift(ifftn(fftshift(Hhat2)));

% Point spread guess based on phase
phat1 = abs(hhat1).^2;
phat2 = abs(hhat2).^2;

Phat1 = ifftshift(fftn(fftshift(phat1)));
Phat2 = ifftshift(fftn(fftshift(phat2)));

% Guess of the object
Ohat = (conj(Phat1) .* Z1 + conj(Phat2) .* Z2) ./ (abs(Phat1).^2 + abs(Phat2).^2 + 1e-20);

% Initial error in guess
error = sum(sum(abs(z1 - ifftshift(ifftn(fftshift(Ohat.*Phat1)))).^2)) + sum(sum(abs(z2 - ifftshift(ifftn(fftshift(Ohat.*Phat2)))).^2));

% Variables to change direction of c vector
q1 = 1;
q2 = 1;
q3 = 1;
q4 = 1;
q5 = 1;
q6 = 1;
cg(1:11) = 0;
counter = 0;
errorlist = [];




%% Algorithm changes the coefficients of the c vector on at a time, calculates a new phase, point spread function, and object estimate,
... and then calculates a new error.
... If that error is smaller than the previous error, we continue to move the coefficients in that direction
... If the error is bigger than the previous error, we move in the opposite direction
... We only search for the fourth (defocus) and fifth zernike polynomial here, but it can be extended to others
... Stop at about 30 iterations, at which the error has most likely converged
while 1,
    % don't need to optimize for first-third (piston, tip, tilt)
%     %FIRST zernike
%     c = cg;
%     c(1) = c(1) +q1*.01;
% 
% %     Thetag = zernikeps(c,r,theta,one);
%     Thetag = ZernikeMap(c,X,Y,[0 0 1]);
%     
%     Hhatg1 = a.*exp(2*pi*1i.*Thetag);
%     Hhatg2 = a.*exp(2*pi*1i.*(Thetag)).*phi;
%     
%     hhatg1 = ifftshift(ifftn(fftshift(Hhatg1)));
%     hhatg2 = ifftshift(ifftn(fftshift(Hhatg2)));
%     
%     phatg1 = conj(hhatg1).*hhatg1;
%     phatg2 = conj(hhatg2).*hhatg2;
%     
%     Phatg1 = ifftshift(fftn(fftshift(phatg1)));
%     Phatg2 = ifftshift(fftn(fftshift(phatg2)));
%     
%     Ohatg = (conj(Phatg1).*Z1 + conj(Phatg2).*Z2) ./ (conj(Phatg1).*Phatg1 + conj(Phatg2).*Phatg2 + 1e-20);
%     
%     ohat = real(ifftshift(ifftn(fftshift(Ohat))));
%     
%     errorg = sum(sum(z1 - ifftshift(ifftn(fftshift(Ohatg.*Phatg1))))) + sum(sum(z2 - ifftshift(ifftn(fftshift(Ohatg.*Phatg2)))));
%     
%     if errorg < error
%         cg(1) = c(1);
%     elseif errorg > error
%         q1 = -1;
%     else
%         q1 = 0;
%         cg(1) = c(1);
%     end
% 
%  
% 
%     %SECOND
%     c = cg;
%     c(2) = c(2) +q2*.01;
% 
% %     Thetag = zernikeps(c,r,theta,one);
%     Thetag = ZernikeMap(c,X,Y,[0 0 1]);
%     
%     Hhatg1 = a.*exp(2*pi*1i.*Thetag);
%     Hhatg2 = a.*exp(2*pi*1i.*(Thetag)).*phi;
%     
%     hhatg1 = ifftshift(ifftn(fftshift(Hhatg1)));
%     hhatg2 = ifftshift(ifftn(fftshift(Hhatg2)));
%     
%     phatg1 = conj(hhatg1).*hhatg1;
%     phatg2 = conj(hhatg2).*hhatg2;
%     
%     Phatg1 = ifftshift(fftn(fftshift(phatg1)));
%     Phatg2 = ifftshift(fftn(fftshift(phatg2)));
%     
%     Ohatg = (conj(Phatg1).*Z1 + conj(Phatg2).*Z2) ./ (conj(Phatg1).*Phatg1 + conj(Phatg2).*Phatg2 +1e-20);
%     
%     ohat = real(ifftshift(ifftn(fftshift(Ohat))));
%     
%     errorg = sum(sum(z1 - ifftshift(ifftn(fftshift(Ohatg.*Phatg1))))) + sum(sum(z2 - ifftshift(ifftn(fftshift(Ohatg.*Phatg2)))));
%     
%     if errorg < error
%         cg(2) = c(2);
%     elseif errorg > error
%         q2 = -1;
%     else
%         q2 = 0;
%         cg(2) = c(2);
%     end
% 
%     
%     %THIRD
%     c = cg;
%     c(3) = c(3) +q3*.01;
% 
% %     Thetag = zernikeps(c,r,theta,one);
%     Thetag = ZernikeMap(c,X,Y,[0 0 1]);
%     
%     Hhatg1 = a.*exp(2*pi*1i.*Thetag);
%     Hhatg2 = a.*exp(2*pi*1i.*(Thetag)).*phi;
%     
%     hhatg1 = ifftshift(ifftn(fftshift(Hhatg1)));
%     hhatg2 = ifftshift(ifftn(fftshift(Hhatg2)));
%     
%     phatg1 = conj(hhatg1).*hhatg1;
%     phatg2 = conj(hhatg2).*hhatg2;
%     
%     Phatg1 = ifftshift(fftn(fftshift(phatg1)));
%     Phatg2 = ifftshift(fftn(fftshift(phatg2)));
%     
%     Ohatg = (conj(Phatg1).*Z1 + conj(Phatg2).*Z2) ./ (conj(Phatg1).*Phatg1 + conj(Phatg2).*Phatg2 + 1e-20);
%     
%     ohat = real(ifftshift(ifftn(fftshift(Ohat))));
%     
%     errorg = sum(sum(z1 - ifftshift(ifftn(fftshift(Ohatg.*Phatg1))))) + sum(sum(z2 - ifftshift(ifftn(fftshift(Ohatg.*Phatg2)))));
%     
%     if errorg < error
%         cg(3) = c(3);
%     elseif errorg > error
%         q3 = -1;
%     else
%         q3 = 0;
%         cg(3) = c(3);
%     end

    
    %FOURTH
    c = cg;
    c(4) = c(4) +q4*0.05;
    
    Thetag = zernikeps(c,r,theta,one);
    %Thetag = ZernikeMap(c,X,Y,[0 0 1]);
    
    Hhatg1 = a.*exp(2*pi*1i.*Thetag);
    Hhatg2 = a.*exp(2*pi*1i.*(Thetag)).*phi;
    
    hhatg1 = ifftshift(ifftn(fftshift(Hhatg1)));
    hhatg2 = ifftshift(ifftn(fftshift(Hhatg2)));
    
    phatg1 = abs(hhatg1).^2;
    phatg2 = abs(hhatg2).^2;
    
    Phatg1 = ifftshift(fftn(fftshift(phatg1)));
    Phatg2 = ifftshift(fftn(fftshift(phatg2)));
    
    Ohatg = (conj(Phatg1).*Z1 + conj(Phatg2).*Z2) ./ (abs(Phatg1).^2 + abs(Phatg2).^2 + 1e-20);
    
    %ohat = real(ifftshift(ifftn(fftshift(Ohat))));
    
    errorg = sum(sum(abs(z1 - ifftshift(ifftn(fftshift(Ohatg.*Phatg1)))).^2)) + sum(sum(abs(z2 - ifftshift(ifftn(fftshift(Ohatg.*Phatg2)))).^2));

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
    c(5) = c(5) +q5*.001;
    
    Thetag = zernikeps(c,r,theta,one);
    %Thetag = ZernikeMap(c,X,Y,[0 0 1]);
    
    Hhatg1 = a.*exp(2*pi*1i.*Thetag);
    Hhatg2 = a.*exp(2*pi*1i.*(Thetag)).*phi;
    
    hhatg1 = ifftshift(ifftn(fftshift(Hhatg1)));
    hhatg2 = ifftshift(ifftn(fftshift(Hhatg2)));
    
    phatg1 = conj(hhatg1).*hhatg1;
    phatg2 = conj(hhatg2).*hhatg2;
    
    Phatg1 = ifftshift(fftn(fftshift(phatg1)));
    Phatg2 = ifftshift(fftn(fftshift(phatg2)));
    
    Ohatg = (conj(Phatg1).*Z1 + conj(Phatg2).*Z2) ./ (conj(Phatg1).*Phatg1 + conj(Phatg2).*Phatg2 + 1e-20);
    
    %ohat = real(ifftshift(ifftn(fftshift(Ohat))));
    
    errorg = sum(sum(abs(z1 - ifftshift(ifftn(fftshift(Ohatg.*Phatg1)))).^2)) + sum(sum(abs(z2 - ifftshift(ifftn(fftshift(Ohatg.*Phatg2)))).^2));
    
    if errorg < error
        cg(5) = c(5);
    elseif errorg > error
        q5 = -1;
    else
        q5 = 0;
        cg(5) = c(5);
    end
    
    error = errorg;
	counter = counter + 1;
    errorlist = [errorlist, error];
    %break after certain number of iterations
    if counter == 30;
        break
    end

end

%% We then take out retrieved c vector, and calculate the phase, point spread function, and recovered object and plot them
... we also plot the error to show its convergence
c

%Thetag = ZernikeMap(c,X,Y,[0 0 1]);
Thetag = zernikeps(cg,r,theta,one);
    
Hhatg1 = a.*exp(2*pi*1i.*Thetag);
Hhatg2 = a.*exp(2*pi*1i.*(Thetag)).*phi;

hhatg1 = ifftshift(ifftn(fftshift(Hhatg1)));
hhatg2 = ifftshift(ifftn(fftshift(Hhatg2)));

phatg1 = conj(hhatg1).*hhatg1;
phatg2 = conj(hhatg2).*hhatg2;

Phatg1 = ifftshift(fftn(fftshift(phatg1)));
Phatg2 = ifftshift(fftn(fftshift(phatg2)));

Ohatg = (conj(Phatg1).*Z1 + conj(Phatg2).*Z2) ./ (conj(Phatg1).*Phatg1 + conj(Phatg2).*Phatg2 + .1e-20);

ohat = real(ifftshift(ifftn(fftshift(Ohatg))));

% Original object
figure
subplot(2,2,1)
imshow(o);
imagesc(o)
title('Simlated Object','FontSize',24)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

% Error convergence
subplot(2,2,2)
scatter([1:length(errorlist)],errorlist,'.')
title('Error Convergence','FontSize', 24)
ylabel('Error','FontSize', 20)
xlabel('Iterations','FontSize', 20)
grid on

% One blurred image
subplot(2,2,3)
imagesc(z1)
title('Blurred Image', 'FontSize',24)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

% Recovered image
subplot(2,2,4)
imagesc(ohat);
title('Recovered Image','FontSize',24)
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

pause(1);