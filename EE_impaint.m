% -----------------------------------------------------------
% -----------------------------------------------------------
% f------The Observed Picture(M*N);
% u------Denoised picture(M*N);
% p------p = (u_x, u_y);
% n------n = (u_x, u_y)/abs((u_x, u_y));
% m------m = n;
% lambda1, lambda2, lambda4-----Lagrangian Multipliers;
% r1, r2, r4------penalty parameters;
% FDx------Forward Difference function in the x direction;
% FDy------Forward Difference function in the y direction;
% BDx------Backward Difference function in the x direction;
% BDy------Backward Difference function in the y direction;
% CDx------Central Difference function in the x direction;
% CDy------Central Difference function in the y direction;
% h------Laplacian Operator Kernel DBxDFx + DByDFy;
% Fh-----Fourier Transformation of h;
% Ac-------Average operator from square nodes to circle nodes;
% As-------Average operator from circle nodes to square nodes;
% ------------------------------------------------------------
% ------------------------------------------------------------


%%
clear all

% Draw the inpainting area on the picture;
[X,Y,u0,M,N] = Interactive;

figure(1)
imshow(u0)
% Modify X,Y; So that 1<=x<=N, 1<=y<=M;
row = (1<=X&X<=N) & (1<=Y&Y<=M);
X(row==0) = [];
Y(row==0) = [];

% Create a D matrix to track the point in the 
% inpainting area D;
index = (X-1)*M + Y;
D = zeros(M,N);
D(index) = 1;
D = logical(D);
tao = 1 - D;
tao = logical(tao);
u0tao = u0(tao);
%%
% Parameters of the model;
r1 =1;
r2 = 800;
r3 = 500;
r4 = 1200;
eta = 500;
a = 1;
b = 20;
eps = 10^(-2);


% Calculate the Fourier Transform of h;
h = zeros(M,N);
h(2,1) = 1;
h(1,2) = 1;
h(M,1) = 1;
h(1,N) = 1;
h(1,1) = -4;
Fh = fft2(h);

% Call the Forward Difference function "FDx(U)", "FDy(U)" and 
% Backward Divergence function "Dive(X,Y);
[FDx,FDy,BDx,BDy,CDx,CDy,avgx,avgy,Ac,As] = fun;

% Initialization
% u and lambda1 are on the stuffed circle nodes;
u = zeros(M,N); 
lambda1 = zeros(M,N);

% p1, n1, lambda21 and lambda41 are on the circle nodes;
% p2, n2, lambda22 and lambda 42 are on the square nodes;
p1 = zeros(M,N);
p2 = p1;
n1 = p1;
n2 = p1;
lambda21 = lambda1;
lambda22 = lambda1;
lambda3 = lambda1;
lambda41 = lambda1;
lambda42 = lambda1;

% Construct g1, g2, g3, g4;
g1 = zeros(M,N);
g2 = g1;
g3 = g1;
g4 = g1;

g1(1,2) = 1;
g1(1,N) = 1;
g1(1,1) = -2;
Fg1 = fft2(g1);

g2(1,N) = 1;
g2(2,1) = 1;
g2(2,N) = -1;
g2(1,1) = -1;
Fg2 = fft2(g2);

g3(M,1) = 1;
g3(1,2) = 1;
g3(M,2) = -1;
g3(1,1) = -1;
Fg3 = fft2(g3);

g4(M,1) = 1;
g4(2,1) = 1;
g4(1,1) = -2;
Fg4 = fft2(g4);

Cal = 1000;

tic
%%
%------------------------------------------------------------------------
%---------------------------Main Loop------------------------------------
%------------------------------------------------------------------------
for k = 1 : Cal
    
    %--------------------------------------------------------------------
    %------------------------ Update v;----------------------------------
    w = u - lambda3/r3;
    wtao = w(tao);
    v = w;
    MM = 1 - eta./(r3*abs(wtao-u0tao));
    MM = max(0,MM);
    v(tao) = u0tao + MM.*(wtao-u0tao);
    
    %--------------------------------------------------------------------
    %------------------------ Update u;----------------------------------
    uold = u;
    F = r2*(BDx(p1)+BDy(p2)) + BDx(lambda21)+BDy(lambda22) - r3*v - lambda3;
    Fu = fft2(F)./(r2*Fh-r3);
    u = real(ifft2(Fu));
    
    %%
    %--------------------------------------------------------------------
    %------------------------ Update m;----------------------------------
    % z at circle;
    zc1 = n1 + 1/r4*(lambda41 + (avgx(lambda1) + r1).*p1);
    zc2 = Ac(n2)+1/r4*(Ac(lambda42) + (avgx(lambda1) + r1).*Ac(p2));
    
    zs1 = As(n1) + 1/r4*(As(lambda41) + (avgy(lambda1) + r1).*As(p1));
    zs2 = n2 + 1/r4*(lambda42 + (avgy(lambda1)+r1).*p2);

    
    zc = sqrt(zc1.^2+zc2.^2);
    zs = sqrt(zs1.^2+zs2.^2);

    m1 = zc1.*(zc<=1);
    m1(zc>1) = zc1(zc>1)./zc(zc>1);
    m2 = zs2.*(zs<=1);
    m2(zs>1) = zs2(zs>1)./zs(zs>1);
    
    %---------------------------------------------------------------------
    %------------------------ Update p------------------------------------
    % at circle nodes;
    cc = a + r1 +avgx(lambda1) + b * (CDx(n1) + BDy(avgx(n2))).^2;
    qc1 = FDx(u) + (r1/r2 + 1/r2*avgx(lambda1)).*m1 - 1/r2*lambda21;
    qc2 = CDy(avgx(u)) + (r1/r2 + 1/r2*avgx(lambda1)).*Ac(m2) - 1/r2*Ac(lambda22);
    p1 = max(0,1-cc./(r2*sqrt(qc1.^2+qc2.^2))).*qc1;
    
    % at square nodes;
    cs = a + r1 + avgy(lambda1) + b * (CDy(n2) + BDx(avgy(n1))).^2;
    qs1 = CDx(avgy(u)) + (r1/r2 + 1/r2*avgy(lambda1)).*As(m1) - 1/r2*As(lambda21);
    qs2 = FDy(u) + (r1/r2 + 1/r2*avgy(lambda1)).*m2 - 1/r2*lambda22;
    p2 = max(0,1 - cs./(r2*sqrt(qs1.^2+qs2.^2))).*qs2;
    %%
    %---------------------------------------------------------------------
    %------------------------ Update n------------------------------------    
    AP1 = [p1(:,1)+p1(:,N), p1(:,2:N)+p1(:,1:N-1)];
    AP1 = AP1/2;
    AP2 = [p2(1,:)+p2(M,:);p2(2:M,:)+p2(1:M-1,:)];
    AP2 = AP2/2;
    AP = sqrt(AP1.^2 + AP2.^2);
    c = 2*b * max(AP(:));
    
    n1old = n1;
    n2old = n2;
    tol = 1;
    
     while tol>eps
        %%
        % Build up the System of Equations:
        % [a11, a12; a21, a22]*[Fn1new; Fn2new] = [Ff1; Ff2];
        a11 = -c * Fg1 + r4;
        a12 = -c * Fg2;
        vv = (c - 2*b*AP).*(BDx(n1old)+BDy(n2old));
        f1 = r4*m1 -lambda41 - FDx(vv);
        Ff1 = fft2(f1);
        
        a21 = -c * Fg3;
        a22 = -c * Fg4 + r4;
        f2 = r4*m2 - lambda42 - FDy(vv);
        Ff2 = fft2(f2);
        
        Det = a11.* a22 - a12.* a21;
        Fn1 = (a22.*Ff1 - a12.*Ff2)./Det;
        Fn2 = (-a21.*Ff1 + a11.*Ff2)./Det;
        n1new = real(ifft2(Fn1));
        n2new = real(ifft2(Fn2));
        
%         tol = sqrt((n1new - n1old).^2 + (n2new - n2old).^2);
        dif1 = -2*FDx(b*AP.*(BDx(n1new)+BDy(n2new)))+r4*(n1new-m1)+lambda41;
        dif2 = -2*FDy(b*AP.*(BDx(n1new)+BDy(n2new)))+r4*(n2new-m2)+lambda42;
        tol = max(abs(dif1(:)))+max(abs(dif2(:)));
        n1old = n1new;
        n2old = n2new;
        
     end
    
    n1 = n1new;
    n2 = n2new;
    %%
    %----------------------------------------------------------------------
    %------------------ Update lambda1 throuth lambda4---------------------
    lambda1old = lambda1;
    lambda21old = lambda21;
    lambda22old = lambda22;
    lambda3old = lambda3;
    lambda41old = lambda41;
    lambda42old = lambda42;
    
    AM1 = [m1(:,1)+m1(:,N), m1(:,2:N)+m1(:,1:N-1)];
    AM1 = AM1/2;
    AM2 = [m2(1,:)+m2(M,:);m2(2:M,:)+m2(1:M-1,:)];
    AM2 = AM2/2;
    
    R = AP - (AM1.*AP1+AM2.*AP2);
    R1(k) = sum(abs(R(:)));
    
    
    lambda1 = lambda1 + r1* R.*(R>=10^(-12));
    lambda21 = lambda21 + r2*( p1 - FDx(u) );
    lambda22 = lambda22 + r2*( p2 - FDy(u) );
    lambda3 = lambda3 + r3*( v - u );
    lambda41 = lambda41 + r4*( n1 - m1 );
    lambda42 = lambda42 + r4*( n2 - m2 );
    
    %----------------------------------------------------------------------
    %-------------------- Relative Errors 1 through 4 ---------------------
    L1(k) = sum(sqrt((lambda1(:)-lambda1old(:)).^2))...
        /sum(sqrt(lambda1old(:).^2));
    
    L2(k) = (sum(sqrt((lambda21(:)-lambda21old(:)).^2)) ...
        + sum(sqrt((lambda22(:)-lambda22old(:)).^2)))...
        /(sum(sqrt(lambda21old(:).^2))+sum(sqrt(lambda22old(:).^2)));
    
    L3(k) = sum(sqrt((lambda3(:)-lambda3old(:)).^2)) ...
        /sum(sqrt(lambda3old(:).^2));
    
    L4(k) = (sum(sqrt((lambda41(:)-lambda41old(:)).^2))...
        + sum(sqrt((lambda42(:)-lambda42old(:)).^2)))...
        /(sum(sqrt(lambda41old(:).^2))+sum(sqrt(lambda42old(:).^2)));
    
    %----------------------------------------------------------------------
    %----------------------- Residues 1 through 4 -------------------------

    R = sqrt((p1 - FDx(u)).^2 +(p2 - FDy(u)).^2);
    R2(k) = sum(R(:));
    R = sqrt((u-v).^2);
    R3(k) = sum(R(:));
    R = sqrt((n1 - m1).^2+(n2 - m2).^2);
    R4(k) = sum(R(:));
    
    %----------------------------------------------------------------------
    %------------------------------ Energy --------------------------------
    Au = sqrt(FDx(u).^2+FDy(u).^2);
    Ra1 = FDx(u)./Au;
    Ra2 = FDy(u)./Au;
    R = (a + (BDx(Ra1)+BDy(Ra2)).^2*b).*Au;
%     R = (BDx(n1)+BDy(n2)).^2;
    E(k) = sum(R(:));
    R = eta * abs(u(tao)-u0tao);
    E(k) = E(k) + sum(R(:));
    disp = E(k)
    
    %----------------------------------------------------------------------
    %-------------------------- Relative Error ----------------------------
    RE(k) = sum(abs((u(:)-uold(:)))) ...
        /sum(abs(uold(:)));
end
toc
%%
figure(3)
imshow(u)

x = 1:Cal;
figure(4)
plot(x,log10(L1),'r', x,log10(L2),'g',x,log10(L3),'y',x,log10(L4),'b')
legend('L1','L2','L3','L4')

figure(5)
plot(x,log10(R1/(M*N)),'r', x,log10(R2/(M*N)),'g',...
        x,log10(R3/(M*N)),'y',x,log10(R4/(M*N)),'b')
legend('Residue1','Residue2','Residue3','Residue4')

figure(6)
plot(x,log10(E))
legend('Energy')

figure(7)
plot(x,log10(RE))
legend('Relative error in u')

figure(8)
imshow(u0)
hold on
x = size(X,1);
for i = 1:x
    rectangle('Position',[X(i),Y(i),1,1],...
        'FaceColor','r','EdgeColor','r')
end