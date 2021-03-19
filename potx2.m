function potx2
% Questa funzione mostra l'andamento dell'energia valutata con i metodi:
% H tridiagonale, H pentadiagonale, H completa.
% confronto con Energia aspettata per V = x^2

N = 600;
x = (-(N-1)/2:(N-1)/2)'; % questa è la grid, passo reticolare 1
V = 0.5*abs(x).^2;

n = floor(N/2); % produco lo spazio di fourier
nn = floor((N-1)/2);
k = (2*pi/N)*(-n:nn)';
maxT = pi^2/2;  % by construction, this is just pi^2/2
maxV = 0.5*(max(x))^2;  % max x perche prendo l'estremo della scatola
dx = (maxT/maxV)^(1/4); % "optimal" lattice spacing
V = V*dx^2; % scalo

% Metodo del laplaciano tridiagonale 
T = (1/dx^2/2)*(diag(ones(N,1))-diag(ones(N-1,1),1));
T=(T+T');
T(1,end) = T(1,2); T(end,1) = T(2,1);
H = T + diag(V);
[~,e] = eig(H);
E(:,1) = diag(e);

% Metodo del laplaciano pentadiagonale
D0 = diag(ones(N,1));
D1 = diag(ones(N-1,1),1);
D2 = diag(ones(N-2,1),2);
T = (1/dx^2/2)*((5/4)*D0 - (4/3)*D1 + (1/12)*D2);
T = T+T';
T(1,end) = T(1,2); T(end,1) = T(2,1);
T(2,end) = T(1,3); T(end,2) = T(3,1);
T(1,end-1) = T(1,3); T(end-1,1) = T(3,1);
H = T + diag(V);
[~,e] = eig(H);
E(:,2) = diag(e);

% Metodo dell'Hamiltoniana completa con trasformata di Fourier
k = fftshift(k/dx);

T = ifft(diag(k.^2/2)*fft(eye(N))); %energia cinetica
T = real(T+T')/2;
H = T + diag(V);
[~,e] = eig(H);
E(:,3) = diag(e);


Ex = 0.5*(pi*2/2/beta(3/2,1/2))^1*(1/2:N-1/2)'.^1;
subplot(2,1,1)

% Plot e grafica
plot(Ex,'k')
hold on
plot(E)
title('Autovalori potenziale armonico')
h1 = legend('Expected (Ex)','Tridiagonal','Pentadiagonal','Complete','location','nw');
set(h1,'FontSize',14);
hold off

subplot(2,1,2)
semilogy(abs((E./Ex)-1))
title('Differenza relativa | (E/Ex) - 1 |')

set(gcf,'numbertitle','off','menu','none', ...
    'name','Discretizzazione di Hamiltoniana')

end