function potxa(alp)
% Questa funzione mostra l'andamento dell'energia valutata con il 
% metodo di Fourier (Laplaciano Completo) e usando la Formula Semiclassica 
% per V = |x|^alpha
% confronta le due stime di livelli energetici al variare di alpha
%
% potxa senza input mostra l'andamento per diversi alpha e limite buca pot.
% potxa(alp) per selezionare il potenziale con l'alpha desiderato
% si suggerisce di non inserire alp grandi (maggiori di 1e2 circa)

%% Parte preliminare

narginchk(0,1)
if nargin == 0 
    shg
    clf reset
    set(gcf,'numbertitle','off','menu','none', ...
    'name','Potenziali |x|^alpha')

    a = 2;
    h.minus = uicontrol('units','norm','pos',[.38 .01 .06 .05], ...
          'fontsize',12,'string','<','callback','potxa(''a--'')');
    h.alp = uicontrol('units','norm','pos',[.46 .01 .12 .05], ...
          'fontsize',12,'userdata',a,'callback','potxa(''a=2'')');
    h.plus = uicontrol('units','norm','pos',[.60 .01 .06 .05], ...
          'fontsize',12,'string','>','callback','potxa(''a++'')');
    h.close = uicontrol('units','norm','pos',[.80 .01 .10 .05], ...
          'fontsize',12,'string','close','callback','close');
    set(gcf,'userdata',h)
    alp = 'a=2';
end

h = get(gcf,'userdata');

if isnumeric(alp), a = alp; 
else, a = get(h.alp,'userdata'); end

switch alp
   case 'a--', a = a-1;
   case 'a++', a = a+1;
   case 'a=2', a = 2;
end

if ~isnumeric(alp) 
    set(h.alp,'string',['a = ' num2str(a)],'userdata',a);
    if a<2, set(h.minus,'enable','off');
    else, set(h.minus,'enable','on');
    end
    if a>49, set(h.plus,'enable','off');
    else, set(h.plus,'enable','on');
    end
end

%% Parte principale

% confronto tra Fourier Method e Formula Semiclassica
N = 600;
E = Vxtoalp(N,a);

b = 2*a/(a+2);
Es = 0.5*(pi*a/2/beta(3/2,1/a))^b*(1/2:N-1/2)'.^b;

% Grafica
subplot(2,1,1)
plot(E)
hold on
plot(Es)
title('Autovalori potenziali  |x| ^a')
h1 = legend('Complete (E)','Semiclassic (Es)','location','nw');
set(h1,'FontSize',14);
hold off

subplot(2,1,2)
semilogy(abs((E./Es)-1))
title('Differenza relativa  | (E/Es) - 1 |')

% limite a grande --> buca di potenziale
if nargin == 0
    figure(2)
    set(gcf,'numbertitle','off','menu','none', ...
    'name','Limite di buca di potenziale')
    q = linspace(0,50,50);
    E0 = 0.5*(q.^2)*(pi/2)^2;
    plot(q,E0)
    hold on
    for a = 30:10:200
        E = Vxtoalp(60,a);
        plot(E(1:50),'r')
    end
    h2 = legend('Expected','| x | ^a','location','nw');
    set(h2,'FontSize',14);
    title('Limite alpha grande')
end

    function E = Vxtoalp(N,alpha)
        
        x = (-(N-1)/2:(N-1)/2)'; % questa è la grid, passo reticolare 1
        V = 0.5*abs(x).^alpha;

        n = floor(N/2); % produco lo spazio di fourier
        nn = floor((N-1)/2);
        k = (2*pi/N)*(-n:nn)';
        maxT = pi^2/2;  % by construction, this is just pi^2/2
        maxV = 0.5*(max(x))^alpha;  % max x perche prendo l'estremo della scatola
        dx = (maxT/maxV)^(1/(alpha+2)); % "optimal" lattice spacing
        V = V*dx^alpha; % scalo
        k = fftshift(k/dx);

        T = ifft(diag(k.^2/2)*fft(eye(N))); %energia cinetica
        T = real(T+T')/2;
        H = T + diag(V);
        [~,e] = eig(H);
        E = diag(e);

    end
end