function schrodes(psi0,V,odesolver)
% schrodes solve schrodinger eq. by using 'Odesolver Method'

% es: schrodes
% es: schrodes(@(x) exp(-x.^2),@(x) x.^2)

%% Plot dell'evoluzione temporale della psi
    function status = schplot(t,psi,~,~)
        if length(psi) ~= 100, return; end % l'ultimo ciclo esce
        set(h,'YData' ,abs(psi).^2);
        str = sprintf('t = %-8.1f',t); 
        E = real((psi'*H*psi)/(psi'*psi));
        deltaE = [deltaE, E-E0];
        str = [str ' ,  ' sprintf('E-E0 = %-+2.2e',deltaE(end))];
        set(titl,'string',str);
        drawnow
        if get(stop,'value'), status = 1; % dopo stop ripeto ultima volta
        else, tau = [tau, t]; status = 0; end % continua a ripeterlo
    end

%% Parte principale

narginchk(0,3)

% Defaults
if nargin == 0
    psi0 = @(x) exp(-(x-4).^2);
end
if nargin < 2, V = @(x)0.*x; end
if nargin < 3, odesolver = @ode113; end

x = linspace(-10,10,1e2);
N = length(psi0(x));
L = 20;
a = L/N;
x = a*(-(N-1)/2:(N-1)/2)';

% Laplaciano completo
k = (pi/a/(N+1))*(1:N)';
T = sinft(diag(k.^2/2)*sinft(eye(N)));
T = real(T+T')/2;
H = T + diag(V(x));

tau = 0; % sarà un vettore con tutti i tempi in cui valuto psi
E0 = real((psi0(x)'*H*psi0(x))/(psi0(x)'*psi0(x))); % energia iniziale
deltaE = 0; % sarà un vettore con tutte le differenze di Energia nel tempo

% Grafica
shg
clf reset
set(gcf,'numbertitle','off','menu','none', ...
    'name','Odesolver Method')
h = plot(x,abs(psi0(x)).^2,'.-');
axis manual
titl = title(sprintf('t = %-8.1f',tau));
stop = uicontrol('style','toggle','string','stop',...
      'units','normalized','position',[.45 .01 .1 .05]);

% Risoluzione
opts = odeset('reltol',1.e-6,'outputfcn', @schplot);
odesolver(@schreqn, [0, inf], psi0(x), opts);

% Plot dell'andamento di E nel tempo
figure(2)
set(gcf,'numbertitle','off','menu','none', ...
    'name','Odesolver Method')
plot(tau,deltaE/E0)
title('(E(t) - E0) / E0')

set(stop,'string','close','value',0,'callback','close(gcf)')

% Eq. di Schrodinger scritta come ODE
function ydot = schreqn(~,y)
    ydot = -1i*H*y;
end

end