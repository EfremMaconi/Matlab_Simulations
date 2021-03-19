function out = schructure(in)
% schructure solve schrodinger eq. by using differents method
% Boundary Condition: 0 (DBC) or 1 (PBC)
% Method can be:  'OpS' (operator splitting) . 'Com' (complete Hamiltonian)
%                 'Tri' (tridiagonal H) . 'Pen' (pentadiagonal H)

% es: a = schructure; schructure(a) (Default)

% es: a.psi0 = @(x)exp(-(x+1).^2); 
%     a.V = @(x)exp(-(2*(x-2)).^100); (Gradino)

% es: a.psi0 = @(x)(2*x.^2 -1).*exp(-0.5*x.^2); 
%     a.V = @(x)0.5*x.^2; (Autostato pot. armonico)

% es: a.psi0 = @(x)exp(-x.^2); a.V = @(x)50/sqrt(2*pi)*exp(-(50*(x-3)).^2);
% a.dt=0.001; (Approssimazione della Delta)

% es: a.psi0 = @(x)exp(-x.^2); 
%     a.V = @(x)-2*exp(-(x/4).^500)+2; a.dt=0.01; (Buca finita)

narginchk(0,1)

%% Set defaults

dflt.psi0 = @(x) exp(-(x-4).^2); 
dflt.BC = 0; 
dflt.V = @(x) 0.*x; 
dflt.method = 'OpS';

dflt.N = 1024;
dflt.dt = 0.05;

if nargin == 0
    out = dflt;
    return;
end

%% Setup

N = in.N;
L = 20;
a = L/N;
x = a * ( -(N-1)/2 : (N-1)/2 )' ;
dt = in.dt;

% Condizioni al Bordo
switch in.BC
   case 0, k = (pi/a/(N+1))*(1:N)'; % DBC
   case 1, n = floor(N/2); nn = floor((N-1)/2);
           k = (2*pi/N)*(-n:nn)'; k = fftshift(k/a); % PBC
   otherwise
        error('Boundary Condition must be: 0 (DBC) or 1 (PBC)')
end

psi = in.psi0(x);
V = in.V(x);

% Parte grafica
h = plot(x,abs(psi).^2,'.-');
hold on 
axis manual
set(gcf,'numbertitle','off','menu','none', ...
    'name','Risoluzione numerica eq. di Schroedinger dipendente dal tempo')
stop = uicontrol('style','toggle','string','stop',...
      'units','normalized','position',[.45 .01 .1 .05]);

% Studio dei diversi metodi
switch in.method
    
  case 'Com' % Metodo esponenziale di matrice Hamiltoniana completa 

    switch in.BC
        case 0, T = sinft( diag(k.^2/2)*sinft(eye(N)) ); % DBC
        case 1, T = ifft( diag(k.^2/2)*fft(eye(N)) ); % PBC
    end
    T = real(T+T')/2;
    H = T + diag(V);
    U = expm(-1i*dt*H); %operatore evoluzione temporale
    
  case 'Tri' % Metodo esponenziale di matrice Hamiltoniana tridiagonale
        
    T = (1/a^2/2)*(diag(ones(N,1))-diag(ones(N-1,1),1)); % DBC
    T = T+T';
    if in.BC == 1 % PBC
        T(1,end) = T(1,2);
        T(end,1) = T(1,2);
    end
    H = T + diag(V);
    U = expm(-1i*H*dt);
    
  case 'Pen' % Metodo esponenziale di matrice Hamiltoniana pentadiagonale
        
    D0 = diag(ones(N,1));
    D1 = diag(ones(N-1,1),1);
    D2 = diag(ones(N-2,1),2);
    T = (1/a^2/2)*((5/4)*D0 - (4/3)*D1 + (1/12)*D2); % DBC
    T = T+T';
    if in.BC == 1 % PBC
        T(1,end) = T(1,2);    T(end,1) = T(1,2);
        T(1,end-1) = T(1,3);  T(2,end) = T(1,3);
        T(end-1,1) = T(1,3);  T(end,2) = T(1,3);
    end
    H = T + diag(V);
    U = expm(-1i*H*dt);
    
  case 'OpS' % Metodo Operator Splitting
      
    switch in.BC
        case 0, U = diag(exp(-1i*V*dt/2))*...
                sinft( diag(exp(-1i*(0.5*k.^2)*dt))*...
                sinft(eye(N)) )* diag(exp(-1i*V*dt/2)); % DBC
        case 1, U = diag(exp(-1i*V*dt/2))*...
                ifft( diag(exp(-1i*(0.5*k.^2)*dt))*...
                fft(eye(N)) )* diag(exp(-1i*V*dt/2));  % PBC
    end
    H = (1i/dt)*logm(U);
    
   otherwise
        error('method must be one of these: Com, Tri, Pen, OpS')

end

E0 = real((psi'*H*psi)/(psi'*psi));
if E0 < min(V), warning('E0 < min(V)'), end

% Scalo graficamente V ed E0 per rendere visivamente intuitivo il sistema
if max(V) > E0 
    pot = V - min(V);
    e0 = E0 - min(V);
    e0 = e0*max(psi)/max(pot);
    pot = pot*max(psi)/max(pot);
else
    pot = V - min(V);
    e0 = 0.8*max(psi);
    pot = 0.8*pot/(E0-min(V));
end

plot(x,pot)
plot(x,(e0*diag(ones(length(x)))'))

j = 0; % indice (t = j*dt)
titl = title(sprintf('t = %-8.1f',j));

% Evoluzione della psi e andamento dell'energia
while ~get(stop,'value')
    j = j+1;
    E(j) = real((psi'*H*psi)/(psi'*psi));
    psi = U*psi;
    set(h,'ydata',abs(psi).^2);
    drawnow
    str = sprintf('t = %-8.1f',j*dt); 
    str = [str ' ,  ' sprintf('(E-E0)/E0 = %-+2.2e',E(end)-E0)];
    set(titl,'string',str);
end

hold off

% Plotto l'andamento di E(t)
figure(2)
set(gcf,'numbertitle','off','menu','none', ...
    'name','Risoluzione numerica eq. di Schroedinger dipendente dal tempo')
plot((1:j).*dt,(E-E0)/E0)
title('(E(t) - E0) / E0')

set(stop,'string','close','value',0,'callback','close(gcf)')

end
