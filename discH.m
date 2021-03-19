function discH(arg)
% Questa funzione mostra differenti metodi di discretizzazione di
% Hamiltoniane plottando l'energia degli autostati della buca di potenziale
%
% discH; senza input mostra l'andamento per diversi N
% discH(arg) fa una discretizzazione con N = arg 
% si suggerisce di non inserire arg grandi (maggiori di 4e3 circa)

%% Parte preliminare

narginchk(0,1)

if nargin == 0
    shg
    clf reset
    set(gcf,'numbertitle','off','menu','none', ...
    'name','Discretizzazione di Hamiltoniana')

    N = 256;
    h.minus = uicontrol('units','norm','pos',[.38 .01 .06 .05], ...
          'fontsize',12,'string','<','callback','discH(''N/2'')');
    h.n = uicontrol('units','norm','pos',[.46 .01 .12 .05], ...
          'fontsize',12,'userdata',N,'callback','discH(''N=256'')');
    h.plus = uicontrol('units','norm','pos',[.60 .01 .06 .05], ...
          'fontsize',12,'string','>','callback','discH(''N*2'')');
    h.close = uicontrol('units','norm','pos',[.80 .01 .10 .05], ...
          'fontsize',12,'string','close','callback','close');
    set(gcf,'userdata',h)
    arg = 'N=512';
end

h = get(gcf,'userdata');
if isnumeric(arg), N = arg; 
else, N = get(h.n,'userdata'); end

switch arg
   case 'N/2', N = N/2;
   case 'N*2', N = N*2;
   case 'N=256', N = 256;
end

if ~isnumeric(arg) 
    set(h.n,'string',['N = ' num2str(N)],'userdata',N);
    if N<257, set(h.minus,'enable','off');
    else, set(h.minus,'enable','on');
    end
    if N>2047, set(h.plus,'enable','off');
    else, set(h.plus,'enable','on');
    end
end

%% Parte principale

L = 2;
a = L/N;

% Andamento teorico dell'energia
E0 = 0.5*(pi*(1:500)/a/(N+1)).^2;
subplot(2,1,1)
plot(E0,'k')
hold on


% Formula approssimante equialente a metodo del laplaciano tridiagonale
% Eq = 0.5*((2/a)*sin((pi*q)/(2*(N+1)))).^2;

% Metodo del laplaciano tridiagonale 
H0 = (1/a^2/2)*(diag(ones(N,1))-diag(ones(N-1,1),1));
H0=(H0+H0');
[~,e] = eig(H0);
E(:,1) = diag(e);

% Metodo del laplaciano pentadiagonale
D0 = diag(ones(N,1));
D1 = diag(ones(N-1,1),1);
D2 = diag(ones(N-2,1),2);
H0 = (1/a^2/2)*((5/4)*D0 - (4/3)*D1 + (1/12)*D2);
H0 = H0+H0';
[~,e] = eig(H0);
E(:,2) = diag(e);

% Metodo della trasformata di Fourier (Completo)
k = (pi/a/(N+1))*(1:N)';
H0 = sinft(diag(k.^2/2)*sinft(eye(N)));
H0 = real(H0+H0')/2;
[~,e] = eig(H0);
E(:,3) = diag(e);
if N < 500
    plot(E)
else
    plot(E(1:500,:))
end

% Grafica
title('Autovalori particella libera')
h1 = legend('Expected','Tridiagonal','Pentadiagonal','Fourier','location','nw');
set(h1,'FontSize',14);
hold off

subplot(2,1,2)
m = min(N,500);
semilogy(1:m,abs(E0(1:m) - E(1:m,:)')./E0(1:m))
title('Differenza Relativa | (E/Ex) - 1 |')
set(h1,'FontSize',14);

end
