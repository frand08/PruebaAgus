%% Calculo de la Transformada Discreta de Karhunen-Loève
sc=senoidal(5,1,0,0,64,128);
sc=sc + senoidal(10,1,0,0,64,128);
% M=128;
n=randn(1,64);
s=sc+n;

    %% Estimacion de matriz de correlacion
   tic
    rx = xcorr(s,'biased');          % Calculo de la funcion de autocorrelacion/N
    rx = rx( length(rx)/2:end );
    
    Rx = toeplitz(rx);               % Matriz circulante tij = t(i+1) (j+1), como convercion de vector a matriz
    
    [V,D] = eig(Rx);                 % Calculo de valores propios
    [D,I] = sort(diag(D),'descend'); % Ordenamiento de autovalores
    V = V(:,I);                      % Ordenamiento en matriz de autovectorestovectores
    
    S=s*V;
toc
    S(:end)=0;
    Si=S*V';

    figure(1);
    % se�al, transformada, transformada inversa
    subplot(411);plot(s);
    title('signal');
    subplot(412);plot(S);
    title('KLT');
    subplot(413);hold on;plot(Si,'g');plot(sc,'r');
    title('signal RED + KLT truncated');
    subplot(414);plot(sc-Si);
    title(['diferencia con varianza ' num2str(var(sc-Si))]);