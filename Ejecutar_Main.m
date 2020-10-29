clc;
clear;
close all;

%% se declaran los niveles de ruido.
Noise=[40 30 20 10];

%Ciclo que recorre los diferentes niveles de ruido
for jj=1:4
%% Aqui se definen cuantas veces se repite el proceso de estimación de
% los parametros.

Iter=100;

%% Parametros exactos o teoricos del problema inverso.
Var_ent.kappa=0.47;     %Adimensional - Rango de [0,1]

        fprintf('Los parametros teoricos a estimar son: \n'); 
        fprintf('kappa     = %6.2f [adimensional]\n',Var_ent.kappa);

%% Nivel de ruido en las mediciones.
Var_ent.NR=Noise(jj);

%% Método de optimización implementado.
% 0 para el método Trust Region Reflective
% 1 para el método PSO
Method=1;

%% Parametros para correr el programa

Var_ent.nder=30;                % nder: número de nodos en la coordenada radial,
Var_ent.ndet=100;               % ndet: número de nodos en la coordenada angular,
Var_ent.ndt=1000;               % ndt : número de nodos de la variable tiempo,
Var_ent.tf=2*pi;                % tf  : tiempo adimensional de simulación,

%% Crear los datos de medición

[vd_Teo,wd_Teo,rd] = Cyli_WGN_LOW(Var_ent);

%%
for xx=1:Iter

 %% Crear los datos de medición
vd_meas=awgn(vd_Teo,Var_ent.NR,'measured');
wd_meas=awgn(wd_Teo,Var_ent.NR,'measured');

save('Meas','vd_meas','wd_meas','rd');
fprintf('Se han creado los datos de medición con un nivel de ruido= %i dB \n \n',Var_ent.NR);

%%
tic

%Se declara el numero de parametros a encontrar en el problema inverso;
N=1;

%Se declaran las restricciones para el valor de los parametros
x_L=zeros(N,1);  
x_H=1;

%se declara el punto inicial
x0=(x_H-x_L).*rand(N,1)+x_L;

%Se declaran los criterios de parada para el algoritmo de optimización.
iter=1000;
tolx=1e-14;
tolF=1e-14;

if Method==0
    options = optimoptions(@lsqnonlin,'MaxIterations',iter,'StepTolerance',tolx,'MaxFunctionEvaluations',iter*3,...
        'FunctionTolerance',tolF,'Algorithm','trust-region-reflective','Display','iter','UseParallel',false); % default algorithm
    [y,resnorm,residual,exitflag,output]=lsqnonlin(@SE_model_bajos,x0,x_L,x_H,options,Var_ent);

else
    options = optimoptions(@particleswarm,'MaxIterations',iter,...
    'FunctionTolerance',tolF,'Display','none','UseParallel',true); % default algorithm
    [y,fval,exitflag,output]=particleswarm(@(x) RMSE_model_bajos(x,Var_ent),N,x_L,x_H,options);
end

    fprintf('Iteración = %i: \n \n',xx);

    error1=abs(y(1)-Var_ent.kappa)*(100/Var_ent.kappa);
%     error2=abs(y(2)-Var_ent.BmT)*(100/Var_ent.BmT);
    
        fprintf('Los parámetros estimados son: \n');
        fprintf('kappa     = %6.2f [adimensional]\n',y(1));
%         fprintf('B         = %6.1f [mT] \n \n',y(2));
        
        fprintf('El error porcentual para los párametros estimados es: \n');
        fprintf('Error porcentual de kappa     = %6.2f%% \n',error1);
%         fprintf('Error porcentual de B         = %6.1f%% \n \n',error2);
    
time=toc;

kappa_est(xx,jj)=y;
error_kappa(xx,jj)=error1;
Time(xx,jj)=time;

if Method==0
    filenm = sprintf( 'Results_TRR_%i.mat',Var_ent.NR);   
    save(filenm,'kappa_est','error_kappa','Time' );  
else
    filenm = sprintf( 'Results_PSO_%i.mat',Var_ent.NR);   
    save(filenm,'kappa_est','error_kappa','Time' ); 
end

end

Prom_kappa(jj)=mean(kappa_est(:,jj));
desv_kappa(jj)=std(kappa_est(:,jj));
Prom_time(jj)=mean(Time(:,jj));

    fprintf('\n Terminaron las simulaciones para el nivel de ruido = %6.2f [adimensional]\n',Noise(jj));
    fprintf('El promedio de kappa estimado es = %2.4f [adimensional]\n',Prom_kappa(jj));
    fprintf('La desviación estandar de kappa estimado es = %e [adimensional]\n',desv_kappa(jj));
    fprintf('El tiempo de simulacion promedio es = %6.2f [adimensional]\n \n',Prom_time(jj));

end