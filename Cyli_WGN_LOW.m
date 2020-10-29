function [vd,wd,rd]=Cyli_WGN_LOW(Var_ent)
%Cyli_WGN: generates profiles from the Internal Angular Momentum Diffusion Theory (IAMDT).
% This program performs the velocity profiles generation of a ferrofluid sample (WBF1) 
% immersed in a cylindrical container of infinite axial length, and under the 
% effect of a low amplitude rotating magnetic field, with the addition of SNR (10, 15 , 20, 30 and 40 dB).
%
% Cyli_WGN( nder , ndet , ndt , tf ,maxIterT )
%
% nder: radial nodes.
% ndet: angular nodes.
% ndt : time nodes.
% tf  : simulation dimensionless time.
% maxIter: maximum iterations.
%
% Author: Cristian Camilo Jiménez Leiva
% Copyright © 2020 UNIVERSIDAD INDUSTRIAL DE SANTANDER 
% Any comments: cristian.jimenez@correo.uis.edu.co


%% Se definene los parametros para correr el programa

nder=Var_ent.nder;                
ndet=Var_ent.ndet;               
ndt=Var_ent.ndt;               
tf=Var_ent.tf;                

%% CARACTERIZACIÓN DE MUESTRA DE FERROFLUIDO WBF1
FF.name = 'WBF1'; 
FF.eta=1.03e-3;                   % Viscosidad de cizalla
FF.eta0=1.02e-3;                  % Viscosidad del liquido portador de las nanopartículas
FF.ji=0.106;                      % susceptibilidad magnética inicial
FF.phi=2.13e-3;                   % Fracción volumétrica
FF.tau=1.67e-5;                   % Tiempo de relajación (Browniano) [s]
FF.kappa=Var_ent.kappa;                     % VARIABLE QUE RELACIONA A VARIOS PARÁMETROS DEL FERROFLUIDO -- ESTA ES UNA DE LAS VARIABLES QUE QUEREMOS OBTENER CON EL PROBLEMA INVERSO.

%% FF.kappa=sqrt( (4*FF.eta*FF.R0^2*FF.zita ) / ( FF.eta_e*FF.eta' ) )   DEFINICIÓN DE LA VARIABLE FF.kappa. FF.eta' se define como el "spin viscosity" o viscosidad de giro.

FF.u0=4*pi*10^-7;                 % Permeabilidad del aire o vacío
FF.om=150;                        % Frecuencia de rotación del campo magnético [Hz]
FF.omTau=2*pi*FF.om*FF.tau;             % Frecuencia adimensional
FF.zita=1.5*FF.phi*FF.eta0;             % Viscosidad de vórtice
FF.eta_e=FF.eta+FF.zita;                % Constante = FF.eta + FF.zita
FF.R0=24.7e-3;                    % Radio del cilindro REPORTADO POR TORRES-DIAZ.

%% VALOR DEL PASO EN LA COORDENADA RADIAL, ANGULAR Y EN EL TIEMPO    
her=1/(nder-1);                % Paso en coordenada radial
het=2*pi/ndet;                 % Paso en coordenada angular
htt=tf/(ndt-1);                % Paso en coordenada temporal

%% CONSTANTES DEL SISTEMA

FF.lz0=FF.omTau/(1+FF.omTau^2);         % Torque análitico de orden cero.
FF.Md=425e3;                      % Magnetización de los dominios magnéticos
FF.T=294;                         % Temperatura del sistema
FF.d=14.3e-9;                     % Diámetro de las partículas magnéticas
FF.Kb=1.38064852e-23;             % Constante de Boltzmann



%% INICIALIZACIÓN DE VARIABLES DEL PROGRAMA
r=zeros(nder,1);  teta=zeros(ndet,1); t=zeros(ndt,1);  % Inicialización de variables
HrP=zeros(nder,ndt);HtP=zeros(nder,ndt);MrP=zeros(nder,ndt);MtP=zeros(nder,ndt);lzP=zeros(nder,ndt); cont=zeros(nder,1);
g1=zeros(ndt,1);g2=zeros(ndt,1);g3=zeros(ndt,1);g4=zeros(ndt,1);g5=zeros(ndt,1); Mr=zeros(nder,ndet,ndt); Mt=zeros(nder,ndet,ndt); lzm=zeros(1,nder);
    
%% ASIGNACIÓN DE VALORES DE LA COORDENADA RADIAL, ANGULAR Y DEL TIEMPO DE SIMULACIÓN EN EL CONTENEDOR CILÍNDRICO    
    for k=1:nder
        r(k)=her*(k-1);%+(R1/FF.R0); % Se establece el dominio de los puntos del espacio en la coordenada radial r=0 a r=1.
    end
    
    for k=1:ndet
        teta(k)=het*(k-1); % Se establece el dominio de los puntos del espacio en la coordenada tetha tetha=0 a tetha=2pi.
    end
    
    for k=1:ndt
        t(k)=htt*(k-1); % Se establece el dominio de los puntos del espacio en la coordenada tiempo t=0 a t=tf.
    end
      
%% VARIABLES QUE SE VAN A USAR EN LA SOLUCIÓN DEL PROBLEMA HIDRODINÁMICO

B=FF.eta_e/(FF.eta*her^2);           % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
Cx=FF.eta_e/(FF.eta*2*her);          % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
D=-FF.eta_e/FF.eta;                  % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
E=-FF.zita/(her*FF.eta);             % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
F=4*FF.eta/(FF.eta_e*FF.kappa^2);       % Constante de la discretización de las ecuaciones diferenciales hidrodinámicas
   
    
%% EN ESTAS LINEAS SE LE ESPECIFICA AL PROGRAMA EL VALOR DE DENSIDAD DE CAMPO MAGNÉTICO QUE SE VAN A REALIZAR LAS SIMULACIONES

BmT=0.01;                             % VALOR DE LA DENSIDAD DE CAMPO --- -- ESTA ES UNA DE LAS VARIABLES QUE QUEREMOS OBTENER CON EL PROBLEMA INVERSO.

alf=zeros(1,length(BmT));                              % INCIALIZACIÓN DEL PARÁMETRO DE LANVEGIN

filenm = sprintf( 'Datos_%s.mat',FF.name);   
save( filenm , 'FF' );                                 % NOMBRA EL ARCHIVO .mat EN EL QUE SE GUARDAN LOS PARÁMETROS DEL FERROFLUIDO EN ESTA SIMULACIÓN.

for campo=1:length(BmT)
    
    K=BmT(campo)*1e-3/FF.u0;                                                   % VALOR DE INTENSIDAD DE CAMPO MAGNETICO DE LA ITETACIÓN 'campo'
    
    alf(campo)=pi/6*(FF.u0*FF.Md*FF.d^3*K)/(FF.Kb*FF.T);                                   % PARÁMETRO DE LANGEVIN
    epsi=FF.u0*FF.ji*K^2*FF.tau/FF.zita;                                                % PARÁMETRO DE PERTURBACIÓN
    

    fprintf('\n')
    disp('PARÁMETROS CON LOS QUE SE ESTÁ CORRIENDO EL PROGRAMA:'), fprintf('\n')
    disp(FF.name);fprintf('\n')
    fprintf('kappa     = %7.2f [adimensional]\n',FF.kappa);
    fprintf('f         = %7.0f [Hz] \n',FF.om);
    fprintf('omTau     = %7.3f [adimensional] \n',FF.omTau);
    fprintf('B         = %7.3f [mT] \n',K*FF.u0/1e-3);
    fprintf('her       = %7.4f [adimensional] \n',her);
    fprintf('het       = %7.4f [adimensional] \n',het);
    fprintf('htt       = %7.4f [adimensional] \n',htt);
    fprintf('alpha     = %7.3f [adimensional] \n',alf(campo));
    fprintf('epsi      = %7.5f [adimensional] \n',epsi);%,pause
    fprintf('lz0       = %7.4f [adimensional] \n\n\n',FF.lz0);%pause
     
    
    [v, w]=cilindrico_analitica_bessel(nder,FF.eta,FF.eta0,FF.phi,FF.tau,FF.kappa,FF.om);                     % CALCULA PERFILES ADIMENSIONALES v Y w
        
    
    vd = v.*FF.u0*FF.ji*K^2*FF.omTau*FF.R0/FF.zita*1000;                                                   % ¡¡¡ESTA ES LA SALIDA vd_OUT!!!
    wd = w.*FF.u0*FF.ji*K^2*FF.omTau/(FF.zita);                                                            % ¡¡¡ESTA ES LA SALIDA wd_OUT!!!                
        rd = r.*FF.R0*1000;                                                                        
                                                   

end

end

  
% % %% GRAFICA: vd_OUT y vd10_OUT
% %     
% %     % SNR: 10dB
% %     
% %     vd10_OUT = awgn(vd_OUT,10,'measured');                                                                    % ¡¡¡ESTA ES LA SALIDA vd10_OUT!!!      
% %     wd10_OUT = awgn(wd_OUT,10,'measured');                                                                    % ¡¡¡ESTA ES LA SALIDA wd10_OUT!!!
% %     
% %     hoja1=figure(1);                                                            
% %     subplot(2,1,1)
% %     h=plot(rd,vd_OUT,'r',rd,vd10_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2); 
% %     hold on, grid on, title([FF.name,' -- B \approx 0 mT ','-- \kappa = ',num2str(FF.kappa),' --  f = ',num2str(FF.om),' Hz -- 10dB']);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');
% %     legend('vd\_OUT','vd10\_OUT','Location','south'), xlim([0 FF.R0*1000]);
% %     
% %     subplot(2,1,2)                                                        
% %     h=plot(rd,wd_OUT,'r',rd,wd10_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);
% %     hold on, grid on, set(h(1),'linewidth',2);set(h(2),'linewidth',2);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold'); xlim([0 FF.R0*1000]);
% %     legend('wd\_OUT','wd10\_OUT');pause(.1)
% %     
% %     set(hoja1,'Units','Inches');
% %     pos=get(hoja1,'position');
% %     set(hoja1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% %     pdfnm = sprintf( 'Perfil_%s_%s_10dB',FF.name,Bnm); % NOMBRA EL ARCHIVO .pdf DE LA GRÁFICA DEL RUIDO PARA ESTE FERROFLUIDO Y PARA ESTE VALOR DE B.
% %     print(hoja1,pdfnm,'-dpdf','-r0')
% %     
% % 
% %     
% % %% GRAFICA: vd_OUT y vd15_OUT
% %     
% %     vd15_OUT = awgn(vd_OUT,15,'measured');                                                                    % ¡¡¡ESTA ES LA SALIDA vd15_OUT!!!      
% %     wd15_OUT = awgn(wd_OUT,15,'measured');                                                                    % ¡¡¡ESTA ES LA SALIDA wd15_OUT!!!
% %     
% %     hoja5=figure(5);                                                            
% %     subplot(2,1,1)
% %     h=plot(rd,vd_OUT,'r',rd,vd15_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2); 
% %     hold on, grid on, title([FF.name,' -- B \approx 0 mT ','-- \kappa = ',num2str(FF.kappa),' --  f = ',num2str(FF.om),' Hz -- 15dB']);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');
% %     legend('vd\_OUT','vd15\_OUT','Location','south'), xlim([0 FF.R0*1000]);
% %     
% %     subplot(2,1,2)                                                        
% %     h=plot(rd,wd_OUT,'r',rd,wd15_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);
% %     hold on, grid on, set(h(1),'linewidth',2);set(h(2),'linewidth',2);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold'); xlim([0 FF.R0*1000]);
% %     legend('wd\_OUT','wd15\_OUT');pause(.1)
% %     
% %     set(hoja5,'Units','Inches');
% %     pos=get(hoja5,'position');
% %     set(hoja5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% %     pdfnm = sprintf( 'Perfil_%s_%s_15dB',FF.name,Bnm); % NOMBRA EL ARCHIVO .pdf DE LA GRÁFICA DEL RUIDO PARA ESTE FERROFLUIDO Y PARA ESTE VALOR DE B.
% %     print(hoja5,pdfnm,'-dpdf','-r0')
% %     
% % %% GRAFICA: vd_OUT y vd20_OUT
% %     
% %     vd20_OUT = awgn(vd_OUT,20,'measured');                                                                    % ¡¡¡ESTA ES LA SALIDA vd20_OUT!!!    
% %     wd20_OUT = awgn(wd_OUT,20,'measured');                                                                    % ¡¡¡ESTA ES LA SALIDA wd20_OUT!!!
% %     
% %     hoja2=figure(2);                                                            
% %     subplot(2,1,1)
% %     h=plot(rd,vd_OUT,'r',rd,vd20_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2); 
% %     hold on, grid on, title([FF.name,' -- B \approx 0 mT ','-- \kappa = ',num2str(FF.kappa),' --  f = ',num2str(FF.om),' Hz -- 20dB']);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');
% %     legend('vd\_OUT','vd20\_OUT','Location','south'), xlim([0 FF.R0*1000]);
% %     
% %     subplot(2,1,2)                                                         
% %     h=plot(rd,wd_OUT,'r',rd,wd20_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);
% %     hold on, grid on, set(h(1),'linewidth',2);set(h(2),'linewidth',2);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold'); xlim([0 FF.R0*1000]);
% %     legend('wd\_OUT','wd20\_OUT');pause(.1)
% %     
% %     set(hoja2,'Units','Inches');
% %     pos=get(hoja2,'position');
% %     set(hoja2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% %     pdfnm = sprintf( 'Perfil_%s_%s_20dB',FF.name,Bnm); % NOMBRA EL ARCHIVO .pdf DE LA GRÁFICA DEL RUIDO PARA ESTE FERROFLUIDO Y PARA ESTE VALOR DE B.
% %     print(hoja2,pdfnm,'-dpdf','-r0')
% %     
% % %% GRAFICA: vd_OUT y vd30_OUT
% %     
% %     vd30_OUT = awgn(vd_OUT,30,'measured');                                                                   % ¡¡¡ESTA ES LA SALIDA vd30_OUT!!!  
% %     wd30_OUT = awgn(wd_OUT,30,'measured');                                                                   % ¡¡¡ESTA ES LA SALIDA wd30_OUT!!!
% %     
% %     hoja3=figure(3);                                                            
% %     subplot(2,1,1)
% %     h=plot(rd,vd_OUT,'r',rd,vd30_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);
% %     hold on, grid on, title([FF.name,' -- B \approx 0 mT ','-- \kappa = ',num2str(FF.kappa),' --  f = ',num2str(FF.om),' Hz -- 30dB']);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');
% %     legend('vd\_OUT','vd30\_OUT','Location','south'), xlim([0 FF.R0*1000]);
% %     
% %     subplot(2,1,2)                                                         
% %     h=plot(rd,wd_OUT,'r',rd,wd30_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);
% %     hold on, grid on, set(h(1),'linewidth',2);set(h(2),'linewidth',2);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold'); xlim([0 FF.R0*1000]);
% %     legend('wd\_OUT','wd30\_OUT');
% %     
% %     set(hoja3,'Units','Inches');
% %     pos=get(hoja3,'position');
% %     set(hoja3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% %     pdfnm = sprintf( 'Perfil_%s_%s_30dB',FF.name,Bnm); % NOMBRA EL ARCHIVO .pdf DE LA GRÁFICA DEL RUIDO PARA ESTE FERROFLUIDO Y PARA ESTE VALOR DE B.
% %     print(hoja3,pdfnm,'-dpdf','-r0')
% %     
% % %% GRAFICA: vd_OUT y vd40_OUT
% %     
% %     vd40_OUT = awgn(vd_OUT,40,'measured');                                                                   % ¡¡¡ESTA ES LA SALIDA vd40_OUT!!!  
% %     wd40_OUT = awgn(wd_OUT,40,'measured');                                                                   % ¡¡¡ESTA ES LA SALIDA wd40_OUT!!!
% %     
% %     hoja4=figure(4);                                                            
% %     subplot(2,1,1)
% %     h=plot(rd,vd_OUT,'r',rd,vd40_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);
% %     hold on, grid on, title([FF.name,' -- B \approx 0 mT ','-- \kappa = ',num2str(FF.kappa),' --  f = ',num2str(FF.om),' Hz -- 40dB']);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('v [mm/s]','FontWeight','bold');
% %     legend('vd\_OUT','vd40\_OUT','Location','south'), xlim([0 FF.R0*1000]);
% %     
% %     subplot(2,1,2)                                                         
% %     h=plot(rd,wd_OUT,'r',rd,wd40_OUT,'k'); set(h(1),'linewidth',2,'MarkerSize',2);set(h(2),'linewidth',2);
% %     hold on, grid on, set(h(1),'linewidth',2);set(h(2),'linewidth',2);
% %     xlabel('r [mm]','FontWeight','bold');ylabel('\omega [rad/s]','FontWeight','bold'); xlim([0 FF.R0*1000]);
% %     legend('wd\_OUT','wd40\_OUT');
% %     
% %     set(hoja4,'Units','Inches');
% %     pos=get(hoja4,'position');
% %     set(hoja4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% %     pdfnm = sprintf( 'Perfil_%s_%s_40dB',FF.name,Bnm); % NOMBRA EL ARCHIVO .pdf DE LA GRÁFICA DEL RUIDO PARA ESTE FERROFLUIDO Y PARA ESTE VALOR DE B.
% %     print(hoja4,pdfnm,'-dpdf','-r0')
% % 
% % 
% % %% EXPORTACIÓN DE DATOS A HOJA DE CÁLCULO (EXCEL)
% % 
% % grupo = [rd vd_OUT' wd_OUT' vd10_OUT' wd10_OUT' vd15_OUT' wd15_OUT' vd20_OUT' wd20_OUT' vd30_OUT' wd30_OUT' vd40_OUT' wd40_OUT'];
% % datos = num2cell(grupo);
% % Encabezados = {'RADIO','V','W','V10dB','W10dB','V15dB','W15dB','V20dB','W20dB','V30dB','W30dB','V40dB','W40dB'};
% % Tabla = [Encabezados; datos]; 
% % xlsxnm = sprintf( 'Tabulacion_%s_%s',FF.name,nombreTabla );
% % xlswrite(xlsxnm,Tabla)
% % 
% % end


