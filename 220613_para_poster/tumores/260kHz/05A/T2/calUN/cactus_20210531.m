%%%% cactus - D.G. Actis 31/05/2021

% Toma uno o varios archivos con formato de nombre
% 'xxxkHz_yyydA_zzzMss_nombre.txt' para un solo archivo
% 'xxxkHz_yyydA_zzzMss_nombre*.txt' para varios, con * cualquier caracter.
% xxx la frecuencia en kHz, yyy el valor de corriente IDC en A, zzz el valor 
% de muestreo (cuántos puntos por segundo se han registrado). El valor de
% corriente es empleado para la amplitud del campo. El valor del muestreo
% es usado para la base temporal y el de frecuencia para verificar esa base
% temporal. Además el programa busca los archivos
% 'xxxkHz_yyydA_zzzMss_nombre_fondo.txt'
% 'xxxkHz_yyydA_zzzMss_nombre_cal.txt'

% Se ajustan funciones senoidales a las tres señales de referencia. Se elige
% el tiempo cero para que las tres referencias comiencen en fase. Se
% verifica que las frecuencias coincidan entre sí, y con el nombre del
% archivo.

% Se resta señal de muestra y señal de fondo, y señal de calibración y
% señal de fondo correspondiente. Para eso se dilata el tiempo del fondo
% multiplicándolo por el cociente de frecuencias (fondo/muestra, y
% fondo/calibración en cada caso) por si hay alguna pequeña discrepancia,
% tanto para la muestra como para el paramagneto de calibración. A partir
% de ahí se trabaja únicamente con estas restas y con las referencias de
% muestra y de calibración que contienen la información del campo.

% Se filtra el ruido aislado de cada resta, las opciones son
% 0) No filtrar
% 1) discriminando puntos donde la derivada (o menos la derivada) sea alta
% en comparación al resto de la señal, y sus entornos. En esas regiones se
% ajusta un polinomio con los puntos sin ruido a ambos lados de la zona
% ruidosa (se hace lo propio en las mismas regiones temporales que su señal
% para las respectivas referencias).
% 2) Con un filtro pasabajos de Fourier.

% Se recortan las señales para tener un número entero de períodos, y se
% omiten tanto el primer medio período como el último. De esta manera se
% evita el ruido que no pudo ser tratado adecuadamente al principio y al
% final de cada medida.

% Se promedia resta y referencia sobre todos los períodos y se integra.

% Se lleva el campo a unidades de A/m normalizando y multiplicando por el
% campo máximo medido para el valor correspondiente de IDC. Se grafica el
% ciclo en unidades de V*s para la calibración y se ajusta una recta. Con
% la pendiente obtenida se lleva a unidades de A/m el eje de magnetización
% para la muestra.

% Se guarda un archivo con la imagen del gráfico del ciclo de histéresis,
% y otro con los puntos del ciclo en ascii.

% Se calculan la coercitividad y la remanencia.

% Se calcula el valor de SAR integrando los ciclos.

% Se imprime en un archivo de salida la siguiente información:
% Nombre del archivo, Frecuencia (kHz), Campo máximo (kA/m), SAR (W/g),
% Coercitividad (kA/m), Magnetización Remanente (kA/m), Peor quita de ruido
% porcentual.



function primavera_owon_20210428()

% Cierro todas las figuras
    close all
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input necesario de parte del usuario
            
    % ¿Qué archivos desea abrir?
        todos=1;
        % todos = 0 -> Abre el explorador para elegir los archivos
        % todos = 1 -> Abre todos los archivos en la carpeta donde esté el
        % script, cuyo nombre de archivo termine con el nombre de muestra:
            nombre='T2';

    % ¿Qué gráficos desea ver? (1 = sí, ~1 = no)
        graficos=[0;    % Referencias y sus ajustes 
                  0;    % Referencias de muestra y fondo sincronizadas. Señales de fondo y muestra antes de restarse.
                  0;    % Referencias de calibración y fondo sincronizadas. Señales de fondo y calibración antes de restarse
                  0;    % Resta de señales: muestra-fondo y calibracion-fondo
                  0;    % Filtrado de calibración
                  0;    % Filtrado de muestra
                  0;    % Recorte de períodos enteros calibración
                  0;    % Recorte de períodos enteros muestra
                  0;    % Campo y magnetización normalizados calibración
                  1;    % Ciclos de calibración
                  0;    % Campo y magnetización normalizados muestra
                  1;    % Ciclos M vs. H de todas las medidas
                  0;    % SAR vs. amplitud de campo
                  0];   % SAR vs. amplitud de campo al cuadrado
              
    % ¿Desea filtrar las señales? (0 = No, 1 = Filtro Actis, 2 = Filtro
    %                              Fourier, 3 = Filtro Fourier+Actis)
        filtrarcal = 2;         % Filtro para la calibración
        filtrarmuestra = 2;     % Filtro para la muestra
              
    % ¿Quiere generar una imagen png con cada ciclo M vs. H obtenido? escriba
    % guarda_imagen_ciclo=1. Si no deje 0 o cualquier otro valor.
        guarda_imagen_ciclo=1;
            
    % Masa de nanoparticulas sobre volumen de FF en g/m^3, se utiliza para
    % el cálculo de SAR
        concentracion = 7725;    
    % Permeabilidad magnética del vacío en N/A^2
        mu0=4*pi*10^-7;        
    % Nombre del archivo de salida
        nombre_archivo_salida = 'T2260kHz05A.dat';
    % Texto para los datos del ciclo de salida
        textociclo= '_ciclo.dat';
    % Texto para la imagen del ciclo de salida
        textoimagen='_ciclo.png';
    % Texto que identifica a los archivos de fondo
        textofondo = '_fondo';
    % Texto que identifica a los archivos de fond
        textocal = '_cal';

    % Calibración bobina: constante que dimensionaliza al campo en A/m a
    % partir de la calibración realizada sobre la bobina del RF
        pendiente_HvsI=43.18*79.77;
        ordenada_HvsI=2.73*79.77;
        
    % Susceptibilidad del patrón de calibración
        %rho_bulk_Gd2O3=7.41*10^3; % kg/m^3
        rho_patron_Gd2O3=2*10^3;  % kg/m^3
        xi_bulk_Gd2O3_masa=1.35*10^(-4)*4*pi/10^3; %emu*m/g/A=m^3/kg
        xi_patron_vol=xi_bulk_Gd2O3_masa*rho_patron_Gd2O3;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Abre los archivos
    [fa, aux_archivos_muestra]=abrir_archivos_esar(todos,nombre);

% Acá va a ir guardando el nombre de cada archivo tratado
    nombres_archivos=cell(1,fa);
% Matriz que saldrá impresa al final en el archivo de salida    
    data_salida=zeros(fa,6);
    

% Barre sobre todos los archivos    
for k=1:fa

    % Nombre del archivo de muestra
        filenombre_muestra=char(aux_archivos_muestra(k));
        nombres_archivos(k)={filenombre_muestra};
    % Nombre del archivo de fondo
        filenombre_fondo=strcat(filenombre_muestra(1:end-4),textofondo,filenombre_muestra(end-3:end));
    % Nombre del archivo de gadolinio
        filenombre_cal=strcat(filenombre_muestra(1:end-4),textocal,filenombre_muestra(end-3:end));
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Parámetros de la medida a partir de nombre del archivo de muestra
            
    % Corriente IDC en el generador de RF. Esto viene del archivo y
    % si está mal nomenclado habrá un error en la escala de campo
        Idc=str2double(filenombre_muestra(8:10))/10;
    
	% Frecuencia del nombre del archivo. Sirve para comparar con la
	% frecuencia ajustada. Si difieren en más de un porcentaje frena
        frecuencia_nombre=str2double(filenombre_muestra(1:3))*1000;
        
    % Base temporal
        deltat=1e-6/str2double(filenombre_muestra(14:16));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Glosario de variables

        % t: tiempo
        % v: voltaje
        % f: fondo
        % m: muestra
        % c: calibración
        % r: referencia
        
    % Lectura de datos
        % Tiempo, señal y referencia de fondo (bobinas vacías)
            [t_f, v_f, v_r_f]=medida_cruda(filenombre_fondo,deltat);
        % Tiempo, señal y referencia de calibración (paramagneto)
            [t_c, v_c, v_r_c]=medida_cruda(filenombre_cal,deltat);
        % Tiempo, señal y referencia de muestra cuyo SAR se calculará
            [t_m, v_m, v_r_m]=medida_cruda(filenombre_muestra,deltat);
                              
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Ajusta las tres referencias con funciones seno:
            [~,frecuencia_f,fase_f,ajuste_f]=ajusta_seno(t_f,v_r_f);
            [valor_medio_c,frecuencia_c,fase_c,ajuste_c]=ajusta_seno(t_c,v_r_c);
            [valor_medio_m,frecuencia_m,fase_m,ajuste_m]=ajusta_seno(t_m,v_r_m);

        % Comparacion de las referencias y sus ajustes
            if graficos(1)==1
                figure('Name','Ajustes')
                plot(t_m,v_r_m,'o',t_f,v_r_f,'o',t_c,v_r_c,'o',t_m,ajuste_m,t_f,ajuste_f,t_c,ajuste_c)
                title('Comparacion de referencias y ajustes')
                legend('Referencia de fondo','Referencia de muestra','Referencia de la calibración','Ajuste de referencia de fondo','Ajuste de referencia de muestra','Ajuste de referencia de la calibración')
                xlabel('Tiempo (s)')
                ylabel('Referencia (V)')            
            end
                        
        % Si la diferencia entre ambas frecuencias es mayor al 2 %, error.
            if abs(frecuencia_m-frecuencia_f)/frecuencia_f>0.1||...
                    abs(frecuencia_c-frecuencia_f)/frecuencia_f>0.1||...
                    abs(frecuencia_c-frecuencia_m)/frecuencia_m>0.1||...
                    abs(frecuencia_m-frecuencia_nombre)/frecuencia_m>0.1
                fprintf('Error: incompatibilidad de frecuencias\n')
                fprintf('La frecuencia de muestra es %.2f kHz\n',frecuencia_m)
                fprintf('La frecuencia de fondo es %.2f kHz\n',frecuencia_f)
                fprintf('La frecuencia de paramagneto es %.2f kHz\n',frecuencia_c)
                fprintf('La frecuencia nomenclada es %.2f kHz\n',frecuencia_nombre)
                break
            end   
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Resta medida y fondo interpolando para que correspondan ambos
    % tiempos. No trabaja más con la medidas individuales, sólo con las
    % restas
    
    % Muestra:
    
        if graficos(2)==1
            [t_m1,v_r_m1,Resta_m]=resta_inter(t_m,v_m,v_r_m,fase_m,frecuencia_m,...
            valor_medio_m,t_f,v_f,v_r_f,fase_f,frecuencia_f,'muestra');
        else
            [t_m1,v_r_m1,Resta_m]=resta_inter(t_m,v_m,v_r_m,fase_m,frecuencia_m,...
            valor_medio_m,t_f,v_f,v_r_f,fase_f,frecuencia_f,0);
        end
        
    % Calibración:
    
        if graficos(3)==1
            [t_c1,v_r_c1,Resta_c]=resta_inter(t_c,v_c,v_r_c,fase_c,frecuencia_c,...
            valor_medio_c,t_f,v_f,v_r_f,fase_f,frecuencia_f,'calibración');
        else
            [t_c1,v_r_c1,Resta_c]=resta_inter(t_c,v_c,v_r_c,fase_c,frecuencia_c,...
            valor_medio_c,t_f,v_f,v_r_f,fase_f,frecuencia_f,0);
        end
        
    % Gráfico de restas de señales de muestra o calibración menos fondo
            if graficos(4)==1
                figure('Name','Resta de señales')
                plot(t_m1,Resta_m,'o-',t_c1,Resta_c,'o-')
                title('Resta de señales')
                legend('Resta de medida de muestra y medida de fondo',...
                    'Resta de medida de calibración y medida de fondo')
                xlabel('Tiempo (s)')
                ylabel('Resta (V)')
            end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Filtrado
        % Calibración
        if graficos(5)==1
            [t_c2,v_r_c2,Resta_c2]=filtrando_ruido(t_c1,v_r_c1,Resta_c,filtrarcal,'calibración');
        else
            [t_c2,v_r_c2,Resta_c2]=filtrando_ruido(t_c1,v_r_c1,Resta_c,filtrarcal,0);
        end
    
        % Muestra
        if graficos(6)==1
            [t_m2,v_r_m2,Resta_m2]=filtrando_ruido(t_m1,v_r_m1,Resta_m,filtrarmuestra,'muestra');
        else
            [t_m2,v_r_m2,Resta_m2]=filtrando_ruido(t_m1,v_r_m1,Resta_m,filtrarmuestra,0);
        end
       length(Resta_m)   
       length(Resta_m2)
        % Diferencia entre señal sin ruido y señal. Guarda el peor valor.
%             dif_resta_m=Resta_m2-interp1(t_m1,Resta_m,t_m2,'spline');
%             dif_resta_c=Resta_c2-interp1(t_c1,Resta_c,t_c2,'spline');
%             peor_diferencia=max([mean(abs(dif_resta_m))/max(Resta_m),mean(abs(dif_resta_c))/max(Resta_c)]);
 
       

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Recorta un número entero de períodos o ciclos
        % Calibración
        if graficos(7)==1
            [t_c3,v_r_c3,Resta_c3,N_ciclos_c]=recorte(t_c2,v_r_c2,Resta_c2,frecuencia_c,'calibracion');
        else
            [t_c3,v_r_c3,Resta_c3,N_ciclos_c]=recorte(t_c2,v_r_c2,Resta_c2,frecuencia_c,0);
        end
    
        % Muestra
        if graficos(8)==1
            [t_m3,v_r_m3,Resta_m3,N_ciclos_m]=recorte(t_m2,v_r_m2,Resta_m2,frecuencia_m,'muestra');
        else
            [t_m3,v_r_m3,Resta_m3,N_ciclos_m]=recorte(t_m2,v_r_m2,Resta_m2,frecuencia_m,0);
        end

            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Campo y magnetización: se integran las funciones de referencia y
    % resta.
    
        % Últimos ajustes de las referencias
            [~,frecuencia_final_m,~,~]=ajusta_seno(t_m3,v_r_m3);
            [~,frecuencia_final_c,~,~]=ajusta_seno(t_c3,v_r_c3);

        % Promedia los ciclos antes de integrar.
            [tiempo_final_m,Rm,fem_campo_m,delta_t_m]=promediado_ciclos(t_m3,Resta_m3,v_r_m3,frecuencia_final_m,N_ciclos_m);
            [tiempo_final_c,Rc,fem_campo_c,delta_t_c]=promediado_ciclos(t_c3,Resta_c3,v_r_c3,frecuencia_final_c,N_ciclos_c);
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Calibración
        % Calcula las sumas acumuladas y convierte a campo y magnetizacion
            % Constante que dimensionaliza al campo en A/m a partir de la
            % calibración realizada sobre la bobina del RF
                C_norm_campo=Idc*pendiente_HvsI+ordenada_HvsI;
        
            % Integral de la fem inducida, es proporcional al campo más una
            % constante
                campo_ua0_c = delta_t_c*cumtrapz(fem_campo_c); 
            % Se le resta una constante para estar centrado en cero
                %campo_ua_c = campo_ua0_c-(max(campo_ua0_c)+min(campo_ua0_c))/2;
                campo_ua_c = campo_ua0_c-mean(campo_ua0_c);
            % Da unidades de A/m al campo
                campo_c=campo_ua_c*C_norm_campo/max(campo_ua_c);

            % Integral de la resta de fems inducidas, restado también el
            % fondo. Es proporcional a la magnetización más una constante
                magnetizacion_ua0_c=delta_t_c*cumtrapz(Rc);      
            % Se le resta una constante para estar centrado en cero
                magnetizacion_ua_c=magnetizacion_ua0_c-mean(magnetizacion_ua0_c);


            
            if graficos(9)==1
                figure('Name','Campo y magnetización normalizados del paramagneto')
                plot(tiempo_final_c,campo_ua_c/max(campo_ua_c),tiempo_final_c,magnetizacion_ua_c/max(magnetizacion_ua_c),'LineWidth',1.5)
                title('Campo y magnetización normalizados calibración')
                legend('Campo','Magnetizacion')
                xlabel('Tiempo (s)')
                ylabel('u.a')
            end
            

            % Ajusta una recta al ciclo de la calibración
                recta=polyfit(campo_c,magnetizacion_ua_c,1);
                pendiente=recta(1);
                
            if graficos(10)==1    
                figure(21)
                hold on
                plot(campo_c,magnetizacion_ua_c,'-.',campo_c,recta(1).*campo_c+recta(2),'LineWidth',1.5,'Color',[k/fa,1-k/fa,1-k/fa])
                grid on
                hold off
                title('Ciclos de calibración en función del campo')
                legend('Ciclos','Recta ajustada')
                xlabel('H (Oe)')
                ylabel('u.a')
            end
            
        % Calibración para pasar la magnetización a A/m
            calibracion=xi_patron_vol/pendiente;
            
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Muestra

        % Calcula las sumas acumuladas y convierte a campo y magnetización
        
            % Integral de la fem inducida, es proporcional al campo más una
            % constante
                campo_ua0_m = delta_t_m*cumtrapz(fem_campo_m); 
            % Centrado en cero
                campo_ua_m = campo_ua0_m-mean(campo_ua0_m);

            % Integral de la resta de fems inducidas, restado también el
            % fondo. Es proporcional a la magnetización más una constante
                magnetizacion_ua0_m=delta_t_m*cumtrapz(Rm);      
            % Se le resta una constante para estar centrado en cero
                %magnetizacion_ua_m=magnetizacion_ua0_m-(max(magnetizacion_ua0_m)+min(magnetizacion_ua0_m))/2;
                magnetizacion_ua_m=magnetizacion_ua0_m-mean(magnetizacion_ua0_m);

            % Da unidades de A/m a la magnetizacion final utilizando la
            % proporcionalidad obtenida en la calibracion. Este paso podría
            % ser realizado directamente sobre un valor de calibracion
            % proveido por el usuario, de calibraciones anteriores.
                magnetizacion_m=calibracion.*magnetizacion_ua_m;

            if graficos(11)==1
                figure('Name','Campo y magnetización normalizados de la muestra')
                plot(tiempo_final_m,campo_ua_m/max(campo_ua_m),tiempo_final_m,magnetizacion_m/max(magnetizacion_m),'LineWidth',1.5)
                title('Campo y magnetización normalizados muestra')
                legend('Campo','Magnetizacion')
                xlabel('Tiempo (s)')
                ylabel('u.a')
            end
            
             
            % Campo y magnetización finales
                
                campo_m=campo_ua_m/max(campo_ua_m)*C_norm_campo;

            % Grafica el ciclo momentáneamente para guardar una imagen
                if guarda_imagen_ciclo==1
                    titulo=['Ciclo de histéresis del archivo ' filenombre_muestra];
                    fig=figure(400);
                    plot(campo_m,magnetizacion_m,'.','LineWidth',1.5,'Color',[k/fa,1-k/fa,1-k/fa])
                    title(titulo, 'Interpreter', 'none')
                    grid on
                    legend('Magnetizacion vs Campo')
                    xlabel('Campo (A/m)')
                    ylabel('Magnetizacion (A/m)')
                    drawnow
                    frame = getframe(fig);
                    im = frame2im(frame);
                    close(fig)

                % Guarda imagen png del ciclo
                     archivo_imagen=strcat(filenombre_muestra(1:end-4),textoimagen);
                     imwrite(im,archivo_imagen)
                end

            % Va graficando todos los ciclos juntos
                if graficos(12)==1
                    figure(20)
                    hold on
                    plot(campo_m,magnetizacion_m,'-.','LineWidth',1.5,'Color',[k/fa,1-k/fa,1-k/fa])
                    grid on
                    title('Ciclo de histéresis')
                    legend('Magnetizacion vs Campo')
                    xlabel('Campo (A/m)')
                    ylabel('Magnetizacion (A/m)')
                    hold off
                end

        % Exporta ciclo en ascii
            ciclo_out=[(tiempo_final_m-tiempo_final_m(1))', campo_m',magnetizacion_m'];
            filenombre_out=strcat(filenombre_muestra(1:end-4),textociclo);
            fid = fopen(filenombre_out,'wt');
            fprintf(fid, '%s\t %s\t %s\n','Tiempo (s)','Campo (A/m)','Magnetizacion (A/m)');
            fclose(fid);
            dlmwrite(filenombre_out,ciclo_out,'delimiter','\t','precision','%6e','-append')
            
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Cálculo de campo coercitivo y remanencia

        m = magnetizacion_m;
        h = campo_m;

        clear Hc
        ihc=1;
        imr=1;
        Hc=zeros(2,1);
        Mr=zeros(2,1);
        for i=1:length(m)-1
            if (  m(i)>0 && m(i+1)<0  || m(i)<0 && m(i+1)>0 )
                 Hc(ihc)=h(i) - m(i)* ( h(i+1)-h(i) ) / (m(i+1)-m(i));
                 ihc = ihc+1;
            end
            if (  h(i)>0 && h(i+1)<0  || h(i)<0 && h(i+1)>0 )
                 Mr(imr)=m(i) - h(i)* ( m(i+1)-m(i) ) / (h(i+1)-h(i));
                 imr = imr+1;
            end
        end

        Hc_mean = mean(abs(Hc));
        %Hc_error = std(abs(Hc));

        Mr_mean = mean(abs(Mr));
        %Mr_error = std(abs(Mr));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Determinación de áreas

    % Arma un vector auxilar desplazado a valores positivos
        magnetizacion_des = magnetizacion_m+2*abs(min(magnetizacion_m));
        
    % Área del lazo en A^2/m^2
        area=-trapz([campo_m, campo_m(1)],[magnetizacion_des, magnetizacion_des(1)]);
 
    % Cálculo de potencia específica disipada SAR
        SAR=mu0*abs(area)*frecuencia_final_m/concentracion;
                
        % Se grafica el SAR en función del campo máximo y en función del campo
        % máximo al cuadrado 
            if graficos(13)==1
                figure('Name','SAR vs. amplitud de campo')
                hold on
                plot(max(campo_m)/1000,SAR,'.-','MarkerSize',20)
                hold off
                title('SAR vs Campo Maximo')
                xlabel('Campo Maximo (kA/m)')
                ylabel('SAR (W/g)')
            end
            
            if graficos(14)==1
                figure('Name','SAR vs. amplitud de campo al cuadrado')
                hold on
                plot((max(campo_m)/1000).^2,SAR,'.','MarkerSize',20)
                hold off
                title('SAR vs Campo Maximo^2')
                xlabel('Campo Maximo^2 (kA/m)^2')
                ylabel('SAR (W/g)')
            end
            
            % Información impresa en el archivo de salida para cada medida
            % tratada
                data_salida(k,:)=[frecuencia_final_m/1000, C_norm_campo/1000, SAR, Hc_mean/1000, Mr_mean/1000, calibracion];

end


% Exporta archivo de salida
   % Genera el archivo 
        salidaID = fopen(nombre_archivo_salida,'wt');
    % Encabezado del archivo de salida
        fprintf(salidaID, '%s\t %s\t %s\t %s\t %s\t %s\t %s\n','Nombre del archivo','Frecuencia (kHz)','Campo maximo (kA/m)','SAR (W/g)', 'Coercitividad (kA/m)','Magnetizacion remanente (kA/m)', 'Pendiente de calibracion');
    % Data a imprimir   
        salida=[nombres_archivos', num2cell(data_salida)];
        formato = '%s\t %s\t %s\t %s\t %s\t %s\t %s\n';
        [nfilas,~] = size(salida);

    % Imprime los datos en el archivo
        for fila = 1:nfilas
            fprintf(salidaID,formato,salida{fila,:});
        end
        fclose(salidaID);


% Función que identifica el ruido
function [t3, marcador]=encuentra_ruido(t,v,ancho,entorno)
    % Toma una señal (t,v) y calcula la derivada de v respecto de t y su
    % valor absoluto medio "ruido_tranqui". Marca los puntos con derivada en
    % valor absoluto mayor que "ancho" veces "ruido_tranqui" y un entorno de
    % estos puntos igual a "entorno" puntos para cada lado.

        %   Suaviza con un promedio leve
            WindowSize = 5;
            be = (1/WindowSize)*ones(1,WindowSize);
            t1=t(WindowSize+1:end)-WindowSize*(t(2)-t(1))/2;
            v_fe=filter(be,1,v);
            v1=v_fe(WindowSize+1:end);

        %   Calcula la derivada de v respecto de t
            derivada=diff(v1)./diff(t1);
            t2=t1(1:end-1)+(t1(2)-t1(1))/2;

        %   Suaviza la derivada
            t3=t2(WindowSize+1:end)-WindowSize*(t2(2)-t2(1))/2;
            derivada0=filter(be,1,derivada);
            derivada2=derivada0(WindowSize+1:end);
    
        %   El ruido característico de la señal es el valor medio del valor
        %   absoluto de la derivada
            ruido_tranqui=mean(abs(derivada2));
            aux1=zeros(1,length(derivada2));

%         %   Ajusta una función armónica a la señal suavizada.
%                 suavep=fftsmooth(v,30);
%                 %   Valor medio de la señal, es decir constante aditiva
%                     parp(1)=mean(suavep);
%                 %   Amplitud
%                     parp(2)=(max(suavep)-min(suavep))/2;
%                 %   Frecuencia. Para evitar problemas por errores de nomenclatura en los archivos mido tiempo entre picos:
%                     [~, indicesp]=findpeaks(suavep);
%                     tiempo_entre_picosp = mean(diff(t(indicesp)));
%                     parp(3)=1/tiempo_entre_picosp;
%                 %   Fase inicial, a partir del tiempo del primer máximo
%                     parp(4)=2*pi*parp(3)*t(indicesp(1))-pi/2;
%                 %   Ajusto y obtengo los coeficientes.
%                     ajustep=fitnlm(t,suavep, @(b,x)(b(1) + b(2)*sin(2*pi*b(3)*x - b(4))), parp);
%                     coefCIp=ajustep.Coefficients.Estimate;
%                     coeficientesp=coefCIp(:,1);
%                     frec=coeficientesp(3);
%                     amp=coeficientesp(2);
%                     fase=coeficientesp(4);
%         %   Resta
%             derivada3=derivada2-amp*2*pi*frec*cos(2*pi*frec*t3-fase);
              derivada3=derivada2;
                  
        %   Marca los puntos que superan en ancho veces el ruido normal
            for jq=1:length(derivada3)
                if abs(derivada3(jq))>ancho*ruido_tranqui
                    aux1(jq)=1;
                else
                    aux1(jq)=0;
                end
            end
            
        %   Prepara el marcador    
            marcador=zeros(1,length(derivada2));
            
        %   Si hay un solo cambio de signo en la derivada, no lo marca, ya
        %   que debe tratarse de un pico natural de la señal

            for jq=entorno+1:length(derivada3)-entorno
                if max(aux1(jq-entorno:jq+entorno))==1
                    marcador(jq+round(entorno/2))=1;
                else
                    marcador(jq+round(entorno/2))=0;
                end
            end
            
        %   Acomodo los extremos
            for jq=1:entorno
                if marcador(entorno+1)==1
                    marcador(jq)=1;
                end
                if marcador(length(derivada3)-entorno)==1
                    marcador(length(derivada3)+1-jq)=1;
                end
            end
    end        
    
    
	function [s_data_v] = fftsmooth(data_v,freq_n)
        % filtra la transformada discreta de Fourier del vector data_v,
        % anula todaslas frecuencias por encima de freq_n y antitransforma
        
        fft_data_v = fft(data_v);
        s_fft_data_v = zeros(1,length(data_v));
        s_fft_data_v(1:freq_n) = fft_data_v(1:freq_n);
        s_fft_data_v(end-freq_n:end) = fft_data_v(end-freq_n:end);
        s_data_v = real(ifft(s_fft_data_v));
    end

    function [fa, aux_archivos_muestra]=abrir_archivos_esar(todos,muestra)
        % Separamos archivo individual vs. tomar todos los archivos
        if not(todos==1)
            % Abrimos los archivos
                texto=strcat('Seleccione los archivos con las medidas de la muestra//',version);
                [filevis,~] = uigetfile('*.*',texto,'MultiSelect', 'on'); % busca el/los archivo(s) a abrir

                if iscell(filevis)            % Si se abrieron varios archivos
                    archivos=char(filevis);   % fa es el numero de archivos abiertos
                    [fa,~]=size(archivos);    
                else
                    fa=1;
                    archivos=filevis;
                end

            % Nombres de los archivos de muestra
                aux_archivos_muestra = cellstr(archivos);
        else
            % Mira el directorio donde está el script y ve cuántos archivos hay
                projectdir = pwd;
                dinfo = dir(fullfile(projectdir));
                dinfo([dinfo.isdir]) = [];     
                nfiles = length(dinfo);

                contador=0;
                
            % Prepara una variable auxiliar con el maximo tamaño posible
                auxiliar_grande=cell(nfiles,1);
            % Para cada archivo mira que contenga un guión bajo
                for arch = 1 : nfiles
                    aux2=dinfo(arch).name;
                    if ~isempty(regexp(aux2,['\d\d\dkHz_\d\d\ddA_\d\d\dMss_',muestra,'\.[dt][ax]t'], 'once'))||...
                          ~isempty(regexp(aux2,['\d\d\dkHz_\d\d\ddA_\d\d\dMss_',muestra,'.\.[dt][ax]t'], 'once'))
                        auxiliar_grande(contador+1,:)=cellstr(aux2);
                        contador=contador+1;
                    end
                end
                aux_archivos_muestra=auxiliar_grande(1:contador,:);
                fa = contador;
        end
    end

    function [t, v, vr]=medida_cruda(nombre_de_archivo,intervalo_temporal)
        % Carga señal y referencia
            archivo=importdata(nombre_de_archivo);
            datos=archivo.data;
        % Identifica las dos señales de cada canal como señal y referencia
        % Recibe datos en mV y pasa a V.
            t   = intervalo_temporal*(datos(:,1)-datos(1,1));
            v   = datos(:,2)*0.001;
            vr  = datos(:,3)*0.001;
    end

    function [valor_medio,frecuencia,fase,ajuste]=ajusta_seno(t,vr)
        % Ajusta los datos vr(tr) con una función seno
        % V(t)=V0+A*sin(2*pi*f*tr-phi)        
            % Estima valores iniciales para los parámetros de ajuste
                p0=zeros(4,1);
            % Valor medio, es decir constante aditiva
                p0(1)=mean(vr);
            % Amplitud
                p0(2)=(max(vr)-min(vr))/2;
            % Para la frecuencia mide tiempo entre picos
                vr_suave=fftsmooth(vr,round(length(vr)*3/1000));
                [~, indices_picos]=findpeaks(vr_suave);
                tiempo_entre_picos = mean(diff(t(indices_picos)));
                p0(3)=1/tiempo_entre_picos;
            % Fase inicial, a partir del tiempo del primer máximo
                p0(4)=2*pi*p0(3)*t(indices_picos(1))-pi/2;
            % Ajusta y obtiene los coeficientes
                ajuste_todo=fitnlm(t,vr, @(b,x)(b(1) + b(2)*sin(2*pi*b(3)*x - b(4))), p0);
                coefCI_cualq=ajuste_todo.Coefficients.Estimate;
                coeficientes=coefCI_cualq(:,1);
            % Recupera valor medio, frecuencia y fase, y vector ajuste
                valor_medio=coeficientes(1);
                frecuencia=coeficientes(3);
                fase=coeficientes(4);
                ajuste=ajuste_todo.Fitted;
    end

    function [t1,vr1,Resta]=resta_inter(t,v,vr,fase,frecuencia,offset,tf,vf,vrf,fasef,frecuenciaf,graf)
        
        % Desplazamiento temporal para poner en fase las referencias, y
        % resta de valores medios de referencias.

        % Saca el modulo 2 pi de las fases y calcula el tiempo de fase 0
        
            t_fase=mod(fase,2*pi)/(frecuencia*2*pi);
            t_fase_f=mod(fasef,2*pi)/(frecuenciaf*2*pi);
            
        % Desplaza en tiempo para que haya coincidencia de fase entre
        % referencias. La amplitud debería ser positiva siempre por el
        % valor inicial del parámetro de ajuste.
 
            t1 = t-t_fase;
            tf_aux= tf - t_fase_f;
            
        % Corrije por posible diferencia de frecuencias, dilatando el
        % tiempo del fondo. Crea dos tiempos de fondo, uno que lidia con
        
            tf_mod=tf_aux*frecuenciaf/frecuencia;

        % Resta el offset de la referencia
        
            vr1 = vr-offset;
            
        % Resta medida y fondo interpolando para que corresponda mejor el
        % tiempo. No trabaja más con la medidas individuales, sólo con la
        % resta. Se toma la precaución para los casos de trigger distinto
        % entre fondo y medida. Comienza en fase 0 o sea t=0.
            tmin=0;
            tmax=tf_mod(end);
        % Recorta el tiempo a mirar    
            t_aux=t1(t1>=tmin & t1<=tmax);
        % Interpola el fondo a ese tiempo
            interpolacion_aux=interp1(tf_mod,vf,t_aux,'linear');
            interpolacion=zeros(size(v));
            for indice_temporal=1:length(t1)
            % Cambia índice viejo por índice nuevo en la base original
                [~, indice_mas_cercano]=min(abs(t_aux-t1(indice_temporal)));
                interpolacion(indice_temporal)=interpolacion_aux(indice_mas_cercano);
            end
        % Resta
            Resta = v-interpolacion;
        
        % Comparacion de las referencias de muestra y fondo desplazadas en
        % tiempo y offset
            if graf~=0
                figure('Name',['Referencias sincronizadas, ',graf])
                plot(tf_mod,vrf,'o',t1,vr1','o')
                title('Comparacion de referencias desplazadas y restados sus offsets')
                legend('Referencia de fondo',strcat('Referencia de ',graf))
                grid on
                xlabel('Tiempo (s)')
                ylabel('Referencia (V)')
                
                figure('Name',strcat('Señales antes de restarse, ',graf))
                plot(tf_mod,vf,'-o',t1,v,'-o',t1,interpolacion,'-o')
                title('Comparación de medidas')
                legend('Medida de fondo',strcat('Medida de ',graf),strcat('Interpolación del fondo al tiempo de la ',graf))
                xlabel('Tiempo (s)')
                ylabel('Señal (V)')
            end
            
    end
    
    function [t2,vr2,v2]=filtrando_ruido(t,vr,v,filtrar,graf)
        if filtrar==0
            t2=t';
            vr2=vr';
            v2=v';
        elseif filtrar==2 || filtrar==3
        % Filtro por Fourier la calibración
        freq=round(length(vr)/5);
        v2=fftsmooth(v,freq)';
        vr2=fftsmooth(vr,freq)';
        t2=t';

                if graf~=0
                    % Control de que el suavizado final sea satisfactorio
                    figure('Name',['Quitando ruido de la ',graf])
                    plot(t,vr/max(vr2),'o',t2,vr2/max(vr2),...
                        t,v/max(v2)+3,'o',t2,v2/max(v2)+3,t,v+3)
                    title('Filtro de Fourier')
                    legend('Referencia de muestra','Sin ruido','Resta de señales','Sin ruido','Zona de ruido')
                    xlabel('Tiempo (s)')
                    ylabel('u.a.')
                    axis([t(1) t(end) -1.1 4.1])
                end
        elseif filtrar==1
            % Identifica el ruido 

            % Factor del ruido natural a partir del cual se considera que
            % hay que filtrar
                ancho=2.5;
            % Puntos a ambos lados del sector de señal ruidosa que serán
            % incluidos en el ruido.
                entorno=5;

            % Marcador de regiones a filtrar
                [t2, marcador]=encuentra_ruido(t',v',ancho,entorno);

            % Ajuste ruido
            % Parámetros del ajuste: puntos a cada lado de la region a
            % filtrar que serán considerados para el ajuste, y grado del
            % polinomio a ajustar.
                puntos_ajuste=80;
                grado_pol=3;

            % Tiempos y señales en las mismas dimensiones que los marcadores.
                vr2=interp1(t,vr,t2,'spline');
                v2=interp1(t,v,t2,'spline');

                % Comienza a filtrar a partir de que tiene suficientes
                % puntos detrás
                    w=puntos_ajuste+1;

                % Barre la señal. No filtra ni el principio ni el
                % final, por eso más adelante se eliminan el primer y
                % el último semiperíodos.
                    while w<length(v2)
                        % Si no hay ruido deja la señal como estaba
                        if marcador(w-1)==0
                            w=w+1;
                        % Si hay ruido
                        elseif marcador(w-1)==1
                            q=w;
                            % Busca hasta donde llega el ruido
                            while marcador(q-1)==1 && q<length(v2)
                                q=q+1;
                            end
                            % Si tiene suficientes puntos del otro lado
                            % realiza el ajuste
                            if q<length(v2)-puntos_ajuste
                                y=[v2(w-puntos_ajuste:w);v2(1+q:1+q+puntos_ajuste)];
                                x=[t2(w-puntos_ajuste:w);t2(1+q:1+q+puntos_ajuste)];
                                p=polyfit(x,y,grado_pol);
                                v2(w:q-1)=polyval(p,t2(w:q-1));
                                y=[vr2(w-puntos_ajuste:w);vr2(1+q:1+q+puntos_ajuste)];
                                x=[t2(w-puntos_ajuste:w);t2(1+q:1+q+puntos_ajuste)];
                                p=polyfit(x,y,grado_pol);
                                vr2(w:q-1)=polyval(p,t2(w:q-1));                                    
                            end
                            w=q;
                        end
                    end
                    
                    if graf~=0
                        % Control de que el suavizado final sea satisfactorio
                        figure('Name',strcat('Quitando ruido de la',{' '},graf))
                        plot(t,vr/max(vr2),'o',t2,vr2/max(vr2),t2,marcador,...
                            t,v/max(v2)+3,'o',t2,v2/max(v2)+3,t,v+3)
                        title('Filtro ad hoc')
                        legend('Referencia de muestra','Sin ruido','Zona de ruido','Resta de señales','Sin ruido','Zona de ruido')
                        xlabel('Tiempo (s)')
                        ylabel('u.a.')
                        axis([t(1) t(end) -1.1 4.1])
                    end
        end
        
        if filtrar==3
            v=v2';
            vr=vr2';
            t=t2';
            % Identifica el ruido 

            % Factor del ruido natural a partir del cual se considera que
            % hay que filtrar
                ancho=2.5;
            % Puntos a ambos lados del sector de señal ruidosa que serán
            % incluidos en el ruido.
                entorno=5;

            % Marcador de regiones a filtrar
                [t2, marcador]=encuentra_ruido(t',v',ancho,entorno);

            % Ajuste ruido
            % Parámetros del ajuste: puntos a cada lado de la region a
            % filtrar que serán considerados para el ajuste, y grado del
            % polinomio a ajustar.
                puntos_ajuste=80;
                grado_pol=3;

            % Tiempos y señales en las mismas dimensiones que los marcadores.
                vr2=interp1(t,vr,t2,'spline');
                v2=interp1(t,v,t2,'spline');

                % Comienza a filtrar a partir de que tiene suficientes
                % puntos detrás
                    w=puntos_ajuste+1;

                % Barre la señal. No filtra ni el principio ni el
                % final, por eso más adelante se eliminan el primer y
                % el último semiperíodos.
                    while w<length(v2)
                        % Si no hay ruido deja la señal como estaba
                        if marcador(w-1)==0
                            w=w+1;
                        % Si hay ruido
                        elseif marcador(w-1)==1
                            q=w;
                            % Busca hasta donde llega el ruido
                            while marcador(q-1)==1 && q<length(v2)
                                q=q+1;
                            end
                            % Si tiene suficientes puntos del otro lado
                            % realiza el ajuste
                            if q<length(v2)-puntos_ajuste
                                y=[v2(w-puntos_ajuste:w);v2(1+q:1+q+puntos_ajuste)];
                                x=[t2(w-puntos_ajuste:w);t2(1+q:1+q+puntos_ajuste)];
                                p=polyfit(x,y,grado_pol);
                                v2(w:q-1)=polyval(p,t2(w:q-1));
                                y=[vr2(w-puntos_ajuste:w);vr2(1+q:1+q+puntos_ajuste)];
                                x=[t2(w-puntos_ajuste:w);t2(1+q:1+q+puntos_ajuste)];
                                p=polyfit(x,y,grado_pol);
                                vr2(w:q-1)=polyval(p,t2(w:q-1));                                    
                            end
                            w=q;
                        end
                    end
                    
                    if graf~=0
                        % Control de que el suavizado final sea satisfactorio
                        figure('Name',strcat('Quitando ruido de la',{' '},graf))
                        plot(t,vr/max(vr2),'o',t2,vr2/max(vr2),t2,marcador,...
                            t,v/max(v2)+3,'o',t2,v2/max(v2)+3,t,v+3)
                        title('Filtro ad hoc')
                        legend('Referencia de muestra','Sin ruido','Zona de ruido','Resta de señales','Sin ruido','Zona de ruido')
                        xlabel('Tiempo (s)')
                        ylabel('u.a.')
                        axis([t(1) t(end) -1.1 4.1])
                    end
        end
        
    end
    
    function [t2,vr2,v2,N_ciclos]=recorte(t,vr,v,frecuencia,graf)
        % Recorta un número entero de períodos o ciclos,arrancando en fase
        % cero (campo máximo o mínimo según polaridad)
    
        % Calculo del numero de ciclos enteros
            N_ciclos = floor(t(end)*frecuencia);
            
        
        % Indices de los ciclos
            indices_recorte=find(0<=t & t<N_ciclos/frecuencia);
            
        % Si quisiera recortar medio ciclo a cada lado
            %largo=indices_recorte(end)-indices_recorte(1);
            %if mod(largo,N_ciclos)==0
            %elseif mod(largo,N_ciclos)<=0.5
            %    largo=largo-mod(largo,N_ciclos);
            %else
            %    largo=largo+N_ciclos-mod(largo,N_ciclos);
            %end
        
            %indices_recorte=ceil(largo/N_ciclos/2):1:ceil(largo-largo/N_ciclos/2);
            %N_ciclos=N_ciclos-1;

        % Se recortan los vectores
            t2   = t(indices_recorte);
            vr2 = vr(indices_recorte);
            v2   = v(indices_recorte);
            
        % Referencia de muestra sin valor medio y su ajuste,y resta de
        % señal de muestra y fondo, recortadas para que tengan un numero
        % entero de ciclos
            if graf~=0
                figure('Name',['Numero entero de períodos, ',graf])
                plot(t2,vr2,'o',t2,400*v2,'o')
                title('Referencia y resta de señales. N entero de periodos')
                legend('Referencia','Resta reescaleada x400')
                xlabel('Tiempo (s)')
                ylabel('Señal (V)') 
            end
   
    end

    function [tf,vf,vrf,delta_t]=promediado_ciclos(t,v,vr,frecuencia,N_ciclos)
            tf=t(t<=t(1)+1/frecuencia);
            vrf=zeros(size(tf));
            vf=zeros(size(tf));
            for ind_prom=1:N_ciclos
                if tf(end)+(ind_prom-1)/frecuencia<t(end)
                    vrf=vrf+interp1(t,vr,tf+(ind_prom-1)/frecuencia)/N_ciclos;
                    vf=vf+interp1(t,v,tf+(ind_prom-1)/frecuencia)/N_ciclos;
                else
                    vrf=vrf+interp1(t,vr,tf+(ind_prom-1)/frecuencia,'spline')/N_ciclos;
                    vf=vf+interp1(t,v,tf+(ind_prom-1)/frecuencia,'spline')/N_ciclos;
                end
            end
         
            
        % Le quita el valor medio a la resta y al fondo
            vf = vf-mean(vf);
            vrf = vrf-mean(vrf);            
        % Paso temporal
            delta_t=(tf(end)-tf(1))/length(tf);

    end
end
