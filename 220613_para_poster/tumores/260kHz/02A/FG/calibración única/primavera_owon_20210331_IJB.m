
%%%% Primavera owon - D.G. Actis 18/03/2020


% Toma un archivo con formato de nombre 'xxxkHz_yydA_zzzMss_nombre.dat' con xxx
% la frecuencia en kHz, yy el valor de corriente IDC en A, zzz el valor 
% de muestreo, cuántos puntos por segundo se han registrado y ccccc la concentración en g/m^3. El valor de 
% frecuencia sólo es utilizado como semilla para ajustar la frecuencia
% real, y para verificar que no haya ocurrido un error al nombrar los
% archivos. El de corriente es empleado para verificar que el campo
% obtenido a partir de la amplitud sea coherente. El valor del muestreo es
% usado para la base temporal. Además el programa busca los archivos
% 'xxxkHz_yydA_zzzMss_fondo.dat' y 'xxxkHz_yydA_zzzMss_cal.dat'


% Se ajustan funciones senoidales a las tres señales de referencia. Se elige
% el tiempo cero para que las tres referencias comiencen en fase. También
% se dilata el tiempo de señal y referencia de fondo multiplicándolo por el
% cociente de frecuencias, para el caso en que haya alguna pequeña
% discrepancia, tanto para la muestra como para el gadolinio.

% Se resta señal de muestra y señal de fondo, y señal de calibración y
% señal de fondo correspondiente. A partir de ese punto se
% trabaja únicamente con estas restas y con las referencias.

% Se filtra el ruido aislado de cada resta, discriminando puntos donde la
% derivada (o menos la derivada) sea alta en comparación al resto de la
% señal, y sus entornos. En esas regiones se ajusta un polinomio con
% los puntos sin ruido a ambos lados de la zona ruidosa. Se hace lo propio
% en las mismas regiones temporales que su señal para las respectivas
% referencias.

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



function ciclosRF=ciclosRF(x)

% Cierro todas las figuras
    close all
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input necesario de parte del usuario

    % Masa de nanoparticulas sobre volumen de FF en g/m^3, se utiliza para
    % el cálculo de SAR
        concentracion = 10000;
    % ¿Quiero quitarle ruido a la muestra? 2=FOURIER 1=SI 0=NO
        FILTRARMUESTRA = 2;
    % Permeabilidad del vacío en N/A^2
        mu0=4*pi*10^-7;        
    % Nombre del archivo de salida
        nombre_archivo_salida = 'prueba_ESAR.dat';
    % Texto que identifica a los archivos de fondo
        textofondo = '_fondo.txt';
    % Texto que identifica a los archivos de fondo
        textocal = '_cal.txt';
    % Texto que identifica a los archivos de gadolinio 
    % Importante: los caracteres 8, 9 y 10 del nombre del archivo deben ser
    % el valor numérico (en deciamperes) de Idc, y los 12, 13 y 14 los
    % correspondientes al muestreo.
    
    % Calibración bobina: constante que dimensionaliza al campo en A/m a
    % partir de la calibración realizada sobre la bobina del RF
        pendiente_HvsI=43.18*79.77;
        ordenada_HvsI=2.73*79.77;
        
    % Susceptibilidad del patrón de calibración
        rho_bulk_Gd2O3=7.41*10^3; % kg/m^3
        rho_patron_Gd2O3=2*10^3;  % kg/m^3
        xi_bulk_Gd2O3_masa=1.35*10^(-4)*4*pi/10^3; %emu*m/g/A=m^3/kg
        xi_patron_vol=xi_bulk_Gd2O3_masa*rho_patron_Gd2O3;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    nombres_archivos={};
% Matriz que saldrá impresa al final en el archivo de salida    
    data_salida=zeros(fa,6);
    

% Barre sobre todos los archivos    
for k=1:fa

    % Nombre del archivo de muestra
        filenombre_muestra=char(aux_archivos_muestra(k));
    % Nombre del archivo de fondo
        filenombre_fondo=strcat(filenombre_muestra(1:end-4),textofondo);
    % Nombre del archivo de gadolinio
        filenombre_cal=strcat(filenombre_muestra(1:end-4),textocal);

    
        nombres_archivos(k)={filenombre_muestra};


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Lectura de datos

        % Carga señal y referencia de muestra
            arch_muestra=importdata(filenombre_muestra);
            muestra=arch_muestra.data;
        % Carga señal y referencia de fondo
            arch_fondo=importdata(filenombre_fondo);
            fondo=arch_fondo.data;
        % Carga señal y referencia de calibración
            arch_cal=importdata(filenombre_cal);
            cal=arch_cal.data;
            
    % Corriente IDC en el generador de RF. Esto viene del archivo y
    % si está mal nomenclado habrá un error
        Idc=str2double(filenombre_muestra(8:9));
        
    % Concentración de nanopartículas en gramo sobre metro cúbico
%         concentracion=str2double(filenombre_muestra(21:25))
    % Base temporal
        deltat=1e-6/str2double(filenombre_muestra(12:14));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Glosario de variables

        % t: tiempo
        % v: voltaje
        % f: fondo
        % m: muestra
        % c: calibración
        % r: referencia

   
    % Identifica las dos señales de cada canal como señal y referencia
    % para fondo, muestra y calibración. Recibe datos en mV y acá pasa a
    % V.
            % Fondo, es decir bobinas en ausencia de muestra
            t_f   = deltat*(fondo(:,1)-fondo(1,1));
            v_f   = fondo(:,2)*0.001;
            v_r_f = fondo(:,3)*0.001;
            
            % Calibración, material paramagnético, sin histéresis
            t_c   = deltat*(cal(:,1)-cal(1,1));
            v_c   = cal(:,2)*0.001;
            v_r_c = cal(:,3)*0.001;

            % Muestra cuyo SAR calculamos
            t_m   = deltat*(muestra(:,1)-muestra(1,1));
            v_m   = muestra(:,2)*0.001;
            v_r_m = muestra(:,3)*0.001;
            

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Aumenta el número de períodos para compensar por los pocos ciclos

        % Ajusta las tres referencias con funciones seno:
        % V(t)=V0+A*sin(2*pi*f*t-phi)        
            % Estima valores iniciales para los ajustes
            % Valor medio de la señal, es decir constante aditiva
                par_m(1)=mean(v_r_m);
                par_f(1)=mean(v_r_f);
                par_c(1)=mean(v_r_c);
            % Amplitud
                par_m(2)=(max(v_r_m)-min(v_r_m))/2;
                par_f(2)=(max(v_r_f)-min(v_r_f))/2;
                par_c(2)=(max(v_r_c)-min(v_r_c))/2;
            % Para evitar problemas por errores de nomenclatura en los archivos mide tiempo entre picos:
                suave_m=fftsmooth(v_r_m,round(length(v_r_m)*6/1000));
                suave_f=fftsmooth(v_r_f,round(length(v_r_f)*6/1000));
                suave_c=fftsmooth(v_r_c,round(length(v_r_c)*6/1000));
                
                [~, indices_m]=findpeaks(suave_m);
                [~, indices_f]=findpeaks(suave_f);
                [~, indices_c]=findpeaks(suave_c);
                tiempo_entre_picos_m = mean(diff(t_m(indices_m)));
                tiempo_entre_picos_f = mean(diff(t_f(indices_f)));
                tiempo_entre_picos_c = mean(diff(t_c(indices_c)));
                par_m(3)=1/tiempo_entre_picos_m;
                par_f(3)=1/tiempo_entre_picos_f;
                par_c(3)=1/tiempo_entre_picos_c;
            % Fase inicial, a partir del tiempo del primer máximo
                par_m(4)=2*pi*par_m(3)*t_m(indices_m(1))-pi/2;
                par_f(4)=2*pi*par_f(3)*t_f(indices_f(1))-pi/2;
                par_c(4)=2*pi*par_c(3)*t_c(indices_c(1))-pi/2;

            % Ajusta y obtiene los coeficientes. Separa la frecuencia para nombrar
            % el archivo de salida
                ajuste_m_todo=fitnlm(t_m,v_r_m, @(b,x)(b(1) + b(2)*sin(2*pi*b(3)*x - b(4))), par_m);
                ajuste_f_todo=fitnlm(t_f,v_r_f, @(b,x)(b(1) + b(2)*sin(2*pi*b(3)*x - b(4))), par_f);
                ajuste_c_todo=fitnlm(t_c,v_r_c, @(b,x)(b(1) + b(2)*sin(2*pi*b(3)*x - b(4))), par_c);
                coefCI_m=ajuste_m_todo.Coefficients.Estimate;
                coeficientes_m=coefCI_m(:,1);
                coefCI_f=ajuste_f_todo.Coefficients.Estimate;
                coeficientes_f=coefCI_f(:,1);
                coefCI_c=ajuste_c_todo.Coefficients.Estimate;
                coeficientes_c=coefCI_c(:,1);
                
                frecuencia_m=coeficientes_m(3);
                fase_m=coeficientes_m(4);
                ajuste_m=ajuste_m_todo.Fitted;
                frecuencia_f=coeficientes_f(3);
                fase_f=coeficientes_f(4);
                ajuste_f=ajuste_f_todo.Fitted;
                frecuencia_c=coeficientes_c(3);
                fase_c=coeficientes_c(4);
                ajuste_c=ajuste_c_todo.Fitted;


%         % Si hay menos de cuatro períodos, repite parte de la señal
%             if length(indices_f)<=4
%                 % Ubica el primer punto que tiene la misma fase que el
%                 % final de la señal
%                     fase_final=2*pi*frecuencia_f*t_r_f(end)-fase_f;
%                     indice_aux=0;
%                     fase_aux=-fase_f;
%                     while abs(mod(fase_final-fase_aux,2*pi))>0.001 && 2*pi-abs(mod(fase_final-fase_aux,2*pi))>0.001
%                         indice_aux=indice_aux+1;
%                         fase_aux=2*pi*frecuencia_f*t_r_f(indice_aux)-fase_f;
%                     end
%                 % Agrega desde ese punto hasta el final, después del final,
%                 % para señal y referencia.
%                     t_r_f0=[t_r_f(1:end-1)',t_r_f(indice_aux:end-1)'+t_r_f(end)-t_r_f(indice_aux),t_r_f(indice_aux:end-1)'+2*(t_r_f(end)-t_r_f(indice_aux)),t_r_f(indice_aux:end-1)'+3*(t_r_f(end)-t_r_f(indice_aux))];
%                     v_r_f0=[v_r_f(1:end-1)',v_r_f(indice_aux:end-1)',v_r_f(indice_aux:end-1)',v_r_f(indice_aux:end-1)'];
%                     t_f0=[t_f(1:end-1)',t_f(indice_aux:end-1)'+t_f(end)-t_f(indice_aux),t_f(indice_aux:end-1)'+2*(t_f(end)-t_f(indice_aux)),t_f(indice_aux:end-1)'+3*(t_f(end)-t_f(indice_aux))];
%                     v_f0=[v_f(1:end-1)',v_f(indice_aux:end-1)',v_f(indice_aux:end-1)',v_f(indice_aux:end-1)'];
%             else
%                 t_r_f0=t_r_f';
%                 v_r_f0=v_r_f';
%                 t_f0=t_f';
%                 v_f0=v_f';
%             end



            
        % Comparacion de las referencias y sus ajustes
            figure(11)
            plot(t_m,v_r_m,'o',t_f,v_r_f,'o',t_c,v_r_c,'o',t_m,ajuste_m,t_f,ajuste_f,t_c,ajuste_c)
            title('Comparacion de referencias y ajustes')
            legend('Referencia de fondo','Referencia de muestra','Referencia de la calibración','Ajuste de referencia de fondo','Ajuste de referencia de muestra','Ajuste de referencia de la calibración')
            xlabel('Tiempo (s)')
            ylabel('Referencia (V)')            
            
            
            
            
        % Si la diferencia entre ambas frecuencias es muy grande, error.
            if abs(frecuencia_m-frecuencia_f)/frecuencia_f>0.02||abs(frecuencia_c-frecuencia_f)/frecuencia_f>0.02||abs(frecuencia_c-frecuencia_m)/frecuencia_f>0.02
                'Error: incompatibilidad de frecuencias'
                frecuencia_m
                frecuencia_f
                frecuencia_c
                break
            end

        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Desplazamiento temporal para poner en fase las referencias

        % Saca el modulo 2 pi de las fases y calcula el tiempo de fase 0
            tiempo_fase_m=mod(fase_m,2*pi)/(frecuencia_m*2*pi);
            tiempo_fase_f=mod(fase_f,2*pi)/(frecuencia_f*2*pi);
            tiempo_fase_c=mod(fase_c,2*pi)/(frecuencia_c*2*pi);
            



        % Desplaza en tiempo para que haya coincidencia de fase entre
        % referencias. La amplitud debería ser positiva siempre por el
        % semilleo.
 

            t_m1 = t_m-tiempo_fase_m;
            t_f1 = t_f-tiempo_fase_f;
            t_c1 = t_c-tiempo_fase_c;

            
        % Corrije por posible diferencia de frecuencias, dilatando el tiempo
        % del fondo
            t_fm=t_f1*frecuencia_f/frecuencia_m;
            t_cm=t_f1*frecuencia_c/frecuencia_m;

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Resta el offset de las referencias

        v_r_f1 = v_r_f-coeficientes_f(1);
        v_r_m1 = v_r_m-coeficientes_m(1);
        v_r_c1 = v_r_c-coeficientes_c(1);
        v_f1=v_f;
        v_m1=v_m;
        v_c1=v_c;

        ajuste_f1 = ajuste_f-coeficientes_f(1);
        ajuste_m1 = ajuste_m-coeficientes_m(1);
        ajuste_c1 = ajuste_c-coeficientes_c(1);

        % Comparacion de las referencias de muestra y fondo desplazadas en
        % tiempo y offset
            figure(12)
            plot(t_fm,v_r_f1,'o',t_m1,v_r_m1','o',...
                 t_fm,ajuste_f1,t_m1,ajuste_m1)
            title('Comparacion de referencias desplazadas y restados sus offsets')
            legend('Referencia de fondo','Referencia de muestra')
            grid on
            xlabel('Tiempo (s)')
            ylabel('Referencia (V)')
            
        % Comparacion de las referencias de calibración y fondo desplazadas en
        % tiempo y offset
            figure(13)
            plot(t_cm,v_r_f1,'o',t_c1,v_r_c1','o',...
                 t_cm,ajuste_f1,t_c1,ajuste_c1)
            title('Comparacion de referencias desplazadas y restados sus offsets')
            legend('Referencia de fondo','Referencia de calibración')
            grid on
            xlabel('Tiempo (s)')
            ylabel('Referencia (V)')

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Recorte para que ambas medidas tengan igual tiempo inicial y número
    % entero de períodos

        % Resta medida y fondo interpolando para que corresponda mejor el
        % tiempo. No trabaja más con la medidas individuales, sólo con la
        % resta. Se toma la precaución para los casos de trigger distinto
        % entre fondo y medida.
            tmin_m=t_fm(1);
            tmax_m=t_fm(end);
            tiempo_coincidente_m=find(t_m1>=tmin_m & t_m1<=tmax_m);
            t_aux_m=t_m1(tiempo_coincidente_m);
            interpolacion_aux_m=interp1(t_fm,v_f1,t_aux_m,'linear');
            interpolacion_m=zeros(size(v_m1));
            for i=1:length(t_m1)
                [~, j]=min(abs(t_aux_m-t_m1(i)));
                interpolacion_m(i)=interpolacion_aux_m(j);
            end

            Resta_m = v_m1-interpolacion_m;
            
            tmin_c=t_cm(1);
            tmax_c=t_cm(end);
            tiempo_coincidente_c=find(t_c1>=tmin_c & t_c1<=tmax_c);
            t_aux_c=t_c1(tiempo_coincidente_c);
            interpolacion_aux_c=interp1(t_cm,v_f1,t_aux_c,'linear');
            interpolacion_c=zeros(size(v_c1));
            for i=1:length(t_c1)
                [~, j]=min(abs(t_aux_c-t_c1(i)));
                interpolacion_c(i)=interpolacion_aux_c(j);
            end
            

            Resta_c = v_c1-interpolacion_c;
            
            
        % Comparacion de las medidas con el tiempo que pone las referencias en
        % fase
            figure(14)
            plot(t_fm,v_f1,'-o',t_m1,v_m1,'-o',t_m1,interpolacion_m,'-o')
            title('Comparación de medidas')
            legend('Medida de fondo','Medida de muestra','Interpolación del fondo al tiempo de la muestra')
            xlabel('Tiempo (s)')
            ylabel('Señal (V)')
            
        % Comparacion de las medidas con el tiempo que pone las referencias en
        % fase
            figure(15)
            plot(t_cm,v_f1,'-o',t_c1,v_c1,'-o',t_c1,interpolacion_c,'-o')
            title('Comparación de medidas')
            legend('Medida de fondo','Medida de calibración','Interpolación del fondo al tiempo de la calibración')
            xlabel('Tiempo (s)')
            ylabel('Señal (V)')            

        % Resta de señales de muestra o calibración menos fondo
            figure(16)
            plot(t_m1,Resta_m,'o-',t_c1,Resta_c,'o-')
            title('Resta de señales')
            legend('Resta de medida de muestra y medida de fondo','Resta de medida de calibración y medida de fondo')
            xlabel('Tiempo (s)')
            ylabel('Resta (V)')

        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Filtrado
    
        % Filtro por Fourier la calibración
        freq=round(length(v_r_c)/5);
        Resta_c3=fftsmooth(Resta_c,freq)';
        v_r_c3=fftsmooth(v_r_c,freq)';
        t_c3=t_c1';
        ajuste_c3=ajuste_c1';
        
        
    if FILTRARMUESTRA~=2 
        % Identifica el ruido 

            % Factor del ruido natural a partir del cual se considera que
            % hay que filtrar
                ancho=2.5;
            % Puntos a ambos lados del sector de señal ruidosa que serán
            % incluidos en el ruido.
                entorno=6;

            % Marcadores de regiones a filtrar
                [t_m2, marcador_resta_m]=encuentra_ruido(t_m1',Resta_m',ancho,entorno);
                [t_c2, marcador_resta_c]=encuentra_ruido(t_c1',Resta_c',ancho,entorno);

        % Ajuste ruido
            % Parámetros del ajuste: puntos a cada lado de la region a
            % filtrar que serán considerados para el ajuste, y grado del
            % polinomio a ajustar.
                puntos_ajuste=80;
                grado_pol=3;

            % Tiempos y señales en las mismas dimensiones que los marcadores.
                v_r_m2=interp1(t_m1,v_r_m1,t_m2,'spline');
                v_r_c2=interp1(t_c1,v_r_c1,t_c2,'spline');
                Resta_m2=interp1(t_m1,Resta_m,t_m2,'spline');
                Resta_c2=interp1(t_c1,Resta_c,t_c2,'spline');

                
            % Filtrado de señales de resta
            % Aquí podría meter mano para filtrar la calibración si no
            % quiero fourier, permitiendo que m vaya de 1 a 2
                    for m=1:1
                        if m==1
                            voltaje=Resta_m2;
                            referencia=v_r_m2;
                            tiempo=t_m2;
                            marcador=marcador_resta_m;
                        else
                            voltaje=Resta_c2;
                            referencia=v_r_c2;
                            tiempo=t_c2;
                            marcador=marcador_resta_c;
                        end
                    % Comienza a filtrar a partir de que tiene suficientes
                    % puntos detrás
                        w=puntos_ajuste+1;
                        voltaje2=voltaje;
                        referencia2=referencia;

                    % Barre la señal. No filtra ni el principio ni el
                    % final, por eso más adelante se eliminan el primer y
                    % el último semiperíodos.
                        while w<length(voltaje)
                            % Si no hay ruido deja la señal como estaba
                            if marcador(w-1)==0
                                w=w+1;
                            % Si hay ruido
                            elseif marcador(w-1)==1
                                t=w;
                                % Busca hasta donde llega el ruido
                                while marcador(t-1)==1 && t<length(voltaje)
                                    t=t+1;
                                end
                                % Si tiene suficientes puntos del otro lado
                                % realiza el ajuste
                                if t<length(voltaje)-puntos_ajuste
                                    y=[voltaje(w-puntos_ajuste:w);voltaje(1+t:1+t+puntos_ajuste)];
                                    x=[tiempo(w-puntos_ajuste:w);tiempo(1+t:1+t+puntos_ajuste)];
                                    p=polyfit(x,y,grado_pol);
                                    voltaje2(w:t-1)=polyval(p,tiempo(w:t-1));
                                    y=[referencia(w-puntos_ajuste:w);referencia(1+t:1+t+puntos_ajuste)];
                                    x=[tiempo(w-puntos_ajuste:w);tiempo(1+t:1+t+puntos_ajuste)];
                                    p=polyfit(x,y,grado_pol);
                                    referencia2(w:t-1)=polyval(p,tiempo(w:t-1));                                    
                                end
                                w=t;
                            end
                        end
                    
                    % Renombra
                        if m==1
                            if FILTRARMUESTRA==1
                                Resta_m3=voltaje2;
                                v_r_m3=referencia2;
                                t_m3=t_m2;
                            else
                                Resta_m3=Resta_m';
                                v_r_m3=v_r_m';
                                t_m3=t_m1';
                            end
                        elseif m==2
                            Resta_c3=voltaje2;
                            v_r_c3=referencia2;
                        end
                    end
            
        % Nada mas para que tiempo y señal tengan la misma nomenclatura
%             t_c3=t_c2;
            ajuste_m3=interp1(t_m1,ajuste_m1,t_m3,'spline');
%             ajuste_c3=interp1(t_c1,ajuste_c1,t_c3,'spline');
            
        % Control de que el suavizado final sea satisfactorio
            figure(17)
            plot(t_m1,v_r_m1/max(v_r_m3),'o',t_m3,v_r_m3/max(v_r_m3),t_m2,marcador_resta_m,...
                t_m1,Resta_m/max(Resta_m3)+3,'o',t_m3,Resta_m3/max(Resta_m3)+3,t_m2,marcador_resta_m+3)
            title('Quitando el ruido de la muestra')
            legend('Referencia de muestra','Sin ruido','Zona de ruido','Resta de señales','Sin ruido','Zona de ruido')
            xlabel('Tiempo (s)')
            ylabel('u.a.')
            axis([t_m1(1) t_m1(end) -1.1 4.1])
            
%             figure(18)
%             plot(t_c1,v_r_c1/max(v_r_c3),'o',t_c3,v_r_c3/max(v_r_c3),t_c2,marcador_resta_c,...
%                 t_c1,Resta_c/max(Resta_c3)+3,'o',t_c3,Resta_c3/max(Resta_c3)+3,t_c2,marcador_resta_c+3)
%             title('Quitando el ruido de la calibración')
%             legend('Referencia de calibración','Sin ruido','Zona de ruido','Resta de señales','Sin ruido','Zona de ruido')
%             xlabel('Tiempo (s)')
%             ylabel('u.a.')
%             axis([t_c1(1) t_c1(end) -1.1 4.1])
            
        % Diferencia entre señal sin ruido y señal. Guarda el peor valor.
            dif_resta_m=Resta_m3-interp1(t_m1,Resta_m,t_m3,'spline');
            dif_resta_c=Resta_c3-Resta_c;
            peor_diferencia=max([mean(abs(dif_resta_m))/max(Resta_m),mean(abs(dif_resta_c))/max(Resta_c)]);
            
    else
        
        % Frecuencia de corte en Hz;
%         freq_corte=2*10^6;
%         freq_natural=1/(t_f(2)-t_f(1));
        freq=round(length(v_r_m)/5);
        Resta_m3=fftsmooth(Resta_m,freq)';
        v_r_m3=fftsmooth(v_r_m,freq)';
        t_m3=t_m1';
        ajuste_m3=ajuste_m1';
        
        % Control de que el suavizado final sea satisfactorio
            figure(17)
            plot(t_m1,v_r_m1/max(v_r_m3),'o',t_m3,v_r_m3/max(v_r_m3),...
                t_m1,Resta_m/max(Resta_m3)+3,'o',t_m3,Resta_m3/max(Resta_m3)+3)
            title('Quitando el ruido de la muestra')
            legend('Referencia de muestra','Sin ruido','Resta de señales','Sin ruido')
            xlabel('Tiempo (s)')
            ylabel('u.a.')
            axis([t_m1(1) t_m1(end) -1.1 4.1])
            
            figure(18)
            plot(t_c1,v_r_c1/max(v_r_c3),'o',t_c3,v_r_c3/max(v_r_c3),...
                t_c1,Resta_c/max(Resta_c3)+3,'o',t_c3,Resta_c3/max(Resta_c3)+3)
            title('Quitando el ruido de la calibración')
            legend('Referencia de calibración','Sin ruido','Resta de señales','Sin ruido')
            xlabel('Tiempo (s)')
            ylabel('u.a.')
            axis([t_c1(1) t_c1(end) -1.1 4.1])
        
        % Diferencia entre señal sin ruido y señal. Guarda el peor valor.
            dif_resta_m=Resta_m3-Resta_m;
            dif_resta_c=Resta_c3-Resta_c;
            peor_diferencia=max([mean(abs(dif_resta_m))/max(Resta_m),mean(abs(dif_resta_c))/max(Resta_c)]);
        
    end
        
        
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Recorta un número entero de períodos o ciclos
    
        % Calculo del numero de ciclos enteros
            N_ciclos_m = floor((t_m3(end)-t_m3(1))*(frecuencia_m));
            N_ciclos_c = floor((t_c3(end)-t_c3(1))*(frecuencia_c));
        
        % Indices ciclo
            indices_ciclos_m=find(t_m3<t_m3(1)+N_ciclos_m*1/frecuencia_m);
            indices_ciclos_c=find(t_c3<t_c3(1)+N_ciclos_c*1/frecuencia_c);
            largo_m=indices_ciclos_m(end);
            largo_c=indices_ciclos_c(end);
            if mod(largo_m,N_ciclos_m)==0
            elseif mod(largo_m,N_ciclos_m)<=0.5
                largo_m=largo_m-mod(largo_m,N_ciclos_m);
            else
                largo_m=largo_m+N_ciclos_m-mod(largo_m,N_ciclos_m);
            end
            if mod(largo_c,N_ciclos_c)==0
            elseif mod(largo_c,N_ciclos_c)<=0.5
                largo_c=largo_c-mod(largo_c,N_ciclos_c);
            else
                largo_c=largo_c+N_ciclos_c-mod(largo_c,N_ciclos_c);
            end
        
        % Se recorta medio período a cada lado
            indices_recortado_m=ceil(largo_m/N_ciclos_m/2):1:ceil(largo_m-largo_m/N_ciclos_m/2);
            indices_recortado_c=ceil(largo_c/N_ciclos_c/2):1:ceil(largo_c-largo_c/N_ciclos_c/2);
            N_ciclos_c=N_ciclos_c-1;
            N_ciclos_m=N_ciclos_m-1;

        % Se recortan los vectores
            t_m4   = t_m3(indices_recortado_m);
            t_c4   = t_c3(indices_recortado_c);            
            v_r_m4 = v_r_m3(indices_recortado_m);
            v_r_c4 = v_r_c3(indices_recortado_c);
            Resta_m4   = Resta_m3(indices_recortado_m);
            Resta_c4   = Resta_c3(indices_recortado_c);
            ajuste_m4 = ajuste_m3(indices_recortado_m);
            ajuste_c4 = ajuste_c3(indices_recortado_c);


        % Referencia de muestra sin valor medio y su ajuste,y resta de
        % señal de muestra y fondo, recortadas para que tengan un numero
        % entero de ciclos

            figure(19)
            plot(t_m4,v_r_m4,'o',t_m4,400*Resta_m4,'o',t_c4,v_r_c4,'o',t_c4,400*Resta_c4,'o')
            hold on
            plot(t_m4,ajuste_m4,'r',t_c4,ajuste_c4,'r')
            hold off
            title('Referencia de muestra desplazada y sin valor medio, y resta de señales. N entero de periodos')
            legend('Referencia de muestra','Resta de muestra reescaleada x400','Referencia de calibración','Resta de calibración reescaleada x400')
            xlabel('Tiempo (s)')
            ylabel('Señal (V)')            
            
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Campo y magnetización: se integran las funciones de referencia y
    % resta.
    
        % Últimos ajustes de las referencias
        % Valor medio de la señal, es decir constante aditiva
            par_mf(1)=mean(v_r_m4);
            par_cf(1)=mean(v_r_c4);
        % Amplitud
            par_mf(2)=(max(v_r_m4)-min(v_r_m4))/2;
            par_cf(2)=(max(v_r_c4)-min(v_r_c4))/2;
        % Frecuencia
            par_mf(3)=frecuencia_m;
            par_cf(3)=frecuencia_c;
        % Fase inicial
            par_mf(4)=2*pi*par_m(3)*t_m4(indices_m(1))-pi/2;
            par_cf(4)=2*pi*par_c(3)*t_c4(indices_c(1))-pi/2;

        % Ajusta y obtiene los coeficientes.
            ajuste_final_m=fitnlm(t_m4,v_r_m4, @(b,x)(b(1) + b(2)*sin(2*pi*b(3)*x - b(4))), par_mf);
            coefCI=ajuste_final_m.Coefficients.Estimate;
            coeficientes_final_m=coefCI(:,1);
            frecuencia_final_m=coeficientes_final_m(3);
            
            ajuste_final_c=fitnlm(t_c4,v_r_c4, @(b,x)(b(1) + b(2)*sin(2*pi*b(3)*x - b(4))), par_cf);
            coefCI=ajuste_final_c.Coefficients.Estimate;
            coeficientes_final_c=coefCI(:,1);
            frecuencia_final_c=coeficientes_final_c(3);
            

        % Promedia los ciclos antes de integrar.
            tiempo_final_m=t_m4(find(t_m4<=t_m4(1)+1/frecuencia_final_m));
            tiempo_final_c=t_c4(find(t_c4<=t_c4(1)+1/frecuencia_final_c));
            fondo1_m=zeros(size(tiempo_final_m));
            fondo1_c=zeros(size(tiempo_final_c));
            muestra_0=zeros(size(tiempo_final_m));
            cal_0=zeros(size(tiempo_final_c));
            ajuste_m=zeros(size(tiempo_final_m));
            ajuste_c=zeros(size(tiempo_final_c));
            for m=1:N_ciclos_m
                if tiempo_final_m(end)+(m-1)/frecuencia_final_m<t_m4(end)
                    fondo1_m=fondo1_m+interp1(t_m4,v_r_m4,tiempo_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m;
                    muestra_0=muestra_0+interp1(t_m4,Resta_m4,tiempo_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m;
                    ajuste_m=ajuste_m+interp1(t_m4,ajuste_m4,tiempo_final_m+(m-1)/frecuencia_final_m)/N_ciclos_m;
                else
                    fondo1_m=fondo1_m+interp1(t_m4,v_r_m4,tiempo_final_m+(m-1)/frecuencia_final_m,'spline')/N_ciclos_m;
                    muestra_0=muestra_0+interp1(t_m4,Resta_m4,tiempo_final_m+(m-1)/frecuencia_final_m,'spline')/N_ciclos_m;
                    ajuste_m=ajuste_m+interp1(t_m4,ajuste_m4,tiempo_final_m+(m-1)/frecuencia_final_m,'spline')/N_ciclos_m;
                end
            end
            for m=1:N_ciclos_c
                if tiempo_final_c(end)+(m-1)/frecuencia_final_c<t_m4(end)
                    fondo1_c=fondo1_c+interp1(t_c4,v_r_c4,tiempo_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c;
                    cal_0=cal_0+interp1(t_c4,Resta_c4,tiempo_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c;
                    ajuste_c=ajuste_c+interp1(t_c4,ajuste_c4,tiempo_final_c+(m-1)/frecuencia_final_c)/N_ciclos_c;
                else
                    fondo1_c=fondo1_c+interp1(t_c4,v_r_c4,tiempo_final_c+(m-1)/frecuencia_final_c,'spline')/N_ciclos_c;
                    cal_0=cal_0+interp1(t_c4,Resta_c4,tiempo_final_c+(m-1)/frecuencia_final_c,'spline')/N_ciclos_c;
                    ajuste_c=ajuste_c+interp1(t_c4,ajuste_c4,tiempo_final_c+(m-1)/frecuencia_final_c,'spline')/N_ciclos_c;
                end
            end

            
        % Le quita el valor medio a la resta y al fondo
            Rm = muestra_0-mean(muestra_0);
            fondo1_m=fondo1_m-mean(fondo1_m);
            Rc = cal_0-mean(cal_0);
            fondo1_c=fondo1_c-mean(fondo1_c);            
        % Paso temporal
            delta_t_m=(tiempo_final_m(end)-tiempo_final_m(1))/length(tiempo_final_m);
            delta_t_c=(tiempo_final_c(end)-tiempo_final_c(1))/length(tiempo_final_c);
            
            figure(99)
            plot(tiempo_final_m,fondo1_m)
            
            
            
            
%% CALIBRACIÓN

        % Calcula las sumas acumuladas y convierte a campo y magnetizacion
            % Constante que dimensionaliza al campo en A/m a partir de la
            % calibración realizada sobre la bobina del RF
                C_norm_campo=Idc*pendiente_HvsI+ordenada_HvsI;
        
            % Campo en volt*segundo, falta llevar a la amplitud conocida.
                campo_ua0_c = delta_t_c*cumtrapz(fondo1_c); 
            % Centrado en cero
                campo_ua_c = campo_ua0_c-(max(campo_ua0_c)+min(campo_ua0_c))/2;
            % Lo mismo con el campo obtenido del ajuste de la referencia
                campo_fit_ua_c = delta_t_c*cumtrapz(ajuste_c-mean(ajuste_c));
                campo_fit_ua_c = campo_fit_ua_c-mean(campo_fit_ua_c); 

            % Magnetización en volt*segundo, falta el factor geométrico
                magnetizacion_ua0_c=delta_t_c*cumtrapz(Rc);      
            % Centrado en cero
                magnetizacion_ua_c=magnetizacion_ua0_c-mean(magnetizacion_ua0_c);

        % Doy unidades al campo
            campo_c=campo_ua_c*C_norm_campo/max(campo_ua_c);
            campo_fit_ua_c=campo_fit_ua_c*C_norm_campo/max(campo_fit_ua_c);

            figure(20)
            plot(tiempo_final_c,campo_ua_c/max(campo_ua_c),tiempo_final_c,magnetizacion_ua_c/max(magnetizacion_ua_c),'LineWidth',1.5)
            title('Campo y magnetización normalizados calibración')
            legend('Campo','Magnetizacion')
            xlabel('Tiempo (s)')
            ylabel('u.a')
            
            figure(21)
            plot(campo_c,magnetizacion_ua_c,'LineWidth',1.5)
            hold all
            %plot(campo_fit_ua,magnetizacion_ua,'LineWidth',1.5)
            title('Ciclo de calibración en función del campo')
            xlabel('H (Oe)')
            ylabel('u.a')
        
        % Ajusta una recta
            recta=polyfit(campo_c,magnetizacion_ua_c,1);
            polaridad=sign(recta(1));
            pendiente=recta(1);
            
        % Calibración ua a A/m
            calibracion=xi_patron_vol/pendiente;
            
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% MUESTRA

        % Calcula las sumas acumuladas y convierte a campo y magnetizacion
        
            % Campo en volt*segundo, falta llevar a la amplitud conocida.
                campo_ua0_m = delta_t_m*cumtrapz(fondo1_m); 
            % Centrado en cero
                campo_ua_m = campo_ua0_m-mean(campo_ua0_m);
            % Lo mismo con el campo obtenido del ajuste de la referencia
                campo_fit_ua_m = delta_t_m*cumtrapz(ajuste_m-mean(ajuste_m));
                campo_fit_ua_m = campo_fit_ua_m-mean(campo_fit_ua_m); 

            % Magnetización en volt*segundo, falta el factor geométrico
                magnetizacion_ua0_m=delta_t_m*cumtrapz(Rm);      
            % Centrado en cero
                magnetizacion_ua_m=magnetizacion_ua0_m-(max(magnetizacion_ua0_m)+min(magnetizacion_ua0_m))/2;

        % Doy unidades al campo
            campo_fit_ua_m=campo_fit_ua_m*C_norm_campo/max(campo_fit_ua_m);
            magnetizacion_m=calibracion.*magnetizacion_ua_m;

            figure(22)
            plot(tiempo_final_m,campo_ua_m/max(campo_ua_m),tiempo_final_m,magnetizacion_m/max(magnetizacion_m),'LineWidth',1.5)
            title('Campo y magnetización normalizados muestra')
            legend('Campo','Magnetizacion')
            xlabel('Tiempo (s)')
            ylabel('u.a')
            
             
            % Campo y magnetización finales. Se mantiene un campo adicional
            % generado a partir del ajuste senoidal de la referencia para
            % hacer comparaciones (campo_fit).
                
                campo_m=campo_ua_m/max(campo_ua_m)*C_norm_campo;
                campo_fit=campo_fit_ua_m/max(campo_fit_ua_m)*C_norm_campo;

            % Grafica el ciclo momentáneamente para guardar una imagen
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
                
            % Guardo imagen del ciclo
                 archivo_imagen=[filenombre_muestra '.png'];
                 imwrite(im,archivo_imagen)

            % Va graficando todos los ciclos juntos
                figure(23)
                hold on
                plot(campo_m,magnetizacion_m,'-.','LineWidth',1.5,'Color',[k/fa,1-k/fa,1-k/fa])
                grid on
                title('Ciclo de histéresis')
                legend('Magnetizacion vs Campo')
                xlabel('Campo (A/m)')
                ylabel('Magnetizacion (A/m)')
                hold off

        % Exporta ciclo en ascii
            ciclo_out=[campo_m',magnetizacion_m',magnetizacion_m'/max(magnetizacion_m)];
            filenombre_out=strcat(filenombre_muestra,'_ciclop0.dat');
            save(filenombre_out,'ciclo_out','-ASCII')
        


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Cálculo de campo coercitivo y remanencia

        m = magnetizacion_m;
        h = campo_m;

        clear Hc
        p=1;
        q=1;
        for i=1:length(m)-1
            if (  m(i)>0 && m(i+1)<0  || m(i)<0 && m(i+1)>0 )
                 Hc(p)=h(i) - m(i)* ( h(i+1)-h(i) ) / (m(i+1)-m(i));
                 p = p+1;
            end
            if (  h(i)>0 && h(i+1)<0  || h(i)<0 && h(i+1)>0 )
                 Mr(p)=m(i) - h(i)* ( m(i+1)-m(i) ) / (h(i+1)-h(i));
                 q = q+1;
            end
        end

        Hc_mean = mean(abs(Hc));
        Hc_error = std(abs(Hc));

        Mr_mean = mean(abs(Mr));
        Mr_error = std(abs(Mr));



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Determinación de áreas

    % Arma un vector auxilar desplazado a valores positivos
        magnetizacion_des = magnetizacion_m+2*abs(min(magnetizacion_m));        
    % Área del lazo
        area=-trapz([campo_m, campo_m(1)],[magnetizacion_des, magnetizacion_des(1)])
 
 
    % Cálculo de potencia disipada SAR
        sar=mu0*abs(area)*frecuencia_final_m/concentracion;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Archivo de salida
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Salidas de interés
            % Corriente del resonador
                Corriente_A(k)=str2num(filenombre_muestra(8:9));
            % Frecuencia de la referencia en la medida de la muestra
                Frecuencia_muestra_kHz(k)=frecuencia_final_m/1000; 
            % Frecuencia de la referencia en la medida del fondo
                Frecuencia_fondo_kHz(k)=frecuencia_f/1000;  
            % Specific Absorption Rate
                SAR(k)=sar; 
            % Campo maximo en A/m
                Campo_Maximo(k)=C_norm_campo;
            % Campo coercitivo en A/m
                Coercitividad(k)=Hc_mean;
            % Magnetizacion remanente en A/m
                Magnetizacion_Remananete(k)=Mr_mean;
                
        % Se grafica el SAR en función del campo máximo y en función del campo
        % máximo al cuadrado 
    
            figure(24)
            plot(Campo_Maximo/1000,SAR,'.-','MarkerSize',20)
            title('SAR vs Campo Maximo')
            xlabel('Campo Maximo (kA/m)')
            ylabel('SAR (W/g)')

            figure(25)
            plot((Campo_Maximo/1000).^2,SAR,'.','MarkerSize',20)
            title('SAR vs Campo Maximo^2')
            xlabel('Campo Maximo^2 (kA/m)^2')
            ylabel('SAR (W/g)')

            data_salida(k,:)=[frecuencia_final_m/1000, C_norm_campo/1000, sar, Hc_mean/1000, Mr_mean/1000, 10^11/calibracion];

end
    
% Encabezado del archivo de salida
    encabezado={'Nombre_del_archivo', 'Frecuencia_kHz', 'Campo_maximo_kA_m', 'SAR_W_g', 'Coercitividad_kA_m', 'Magnetizacion_Remanente_kA_m', 'Pendiente_de_calibracion'};
    salida=[encabezado; nombres_archivos', num2cell(data_salida)];
    
    %tabla_salida=table(nombres_archivos',data_salida(:,1),data_salida(:,2),data_salida(:,3),data_salida(:,4),data_salida(:,5),data_salida(:,6),'VariableNames',encabezado);
    
    %writetable(tabla_salida,nombre_archivo_salida,'Delimiter','\t')
    
% % Prepara el archivo para escribir la salida
    IDarchivo = fopen(nombre_archivo_salida,'w');
    formato_encabezado = '%-21s %-21s %-21s %-21s %-21s %-33s %-21s %-21s\r\n';
    formato = '%-21s %-21f %-21f %-21f %-21f %-33f %-21f %-21f\r\n';   
    [nfilas,~] = size(salida);
% 
% % Imprime los datos en el archivo
    fprintf(IDarchivo,formato_encabezado,salida{1,:});
    for fila = 2:nfilas
        fprintf(IDarchivo,formato,salida{fila,:});
    end
    fclose(IDarchivo);


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
        % use fft low pass filter to smoothen out a signal
        
        fft_data_v = fft(data_v);
        s_fft_data_v = zeros(1,length(data_v));
        s_fft_data_v(1:freq_n) = fft_data_v(1:freq_n);
        s_fft_data_v(end-freq_n:end) = fft_data_v(end-freq_n:end);
        s_data_v = real(ifft(s_fft_data_v));
    end
end
