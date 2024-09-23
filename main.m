L = input("\nLongitud de la superficie de difusión: ");
t_max = input("\nTiempo total de la simulación: ");
n_x = input("\nNúmero de nodos en el espacio: ");
n_t = input("\nNúmero de nodos en el tiempo: ");
alpha = input("\nConstante del material: ");
T_1j = input("\nTemperatura del contorno: ");
T_i1 = input("\nTemperatura inicial: ");
T_nj = input("\nTemperatura fija del último nodo: ");
caso = input("\nCaso 1, 2 o 3: ");

d_x = L/(n_x - 1); % Se calcula delta de x
d_t = t_max/(n_t - 1); % Se calcula delta de t

a = d_t;
b = -(2*d_t + d_x^2*d_t*alpha^2 + d_x^2);
d = -d_x^2;

if caso == 1
    A_1 = [zeros(1, n_x - 2); a*eye(n_x - 3, n_x - 2)];
    A_2 = b*eye(n_x - 2);
    A_3 = [zeros(n_x - 2, 1) a*eye(n_x - 2, n_x - 3)];
    A = A_1 + A_2 + A_3;

    temp = T_i1*ones(n_x - 2, 1); % Temperatura inicial
    T_A = [T_1j; temp; T_nj]; % Temperatura de contorno

    sol_A = d*T_i1*ones(n_x - 2, 1);
    sol_A(1) = sol_A(1) - a*T_1j;
    sol_A(n_x - 2) = sol_A(n_x - 2) - a*T_nj;

    % disp(T_A)

elseif caso == 2
    A_1 = [zeros(1, n_x - 1); a*eye(n_x - 2, n_x - 1)];
    A_2 = b*eye(n_x - 1);
    A_3 = [zeros(n_x - 1, 1) a*eye(n_x - 1, n_x - 2)];
    A = A_1 + A_2 + A_3;
    A(n_x - 1, n_x - 2) = A(n_x - 1, n_x - 2)*2;

    temp = T_i1*ones(n_x - 1, 1);
    T_A = [T_1j; temp];

    sol_A = d*T_i1*ones(n_x - 1, 1);
    sol_A(1) = sol_A(1) - a*T_1j;

elseif caso == 3
    A_1 = [zeros(1, n_x - 1); a*eye(n_x - 2, n_x - 1)];
    A_2 = b*eye(n_x - 1);
    A_3 = [zeros(n_x - 1, 1) a*eye(n_x - 1, n_x - 2)];
    A = A_1 + A_2 + A_3;
    A(1, 2) = A(1, 2)*2;

    temp = T_i1*ones(n_x - 1, 1); % Temperatura inicial
    T_A = [temp; T_nj]; % Temperatura de contorno

    sol_A = d*T_i1*ones(n_x - 1, 1);
    sol_A(n_x - 1) = sol_A(n_x - 1) - a*T_nj;
    
else
    disp("Número de caso no valido. Inicie otra vez.")
end

% Soluciones

tic

figure(1)
p = plot(0:d_x:L, T_A);
p.LineWidth = 2;
hold on

for i = 2:n_t
    % Se halla la nueva solución
    temp = pivoteoParcial(A, sol_A);
    if caso == 1
        T_A = [T_1j; temp; T_nj];
        % disp(T_A)
    elseif caso == 2
        T_A = [T_1j; temp];
    elseif caso == 3
         T_A = [temp; T_nj];
    else
        disp("Error");
    end
    
    % disp(T_A)
    p = plot(0:d_x:L, T_A);
    p.LineWidth = 2;

    % Se actualiza el vector b para usar en la siguiente iteración
    sol_A = d*temp;
    if caso == 1
        sol_A(1) = sol_A(1) - a*T_1j;
        sol_A(n_x - 2) = sol_A(n_x - 2) - a*T_nj;
    elseif caso == 2
         sol_A(1) = sol_A(1) - a*T_1j;
    elseif caso == 3
         sol_A(n_x - 1) = sol_A(n_x -1) - a*T_nj;
    else
        disp("Error");
    end

end

hold off
toc

%% Estacionario

L = 1;
n_x = 11;
alpha = 0.2;
T_1j = 70;
T_i1 = 20;
T_nj = 70;
caso = input("\nCaso 1, 2 o 3: ");

d_x = L/(n_x - 1); % Se calcula delta de x
a = -(2 + d_x^2*alpha^2);

if caso == 1
    B_1 = [zeros(1, n_x - 2); eye(n_x - 3, n_x - 2)];
    B_2 = a*eye(n_x - 2);
    B_3 = [zeros(n_x - 2, 1) eye(n_x - 2, n_x - 3)];
    B = B_1 + B_2 + B_3;
  
    sol_B = zeros(n_x - 2, 1);
    sol_B(1) = - T_1j;
    sol_B(n_x - 2) = -1*T_nj;

    T_B = pivoteoParcial(B, sol_B);
    T_B = [T_1j; T_B; T_nj];
    
    figure(2)
    p = plot(0:d_x:L, T_B);
    p.Color = "#4DBEEE";
    p.LineWidth = 2;

elseif caso == 2
    B_1 = [zeros(1, n_x - 1); eye(n_x - 2, n_x - 1)];
    B_2 = a*eye(n_x - 1);
    B_3 = [zeros(n_x - 1, 1) eye(n_x - 1, n_x - 2)];
    B = B_1 + B_2 + B_3;
    B(n_x - 1, n_x - 2) = 2;
    
    sol_B = zeros(n_x - 1, 1);
    sol_B(1) = -T_1j;
  
    T_B = pivoteoParcial(B, sol_B);
    T_B = [T_1j; T_B];
    
    figure(2)
    p = plot(0:d_x:L, T_B);
    p.Color = "#4DBEEE";
    p.LineWidth = 2;

elseif caso == 3
    B_1 = [zeros(1, n_x - 1); eye(n_x - 2, n_x - 1)];
    B_2 = a*eye(n_x - 1);
    B_3 = [zeros(n_x - 1, 1) eye(n_x - 1, n_x - 2)];
    B = B_1 + B_2 + B_3;
    B(1, 2) = 2;
    
    sol_B = zeros(n_x - 1, 1);
    sol_B(n_x - 1) = -T_nj;
  
    T_B = pivoteoParcial(B, sol_B);
    T_B = [T_B; T_nj];
    
    figure(2)
    p = plot(0:d_x:L, T_B);
    p.Color = "#4DBEEE";
    p.LineWidth = 2;
end


%% Funciones 

function x = pivoteoParcial(A, b)
    [filas, columnas] = size(A);
    if filas ~= columnas
        error('La matriz A no es cuadrada');
    end

    n = filas;
    Ab = [A, b];
    for k = 1:n-1
        [~, i_max] = max(abs(Ab(k:n, k)));
        i_max = i_max + k-1;
        if i_max ~= k
            Ab([k, i_max], :) = Ab([i_max, k], :);
        end
    
        for i = k+1:n
            factor = Ab(i, k)/Ab(k, k);
            Ab(i, k:n+1) = Ab(i, k:n+1)- factor*Ab(k, k:n+1);
        end
    end

    %Sustitución hacia atrás
    x = zeros(n, 1);
    x(n) = Ab(n, n+1) / Ab(n, n);
    for i = n-1:-1:1
        x(i) = (Ab(i, n+1) - Ab(i, i+1:n) * x(i+1:n)) / Ab(i, i);
    end
end