%% Solution of the Heat Equation with non homogeneous, time varying Dirichlet boundary conditions using FEM

format long g
%solution on a very fine grid to compute error
sol_ref = load('data/solRef.dat.txt'); 
xRef=sol_ref(:,1);
yRef=sol_ref(:,2);
uRef=sol_ref(:,3);
interp = scatteredInterpolant(xRef, yRef, uRef);

mesh_num = 0:4; 
for mesh = 1 : length(mesh_num) 
ind=mesh_num(mesh);
topology = load(sprintf('data/mesh%i.topol', ind)); 
coordinates = load(sprintf('data/mesh%i.coord', ind));
boundary_conditions = load(sprintf('data/mesh%i.bound', ind));
track= load(sprintf('data/mesh%i.track',  ind));
trace= load(sprintf('data/mesh%i.trace',  ind));

x= coordinates(:,1);
y=coordinates(:,2);
sol = zeros(size(coordinates, 1), 1);
for i=1:size(coordinates, 1)
   sol(i) = interp(x(i), y(i));
end

    num_track = length(track);
    node_indices = trace(:, 1);
    arc_lengths = trace(:, 2);
    solutions_at_times = zeros(length(node_indices), 4);
    
    [H, M] = compute_local_matrices(topology, coordinates);
    H=sparse(H);
    M=sparse(M);
    % Initialize  solution vector (u) to zeros
    num_nodes = size(coordinates, 1);
    
    %u_chol= zeros(num_nodes, 1);
    delta_t = 0.02;
    t_max = 10;
    theta = 0.5;
    num_steps = t_max / delta_t;
    b_chol = zeros(num_nodes,num_steps+1 );
    b_j = zeros(num_nodes,num_steps+1 );
    u_chol= zeros(num_nodes, num_steps+1);
    u_j= zeros(num_nodes, num_steps+1);
    track_solutions = zeros(num_track, num_steps + 1);  % +1 perché includiamo lo stato iniziale

    % Precompute  matrices for the Crank-Nicolson method
    K1 = (M / delta_t) + (theta * H);
    K2 = (M / delta_t) - ((1 - theta) * H);
   % J= sparse(zeros(size(K1)));
    
 for step = 1:num_steps     % Update the time-dependent boundary conditions
        t_current = step * delta_t;
        for c=1:size(boundary_conditions,1)
            if (step*delta_t<=5 && boundary_conditions(c,2)~=0)
                boundary_conditions(c,2)=((2*step*delta_t)/10)*1;
            elseif (step*delta_t>5 && boundary_conditions(c,2)~=0)
                boundary_conditions(c,2)=1;
            end
        end
         b_chol(:,step)= K2*u_chol(:,step)+ theta*(b_chol(:, step+1)+ b_chol(:,step));
       %  b_j(:,step)= K2*u_j(:,step)+ theta*(b_j(:, step+1)+ b_j(:,step));
         K1_new= K1;
         K1_new(boundary_conditions(:,1),:)=0;
         b_chol(:,step)=b_chol(:,step)-K1(:, boundary_conditions(:,1))*boundary_conditions(:,2);
        % b_j(:,step)=b_j(:,step)-K1(:, boundary_conditions(:,1))*boundary_conditions(:,2);
         b_chol(boundary_conditions(:,1),step)=boundary_conditions(:,2);
        % b_j(boundary_conditions(:,1),step)=boundary_conditions(:,2);
         K1_new(:,boundary_conditions(:,1))=0;
         for c=1: size(boundary_conditions,1)
             K1_new(boundary_conditions(c,1),boundary_conditions(c,1))=1;
          
         end
        
    L = ichol(sparse(K1_new));
    J=diag(diag(K1_new)); 
 [u_chol(:,step+1), flag_chol, relres_chol, iter, resvec_chol]=pcg(K1_new, b_chol(:,step), 1e-8, 500, L, L');
%[u_j(:,step+1), flag_j, relres_j, iter_j, resvec_j]=pcg(K1_new, b_j(:,step), 1e-8, 500, J, J');   
  
 track_solutions(:, step + 1) = u_chol(track, step + 1);
      
  if abs(t_current - t_max/4) < eps || abs(t_current - t_max/2) < eps || ...
       abs(t_current - 3*t_max/4) < eps || abs(t_current - t_max) < eps
      if abs(t_current - t_max/4) < eps
            col = 1;
        elseif abs(t_current - t_max/2) < eps
            col = 2;
        elseif abs(t_current - 3*t_max/4) < eps
            col = 3;
        else
            col = 4;
        end
        solutions_at_times(:, col) = u_chol(node_indices, step + 1);
    end
      % Store u_new in an array  to visualize the solution later
       % solutions(:, step) = u_new;
end
      
%disp(track_solutions(:,127));
%disp(track_solutions(:,251));
%disp(track_solutions(:,376));
%disp(track_solutions(:,end));
time_steps = [127, 251, 376, size(track_solutions, 2)]; 
solutions_times = track_solutions(:, time_steps);  % Extract the solutions at the specified time steps


row_names = {'u(p1)', 'u(p2)', 'u(p3)'};      
T = table(solutions_times(:,1), solutions_times(:,2), solutions_times(:,3), solutions_times(:,4), ...
    'VariableNames', {'t=2.5', 't=5', 't=7.5', 't_max'}, 'RowNames', row_names);
disp('Solution at tracking points at t=2.5, 5, 7.5, t = t_max:');
disp(T);

figure;
for i = 1:num_track
    plot(linspace(0, t_max, num_steps + 1), track_solutions(i, :));
    
    hold on;
end
xlabel('Time');
ylabel('Solution');
title('Solution at Tracking Points Over Time');
legend('P1', 'P2', 'P3');
hold off;


% Graph solutions at given times along the boundary ΓN
figure;
plot(arc_lengths, solutions_at_times(:, 1), 'DisplayName', 't = tmax/4');
hold on;
plot(arc_lengths, solutions_at_times(:, 2), 'DisplayName', 't = tmax/2');
plot(arc_lengths, solutions_at_times(:, 3), 'DisplayName', 't = 3/4tmax');
plot(arc_lengths, solutions_at_times(:, 4), 'DisplayName', 't = tmax');
hold off;
xlabel('Arc Length');
ylabel('Solution');
title('Solution along the boundary  ΓN at given times');
legend('show');
grid on;


    figure;
    tr = triangulation(topology, coordinates(:,1), coordinates(:,2), u_chol(:, end));
    trisurf(tr, 'FaceColor', 'interp', 'EdgeColor', 'none');
    axis equal tight;    
    view(3);             % Set the view to 3D
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    zlabel('Solution u');
    title('FEM Simulation Results');
    colorbar;    

%% steady state solution
b2_chol = zeros(num_nodes,1);
b2_j = zeros(num_nodes,1);
H_new= H;
u_chol2= zeros(num_nodes,1);
u_j2= zeros(num_nodes,1);
H_new(boundary_conditions(:,1),:)=0;
b2_chol = b2_chol- H_new(:,boundary_conditions(:,1))*boundary_conditions(:,2);
b2_j=b2_j- H_new(:,boundary_conditions(:,1))*boundary_conditions(:,2);
b2_chol(boundary_conditions(:,1))=boundary_conditions(:,2);
b2_j(boundary_conditions(:,1))=boundary_conditions(:,2);
H_new(:,boundary_conditions(:,1))=0;
for c=1: size(boundary_conditions,1)
      H_new(boundary_conditions(c,1),boundary_conditions(c,1))=1;
          
end
L = ichol(H_new);
%J=diag(diag(H_new));
for d =1:size(H,1)
    J(d,d)=H(d,d);
end
[u_steady, flag_chol, relres_chol, iter, resvec_chol]=pcg(H_new,b2_chol, 1e-9, 1000, L, L');
[u_steady_j, flag_j, relres_j, iter_j, resvec_j]=pcg(H_new,b2_j, 1e-9, 2000, J);

figure;
cr = triangulation(topology, coordinates(:,1), coordinates(:,2), u_steady);
trisurf(cr, 'FaceColor', 'interp', 'EdgeColor', [0.8 0.8 0.8], 'EdgeAlpha', 0.4); % Colore grigio chiaro e semi-trasparente per i bordi
%trisurf(cr, 'FaceColor', 'interp', 'EdgeColor', [0.9 0.9 0.9], 'LineWidth', 0.1, 'EdgeAlpha', 0.1);
axis equal tight; 
view(3); 
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Solution u');
title('FEM Simulation Steady State, Cholesky ');
colorbar;


%ERROR ANALYSIS
epsilon=0;
[ f] = nodal_areas(topology, coordinates);
for i=1:size(coordinates, 1)
     epsilon=epsilon+((u_steady(i)-sol(i))^2)*(f(i));
end
Errore(mesh) = (epsilon)^(0.5);


data_displayed = [
    0.2434390, 0.6046775, 0.7454968, 0.7751273;
    0.0772287, 0.2751718, 0.4328630, 0.4805514;
    0.0183037,  0.0928241, 0.1716526, 0.2008722

];

% Difference between given data and the ones obtained stored in T
differenza = abs(T{:,:} - data_displayed);
Differenza_T = array2table(differenza, 'VariableNames', T.Properties.VariableNames, 'RowNames', row_names);
disp('Difference :');
disp(Differenza_T);

h_current = 0;

% Iterating on all the elements K in the triangulation
for k = 1:size(topology, 1)
    % Find indices of the nodes composing element K
    nodes_indices = topology(k, :) - 1; 
    x_K = coordinates(nodes_indices + 1, 1); % +1 for one-based indexing in MATLAB
    y_K = coordinates(nodes_indices + 1, 2);
    distances_K = zeros(length(nodes_indices), length(nodes_indices));
   % Calculate all pairwise distances within the element K
    for i = 1:length(nodes_indices)
        for j = 1:length(nodes_indices)
            if i ~= j % No need to calculate distance from a point to itself
                distances_K(i, j) = sqrt((x_K(i) - x_K(j))^2 + (y_K(i) - y_K(j))^2);
            end
        end
    end
    hK = max(distances_K(:)); % The diameter of the element is the max distance
    h_current = max(h_current, hK); 
end

h(mesh) = h_current;
disp(['Max diam h = ', num2str(h_current)]);

   ratio_c(1)=0;
        if mesh>=2
            ratio_c(mesh)= (Errore(mesh-1)/Errore(mesh))*((h(mesh))/(h(mesh-1)))^2;
        end

   %convergence plot for pcg, both preconditioners
   figure;
   semilogy(0:iter,resvec_chol,'-o');
   xlabel('Number of iterations N')
   ylabel('Residual norm')
   title('Cholesky preconditioner'); 
   hold on;
   semilogy(0:iter_j,resvec_j,'-*');
   xlabel('Number of iterations N')
   ylabel('Residual norm')
   title('Jacobi vs Cholesky preconditioner'); 
   hold off;
end
mesh = [0 1 2 3 4]';
NT = table(mesh,Errore',h', ratio_c');
NT.Properties.VariableNames = {'mesh','Error','h', 'ratio'};
disp(NT);

