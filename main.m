verbose = false;
%%% directory
dir = "./data/"; 
[status, msg, msgID] = mkdir(dir); 

%% system init
%%% 0 for new simulation, otherwise continue from geo"start".mat
start = 2; 
if start == 0
    %%% geometry
    load("./spheroid.mat")
    geo = Geometry(M, P);
    %%% parameters
    p.dt = 5e-2; % time 
    p.kappa = 5e-2; % bending
    p.T = 200; % total time (frames)
    %%% optimizer
    o.h = 0.2; o.eta = 1e2;
    o.tol_f = 1e-5; o.tol_d = 1e-5; o.max_iter = 10000;
    %%% initialize
    velocity = zeros(size(P, 1) * 3, 1);
    pressure = zeros(size(M, 1), 1);
else
    load(dir + sprintf("geo%d.mat", start), ...
     "M", "P", "velocity", "pressure", "p", "o", "r"); 
    geo = Geometry(M, P);
end

%% main loop (not supposed to be modified)
for t = (start + 1):p.T
    %%% incremental potential minimization
    [~, K, ~, div, KTK, DTD] = geo.evolving_operators();
    P0 = P(:); P(:) = P0 + p.dt * velocity;
    eps_f = Inf; eps_d = Inf; j = 0;
    while ((eps_f > o.tol_f) || (eps_d > o.tol_d)) && (j < o.max_iter)
        %%% update bending force
        fb = geo.bending_force(p.kappa);
        b = - 2 * KTK * (P(:) - P0) - p.dt * (- fb(:) - div' * pressure);
        expan = div * (P(:) - P0);
        %%% measure residual
        eps_f = norm_f(b, geo.v_area); eps_d = norm_d(expan, geo.f_area);
        if verbose fprintf("t = %d, j = %d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d); end
        %%% preconditioner
        H = blkdiag(geo.lap, geo.lap, geo.lap);
        %%% gradient descent/ascent
        P(:) = P(:) + o.h * (H \ b);
        pressure = pressure - o.eta * o.h * expan ./ geo.f_area;
        %%% update geometry
        geo = Geometry(M, P);
        j = j + 1;
    end
    %%% save data
    [P, velocity] = rm_rigid(P, (P(:) - P0) / p.dt, geo.v_area);
    save(dir + sprintf("geo%d.mat", t), "M", "P", "velocity", "pressure", "fb", "p", "o", "r");
    fprintf("Save geo%d.mat at j =%d, eps_f = %0.4g, eps_d = %0.4g \n", t, j, eps_f, eps_d);
    %%% update geometry, remesh if needed
    geo = Geometry(M, P);
end

%% helper functions
function e = norm_f(b, area)
    % normalized norm of vector field
    e = sqrt(sum(b.^2 ./ [area; area; area])) / sum(area);
end

function e = norm_d(expan, area)
    % normalized norm of scalar field
    e = sqrt(sum(expan.^2 ./ area)) / sum(area);
end
