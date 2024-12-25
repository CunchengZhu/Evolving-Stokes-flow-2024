function [vertex, velocity] = rm_rigid(vertex, velocity, v_area)
    % remove rigid motion (translation and rotation) from velocity
    % and shift center of mass to origin
    % inputs:
    %   vertex: n_v by 3, vertex coordinate
    %   velocity: n_v by 3, velocity
    %   v_area: n_v by 1, vertex area
    % outputs:
    %   vertex: n_v by 3, vertex coordinate
    %   velocity: n_v by 3, velocity

    % format velocity 
    n_v = size(vertex, 1);
    velocity = reshape(velocity, n_v, 3);

    % remove translation
    mass = sum(v_area); % total mass
    com = sum(vertex .* v_area) / mass; % center of mass
    u_L = sum(velocity .* v_area) / mass; % trans vel

    % shift com and remove translation
    velocity = velocity - repmat(u_L, n_v, 1);
    vertex = vertex - repmat(com, n_v, 1);
    ang_momentum = sum(cross(vertex, velocity, 2).* v_area, 1);
    
    % construct moment of inertia tensor
    Rsq_id = eye(3);
    Rsq_id = permute(Rsq_id(:, :, ones(n_v, 1)), [3, 1, 2]);
    Rsq_id = Rsq_id .* sum(vertex .* vertex, 2);
    VVT = vertex(:, :, ones(3, 1)) .* permute(vertex(:, :, ones(3, 1)), [1, 3, 2]);
    moment = sum((Rsq_id - VVT) .* v_area, 1);

    % remove rotation
    w = (squeeze(moment) \ ang_momentum')'; % angular velocity
    w = repmat(w, n_v, 1);
    velocity = velocity - cross(w, vertex, 2);
    velocity = velocity(:);
end