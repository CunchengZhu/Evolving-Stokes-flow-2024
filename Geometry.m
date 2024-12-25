classdef Geometry
    properties
        V {mustBeNumeric};
        F {mustBeNumeric};
        mesh
        volume {mustBeNumeric};
        area {mustBeNumeric};
        f_area {mustBeNumeric};
        f_normal {mustBeNumeric};
        f_basis_u {mustBeNumeric};
        f_basis_v {mustBeNumeric};
        f_center {mustBeNumeric};
        v_area {mustBeNumeric};
        he_length {mustBeNumeric};
        he_dihedral {mustBeNumeric};
        he_corner {mustBeNumeric};
        he_cotan_weight {mustBeNumeric};
        he_bary_weight {mustBeNumeric};
        he_mean_curvature {mustBeNumeric};
        v_mean_curvature {mustBeNumeric};
        v_gaussian_curvature {mustBeNumeric};
        f_gaussian_curvature {mustBeNumeric};
        he_mean_curvature_vec {mustBeNumeric};
        he_gaussian_curvature_vec {mustBeNumeric};
        he_schlafli_vec1 {mustBeNumeric};
        he_schlafli_vec2 {mustBeNumeric};
        v_gaussian_curvature_vec {mustBeNumeric};
        v_mean_curvature_vec {mustBeNumeric};
        
        f_f_dihedral;
        grad;
        mass2;
        mass0;
        lap;
    end
    methods
        function obj = Geometry(F, V)
            obj.mesh = Mesh(F, V);
            obj.V = V;
            obj.F = F;

            obj.volume = obj.signed_volume();
            obj.f_center = obj.center();
            [obj.f_area, obj.f_normal, obj.f_basis_u, obj.f_basis_v] = obj.face();
            obj.area = sum(obj.f_area);

            [v_area_sum, ~] = obj.mesh.face_to_vertex(obj.f_area);
            obj.v_area = v_area_sum ./ 3;
            
            [obj.he_length, obj.he_dihedral, obj.f_f_dihedral] = obj.dihedral();
            obj.he_corner = obj.corner_angle();
            [obj.he_cotan_weight, obj.he_bary_weight] = obj.weight();
            [obj.v_gaussian_curvature, obj.f_gaussian_curvature] = obj.gaussian_curvature();
            [obj.he_mean_curvature, obj.v_mean_curvature] = obj.mean_curvature();
            [obj.mass0, obj.mass2, obj.grad, obj.lap] = obj.operators();

            obj.he_mean_curvature_vec = obj.mean_curvature_vec();
            obj.he_gaussian_curvature_vec = obj.gaussian_curvature_vec();
            [obj.he_schlafli_vec1, obj.he_schlafli_vec2] = obj.schlafli_vec();

            [obj.v_gaussian_curvature_vec,~] = obj.mesh.halfedge_to_vertex(obj.he_gaussian_curvature_vec);
            [obj.v_mean_curvature_vec,~] = obj.mesh.halfedge_to_vertex(obj.he_mean_curvature_vec);

        end
        
        function [volume] = signed_volume(obj)
            v = @(i) obj.V(obj.F(:, i), :);
            volume = sum(dot(v(1), cross(v(2), v(3), 2), 2)) / 6;
        end

        function [center] = center(obj)
            center = (obj.V(obj.F(:,1),:) + obj.V(obj.F(:,2),:) + obj.V(obj.F(:,3),:)) ./ 3;
        end

        function [area, normal, f_basis_u, f_basis_v] = face(obj)
            one = obj.V(obj.F(:,2),:)-obj.V(obj.F(:,1),:);
            second = obj.V(obj.F(:,3),:)-obj.V(obj.F(:,2),:);
            cross_product = cross(one, second, 2);
            area = 0.5 * sqrt(sum(cross_product.^2, 2));
            normal = 0.5 * cross_product ./ area; % ccw oriented
            f_basis_u = one ./ vecnorm(one, 2, 2);
            f_basis_v = - cross(normal, f_basis_u, 2);
        end
        
        function [he_length, he_dihedral, f_f_dihedral] = dihedral(obj)
            he_face2 = obj.mesh.he_face(obj.mesh.he_flip);
            n1 = obj.f_normal(obj.mesh.he_face, :);
            n2 = obj.f_normal(he_face2, :);
            e = obj.V(obj.mesh.he_dst, :) - obj.V(obj.mesh.he_src, :);
            he_length = vecnorm(e, 2, 2);
            e = e ./ he_length;
            he_dihedral = atan2(dot(cross(n1,n2,2),e,2),dot(n1,n2,2));
            f_f_dihedral = sparse(obj.mesh.he_face, he_face2, he_dihedral, obj.mesh.n_f, obj.mesh.n_f);
        end
        
        function [he_corner] = corner_angle(obj)
            length = reshape(obj.he_length, [], 3);
            l1 = length(:, 1);
            l2 = length(:, 2);
            l3 = length(:, 3);
            c1 =  acos((l3.^2 + l1.^2 - l2.^2)./(2.*l3.*l1));
            c2 = acos((l1.^2 + l2.^2 - l3.^2)./(2.*l1.*l2));
            c3 = acos((l2.^2 + l3.^2 - l1.^2)./(2.*l2.*l3));
            he_corner = [c1; c2; c3];
        end

        
        function [he_weight, he_bary_weight] = weight(obj)
            he_weight = 0.5 * cot(obj.he_corner(obj.mesh.he_prev));
            l1 = obj.he_length;
            a = obj.f_area(obj.mesh.he_face);
            % he_weight = (l2.^2 + l3.^2 - l1.^2)./(8.* area); % cotan
            he_bary_weight = 2./3.* a./(l1.^2);
        end

        function [v_K, f_K] = gaussian_curvature(obj)
            v_angle_sum = accumarray(obj.mesh.he_src, obj.he_corner, [obj.mesh.n_v, 1]);
            v_alpha = 2 * pi ./ v_angle_sum; % polar_turning_factor
            v_K = 2 * pi - v_angle_sum;
            
            he_alpha = v_alpha(obj.mesh.he_src);
            f_angle_sum = accumarray(obj.mesh.he_face, obj.he_corner .* he_alpha, [obj.mesh.n_f, 1]);
            f_K = f_angle_sum - pi;
        end
        
        function [he_mean_curvature, v_mean_curvature] = mean_curvature(obj)
            he_mean_curvature = 0.25 * obj.he_length .* obj.he_dihedral;
            v_mean_curvature = accumarray(obj.mesh.he_src, he_mean_curvature, [obj.mesh.n_v, 1]);
        end
        
        function [mass0, mass2, grad, lap]  = operators(obj)
            % compute DEC operators 
            % Outputs:
            %   mass0: n_v by n_v mass matrix
            %   mass2: n_f by n_f mass matrix
            %   grad: n_f by n_v by 3 gradient sparse tensor
            %   lap: n_v by n_v Laplacian matrix
            edge = obj.V(obj.mesh.he_dst, :) - obj.V(obj.mesh.he_src, :);
            normal = [obj.f_normal; obj.f_normal; obj.f_normal];
            edge_perp = 0.5 * reshape(cross(normal, edge, 2), [3 * obj.mesh.n_he, 1]);
            point_ind = [obj.F(:, 3); obj.F(:, 1); obj.F(:, 2)];
            point_ind = [point_ind; point_ind; point_ind];
            vec_ind = [ones(obj.mesh.n_he, 1); 2 * ones(obj.mesh.n_he, 1); 3 * ones(obj.mesh.n_he, 1)];
            face_ind = [obj.mesh.he_face; obj.mesh.he_face; obj.mesh.he_face];
            coords = [face_ind, point_ind, vec_ind];
            I_f = (1:obj.mesh.n_f)'; 
            
            % mass is (F, F)
            mass2_inv = sptensor([I_f, I_f], 1./obj.f_area, [obj.mesh.n_f, obj.mesh.n_f]);
            mass2 = sparse(I_f, I_f, obj.f_area, obj.mesh.n_f, obj.mesh.n_f);

            % mass0 is (V, V)
            mass0 = sparse(1:obj.mesh.n_v, 1:obj.mesh.n_v, obj.v_area, obj.mesh.n_v, obj.mesh.n_v);
            
            % grad is (F, V, T*M)
            grad = sptensor(coords, edge_perp, [obj.mesh.n_f, obj.mesh.n_v, 3]);
            
            % lap is (V, V)
            lap = ttt(grad, ttt(mass2_inv, grad, 2, 1), [1, 3], [1, 3]);
            coords = [lap.subs(:, 1), lap.subs(:, 2)];
            lap = sparse(coords(:, 1), coords(:, 2), lap.vals, obj.mesh.n_v, obj.mesh.n_v);
        end
        
        function he_mean_curvature_vec = mean_curvature_vec(obj)
            % Compute the mean curvature vector
            % Outputs:
            %   he_mean_curvature_vec: n_he by 3 mean curvature vector
            h1 = obj.mesh.he_next;
            h2 = obj.mesh.he_next(obj.mesh.he_next(obj.mesh.he_flip));
            v1 = obj.V(obj.mesh.he_dst(h1), :) - obj.V(obj.mesh.he_src(h1), :);
            v2 = obj.V(obj.mesh.he_dst(h2), :) - obj.V(obj.mesh.he_src(h2), :);
            n1 = obj.f_normal(obj.mesh.he_face, :);
            n2 = obj.f_normal(obj.mesh.he_face(obj.mesh.he_flip), :);
            he_mean_curvature_vec = -0.125 * (cross(v1, n1, 2) + cross(v2, n2, 2));
        end

        function he_gaussian_curvature_vec = gaussian_curvature_vec(obj)
            % Compute the Gaussian curvature vector
            % Outputs:
            %   he_gaussian_curvature_vec: n_he by 3 Gaussian curvature vector
            v1 = obj.V(obj.mesh.he_dst, :) - obj.V(obj.mesh.he_src, :);
            v1_prime = v1 ./ vecnorm(v1, 2, 2);
            he_gaussian_curvature_vec = -0.5 * v1_prime .* obj.he_dihedral;
        end

        function [he_schlafli_vec1, he_schlafli_vec2] = schlafli_vec(obj)
            % Compute the Schlafli vector, cf 
            % https://www.sciencedirect.com/science/article/pii/S2667074722000192
            % Outputs:
            %   he_schlafli_vec1: n_he by 3 Schlafli vector 1
            %   he_schlafli_vec2: n_he by 3 Schlafli vector 2
            h1a = obj.mesh.he_next;
            h1b = obj.mesh.he_flip;
            h2a = obj.mesh.he_next(h1a);
            h2b = obj.mesh.he_prev(h1b);
            na = obj.f_normal(obj.mesh.he_face, :);
            nb = obj.f_normal(obj.mesh.he_face(obj.mesh.he_flip), :);
            w1a = 1 ./ tan(obj.he_corner(h1a));
            w1b = 1 ./ tan(obj.he_corner(h1b));
            w2a = 1 ./ tan(obj.he_corner(h2a));
            w2b = 1 ./ tan(obj.he_corner(h2b));
            he_schlafli_vec1 = 0.5 * (na .* w1a + nb .* w1b);
            he_schlafli_vec2 = -0.5 * (na .* w2a + nb .* w2b);
        end 

        function v_force = bending_force(obj, Kb)
            % Compute the bending force
            % Inputs:
            %   Kb: bending modulus
            % Outputs:
            %   v_force: n_v by 3 bending force

            Hi = obj.v_mean_curvature(obj.mesh.he_src) ./ obj.v_area(obj.mesh.he_src);
            Hj = obj.v_mean_curvature(obj.mesh.he_dst) ./ obj.v_area(obj.mesh.he_dst);
            he_force = Kb * ( - obj.he_gaussian_curvature_vec .* (Hi + Hj) ...
                + 2 * obj.he_mean_curvature_vec .* (Hi.^2 / 3 + Hj.^2 / 3 * 2) ...
                - (obj.he_schlafli_vec1 .* Hi + obj.he_schlafli_vec2 .* Hj));
            v_force = zeros(obj.mesh.n_v, 3);
            for i = 1:3
                v_force(:, i) = accumarray(obj.mesh.he_src, he_force(:, i), [obj.mesh.n_v, 1]);
            end
        end

        function energy = willmore_energy(obj, Kb)
            H = obj.v_mean_curvature ./ obj.v_area; 
            energy = Kb * (H .* obj.v_area)' * H;
        end

        function [pgrad, K, W, div, KTK, DTD] = evolving_operators(obj)
            % Compute the differential operators of evolving surface
            % Outputs:
            %   pgrad: 5d sparse tensor, (F, V, T*M, R3*, T*M) extrinsic gradient 
            %   K: 5d sparse tensor, (F, V, [T*M, R3*, T*M]) Killing
            %   operator (symmetric part of pgrad)
            %   W: 5d sparse tensor, (F, V, \T*M, R3*, T*M\) antisymmetric part of pgrad
            %   div: F by 3*V sparse matrix, (F, V, R3*) divergence operator
            %   KTK: 3*V by 3*V sparse matrix, (V, R3, V, R3) viscosity Laplacian
            %   DTD: 3*V by 3*V sparse matrix, gradient of divergence operator

            % n_outer_n_grad is (F, V, T*M, R3*, TN)
            grad_ = squeeze(ttt(ttt(obj.grad, sptensor(ones(3, 1))), sptensor(ones(3,1))));
            n_outer_n = (obj.f_normal(:, :, ones(1,3)) .* permute(obj.f_normal(:,:,ones(1,3)), [1,3,2]));            
            nonzero = sub2ind(size(n_outer_n), grad_.subs(:, 1), grad_.subs(:, 4), grad_.subs(:, 5));
            n_outer_n_grad = sptensor(grad_.subs, n_outer_n(nonzero), size(grad_)) .* grad_;

            % pgrad is (F, V, T*M, R3*, TM) = (F, V, T*M, R3*, R3) - (F, V, T*M, R3*, TN)
            id = sptensor(eye(3));
            pgrad = ttt(obj.grad, id) - n_outer_n_grad;
  
            % K is (F, V, [T*M, R3*, T*M])
            K = (pgrad + permute(pgrad, [1, 2, 5, 4, 3])) / 2;

            % W is (F, V, \T*M, R3*, T*M\)
            W = pgrad - K;

            % div is (F, V, R3*)
            div = ttt(K, id, [3, 5],[1, 2]);

            % L0 = KTK is (V, R3, V, R3)
            I_f = (1:obj.mesh.n_f)'; 
            mass2_inv = sptensor([I_f, I_f], 1./obj.f_area, [obj.mesh.n_f, obj.mesh.n_f]);
            KTK = ttt(K, ttt(mass2_inv, K, 2, 1), [1, 3, 5], [1, 3, 5]);

            % see https://www.tensortoolbox.org/sptenmat_doc.html#20 for
            % sptensor to sparse conversion
            coords = [KTK.subs(:, 1) + (KTK.subs(:, 2) - 1) * obj.mesh.n_v, KTK.subs(:, 3) + (KTK.subs(:, 4) - 1) * obj.mesh.n_v];
            KTK = sparse(coords(:, 1), coords(:, 2), KTK.vals, 3 * obj.mesh.n_v, 3 * obj.mesh.n_v);
            coords = [div.subs(:, 1), div.subs(:, 2) + (div.subs(:, 3) - 1) * obj.mesh.n_v];
            div = sparse(coords(:, 1), coords(:, 2), div.vals, obj.mesh.n_f, 3 * obj.mesh.n_v);
            mass2_inv = sparse(I_f, I_f, 1./obj.f_area, obj.mesh.n_f, obj.mesh.n_f);
            DTD = div' * mass2_inv * div;
        end

    end
end