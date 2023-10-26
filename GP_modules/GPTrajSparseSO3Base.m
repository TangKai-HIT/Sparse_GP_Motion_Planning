classdef GPTrajSparseSO3Base < handle
    %GPTRAJSPARSESO3BASE base class with abstract methods for sparse GP trajectory on SO(3) , dim = 6
    %   sparse GP trajectory with O(1) query time complexity, derived from local Guass-Markov SDE with sparse kernel on SO(3)
    
    properties
        support_SO3; %support/data SO3, 3 X 3 X N 
        support_so3; %support/data so3 in local frame, 3 X N 
        supportId; %supportId index
        Q_c; %hyperparameters, white noise covariance
        kappa0; %kappa of initial state, dim X dim
    end

    properties(SetAccess=protected, GetAccess=public)
        numsupports; %number of support/data states
        dim = 6; %dimension of support/data states

        %lifted form parameters
        supportKappa; %kernel matrix of support states, dim X dim
        supportKappa_inv; %inverse of kernel matrix

        % util variables
        supportA; %(numsupports*dim) X (numsupports*dim), lower triangular mat
        supportA_inv; %inverse of supportA
        supportQ; %(numsupports*dim) X (numsupports*dim), block diagonal mat
        supportQ_inv; %inverse of supportQ
    end
    
    methods(Abstract)
        transMat(obj, t, s); %state transition matrix (s->t)
        QMat(obj, a, b); %Q matrix, Q = intergrate Phi(b,s)*F(s)*Q_c*F(s)'*Phi(b,s)'*ds from a to b, self-defined
         VVec(obj, a, b); %V vector, V = intergrate Phi(b,s)*v(s)*ds from a to b, self-defined
    end

    %public methods
    methods
        function obj = GPTrajSparseSO3Base(support_SO3, support_so3, supportId, Q_c, kappa0)
            %GPTRAJSPARSE 构造此类的实例
            %   support_SO3: 3 X 3 X N 
            %   support_so3: 3 X N 

                obj.support_SO3 = support_SO3;
                obj.support_so3 = support_so3;
                obj.supportId = supportId;
                obj.Q_c = Q_c;
                obj.numsupports = size(support_SO3, 3);
                obj.kappa0 = kappa0;
        end
       
        function [query_SO3, query_so3, queryKappa] = query(obj, queryId)
            %QUERY query function
            %   fast O(1) query interpolation, note: do not support extrapolation yet  
            [l, u, l_id, u_id]  = findInterval(queryId, obj.supportId); %find t_i, t_{i+1}

            Q_query = obj.QMat(l, queryId);

            rowId_u = l_id*obj.dim+1 : u_id*obj.dim; %term in i+1 index
            rowId_l = (l_id-1)*obj.dim+1 : l_id*obj.dim; %term in i index
            Q_u_inv = obj.supportQ_inv(rowId_u, rowId_u); %Q_{i+1} inverse
        
            Phi_tau_u = obj.transMat(u, queryId);
            Psi_query = Q_query * Phi_tau_u' * Q_u_inv;

            Phi_l_tau = obj.transMat(queryId, l);
            Phi_l_u = obj.transMat(u, l);
            Lambda_query = Phi_l_tau - Psi_query * Phi_l_u;
            
            %interpolation
            [gamma_l, gamma_u] = obj.gammaVec(l_id, u_id);
            gamma_query = Lambda_query*gamma_l + Psi_query*gamma_u;

            R_l = obj.support_SO3(:, :, l_id);
            query_SO3 = R_l * expMapSO3(gamma_query(1:3));

            eps_l_query = logMapSO3(R_l' * query_SO3); % eps_i(tau)
            jacobR_query = jacobRSO3(eps_l_query);
            query_so3 = jacobR_query * gamma_query(4:6);
            
            %if require kappa query
            if nargout>1
                % Phi_0_tau = obj.transMat(queryId, obj.supportId(1));
                % kappa_tau_tau = Phi_0_tau * obj.kappa0 * Phi_0_tau' + obj.QMat(obj.supportId(1), queryId);
                kappa_tau_tau = Phi_l_tau * obj.supportKappa(rowId_l, rowId_l) * Phi_l_tau' + Q_query;

                Phi_l_vec = obj.supportA(rowId_l, 1:l_id*obj.dim); %block row mat: Phi(ti, t0~t_{i-1})
                Q_l_diag = obj.supportQ(1:l_id*obj.dim, 1:l_id*obj.dim); %block diag mat: Q_0 (kappa0) ~ Q_t_{i-1}
                Q_l = obj.supportQ(rowId_l, rowId_l); %Q_t_i

                H = Phi_l_vec * Q_l_diag * Phi_l_vec';
                kappa_post = Phi_l_tau * H * Phi_l_tau' + Psi_query * Phi_tau_u * Q_query';

                queryKappa = kappa_tau_tau - kappa_post;
            end
        end
        
        function [gamma_l, gamma_u] = gammaVec(obj, l_id, u_id) 
            %GAMMAVEC local support vectors (dim X 1) on TS of SO3 at support index l_id, i.e. gmma_i, gamma_{i+1}
            omega_l = obj.support_so3(:, l_id);
            omega_u = obj.support_so3(:, u_id);
            R_l = obj.support_SO3(:, :, l_id);
            R_u = obj.support_SO3(:, :, u_id);
            eps_l_u = logMapSO3(R_l' * R_u); % eps_i(t_{i+1})

            gamma_l = [zeros(3,1); omega_l];
            gamma_u = [eps_l_u; jacobRInvSO3(eps_l_u) * omega_u];
        end

        function updateLifted(obj)
            %UPDATELIFTED update lifted form parameters by support states
            obj.supportA = eye(obj.numsupports*obj.dim);
            obj.supportA_inv = eye(obj.numsupports*obj.dim);
            obj.supportQ = zeros(obj.numsupports*obj.dim, obj.numsupports*obj.dim);
            obj.supportQ_inv = zeros(obj.numsupports*obj.dim, obj.numsupports*obj.dim);

            obj.supportQ(1:obj.dim, 1:obj.dim) = obj.kappa0;
            obj.supportQ_inv(1:obj.dim, 1:obj.dim) = inv(obj.kappa0);

            M = obj.numsupports - 1;

            for i=1:M
                Q_i = obj.QMat(obj.supportId(i), obj.supportId(i+1));
                rowId = i*obj.dim+1:(i+1)*obj.dim;

                obj.supportQ(rowId, rowId) = Q_i;
                obj.supportQ_inv(rowId, rowId) = inv(Q_i);

                for j = i:-1:1
                    colunmId = (j-1)*obj.dim+1:j*obj.dim;
                    obj.supportA(rowId, colunmId) = obj.transMat(obj.supportId(i+1), obj.supportId(j));

                    if j==i
                        obj.supportA_inv(rowId, colunmId) = -obj.supportA(rowId, colunmId); %inverse block
                    end
                end
            end

            %support kappa
            obj.supportKappa = obj.supportA * obj.supportQ * obj.supportA';
            obj.supportKappa_inv = obj.supportA_inv' * obj.supportQ_inv * obj.supportA_inv;
        end

    end
    
end

