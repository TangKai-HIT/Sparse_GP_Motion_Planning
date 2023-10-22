classdef GPTrajSparseBase < handle
    %GPTRAJSPARSEBASE base class with abstract methods for sparse GP trajectory
    %   sparse GP trajectory with O(1) query time complexity, derived from Guass-Markov SDE with sparse kernel
    
    properties
        supportPts; %support/data points
        supportId; %supportId index
        Q_c; %hyperparameters, white noise covariance
        mu0; %mu of initial state, dim X 1
        kappa0; %kappa of initial state, dim X dim
    end

    properties(SetAccess=protected, GetAccess=public)
        numsupports; %number of support/data states
        dim; %dimension of support/data states

        %lifted form parameters
        supportKappa; %kernel matrix of support states, dim X dim
        supportKappa_inv; %inverse of kernel matrix
        supportMu; %Mu of support states, dim X 1

        % util variables
        supportA; %(numsupports*dim) X (numsupports*dim), lower triangular mat
        supportA_inv; %inverse of supportA
        supportV; %(numsupports*dim) X 1, vec
        supportQ; %(numsupports*dim) X (numsupports*dim), block diagonal mat
        supportQ_inv; %inverse of supportQ
    end
    
    %abstract methods
    methods(Abstract)
        % query(obj, queryId); %query function
        transMat(obj, t, s); %state transition matrix (s->t)
        QMat(obj, a, b); %Q matrix, Q = intergrate Phi(b,s)*F(s)*Q_c*F(s)'*Phi(b,s)'*ds from a to b, self-defined
         VVec(obj, a, b); %V vector, V = intergrate Phi(b,s)*v(s)*ds from a to b, self-defined
    end

    %public methods
    methods
        function obj = GPTrajSparseBase(supportPts, supportId, Q_c, mu0, kappa0)
            %GPTRAJSPARSE 构造此类的实例
            %   supportPts: N X dim

                obj.supportPts = supportPts;
                obj.supportId = supportId;
                obj.Q_c = Q_c;
                [obj.numsupports, obj.dim] = size(supportPts);
                obj.mu0 = mu0;
                obj.kappa0 = kappa0;
        end
       
        function [queryPt, queryKappa] = query(obj, queryId)
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
            mu_l = obj.supportMu(rowId_l, 1);
            mu_u = obj.supportMu(rowId_u, 1);
            xbar_l = obj.supportPts(l_id, :)';
            xbar_u = obj.supportPts(u_id, :)';
            mu_query = Phi_l_tau * mu_l + obj.VVec(l, queryId);
            queryPt = (mu_query + Lambda_query*(xbar_l-mu_l) + Psi_query*(xbar_u-mu_u))';
            
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

        function updateLifted(obj)
            %UPDATELIFTED update lifted form parameters by support states
            obj.supportA = eye(obj.numsupports*obj.dim);
            obj.supportA_inv = eye(obj.numsupports*obj.dim);
            obj.supportQ = zeros(obj.numsupports*obj.dim, obj.numsupports*obj.dim);
            obj.supportQ_inv = zeros(obj.numsupports*obj.dim, obj.numsupports*obj.dim);
            obj.supportV = zeros(obj.numsupports*obj.dim, 1);

            obj.supportQ(1:obj.dim, 1:obj.dim) = obj.kappa0;
            obj.supportQ_inv(1:obj.dim, 1:obj.dim) = inv(obj.kappa0);
            obj.supportV(1:obj.dim, 1) = obj.mu0;

            M = obj.numsupports - 1;

            for i=1:M
                Q_i = obj.QMat(obj.supportId(i), obj.supportId(i+1));
                rowId = i*obj.dim+1:(i+1)*obj.dim;

                obj.supportQ(rowId, rowId) = Q_i;
                obj.supportQ_inv(rowId, rowId) = inv(Q_i);
                
                obj.supportV(rowId, 1) = obj.VVec(obj.supportId(i), obj.supportId(i+1));

                for j = i:-1:1
                    colunmId = (j-1)*obj.dim+1:j*obj.dim;
                    obj.supportA(rowId, colunmId) = obj.transMat(obj.supportId(i+1), obj.supportId(j));

                    if j==i
                        obj.supportA_inv(rowId, colunmId) = -obj.supportA(rowId, colunmId); %inverse block
                    end
                end
            end

            %support mu & kappa
            obj.supportMu = obj.supportA * obj.supportV;
            obj.supportKappa = obj.supportA * obj.supportQ * obj.supportA';
            obj.supportKappa_inv = obj.supportA_inv' * obj.supportQ_inv * obj.supportA_inv;
        end

    end
    
end

