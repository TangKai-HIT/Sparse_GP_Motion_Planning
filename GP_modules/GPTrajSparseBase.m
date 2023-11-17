classdef GPTrajSparseBase < handle
    %GPTRAJSPARSEBASE base class with abstract methods for sparse GP trajectory in Euclidean space
    %   sparse GP trajectory with O(1) query time complexity, derived from Guass-Markov SDE with sparse kernel
    
    properties
        supportPts; %support/data points, N X dim
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
        supportMu; %Mu of support states, (numsupports*dim) X 1

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
                queryKappa = Q_query - Psi_query * Phi_tau_u * Q_query';
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
    
        function [tau_samp, Lambda_samp, Psi_samp, mu_samp, sampPts, kappa_samp] = upSampleInterval_uniform(obj, l_id, u_id, n_sp)
            %UPSAMPLEINTERVAL_UNIFORM uniformly upsample in a support state interval  
            %   Inputs:
            %       l_id, u_id: lower id and upper id of the interval 
            %       n_sp: number of samples
            %   Outputs:
            %       tau_samp: 1 X n_sp
            %       Lambda_samp, Psi_samp, mu_samp, sampPts(optional): array, (n_sp*dim) X 1
            %       mu_samp: (n_sp*dim) X 1, mu of sampled points 
            %       sampPts (optional): (n_sp*dim) X 1, sampled points 
            %       kappa_samp (optional): (n_sp*dim) X (n_sp*dim)
            
            l = obj.supportId(l_id);
            u = obj.supportId(u_id);

            tau_samp = linspace(l, u, n_sp+2);
            tau_samp = tau_samp(2:end-1);

            if nargout>5 %optional output
                [Lambda_samp, Psi_samp, mu_samp, sampPts, kappa_samp] = obj.upSampleInterval(l_id, u_id, tau_samp);
            elseif nargout>4 %optional output
                [Lambda_samp, Psi_samp, mu_samp, sampPts] = obj.upSampleInterval(l_id, u_id, tau_samp);
            else
                [Lambda_samp, Psi_samp, mu_samp] = obj.upSampleInterval(l_id, u_id, tau_samp);
            end
        end

        function [Lambda_samp, Psi_samp, mu_samp, sampPts, kappa_samp] = upSampleInterval(obj, l_id, u_id, tau_samp)
            %UPSAMPLEINTERVAL upsample in a support state interval given samples of index 
            %   Inputs:
            %       l_id, u_id: lower id and upper id of the interval 
            %       tau_samp: interpolation index samples
            %   Outputs:
            %       Lambda_samp, Psi_samp, mu_samp, sampPts(optional): array, (n_sp*dim) X 1
            %       mu_samp: (n_sp*dim) X 1, mu of sampled points 
            %       sampPts (optional): (n_sp*dim) X 1, sampled points 
            %       kappa_samp (optional): (n_sp*dim) X (n_sp*dim)
            
            l = obj.supportId(l_id);
            u = obj.supportId(u_id);
            n_sp = length(tau_samp);

            %compute intermedian terms
            Q_l_samp = zeros(n_sp*obj.dim, n_sp*obj.dim);
            Phi_samp_u_T = zeros(n_sp*obj.dim, obj.dim);
            Phi_l_samp = zeros(n_sp*obj.dim, obj.dim);
            VVec_samp = zeros(n_sp*obj.dim, 1);

            for i=1:n_sp
                cur_id = (i-1)*obj.dim+1 : i*obj.dim;
                Q_l_samp(cur_id, cur_id) = obj.QMat(l, tau_samp(i));
                Phi_samp_u_T(cur_id, :) = obj.transMat(u, tau_samp(i))';
                Phi_l_samp(cur_id, :) = obj.transMat(tau_samp(i), l);
                VVec_samp(cur_id, :) = obj.VVec(l, tau_samp(i));
            end

            rowId_u = (u_id-1)*obj.dim+1 : u_id*obj.dim; %term in i+1 index
            rowId_l = (l_id-1)*obj.dim+1 : l_id*obj.dim; %term in i index

            if u_id>l_id+1 %support state not adjancent (vital)
                Q_u_inv = inv(obj.QMat(l, u));
            else %support state adjacent
                Q_u_inv = obj.supportQ_inv(rowId_u, rowId_u); %Q_{i+1} inverse
            end
            
            %compute interpolator
            gamma_samp = Q_l_samp * Phi_samp_u_T; %get gamma_{1~n_sp} sampled vector
            Psi_samp = gamma_samp * Q_u_inv; %get Psi_{1~n_sp} sampled vector

            Phi_l_u = obj.transMat(u, l);
            Lambda_samp = Phi_l_samp - Psi_samp * Phi_l_u; %get Lambda_{1~n_sp} sampled vector

            %interpolation
            mu_l = obj.supportMu(rowId_l, 1);
            mu_u = obj.supportMu(rowId_u, 1);
            xbar_l = obj.supportPts(l_id, :)';
            xbar_u = obj.supportPts(u_id, :)';
            mu_samp = Phi_l_samp * mu_l + VVec_samp;
            
            if nargout>3 %optional output
                sampPts = mu_samp + Lambda_samp*(xbar_l-mu_l) + Psi_samp*(xbar_u-mu_u);
            end

            if nargout>4 %optional output
                %A--trans mat
                A = eye(n_sp*obj.dim);
                for i = 2:n_sp
                    cur_rowId = (i-1)*obj.dim+1:i*obj.dim;
                    for j = 1:i-1
                        cur_colId = (j-1)*obj.dim+1:j*obj.dim; 
                        A(cur_rowId, cur_colId) = obj.transMat(tau_samp(i), tau_samp(j));
                    end
                end
                
                H = A*Q_l_samp;
                kappa_samp = (H + H' - Q_l_samp) - Psi_samp*gamma_samp';
            end
        end

    end
    
end

