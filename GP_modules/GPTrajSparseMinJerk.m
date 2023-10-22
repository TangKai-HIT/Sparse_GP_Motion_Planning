classdef GPTrajSparseMinJerk < GPTrajSparseBase
    %GPTRAJSPARSEMINJERK sparse GP trajectory with zeros jerk prior
    %   sparse GP trajectory with O(1) query time complexity, derived from Guass-Markov SDE with sparse kernel
    
    properties(SetAccess=private, GetAccess=public)
        % end state variables (for whole trajectory conditioning)
        % (to do ... , note that conditioning on the observation at the end point won't affect the query of mid local points)
    end
    
    methods
        function obj = GPTrajSparseMinJerk(supportPts, supportId, Q_c, mu0, kappa0)
            %GPTRAJSPARSEMINJERK 构造此类的实例
            %   supportId must include initial condition
            if nargin == 0
                supportPts=1; supportId=1; Q_c=1;
                mu0 = 0; kappa0 = 1;
            end
            obj@GPTrajSparseBase(supportPts, supportId, Q_c, mu0, kappa0); %father class init
        end

        function Phi = transMat(obj, t, s) 
            %TRANSMAT state transition matrix (s->t)
            Phi = transMatMinJerk(t, s, obj.dim/3);
        end

        function Q = QMat(obj, a, b) 
            %QMAT Q matrix, Q = intergrate Phi(b,s)*F(s)*Q_c*F(s)'*Phi(b,s)'*ds from a to b
            Q = MatQ_minJerk(a, b, obj.Q_c);
        end

        function V = VVec(obj, a, b)
            %VVec, V vector, V = intergrate Phi(b,s)*v(s)*ds from a to b, self-defined
            V = zeros(obj.dim, 1);
        end

    end
end

