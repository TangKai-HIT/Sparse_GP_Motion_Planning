classdef GPTrajSparseSO3MinAcc < GPTrajSparseSO3Base
    %GPTRAJSPARSESO3MINACC sparse GP trajectory with zeros acceleration prior
    %   sparse GP trajectory with O(1) query time complexity, derived from Guass-Markov SDE with sparse kernel
    
    properties(SetAccess=private, GetAccess=public)
        % (reserved)
    end
    
    methods
        function obj = GPTrajSparseSO3MinAcc(support_SO3, support_so3, supportId, Q_c, kappa0)
            %GPTRAJSPARSEMINJERK 构造此类的实例
            %   supportId must include initial condition
            %   support_SO3: 3 X 3 X N 
            %   support_so3: 3 X N 

            if nargin == 0
                support_SO3=eye(3); supportId=1; Q_c=eye(3);
                support_so3 = [0;0;0];
                kappa0 = 1e-2*eye(6);
            end
            obj@GPTrajSparseSO3Base(support_SO3, support_so3, supportId, Q_c, kappa0); %father class init
        end

        function Phi = transMat(obj, t, s) 
            %TRANSMAT state transition matrix (s->t)
            Phi = transMatMinAcc(t, s, obj.dim/2);
        end

        function Q = QMat(obj, a, b) 
            %QMAT Q matrix, Q = intergrate Phi(b,s)*F(s)*Q_c*F(s)'*Phi(b,s)'*ds from a to b
            Q = MatQ_minAcc(a, b, obj.Q_c);
        end

        function V = VVec(obj, a, b)
            %VVec, V vector, V = intergrate Phi(b,s)*v(s)*ds from a to b, self-defined
            V = zeros(obj.dim, 1);
        end

    end
end

