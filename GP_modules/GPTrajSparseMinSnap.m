classdef GPTrajSparseMinSnap < GPTrajSparseBase
    %GPTRAJSPARSEMINJERK sparse GP trajectory with zeros snap prior
    %   sparse GP trajectory with O(1) query time complexity, derived from Guass-Markov SDE with sparse kernel
    
    properties(SetAccess=private, GetAccess=public)
        % end state variables (for whole trajectory conditioning)
        % (to do ... , note that conditioning on the observation at the end point won't affect the query of mid local points)
    end
    
    methods
        function obj = GPTrajSparseMinSnap(supportPts, supportId, Q_c, mu0, kappa0)
            %GPTRAJSPARSEMINJERK 构造此类的实例
            %   supportId must include initial condition
            if nargin == 0
                supportPts=1; 
                supportId=1; Q_c=1;
                mu0 = 0; kappa0 = 1;
            end
            obj@GPTrajSparseBase(supportPts, supportId, Q_c, mu0, kappa0); %father class init
        end
        
        function copyObj = copy(obj)
            %COPY return the handle of a newly copied obj
            copyObj = GPTrajSparseMinSnap(obj.supportPts, obj.supportId, obj.Q_c, obj.mu0, obj.kappa0);
        end

        function Phi = transMat(obj, t, s) 
            %TRANSMAT state transition matrix (s->t)
            stateDim = obj.dim/4;
            I = eye(stateDim);
            zero = zeros(stateDim,stateDim);

            delta = t-s; 
            Phi = [I,  delta.*I,  0.5*delta^2.*I,   1/6*delta^3.*I;
                    zero,  I,     delta.*I,   0.5*delta^2.*I;
                    zero, zero,    I,         delta.*I;
                    zero, zero,   zero,          I];
        end

        function Q = QMat(obj, a, b) 
            %QMAT Q matrix, Q = intergrate Phi(b,s)*F(s)*Q_c*F(s)'*Phi(b,s)'*ds from a to b
            Delta = b-a;
            Q = [(Delta^7*obj.Q_c)/252, (Delta^6*obj.Q_c)/72, (Delta^5*obj.Q_c)/30, (Delta^4*obj.Q_c)/24; 
            (Delta^6*obj.Q_c)/72, (Delta^5*obj.Q_c)/20, (Delta^4*obj.Q_c)/8, (Delta^3*obj.Q_c)/6; 
            (Delta^5*obj.Q_c)/30, (Delta^4*obj.Q_c)/8, (Delta^3*obj.Q_c)/3, (Delta^2*obj.Q_c)/2; 
            (Delta^4*obj.Q_c)/24, (Delta^3*obj.Q_c)/6, (Delta^2*obj.Q_c)/2, Delta*obj.Q_c];
        end

        function V = VVec(obj, a, b)
            %VVec, V vector, V = intergrate Phi(b,s)*v(s)*ds from a to b, self-defined
            V = zeros(obj.dim, 1);
        end

    end
end

