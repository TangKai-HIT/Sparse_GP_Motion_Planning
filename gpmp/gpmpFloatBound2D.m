classdef gpmpFloatBound2D < gpmpBase
    %GPMPFLOATBOUND2D gpmp with floating bound(defined by observed variance) in 2D Euclidean space, with 2D signed distance field (SDF)
    %   using conditioned distribution between 2 near states

    properties
        disFieldfn;
        disGradfn;
        % jacobFcn;
        % bodyPts;
    end

    properties(SetAccess=private, GetAccess=public)
        kappa_bound; %kappa of variable states, conditioned by fixed states
        kappa_bound_inv;
        mu; %mu of variable states
    end

    methods
        function obj = gpmpFloatBound2D(gpSparseSets, varDim, order, eta, omega, lambda, kappa_goal, options, constraints)
            %GPMPBOUND2D suppose bounded by 2 fixed ends
            %   gpSparseSets: init struct generate by gpSparse_init_set() function, min-jerk/min-snap traj is required
            %   order: 3~4
            
            if order==3
                gpTrajSparse = GPTrajSparseMinJerk(gpSparseSets.supportPts, gpSparseSets.supportId, gpSparseSets.Q_c, gpSparseSets.mu0, gpSparseSets.kappa0);
            elseif order==4
                gpTrajSparse = GPTrajSparseMinSnap(gpSparseSets.supportPts, gpSparseSets.supportId, gpSparseSets.Q_c, gpSparseSets.mu0, gpSparseSets.kappa0);
            end

            obj@gpmpBase(gpTrajSparse, varDim, order, eta, omega, lambda, options, constraints);

            %update conditioned kappa on middle states
            midState_Id = obj.gpTrajSparse.supportId(2:end-1);
            [~, ~, ~, ~, obj.mu, obj.kappa_bound, obj.kappa_bound_inv] = ...
                                            obj.gpTrajSparse.upSampleInterval_observeGoal(1, obj.numSupports, midState_Id, gpSparseSets.kappa0, kappa_goal);
        end
        
        %% Prior cost
        function cost = gpPriorCost(obj, states)
            %GPPRIORCOST sparse GP prior distribution cost (negative log-likelihood)  of upsampled states
            meanDiff = states - obj.mu;
            cost = 0.5* meanDiff' * obj.kappa_bound_inv * meanDiff;
        end

        %% Obstacle cost & grad
        function cost = obsCost(obj, statesVec)
            %OBSCOST obstacle cost
            %   此处显示详细说明
            
            sampStatesVec = obj.getSampleStates(statesVec);
            states = obj.getStateRowArrange(sampStatesVec);
            
            sdfCost = obj.disFieldfn(states(:, 1:2));
            velNorm = sqrt(sum(states(:, 3:4) .* states(:, 3:4), 2));

            cost = trapz(obj.tau_samp, sdfCost .* velNorm);
        end

        function grad = obsGrad(obj, statesVec) 
            %OBSGRAD gradient of obstacle cost 
            
            sampStatesVec = obj.getSampleStates(statesVec);
            states = obj.getStateRowArrange(sampStatesVec);
            grad = zeros(size(states));

            c = obj.disFieldfn(states(:, 1:2)); % n X 1
            delta_c = obj.disGradfn(states(:, 1:2)); % n X Dim
            
            vel = states(:, 3:4);
            acc = states(:, 5:6);
            velNorm = sqrt(sum(vel .* vel, 2)); % n X 1
            
            unitVel = vel ./ velNorm;

            %position gradient
            directGrad = velNorm .* (delta_c - c .* acc ./ (velNorm.^2)); %direct part
            grad(:, 1:2) = directGrad - unitVel .* sum(unitVel .* directGrad, 2); %after projection 

            %velocity gradient
            grad(:, 3:4) = c .* unitVel;

            grad(isnan(grad) | isinf(grad)) = 0;

            grad = obj.getStateVecArrange(grad);
        end
        
        %% Optimization:
        function results = solve(obj)
            %SOLVE solve optimization problem with fixed support states at 2 ends, use interpolation
            %   exit_flag: 1-converged; 0-reached max iterations; -1-failed

            results = struct("x0", obj.gpTrajSparse.supportPts, "x_solve", [], "cost", [], "costHistory", [], "stateHistory", [], "exitFlag", -1, "steps", []);

            %first step
            x_next = reshape(obj.gpTrajSparse.supportPts', obj.numSupports*obj.stateDim, 1);
            obs_cost_next = obj.obsCost(x_next);
            prior_cost_next = obj.gpPriorCost(x_next);
            total_cost_next = obj.lambda * prior_cost_next + obs_cost_next;
            
            if obj.options.RecCostHis
                obj.costHistory.obsCost = obs_cost_next;
                obj.costHistory.gpCost = prior_cost_next;
                obj.costHistory.totalCost = total_cost_next;
            end

            if obj.options.RecStateHis
                obj.stateHistory = x_next;
            end
            
            %constrained problem settings
            if ~isempty(obj.constraints)
                %init qp params
                H = obj.eta*obj.kappa_bound_inv; %Hessian in QP
                H=(H+H')/2;
                
                %init warm-start obj
                if obj.options.UseOSQP %use osqp solver
                    %to do...
                else
                    qp_options = optimoptions('quadprog','Algorithm','active-set');
                    qp_ws = optimwarmstart(x_next, qp_options); %qp warm-start
                end

                %init constraints
                M_bar = obj.M_samp;

                A = [M_bar; -M_bar];
                numActiveSamp = obj.numTotal*obj.stateDim;
                b = [inf*ones(numActiveSamp, 1); -inf*ones(numActiveSamp, 1)];                
                
                if ~isempty(obj.constraints.posRange) %pos constraints
                    posId = 1 : obj.stateDim : numActiveSamp;
                    for i=1:obj.varDim
                        posId = posId + (i-1);
                        %upper bound
                        b(posId) = obj.constraints.posRange(i, 2);
                        %lower bound
                        b(posId + numActiveSamp) = -obj.constraints.posRange(i, 1);
                    end
                end

                if ~isempty(obj.constraints.velRange) %vel constraints
                    velId = obj.varDim+1 : obj.stateDim : numActiveSamp;
                    for i=1:obj.varDim
                        velId = velId + (i-1);
                        %upper bound
                        b(velId) = obj.constraints.velRange(i, 2);
                        %lower bound
                        b(velId + numActiveSamp) = -obj.constraints.velRange(i, 1);
                    end
                end

                if ~isempty(obj.constraints.accRange) %acc constraints
                    accId = 2*obj.varDim+1 : obj.stateDim : numActiveSamp;
                    for i=1:obj.varDim
                        accId = accId + (i-1);
                        %upper bound
                        b(accId) = obj.constraints.accRange(i, 2);
                        %lower bound
                        b(accId + numActiveSamp) = -obj.constraints.accRange(i, 1);
                    end
                end

                if ~isempty(obj.constraints.jerkRange) %jerk constraints
                    jerkId = 3*obj.varDim+1 : obj.stateDim : numActiveSamp;
                    for i=1:obj.varDim
                        jerkId = jerkId + (i-1);
                        %upper bound
                        b(jerkId) = obj.constraints.jerkRange(i, 2);
                        %lower bound
                        b(jerkId + numActiveSamp) = -obj.constraints.jerkRange(i, 1);
                    end
                end
                
                rmId = find(b==inf | b==-inf);
                b(rmId) = []; %remove unbounded constraints
                A(rmId, :) = [];
            end

            %iteration
            results.exitFlag = 0;
            for i = 0:obj.options.MaxIter
                %end condition
                if i>1 && abs((total_cost_old - total_cost_next)/total_cost_old) < obj.options.TolFun
                    results.exitFlag = 1;
                    break;
                end
                
                x_old = x_next;
                total_cost_old = total_cost_next;
                
                %one-step decent
                if isempty(obj.constraints) %use covariant gradient decent
                    obs_grad = obj.M_samp' * obj.obsGrad(x_old);
                    step = 1/obj.eta * (obj.lambda * (x_old - obj.mu) + obj.omega * obj.kappa_bound * obs_grad);
                    x_next = x_old - step;
                else %use SQP with trust region/regulation
                    obs_grad = obj.M_samp' * obj.obsGrad(x_old);
                    f = obj.omega*obs_grad + obj.kappa_bound_inv*(obj.lambda * (x_old - obj.mu) - obj.eta*x_old);
                    [qp_ws, ~, qp_exitflag, ~, ~] = quadprog(H,f,A,b,[],[],[],[],qp_ws);
                    x_next = qp_ws.X;
                end

                obs_cost_next = obj.obsCost(x_next);
                prior_cost_next = obj.gpPriorCost(x_next);
                total_cost_next = obj.lambda * prior_cost_next + obj.omega * obs_cost_next;
                
                %record history
                if obj.options.RecCostHis
                    obj.costHistory.obsCost(end+1) = obs_cost_next;
                    obj.costHistory.gpCost(end+1) = prior_cost_next;
                    obj.costHistory.totalCost(end+1) = total_cost_next;
                end
    
                if obj.options.RecStateHis
                    obj.stateHistory(:, end+1) = x_next;
                end
            
            end
            
            %update results
            if obj.options.RecCostHis
                results.costHistory = obj.costHistory;
            end
            if obj.options.RecStateHis
                results.stateHistory = obj.stateHistory;
            end

            results.x_solve = obj.getStateRowArrange(x_next);
            results.cost = total_cost_next;
            results.steps = i;

            % obj.gpTrajSparse.supportPts = results.x_solve;
        end       

        %% Utils
        function plotDeformHis(obj, ax, cmapHandle, numPlots)
            %PLOTDEFORMHIS 
            numHistory = size(obj.stateHistory, 2);
            
            Id = ceil(linspace(1, numHistory, numPlots));
            
            hold(ax, "on");
            cc=cmapHandle(numPlots);

            N = 100;
            plotId = linspace(obj.gpTrajSparse.supportId(1), obj.gpTrajSparse.supportId(end), N);
            plotStates = zeros(N, obj.gpTrajSparse.dim);

            for i = 1:numPlots
                rowStates = obj.getStateRowArrange(obj.stateHistory(:, Id(i)));
                sampStates = obj.getStateRowArrange(obj.getSampleStates(obj.stateHistory(:, Id(i))));
                
                %plot traj
                obj.gpTrajSparse.supportPts = rowStates;
                for k = 1:N
                    plotStates(k, :) = obj.gpTrajSparse.query(plotId(k));
                end
                plot(ax, plotStates(:, 1), plotStates(:,2),'color',cc(i,:), 'LineStyle', '-', 'LineWidth', 1);
                
                %plot sampled states
                if Id(i)==1 || Id(i)==numHistory
                    scatter(ax, sampStates(:, 1), sampStates(:,2), 20, 'filled', "MarkerFaceColor", cc(i,:));
                end
                
                %plot support states
                scatter(ax, rowStates(:, 1), rowStates(:,2), 50, 'filled', "MarkerFaceColor", cc(i,:));
            end
            
            %recover initial states
            obj.gpTrajSparse.supportPts = obj.getStateRowArrange(obj.stateHistory(:, 1));
        end

    end
end