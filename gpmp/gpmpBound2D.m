classdef gpmpBound2D < gpmpBase
    %GPMPBOUND2D gpmp with fixed bound in 2D Euclidean space, with 2D cost map (SDF)
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
        mu_bound; %mu of variable states, conditioned by fixed states
        activeId; %active states index between fixed bounds
    end

    methods
        function obj = gpmpBound2D(gpSparseSets, varDim, eta, omega, lambda, options)
            %GPMPBOUND2D suppose bounded by 2 fixed ends
            %   gpSparseSets: init struct generate by gpSparse_init_set() function, min-jerk traj is required

            gpTrajSparse = GPTrajSparseMinJerk(gpSparseSets.supportPts, gpSparseSets.supportId, gpSparseSets.Q_c, gpSparseSets.mu0, gpSparseSets.kappa0);
            
            order = 3;

            obj@gpmpBase(gpTrajSparse, varDim, order, eta, omega, lambda, options, []);

            %update conditioned kappa on middle states
            midState_Id = obj.gpTrajSparse.supportId(2:end-1);
            [~, ~, ~, obj.mu_bound, obj.kappa_bound] = obj.gpTrajSparse.upSampleInterval(1, obj.numSupports, midState_Id);
            
            obj.kappa_bound_inv = inv(obj.kappa_bound); %(to do: find simpler formula to compute kappa_bound_inv)

            obj.activeId = obj.stateDim+1:(obj.numSupports-1)*obj.stateDim;
        end
        
        %% Prior cost
        function cost = gpPriorCost(obj, states)
            %GPPRIORCOST sparse GP prior distribution cost (negative log-likelihood)  of upsampled states
            meanDiff = states(obj.activeId) - obj.mu_bound;
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
            %SOLVE solve unconstrained optimization problem with fixed support states at 2 ends, take interpolating points as optimization variables
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

                obs_grad = obj.M_samp' * obj.obsGrad(x_old);
                step = 1/obj.eta * (obj.lambda * (x_old(obj.activeId) - obj.mu_bound) + obj.omega * obj.kappa_bound * obs_grad(obj.activeId));
                x_next(obj.activeId) = x_old(obj.activeId) - step;

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