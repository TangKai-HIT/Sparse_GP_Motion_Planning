classdef gpmpBase < handle
    %GPMP gaussian process motion planner base class
    %   此处显示详细说明
    
    properties
        gpTrajSparse; %GPTrajSparseBase obj handle
        options; % gpmp optim options struct
        constraints; % gpmp constraints struct
        
        %coeffients
        omega; %obs cost coeffs
        eta; %regularization cost coef
        lambda; %gp smooth cost coef
    end
    
    properties(SetAccess=protected, GetAccess=public)
        varDim; %dimension of variable (position)
        stateDim; %dimension of state (dimension in GP sparse)
        order; %order of state (pos+vel--2, pos+vel+acc--3, pos+vel+acc+jerk--4)
        
        numIntervals; %number of support states' intervals
        resolutions; %1 X numIntervals, define resolutions on each support state intervals
        numSamples; %1 X numIntervals, number of samples (determined by resolutions) on each support state intervals, (not including lower & upper bound of the interval)
        numSupports; %number of support states
        numTotal; %total number of points
        
        M_samp = 1; %sampling matrix
        mu_samp; %sampled state mu vector, ((numSamples+numSupports)*stateDim X 1)
        tau_samp; %sampled state index

        costHistory = struct("obsCost",[], "gpCost", [], "totalCost", []);
        stateHistory;
    end

    methods
        function obj = gpmpBase(gpTrajSparse, varDim, order, eta, omega, lambda, options, constraints)
            %GPMP 构造此类的实例
            %   此处显示详细说明
            gpTrajSparse.updateLifted();
            obj.gpTrajSparse = gpTrajSparse;

            obj.stateDim = gpTrajSparse.dim;
            obj.varDim = varDim;
            obj.order = order;
            
            obj.numSupports = gpTrajSparse.numsupports;
            obj.numIntervals = obj.numSupports - 1;
            
            obj.eta = eta;
            obj.lambda = lambda;
            obj.omega = omega;

            obj.options = options;
            obj.constraints = constraints;

            obj.upSampleByNums(zeros(1, obj.numIntervals)); %upsample
        end
        
        %% Upsample functions:
        function upSampleByNums(obj, numSamples)
            %UPSAMPLEBYNUMS upSample By defining sample Numbers at each intervals (not including lower & upper bound of the interval)
            %   numSamples: 1 X numIntervals

            obj.numSamples = numSamples;
            intervalDiff = diff(obj.gpTrajSparse.supportPts(:, 1:obj.varDim), 1, 1);
            intervalLen = sqrt(sum(intervalDiff.*intervalDiff, 2))';

            obj.resolutions = intervalLen ./ (numSamples+1);
            
            obj.numTotal = sum(obj.numSamples) + obj.numSupports;

            %update M matrix & mu vector
            obj.M_samp = zeros(obj.numTotal*obj.stateDim, obj.numSupports*obj.stateDim);
            obj.mu_samp = zeros(obj.numTotal*obj.stateDim, 1);
            obj.tau_samp = zeros(1, obj.numTotal);
            
            id_l = 1;
            for i=1:obj.numIntervals
                id_u = id_l + obj.numSamples(i) + 1;
                rows_i_eye = (id_l-1)*obj.stateDim+1 : id_l*obj.stateDim; %row indexes of identity matrix in M and Mu
                rows_i_samp = id_l*obj.stateDim+1 : (id_u-1)*obj.stateDim; %row indexes of sampling matrix in M and Mu
                cols_i_eye =  (i-1)*obj.stateDim+1 : i*obj.stateDim; %column indexes of identity matrix in M and Mu
                cols_i_samp =  (i-1)*obj.stateDim+1 : (i+1)*obj.stateDim; %column indexes of sampling matrix in M and Mu

                [tau_samp_i, Lambda_samp_i, Psi_samp_i, mu_samp_i] = obj.gpTrajSparse.upSampleInterval_uniform(i, i+1, obj.numSamples(i));
                
                obj.tau_samp(id_l) = obj.gpTrajSparse.supportId(i);
                obj.tau_samp(id_l+1 : id_l + obj.numSamples(i)) = tau_samp_i;

                obj.M_samp(rows_i_eye, cols_i_eye) = eye(obj.stateDim);
                obj.M_samp(rows_i_samp, cols_i_samp) = [Lambda_samp_i, Psi_samp_i];

                obj.mu_samp(rows_i_eye) = obj.gpTrajSparse.supportMu(cols_i_eye);
                obj.mu_samp(rows_i_samp) = mu_samp_i;

                id_l = id_u;
            end

            rows_i_eye = (id_l-1)*obj.stateDim+1 : id_l*obj.stateDim; %row indexes of identity matrix in M and Mu
            cols_i_eye =  i*obj.stateDim+1 : (i+1)*obj.stateDim; %column indexes of identity matrix in M and Mu
            obj.M_samp(rows_i_eye, cols_i_eye) = eye(obj.stateDim);
            obj.tau_samp(id_l) = obj.gpTrajSparse.supportId(end);
            obj.mu_samp(rows_i_eye) = obj.gpTrajSparse.supportMu(cols_i_eye);

        end

        function upSampleByResol(obj, resolutions)
            %UPSAMPLEBYRESOL upSample By defining resolutions at each intervals
            %   resolutions: 1 X numIntervals

            intervalDiff = diff(obj.gpTrajSparse.supportPts(:, 1:obj.varDim), 1, 1);
            intervalLen = sqrt(sum(intervalDiff.*intervalDiff, 2))';

            obj.numSamples =  ceil(intervalLen ./ resolutions) - 1; 
            upSampleByNums(obj, obj.numSamples);
        end
        
        function sampleStates = getSampleStates(obj, supportStatesVec)
            %GETSAMPLESTATES return sampled states, (N * stateDim) X 1
            %   supportStatesVec: (numSupports * stateDim) X 1 vec

             % supportStates = reshape(obj.gpTrajSparse.supportPts', obj.numSupports*obj.stateDim, 1);
             sampleStates = obj.M_samp * (supportStatesVec - obj.gpTrajSparse.supportMu) + obj.mu_samp;
        end
        
        %% Cost functions:
        function cost = gpPriorCost(obj, states)
            %GPPRIORCOST sparse GP prior distribution cost (negative log-likelihood)  of upsampled states
            meanDiff = states - obj.gpTrajSparse.supportMu;
            cost = 0.5* meanDiff' * obj.gpTrajSparse.supportKappa_inv * meanDiff;
        end

        function cost = totalCost(obj, states)
            %TOTALCOST total cost functions  of upsampled states
            cost = obj.lambda * gpPriorCost(obj, states) + obsCost(obj, states);
        end
        
        %% Optimization:
        function results = solve(obj)
            %SOLVE solve unconstrained optimization problem with all support states (default, overload this function in derived class if any constraint involved)
            %   exit_flag: 1-converged; 0-reached max iterations; -1-failed

            results = struct("x0", obj.gpTrajSparse.supportPts, "x_solve", [], "cost", [], "costHistory", [], "stateHistory", [], "exitFlag", -1, "steps", []);
            
            % %get active states indexes (remained consideration)
            % activeStateId = setdiff(1:obj.numIntervals+1, obj.options.FixedStateId);
            % activeIndex = kron((activeStateId-1).*obj.stateDim, ones(1,obj.stateDim)) +...
            %                                                 kron(ones(1,length(activeStateId)), 1:obj.stateDim); %1:dim + (n-1)d

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

                % samp_old = obj.getSampleStates(x_old);
                step = 1/obj.eta * (obj.lambda * (x_old - obj.gpTrajSparse.supportMu) + obj.omega * obj.gpTrajSparse.supportKappa * obj.M_samp' * obj.obsGrad(x_old));
                x_next = x_old - step;

                obs_cost_next = obj.obsCost(x_next);
                prior_cost_next = obj.gpPriorCost(x_next);
                total_cost_next = obj.lambda * prior_cost_next + obs_cost_next;
                
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
        
        %% Utils:
        function statesRow = getStateRowArrange(obj, statesVec)
            %GETSTATEROWARRANGE turn states vector (N*stateDim) X 1 into numSupports X stateDim row array matrix
            N = length(statesVec)/obj.stateDim;
            statesRow = reshape(statesVec, obj.stateDim, N)';
        end
    
        function statesVec = getStateVecArrange(obj, statesRow)
            %GETSTATEVECARRANGE turn N X stateDim row array matrix into states vector (N*stateDim) X 1
            statesVec = reshape(statesRow', size(statesRow, 1)*obj.stateDim, 1);
        end
        
        function plotCostHistory(obj)
            %PLOTCOSTHISTORY
            if obj.options.RecCostHis
                figure("Position",[300,300,800,600]);
    
                subplot(3,1,1);
                plot(obj.costHistory.obsCost, '-b', 'LineWidth', 0.8);
                title("obstacle cost");
                xlabel('Iterations');
                ylabel('Cost');
                
                subplot(3,1,2);
                plot(obj.costHistory.gpCost, '-b', 'LineWidth', 0.8);
                title("gp prior cost");
                xlabel('Iterations');
                ylabel('Cost');
    
                subplot(3,1,3);
                plot(obj.costHistory.totalCost, '-b', 'LineWidth', 0.8);
                title("total cost");
                xlabel('Iterations');
                ylabel('Cost');
    
                sgtitle('Cost History');
            end
        end

    end

    %% Abstract interfaces:
    methods(Abstract)
        obsCost(obj, states); %obstacle cost of upsampled states
        obsGrad(obj, states); %obstacle gradient of upsampled states
    end
end

