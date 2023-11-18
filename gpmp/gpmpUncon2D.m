classdef gpmpUncon2D < gpmpBase
    %GPMPUNCON (deprecated) unconstrained gpmp in 2D Euclidean space, with 2D cost map (SDF)
    %   not fixed or bounded anywhere, would just flow away by gradient decent

    properties
        disFieldfn;
        disGradfn;
        % jacobFcn;
        % bodyPts;
    end

    methods
        function obj = gpmpUncon2D(gpTrajSparse, varDim, order, eta, omega, lambda, options)
            %GPMPUNCON 构造此类的实例
            %   此处显示详细说明

            obj@gpmpBase(gpTrajSparse, varDim, order, eta, omega, lambda, options, []);
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
        
        %% Utils
        function plotDeformHis(obj, ax, cmapHandle)
            %PLOTDEFORMHIS 
            numHistory = size(obj.stateHistory, 2);

            hold(ax, "on");
            cc=cmapHandle(numHistory);
            for i = 1:numHistory
                rowStates = obj.getStateRowArrange(obj.stateHistory(:, i));
                sampStates = obj.getStateRowArrange(obj.getSampleStates(obj.stateHistory(:, i)));

                plot(ax, sampStates(:, 1), sampStates(:,2),'color',cc(i,:), 'LineStyle', '--', 'LineWidth', 0.5);

                if i==1 || i==numHistory
                    scatter(ax, sampStates(:, 1), sampStates(:,2), 25, 'filled', "MarkerFaceColor", cc(i,:));
                end

                scatter(ax, rowStates(:, 1), rowStates(:,2), 50, 'filled', "MarkerFaceColor", cc(i,:));
            end
        end

    end
end