function [D0, X0] = LRSDL_initD_barZ_bar(Y, opts, seed)
    % cost = \arg\min .5*\|Y - D0X0\|_F^2 + \lambda1\|X0\|_1 + .5*lambda2*\|X0 - M0\|_F^2 + eta*\|D0\|_F^2
    lambda1 = opts.lambda1;
    lambda2 = opts.lambda2;
    lambda3 = 0.3;
    k0 = opts.k0;
    rng(seed, "twister");
    D0 = rand(size(Y,1),k0); %KKD %12/13/18
    X0 = PickDfromY(Y', [0, size(Y,1)],k0);%%KKD 12/13/18
    X0=X0';
    %% cost
    function cost = calc_cost(D0, X0)
%         cost = .5*normF2(Y - D0*X0) + lambda3*nuclearnorm(D0) + ...
%             5*lambda2*normF2(X0 - buildMean(X0)) + lambda1*norm1(X0);
     cost = .5*normF2(Y - D0*X0) + lambda3*nuclearnorm(D0) + ...
         lambda1*norm1(X0);
    end
    %%
    iter = 0;
    cost_old = calc_cost(D0, X0);
    opts.verboseX=0;
    opts.verboseD=0;
    optsX = opts;
    optsX.max_iter =10000;
    optsD = opts;
    optsD.max_iter = 10000;
    while iter < opts.max_iter 
        iter = iter + 1;
        %% ========= update X0 ==============================
        %X0 = \argmin_X0 .5*\|Y - D0X0\|_F^2 + \lambda1\|X0\|_1 + .5*lambda2*\|X0 - M0\|_F^2
        X0 = myLassoWIntrasmall_fista(Y, D0, lambda1, lambda2, X0, optsX);
%         imagesc(X0);
%         if opts.verbose
%             disp(iter);
%             costX0 = calc_cost(D0, X0);
%             disp(costX0);
%         end 

        %% ========= update D ==============================
        % \arg\min .5*\|Y - D0X0\|_F^2 + eta*\|D0\|_F^2
%         D0 = minRankDict0(Y, X0, lambda3, D0, optsD);
        D0 = min_rank_dict(D0, Y*X0', X0*X0', 2*lambda3, optsD);
        costD0(iter) = calc_cost(D0, X0);%Rui:check convergence
        if opts.verbose
            disp(costD0(iter));
        end
        if(rem(iter,5)==0)
            if(abs(costD0(iter) - cost_old) < 1e-4)
               break;
            end
        end
        cost_old = costD0(iter);
    end
end 