function [D, Z] = DL_smallIntraCoeff(Y, k, lambda1, lambda2, opts, seed)
    rng(seed);
    D = rand(size(Y,2),k); 
    Z = PickDfromY(Y, [0, size(Y,2)],k);

    Y=Y';Z=Z';
	function cost = calc(D, Z)
% 		cost = 0.5*normF2(Y - D*Z) + .5*lambda2*normF2(Z - buildMean(Z)) + lambda1*norm1(Z);		
        cost = 0.5*normF2(Y - D*Z) + .5*lambda2*normF2(D' - buildMean(D')) + lambda1*norm1(Z);	%KKD 12/12/18	
    end
    %%
	it = 0;
    cost_old = calc(D, Z);    
    if opts.verbose
        fprintf('%f', cost_old);
        fprintf('.');
    end 
    %%
    opts.verboseZ   = 0;
    opts.verboseD   = 0;
    optsZ = opts;
    optsZ.max_iter = 10000;
    optsD = opts;
    optsD.max_iter = 10000;
    Z_old=Z;
    D_old=D;
    sizeD = numel(D);
    sizeZ = numel(Z);
    while it < opts.max_iter
		it = it + 1;
		%% ========= update Z ==============================
        % Z = \arg\min_Z 0.5*normF2(Y - D*Z) + 
        %     .5*lambda2*normF2(Z - buildMean(Z)) + lambda1*norm1(Z);
		Z = myLassoWIntrasmall_fista(Y, D, lambda1, lambda2, Z, optsZ);
        %% ========= update D ==============================
		E = Y*Z';
		F = Z*Z';
		D = ODL_updateD(D, E, F, optsD);
		%% ========= For debugging ==============================
%       cost_new(it) = calc(D, Z);
        if 1
            cost_new(it) = calc(D, Z);
            if(rem(it,5)==0)
                if(abs(cost_new(it) - cost_old) < 1e-4)
                   break;
                end
            end
            cost_old = cost_new(it);
        end
        if 0
            if rem(it,3) == 0
                ErrZ(it) = norm(D - D_old, 'fro')/sizeD;
                ErrD(it) = norm1(Z - Z_old)/sizeZ;
                if (abs(ErrZ(it))<=opts.tol) && (abs(ErrD(it))<=opts.tol*1e-2)
                    break;
                end
            end
        end
        Z_old=Z;
        D_old=D;
        if opts.verbose
        fprintf('%f', cost_new(it));
        end
    end 
     if opts.verbose
        fprintf('%f', cost_new(it));
    end
    if nargin == 0
        Z = [];
        D = [];
    end
 end 

