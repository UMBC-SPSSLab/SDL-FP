function [X, X0] = LRSDL_updateXX0_mod(Y, Y_range, D, D_range, D0, X, X0, W, opts)
    if nargin == 0
        error('nargin ==0 not supported');
        addpath('../utils');
        addpath('../sparse_coding');
        tic
        d       = 30;
        N       = 7;
        k       = 5;
        k0      = 5;
        C       = 3 ;
        
        Y       = normc(rand(d,N*C));
        D       = normc(rand(d,k*C));
        D0      = normc(rand(d, k0));
        Y_range = N* (0:C);
        D_range = k* (0:C);
        X       = randn(size(D,2), size(Y,2));
        X0      = randn(size(D0, 2), size(Y, 2));
        
        opts.k0       = k0;
        opts.lambda1  = 0.01;
        opts.lambda2  = 0.002;
        opts.lambda3  = 0.1;
        opts.max_iter = 250;
        opts.show     = true;    
        opts.verbose   = true;
        opts          = initOpts(opts); % other attributes
    end
    %% 
    function cost = calc_f(X1)
        [X, X0] = extractFromX1(X1);
        Ybar = Y - D0*X0;   %SK: added   
        cost = 0.5*normF2(Ybar - D*X); % + ...   % SK: modified
               %0.5*opts.lambda2* (FDDL_discriminative(W'*X, Y_range) + ...
               %                     normF2(W'*(X0 - buildMean(X0))));
    end 
    %% Total cost 
    function cost = calc_F(X1)
        cost = calc_f(X1) + lambda1*norm1(X1);
    end 
    %% [X, X0] = [X; X0]
    function [X, X0] = extractFromX1(X1)
        X = X1(1:D_range(end), :);
        X0 = X1(D_range(end) + 1: end, :);
    end    
    %% Gradient for FISTA 
    function g1 = grad(X1)
        [X, X0] = extractFromX1(X1);
        %DtY     = DtY0 - DtD0*X0;
        %Y_0     = buildMhat(DtY, D_range, Y_range);
        %g       = Dhat*X - Y_0 + buildM_2Mbar(X, Y_range, lambda2);
        %g0      = A*X0 - D0tY2 + D0tD*buildMhat(X, D_range, Y_range)...
        %             - lambda2*buildMean(X0);        
        %g1      = [g; g0];
        g1 = -DbartY + DbartDbar*X1 + ...
            [W*W'*(lambda2*2*X + buildM_2Mbar(X, Y_range, lambda2));
             lambda2*(X0 - buildMean(X0))]; % SK: refer to my note
    end
    %% 
    if opts.k0 == 0
        error('k0 == 0 not allowed');
        X = FDDL_updateX(Y, Y_range, D, D_range, X, opts);
        X0 = [];
    else 
        %% Prepare data 
        lambda1  = opts.lambda1;
        lambda2  = opts.lambda2;
        %DtD      = D'*D;
        %D_0      = buildMhat(DtD, D_range, D_range);  
        %Dhat     = D_0 + 2*opts.lambda2*eye(size(D_0,1));
        %D0tD0    = D0'*D0;
        %SK: A        = 2*D0tD0 + opts.lambda2*eye(size(D0,2));
        A        = opts.lambda2*eye(opts.k0);
%         DtY0     = D'*Y;
%         DtD0     = D'*D0;
%         D0tD     = DtD0';
%         D0tY2    = 2*D0'*Y;
        %SK: added
        Dbar = [D D0];
        DbartY = Dbar'*Y;
        DbartDbar = Dbar'*Dbar;
        %% check gradient
        if opts.check_grad &&~check_grad(@calc_f, @grad, [X; X0])
            fprintf('Check gradient or cost again!\n')
            pause
        end       
        %% ========= Main FISTA ==============================
        optsXX0          = opts;
        optsXX0.verbose = false;
        optsXX0.max_iter = 300;
        %SK: L = max(eig(Dhat)) + max(eig(A)) + 4*lambda2 + 1;
        L = max(eig(DbartDbar)) + max(eig(A)) + 4*lambda2 + 1;
        X1 = [X; X0];
        X1      = fista(@grad, X1, L, opts.lambda1, optsXX0, @calc_F);
        [X, X0] = extractFromX1(X1);
    end 
    %%
    if nargin == 0   
        fprintf('done, press any key to see results\n');
        pause;
    end
end 

