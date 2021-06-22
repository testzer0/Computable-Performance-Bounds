rng('default');
vals = 1:40;
mbyns = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
colors2 = ['--b','--g','--r','--c','--m','--y','--k'];
colors1 = ['b','g','r','c','m','y','k'];
mbyn = mbyns(4);
[mcbounds,omegabounds] = getvals_graph_omega_mc(mbyn);
save('outputs_2_4.mat');

function [mcbounds,omegabounds] = getvals_graph_omega_mc(mbyn)
    % Evaluates the MC and Omega bounds at [1,2,...,40] and returns them
    % for a random Hadamard sensing matrix with n = 2048 and m = mbyn*n.
    n = 2048;
    m = round(n*mbyn);
    H = hadamard(n);
    Q = H(randperm(size(H,1)),:);
    Q = Q(1:m,:);
    mcbounds = zeros(40,1);
    omegabounds = zeros(40,1);
    mc = mutual_coherence(Q);
    for i = 1:1:40
        mcbounds(i) = 2/(1 - mc*(4*i-1));
        if (mcbounds(i) < 0)
            mcbounds(i) = Inf;
        end
        omegabounds(i) = bp_omega_bound(Q,i);
        fprintf("Omega bounds for %d : %f\n", i, omegabounds(i));
        fprintf("MC bounds for %d : %f\n", i, mcbounds(i));
    end
end

function disp_bp_bounds(mattype)
% mattype 1 -> Bernoulli, 2 -> Guassian, 3 -> Hadamard
    n = 256;
    for m = [51,77,102,128,154,179,205]
        fprintf("For m = %d\n\n", m);
        if (mattype == 1)
            % Bernoulli
            Q = rand(m,n);
            Q = (Q < 0.5);
            Q = 2*Q - 1;
        elseif (mattype == 2)
            % Guassian
            Q = randn(m,n);
        else
            % Hadamard
            H = hadamard(n);
            Q = H(randperm(size(H,1)),:);
            Q = Q(1:m,:);
        end
        
        for i = 1:size(Q,2)
            Q(:,i) = Q(:,i)/norm(Q(:,i));
        end
        s_star = omega(Q,0,3);
        k_star = floor(s_star/2);
        fprintf("s_star = %f\t\tk_star = %f\n", s_star, k_star);
        for k = 1:k_star
            bd1 = bp_omega_bound(Q,k);
            bd2 = bp_ric_bound(Q,k);
            fprintf("\t\tk = %d\tOmega\t=>\t%f\n",k,bd1);
            fprintf("\t\t      \tRIC\t=>\t%f\n",k,bd2);
        end
        fprintf("-------------------------------\n");
    end
end

function disp_dantzig_bounds(mattype)
% mattype 1 -> Bernoulli, 2 -> Guassian, 3 -> Hadamard
    n = 256;
    for m = [51,77,102,128,154,179,205]
        fprintf("For m = %d\n\n", m);
        if (mattype == 1)
            % Bernoulli
            Q = rand(m,n);
            Q = (Q < 0.5);
            Q = 2*Q - 1;
        elseif (mattype == 2)
            % Guassian
            Q = randn(m,n);
        else
            % Hadamard
            H = hadamard(n);
            Q = H(randperm(size(H,1)),:);
            Q = Q(1:m,:);
        end
        
        for i = 1:size(Q,2)
            Q(:,i) = Q(:,i)/norm(Q(:,i));
        end
        s_star = omega(Q,0,3);
        k_star = floor(s_star/2);
        fprintf("s_star = %f\t\tk_star = %f\n", s_star, k_star);
        for k = 1:k_star
            bd1 = dantzig_omega_bound(Q,k);
            bd2 = dantzig_ric_bound(Q,k);
            fprintf("\t\tk = %d\tOmega\t=>\t%f\n",k,bd1);
            fprintf("\t\t      \tRIC\t=>\t%f\n",k,bd2);
        end
        fprintf("-------------------------------\n");
    end
end

function val = bp_omega_bound(A,k)
    omval = omega(A,2*k,2);
    val = 2*sqrt(2*k)/omval;
end

function val = dantzig_omega_bound(A,k)
    omval = omega(A'*A,2*k,inf);
    val = 2*sqrt(2*k)/omval;
end

function val = bp_ric_bound(A,k)
    r = ric(A,2*k,5000);
    if (r < (sqrt(2)-1))
        val = (4*sqrt(1+r))/(1-(1+sqrt(2))*r);
    else
        val = Inf;
    end
end

function val = dantzig_ric_bound(A,k)
    r2k = ric(A,2*k,5000);
    r3k = ric(A,3*k);
    if (r2k+r3k < 1)
        val = (4*sqrt(k))/(1-r2k-r3k);
    else
        val = Inf;
    end
end

function val = bp_cmsv_bound(A,k)
    % Memory requirements are ~35 GB !!
    r = cmsv(A,4*k);
    val = 2/r;
end

function val = dantzig_cmsv_bound(A,k)
    % Memory requirements are ~35 GB !!
    r = cmsv(A,4*k);
    val = 4*sqrt(k)/(r^2);
end

function val = mc_bound(A,k)
    % This function is not called directly since it's better to precompute
    % mu for various values of k.
    mu = mutual_coherence(A);
    val = 2/(1-mu*(4*k-1));
end

function val = mutual_coherence(A)
    % Only one function so best implemented here itself.
    % First we unit normalize A's columns.
    n = size(A,2);
    for i = 1:n
        A(:,i) = A(:,i)/norm(A(:,i));
    end
    % Now the answer is the max entry in the absolute value of A'*A
    M = A'*A;
    M = M.*~eye(size(M));
    val = max(max(abs(M)));
end