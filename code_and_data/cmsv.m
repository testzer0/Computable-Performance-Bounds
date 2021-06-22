function val = cmsv(A,s,mode)
    if (nargin == 2 || mode == 0)
        val = cmsv_montecarlo(A,s,5000);
    else
        n = size(A,2);
        cvx_begin quiet
            variable Z(n,n);
            minimize trace((A'*A)*Z);
            subject to
                Z == semidefinite(n);
                trace(Z) == 1;
                norm(Z,1) <= s;
        cvx_end
        val = trace((A'*A)*Z);
    end
end

function val = cmsv_montecarlo(A,s,numsamples)
    if nargin > 2
        ns = numsamples;
    else
        ns = 1000;
    end
    val = Inf;
    for i = 1:ns
        val = min(val,get_sub_cmsv(A,s));
    end
end

function c = get_sub_cmsv(A,s)
    Asub = get_random_submatrix(A,size(A,1),s);
    smin = svds(Asub,1,'smallest');
    c = max(0,1 - smin);
end

function Asub = get_random_submatrix(A,a,b)
    Atemp = A(randperm(size(A,1)),:);
    Atemp = Atemp(1:a,:);
    Atemp = Atemp(:,randperm(size(Atemp,2)));
    Asub = Atemp(:,1:b);
end