function r = ric(A,kmult,numsamples)
    if nargin > 2
        ns = numsamples;
    else
        ns = 1000;
    end
    r = 0;
    for i = 1:ns
        r = max(r,get_sub_ric(A,kmult));
    end
end

function r = get_sub_ric(A,kmult)
    Asub = get_random_submatrix(A,size(A,1),kmult);
    smax = svds(Asub,1);
    smin = svds(Asub,1,'smallest');
    r = max(smax^2 - 1, 1 - smin^2);
end

function Asub = get_random_submatrix(A,a,b)
    Atemp = A(randperm(size(A,1)),:);
    Atemp = Atemp(1:a,:);
    Atemp = Atemp(:,randperm(size(Atemp,2)));
    Asub = Atemp(:,1:b);
end