function val = omega(Q,s,diamond)
    if diamond == 2
        val = omega2(Q,s);
    elseif diamond == 3
        val = s_star(Q);
    else
        val = omega1inf(Q,s,diamond);
    end
end

function val = omega2(Q,s)
    val = 100000000;
    for i = 1:size(Q,2)
        val = min(val,omega2_i(Q,i,s));
    end
end

function val = omega2_i(Q,i,s)
    Qminusi = Q;
    Qminusi(:,i) = [];
    [x,r,g,info] = spg_lasso(Qminusi,Q(:,i),s-1,'verbosity',0);
    % Remove the verbosity=0 part if you want to see the iterations.
    val = norm(r,2);
end

function val = omega2_i_alt(Q,i,s)
    n = size(Q,2);
    Qminusi = Q;
    Qminusi(:,i) = [];
    cvx_begin quiet
        variable lambda(n-1);
        minimize norm(Q(:,i) - Qminusi*lambda, 2);
        subject to
            norm(lambda,1) <= s-1;
    cvx_end
    % Remove the verbosity=0 part if you want to see the iterations.
    val = norm(Q(:,i)- Qminusi*lambda, 2);
end

function val = omega1inf(Q,s,lval)
    val = 100000000;
    for i = 1:size(Q,2)
        val = min(val,omega1inf_i(Q,i,s,lval));
    end
end

function val = omega1inf_i(Q,i,s,lval)
    n = size(Q,2);
    Qminusi = Q;
    Qminusi(:,i) = [];
    cvx_begin quiet
        variable lambda(n-1);
        minimize norm(Q(:,i) - Qminusi*lambda, lval);
        subject to
            norm(lambda,1) <= s-1;
    cvx_end
    val = norm(Q(:,i)- Qminusi*lambda, lval);
end

function s = s_star(Q)
    n = size(Q,2);
    s = 0;
    for i = 1:n
        cvx_begin quiet
            variable z(n);
            maximize z(i);
            subject to
                Q*z == 0;
                norm(z,1) <= 1;
        cvx_end
        s = max(s,z(i));
    end
    s = 1/s;
end
