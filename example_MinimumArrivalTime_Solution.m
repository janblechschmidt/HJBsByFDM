function solfunc = solution(mu)
    % TODO not working in instationary case

        % Set parameters
        c1 = -exp(mu-0.5)
        vminus = @(x) c1 * (exp(x) - 1) + x *(1 + mu);
        vplus = @(x) vminus(1-x);
        % vplus = @(x) (1+mu) * (1-x) + c3 * (exp(-x) - exp(-1));
        xminus = -log(-c1);
        a = vminus(xminus)+0.5*(xminus-0.5)^2;
        b = 0;
        c = -0.5;
        vzero = @(x) a + b*(x-0.5) + c*(x-0.5).^2
        % x = linspace(0,1,501);
        % iminus = (x < xminus);
        % izero = (x>=xminus & x<= 1-xminus);
        % iplus = (x > 1-xminus);

        solfunc = @(t,x) (x < xminus) .* vminus(x) + ...
            (x>=xminus & x<= 1-xminus) .* vzero(x) + ...
            (x > 1-xminus) .* vplus(x);
end
