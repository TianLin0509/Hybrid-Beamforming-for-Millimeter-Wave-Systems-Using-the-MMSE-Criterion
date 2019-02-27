function [cost, grad] = getCostGrad(problem, x)
cost = problem.cost(x);
egrad = problem.egrad(x);
grad = problem.M.egrad2rgrad(x, egrad);
end

