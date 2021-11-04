function [x] = kalman_filter(y,A,C,Qw,Qv)
%   Input y  : noisy series
%   Input A  : State dynamics
%   Input C  : Observation dynamics
%   Input Qw : 
%   Input Qv : 

%   Output x : series after kalman filtering



order = size(A,1);
T = length(y);
x = zeros(T+order-1,1); %Make empty x vector

%initialisation => first iteration outside of loop
x_predict(1:1+order-1) = A*y(1:1+order-1);
p_predict = A*y(1:order)*y(1:order)'*A + Qw*eye(order);

gain = (p_predict*C')/(C*p_predict*C'+Qv);
p = (eye(order)-gain*C)*p_predict;
x(1:1+order-1) = x_predict(1:1+order-1)' + gain*(y(1) - C*x_predict(1:1+order-1)');
p_previous = p;
for i = 2:T
 
x_predict(i:i+order-1) = A*x((i-1):(i-1)+order-1); %predicted x is A*previous
p_predict = A*p_previous*A' + Qw*eye(order);

gain = (p_predict*C')/(C*p_predict*C'+Qv);
p = (eye(order)-gain*C)*p_predict;
x(i:i+order-1) = x_predict(i:i+order-1)' + gain*(y(i) - C*x_predict(i:i+order-1)');

p_previous = p;
end

x = x(order:T+order-1); %do not return initilisation.

end

