function [Ea,Eb,Ec,Ed,Ee,Ef,Eg,semimajor_axis, semiminor_axis,phi, X, Y] = ellipse_fit(x,y,steps,x0,y0)
%
% ellipse_fit - Given a set of points (x,y), ellipse_fit returns the
% best-fit ellipse (in the Least Squares sense)
%
% Input:
%                       x - a vector of x measurements
%                       y - a vector of y measurements
%
% Output:
%
%                   semimajor_axis - Magnitude of ellipse longer axis
%                   semiminor_axis - Magnitude of ellipse shorter axis
%                   x0 - x coordinate of ellipse center
%                   y0-  y coordinate of ellipse center
%                   phi - Angle of rotation in radians with respect to
%                   the x-axis
%--------------------------------------------------------------------------

x = x(:);
y = y(:);

%Construct M
M = [2*x.*y y.^2 2*x 2*y ones(size(x))];

% Multiply (-X.^2) by pseudoinverse(M)
e = M\(-x.^2);

%Extract parameters from vector e
a = 1;
b = e(1);
c = e(2);
d = e(3);
f = e(4);
g = e(5);
Ea = a;
Eb = b;
Ec = c;
Ed = d;
Ee = e;
Ef = f;
Eg = g;

%Use Formulas from Mathworld to find semimajor_axis, semiminor_axis, x0, y0
%, and phi

delta = b^2-a*c;

% x0 = (c*d - b*f)/delta;
% y0 = (a*f - b*d)/delta;

phi = 0.5 * acot((c-a)/(2*b));

nom = 2 * (a*f^2 + c*d^2 + g*b^2 - 2*b*d*f - a*c*g);
s = sqrt(1 + (4*b^2)/(a-c)^2);

a_prime = sqrt(nom/(delta* ( (c-a)*s -(c+a))));

b_prime = sqrt(nom/(delta* ( (a-c)*s -(c+a))));


semimajor_axis = max(a_prime, b_prime);
semiminor_axis = min(a_prime, b_prime);


if (a > c)
    phi = phi + pi/2;
end

phi = phi + 0.5*pi;

sinbeta = sin(phi);
cosbeta = cos(phi);

alpha = linspace(0, 360, steps)' .* (pi / 180);
sinalpha = sin(alpha);
cosalpha = cos(alpha);

X = x0 + (semimajor_axis * cosalpha * cosbeta - semiminor_axis * sinalpha * sinbeta);
Y = y0 + (semimajor_axis * cosalpha * sinbeta + semiminor_axis * sinalpha * cosbeta);

if nargout==1, X = [X Y]; end