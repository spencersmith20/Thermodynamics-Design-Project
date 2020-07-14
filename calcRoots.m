%calculates the roots of the cubic Peng-Robinson equation of state
%inputs: amix, bmix (total mixture coefficients), T, P (system pressure and
%temperature

function [vL, vV] = calcRoots(A, B, T, P)

    r = 8.3144598; %gas constant
    
    roots = zeros(1, 3); %coefficients
    b = B - r*T/P;
    c = A/P - 3*B^2 - 2*r*T*B/P;
    d = B^3 + r*T*B^2/P - A*B/P;
    
    Q = (b^2 - 3*c)/9;
    R = (2*b^3 - 9*b*c + 27*d)/54;
        
    if Q^3 - R^2 >= 0 %three root case
        th = acos(R/(Q^(3/2)));
        roots(1) = -2*sqrt(Q)*cos(th/3) - b/3;
        roots(2) = -2*sqrt(Q)*cos((th + 2*pi)/3) - b/3;
        roots(3) = -2*sqrt(Q)*cos((th + 4*pi)/3) - b/3;
        
        vL = min(roots);
        vV = max(roots);
       
    else %one root case
        vV = -sign(R)*((sqrt(R^2 - Q^3) + abs(R))^(1/3) + Q/((R^2 - Q^3)^(1/2) + abs(R))^(1/3)) - b/3;
        vL = 0;
    end
end