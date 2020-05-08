function [airy] = airyjianmingjin(x)

c1 = 3^(-2/3) / gamma(2/3); % notation same as A&S p. 446, eqs. 10.4.4 and 10.4.5jjjj
c2 = 3^(-1/3) / gamma(1/3); 

v = 1.0/3.0;
gammaoneplusv = gamma(1 + 1.0/3.0); % 0.892979511569249;
gammaoneminusv = gamma(1.0 - 1.0/3.0); % 1.3541179394264;

nx = length(x);
for ix = 1:nx
    if x(ix) < 0.0
        airy(ix) = 0.0;
    elseif x(ix) == 0.0
        airy(ix) = c1;
    else
        
        z = 2/3 * x(ix)^(3/2);
        z2 = z * z;
          
        if (z <= 9.0)
        % use 9.6.10 on p. 375 of A&S to compute I_{-1/3} (z)
            sumasymptotic = 1.0;
            ratio = 1.0;
            for k = 1:60
                ratio = ratio * z2 / (4 * k * (k - v));
                sumasymptotic = sumasymptotic + ratio;
            end
            Iminusonethird = sumasymptotic * (2.0 / z)^v / gammaoneminusv;
                     
            % I_{1/3} (z) - calculate the modified Bessel function of the first kind of order 1/3       
            % with eq. 9.6.10 on p. 375 in A&S
            sumasymptotic = 1.0;
            ratio = 1.0;
            for k = 1:40
                ratio = ratio * z2 / ( 4.0 * k * ( k + v) ); % divide by k is the factorial term, divide by (k + vl) is the gamma function (without the gp1 term)
                sumasymptotic = sumasymptotic + ratio;
            end
            Iplusonethird = (0.5 * z)^v / gammaoneplusv * sumasymptotic;            
            
            Konethird = 0.5 * pi * (Iminusonethird - Iplusonethird) / sin(pi / 3);
        else % use the asymptotic formula for K_{1/3} from 9.7.2 on p. 378 in A&S

            if (z < 35.0) 
                k0 = 12;
            elseif ( z < 50.0 ) 
                k0 = 10;
            else % z > 50
                k0 = 8;
            end    
            
            sumasymptotic = 1.0;
            ratio = 1.0;
            for k = 1:k0
                % 'k' in the denominator is the factorial in 9.7.2
                ratio = ratio * (4.0 / 9.0 - ( 2.0 * k - 1.0 )^2) / (8.0 * k * z);
                sumasymptotic = sumasymptotic + ratio;
            end
            Konethird = exp(-z) * sqrt(0.5 * pi / z) * sumasymptotic;
            
%             Kwhittaker = sqrt(pi / 2 / z) * whittakerW(0, 1/3, 2 * z);
%             Kkummmer = sqrt(pi / 2 / z) * exp(- z) * (2 * z)^(5/6) * kummerU(5/6, 5/3, 2 * z);


        end
     
        airy(ix) = sqrt(x(ix)) / sqrt(3) * Konethird / pi;    
    end
end