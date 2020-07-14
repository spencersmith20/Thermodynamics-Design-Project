clear
clc

R = 8.3144598;

%i. import file contents
fn = input('Please enter the name of your file: ', 's'); %get filename

if ~exist(fn, 'file') %check for file existence
    fn = input('Try again! Enter file name: ', 's');
end

file = fopen(fn, 'r'); %get file Id

raw = fscanf(file, '%f'); %array of all file data

n = raw(1); %number of components

%format: Tc Pc w x
compData = zeros(n, 4); %initializes pure component properties matrix

for i = 1:n %each component's props
    current = (i-1)*4 + 1; %sets "currrent var", initializes for each row
    for k = 1:4
        if k == 2 %converts bar to Pa
            compData(i,k) = raw(current + k)*10^5;
        else %all other data stays the same
            compData(i,k) = raw(current + k);
        end
    end
end

fclose(file);

current = n*4 + 2; %sets new current
T = raw(current); %temp and pressure of system
P = raw(current + 1)*10^5;
    
badboy = badData(compData, n);

while (badboy == 1) %fix the mole fractions
    s = 0;
    fprintf('you have some issues in your data. let us fix them.\n');
    
    for i=1:n %replace mole fractions
       fprintf('\nenter a mole fraction for component %d of %d (current sum is %1.1f)\n', i, n, s);
       compData(i,4) = input('enter here ');
       fprintf('\n');
       
       s = s + compData(i,4); %increment mole fraction sum
    end
    
    badboy = badData(compData, n); %check for error again
    
    if badboy == 1 %if there is still an issue, then stay in while loop
        fprintf('\nyou still did not do it right. try again!\n');
    end
end

if badboy == 0 %if there are no mole fraction discrepencies

    % ii. calculate mixing coefficients 

    a = .45724*R^2*(compData(:, 1).^2)./(compData(:, 2));
    b = .0778*R*compData(:,1)./(compData(:, 2));
    k = .37464 + 1.54226*compData(:, 3) - .26992*(compData(:,3).^2);
    alpha = (1 + k.*(1 - sqrt(T./compData(:, 1)))).^2;
    aAlpha = a.*alpha;

    aij = zeros(n);
    amix = 0;
    bmix = 0;

    for i=1:n %calculate aij, binary mixing coefficients
        for j=1:n
            if n > 1
                aij(i,j) = sqrt(aAlpha(i).*aAlpha(j));
            else
                aij(i,j) = aAlpha(i,j);
            end
        end
    end

    for i=1:n %calculate amix and bmix
        bmix = bmix + compData(i, 4)*b(i);

        for j=1:n
            if n > 1
               amix = amix + compData(i, 4)*compData(j,4)*aij(i,j);
            else
               amix = amix + compData(i,4)*compData(i,4)*aij(i,j); 
            end
        end
    end

    vL = zeros(1, n+1);
    vV = vL;

    for i=1:n+1 %calculate roots
        if i ~= n+1 %calculate pure roots
           [vL(i), vV(i)] = calcRoots(aAlpha(i), b(i), T, P);
        else %calculate mix roots (index n+1)
           [vL(i), vV(i)] = calcRoots(amix, bmix, T, P);
        end
    end

    h = zeros(n, 2); %liq, vap; liq, vap; liq, vap

    for i=1:n %departure function for each component in mixture
        v = [vL(i), vV(i)];
        Tr = T/compData(:,1); %reduced temperature

        for j=1:2 %departure function for each volume root calculated
            if v(j) == 0 %for one root case
                h(i,j) = 0;
            else %for two-root case
                h(i,j) = P*v(j) - R*T - (aAlpha(i)/(2*sqrt(2)*b(i)))*(1 + k(i)*sqrt(Tr(i)/aAlpha(i)))*log((v(j) + (1 + sqrt(2))*b(i))/(v(j) + (1 - sqrt(2))*b(i)));
            end
        end
    end

    logPhi = zeros(2, n+1); %initialize results matrices
    phi = logPhi;

    v = [vL(n+1), vV(n+1)]; %mixture volume arrays
    z = P*v./(R*T);

    for i=1:n+1 %for each componnet    

        temp = 0;
        if i ~= n+1 %don't do this for mixing
            for j = 1:n 
                temp = temp + compData(j,4)*aij(i,j);
            end
        end

        for k = 1:2 %for each volume root
            if v(k) == 0 %if volume root does not exist
                logPhi(k,i) = 0;

            else %for existant volume root
                if i == n+1 %use mixture!
                    logPhi(k,i) = z(k) - 1 - log((v(k) - bmix)*P/(R*T)) - (amix/(2*sqrt(2)*bmix*R*T))*log((v(k)+(1+sqrt(2))*bmix)/(v(k)+(1-sqrt(2))*bmix));

                else %use individual
                    logPhi(k,i) = (b(i)/bmix)*(z(k)-1) - log((v(k)-bmix)*P/(R*T));

                    %calculate final value of logPhi
                    logPhi(k,i) = logPhi(k,i) + (amix/(2*sqrt(2)*bmix*R*T))*(b(i)/bmix - 2*temp/amix)*log((v(k)+(1+sqrt(2))*bmix)/((v(k)+(1-sqrt(2))*bmix)));
                end

                phi(k,i) = exp(logPhi(k,i)); %calculate phi
            end
        end
    end

    fug = phi.*P/(10^5); %put fugacities into bar
    h = h./1000; %put hdep into kJ/mol

    %prompts user to enter name of output file:
    fn = input('Please enter the name of your output file: ', 's');
    file = fopen(fn, 'wt'); %get file ID

    %header for file and output
    fprintf(file, 'peng-robinson equation of state for %d-component system @ %.2f bar and %.2f K\n\n', n, P/10^5, T);
    fprintf('\n\npeng-robinson equation of state for %d-component system @ %.2f bar and %.2f K\n\n', n, P/10^5, T);

    %print critical data and mole fractions for each 
    for i = 1:n
       fprintf('COMPONENT %d: Tc = %.2f K, Pc = %.2f bar, w = %.3f, x = %.2f\n', i, compData(i,1), compData(i,2)/10^5, compData(i,3), compData(i,4));
       fprintf(file, 'COMPONENT %d: Tc = %.2f K, Pc = %.2f bar, w = %.3f, x = %.2f\n', i, compData(i,1), compData(i,2)/10^5, compData(i,3), compData(i,4)); 
    end

    fprintf('\nPURE PROPERTIES\n\n');
    fprintf(file, '\nPURE PROPERTIES\n\n');    

    for i=1:n %display PURE species volumes and fugacities
        fprintf('COMPONENT %d\n', i);
        fprintf(file, 'COMPONENT %d\n', i);

        switch(vL(i))
            case 0 %if second element == 0, there is only one root

                %print volume for one root
                fprintf('volume: %.8f m^3/mol\n', vV(i));
                fprintf(file, 'volume: %.8f m^3/mol\n', vV(i));

                %departure function
                fprintf('enthalpy departure function: %.3f kJ/mol\n\n', h(i,2));
                fprintf(file, 'enthalpy departure function: %.3f kJ/mol\n\n', h(i,2));

                %fugacity coefficient
                fprintf('mixture fugacity coefficient: %.4f\n\n', phi(2,i));
                fprintf(file, 'mixture fugacity coefficient: %.4f\n\n', phi(2,i));

                %fugacities 
                fprintf('fugacity: %.4f bar\n\n', fug(2,i));
                fprintf(file, 'fugacity: %.4f bar\n\n', fug(2,i));

            otherwise %three real volume roots

                %volumes
                fprintf('liquid volume: %.8f m^3/mol\nvapor volume: %.8f m^3/mol\n\n', vL(i), vV(i));
                fprintf(file, 'liquid volume: %.8f m^3/mol\nvapor volume: %.8f m^3/mol\n\n', vL(i), vV(i));

                %departure functions
                fprintf('liquid enthalpy departure: %.3f kJ/mol\nvapor enthalpy departure: %.3f kJ/mol\n\n', h(i,1), h(i,2));
                fprintf(file, 'liquid enthalpy departure: %.3f kJ/mol\nvapor enthalpy departure: %.3f kJ/mol\n\n', h(i,1), h(i,2));

                if vL(n+1) ~= 0 %mixture exists as vapor and liquid
                    %fugacity coefficients
                    fprintf('liquid fugacity coefficient: %.4f\nvapor fugacity coefficient: %.4f\n\n', phi(1,i), phi(2,i));
                    fprintf(file, 'liquid fugacity coefficient: %.4f\nvapor fugacity coefficient: %.4f\n\n', phi(1,i), phi(2,i));

                    %fugacities
                    fprintf('liquid fugacity: %.4f bar\nvapor fugacity: %.4f bar\n\n', fug(1,i), fug(2,i));
                    fprintf(file, 'liquid fugacity: %.4f bar\nvapor fugacity: %.4f bar\n\n', fug(1,i), fug(2,i));

                else %mixture does not exist in two phases
                    fprintf('mixture fugacity coefficient: %.4f\n\n', phi(2,i));
                    fprintf(file, 'mixture fugacity coefficient: %.4f\n\n', phi(2,i));

                    %fugacities 
                    fprintf('fugacity: %.4f bar\n\n', fug(2,i));
                    fprintf(file, 'fugacity: %.4f bar\n\n', fug(2,i)); 
                end
        end
    end

    if n > 1 %for non-pure component systems

        fprintf('MIXTURE PROPERTIES\n\n');
        fprintf(file, 'MIXTURE PROPERTIES\n\n');

        %print volume roots and total fugacity coefficients
        if vL(n+1) == 0 %only one mixture root

            %volume
            fprintf('volume: %.8f m^3/mol\n\n', vV(n+1));
            fprintf(file, 'volume: %.8f m^3/mol\n\n', vV(n+1));

            %fugacity coefficient
            fprintf('fugacity coefficient: %.3f\n\n', phi(2, n+1));
            fprintf(file, 'fugacity coefficient: %.3f\n\n', phi(2, n+1));

            %fugacity
            fprintf('fugacity: %.4f bar', fug(2, n+1));
            fprintf(file, 'fugacity: %.4f bar', fug(2, n+1));

        else %three mixture roots

            %volumes
            fprintf('liquid volume: %.8f m^3/mol\nvapor volume: %.8f m^3/mol\n\n', vL(n+1), vV(n+1));
            fprintf(file, 'liquid volume: %.8f m^3/mol\nvapor volume: %.8f m^3/mol\n\n', vL(n+1), vV(n+1));

            %fugacity coefficients
            fprintf('liquid fugacity coefficient: %.4f\nvapor fugacity coefficient: %.4f\n\n', phi(1, n+1), phi(2, n+1));
            fprintf(file, 'liquid fugacity coefficient: %.4f\nvapor fugacity coefficient: %.4f\n\n', phi(1, n+1), phi(2, n+1));

            %fugacities
            fprintf('liquid fugacity: %.4f bar\nvapor fugacity: %.4f bar\n\n', fug(1, n+1), fug(2, n+1));
            fprintf(file, 'liquid fugacity: %.4f bar\nvapor fugacity: %.4f bar\n\n', fug(1, n+1), fug(2, n+1));
        end
    end

    fclose(file);
end