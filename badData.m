function badboy = badData(compData, n)
    
    badboy = 0;
    
    for i=1:n %check for negative mole fractions
        if compData(i,4) < 0
            badboy = 1;
            fprintf('\ncomponent %d has a negative mole fraction!\n', i);
        end
    end

    if(sum(compData(:,4)) > 1) %check for sum greater than 1
        badboy = 1;
        fprintf('\nthe sum of your mole fractions is greater than one!\n');
    end
    
    if (sum(compData(:, 4)) < 1)
        badboy = 1;
        fprintf('\nthe sum of your mole fractions is less than one!\n');
    end
end