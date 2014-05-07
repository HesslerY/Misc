for i=1:101
    for j=1:1000
        if tempsessaisansNAN(i,j)==NaN
            tempsessaisansNAN(i,j)=10^12;
        end
    end
end


B=isnan(tempsessaisansNAN);

tempsessaisansNAN(B)=10^12