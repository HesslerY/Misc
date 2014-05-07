

function fe = felem(theta, phi, S)

switch S
    case 1 %Antenne isotrope
        fe = 1;
    case 2  %Dipôle court suivant axe z
        fe = sin(theta);
    case 3  %Dipôle court suivant axe x
        fe = sin(acos(cos(phi).*sin(theta)));
    case 4 %Dipôle court suivant axe y
        fe = sin(acos(sin(phi).*sin(theta)));
end        