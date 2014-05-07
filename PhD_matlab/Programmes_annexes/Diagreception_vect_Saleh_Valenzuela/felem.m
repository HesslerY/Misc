

function fe = felem(theta, phi, S)

switch S
    case 1 %Antenne isotrope
        fe = 1;
    case 2  %Dip�le court suivant axe z
        fe = sin(theta);
    case 3  %Dip�le court suivant axe x
        fe = sin(acos(cos(phi).*sin(theta)));
    case 4 %Dip�le court suivant axe y
        fe = sin(acos(sin(phi).*sin(theta)));
end        