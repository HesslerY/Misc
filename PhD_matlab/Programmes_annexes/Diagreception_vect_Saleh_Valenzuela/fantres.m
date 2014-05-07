

function fa = fantres(fr, theta, phi, S)

fe = felem(theta, phi, S);

fa = abs(fr) .* fe;

fa = fa / max(max(fa));
%fa = fa / mean(mean(fa));