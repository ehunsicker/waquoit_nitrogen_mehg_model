function [meani, vari] = welford( datum, meani, vari, n )
% calculate variance and mean on the fly during iterations
S = vari*(n-2);
del = datum - meani;
meani = meani + del/n;
del2 = datum - meani;
S = S + del.*del2;
vari = S/(n-1);

end

