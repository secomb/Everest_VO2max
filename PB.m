function Patm = PB(alt)
    coeff = [-0.00149 -0.1112 6.63268];
    Patm = exp(polyval(coeff,alt)); % West 1999
end