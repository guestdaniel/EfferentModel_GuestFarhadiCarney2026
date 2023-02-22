tdres = 100e3
TWOPI = 2π 
C = 2.0/tdres
Fc = 1000.0
c1LP = ( C - TWOPI*Fc ) / ( C + TWOPI*Fc )
c2LP = TWOPI*Fc / (TWOPI*Fc + C)

ihc = zeros(8)
ihcl = zeros(8)
gain = 1.0
x = 3.7

ihc[1] = x*gain;

for i in 1:7
    ihc[i+1] = c1LP*ihcl[i+1] + c2LP*(ihc[i]+ihcl[i]);
end
for j in 1:7
    ihcl[j] = ihc[j]
end