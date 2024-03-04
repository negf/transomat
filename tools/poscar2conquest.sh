cat $1 | awk 'BEGIN{bohr=.52917720859; loutput=0};\
function atom2number(species, atom) {for (i in species) {if (species[i]==atom) {return i} }} \
function atom2number2(nspecies, j) {for (i=1; i<=length(nspecies);i++) { if (j<=nspecies[i]) {return i} }} \
function mod(a, b) {m = a % b ; if (m<0) {m = m + 1} ; return m} \
{if (NR==3) {split($0,a," ") ; printf ("%24.12e %24.12e %24.12e \n", a[1]/bohr, a[2]/bohr, a[3]/bohr)} \
else if (NR==4) {split($0,b," "); printf ("%24.12e %24.12e %24.12e \n", b[1]/bohr, b[2]/bohr, b[3]/bohr)} \
else if (NR==5) {split($0,c," "); printf ("%24.12e %24.12e %24.12e \n", c[1]/bohr, c[2]/bohr, c[3]/bohr)} \
else if (NR==6) {split($0,species," ")} \
else if (NR==7) {split($0,nspecies," "); for(i=1; i<=NF;i++) {nat+=$i ; nspecies[i]=nat}; print nat} \
else if (tolower($1)=="direct") {loutput=1; start_line=NR; a[1]=1.0;b[2]=1.0;c[3]=1.0} \
else if (tolower($1)=="cartesian") {loutput=1; start_line=NR} \
else if (loutput==1) {if ($1=="") {exit};printf ("%24.12e %24.12e %24.12e %3i T T T\n", \
  mod($1/a[1],1),mod($2/b[2],1),mod($3/c[3],1),atom2number2(nspecies,NR-start_line))}}'
