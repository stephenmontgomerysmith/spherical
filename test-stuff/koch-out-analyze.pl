@lines = split /(?=C1)/,join '',<>;

print "$#lines\n";

while (1) {
  open(M,"| math");
  $r = int rand $#lines;
  print "$r\n";
  print M $lines[$r];
  print M "CC = C1*D1*IdentityMatrix[3]+C2*D2;\n";
  print M "CC2 = C3*D2;\n";
  print M "a4xxCC = Table[Sum[CC[[i,j]]a4[[i,j,k,l]],{i,1,3},{j,1,3}],{k,1,3},{l,1,3}];\n";
  print M "adot22 = 2CC-2Tr[CC]a2-5CC.a2-5a2.CC+10a4xxCC;\n";
  print M "adot23 = 2 Tr[CC2] IdentityMatrix[3] - 2 Tr[a2.CC2] IdentityMatrix[3] + 3(a2.CC2+CC2.a2) - 4 Tr[CC2] a2 - 2CC2;\n";
  print M "m=Chop[adot2-adot22-adot23];\n";
  print M "Tr[m.m]\n";
}
close M;
