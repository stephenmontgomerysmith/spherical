#/usr/bin/perl -w

use strict;

use Expect::Simple;

my $rank = $ARGV[0];
die if $rank =~ /\D/ || $rank<=0; 

my %attr = (Cmd => "env TENSOR=1 OUTPUT_REAL=1 perl expand-method.pl",
  Prompt => '> ',
  DisconnectCmd => "q"
);

my $methodtensor = Expect::Simple->new(\%attr);

sub comp {
  my $c = $_[0];
  $c =~ s/0/x/g;
  $c =~ s/1/y/g;
  $c =~ s/2/z/g;
  $methodtensor->send(join "*",split"",$c);
  my $response = $methodtensor->before;
  $response =~ s/\r//g;
  $response =~ s/.*?\n//;
  chomp($response);
  return $response;
}

print "/* Created by \"perl $0 $ARGV[0]\". */\n";
print "/* Returns 1/sqrt(4 pi) times moment tensor of rank $rank. */\n";
print "void tensor${rank}(REAL *psi, REAL a", "[3]" x $rank, ") {\n";

my @index;
for (my $i=0;$i<$rank;$i++) {
  $index[$i] = 0;
}
my %comp;
while (1) {
  my $sort = join('',sort @index);
  if (!defined($comp{$sort})) {
    my $response = comp($sort);
    foreach (split "\n", $response) {
      my ($l,$m,$r,$i) = split ";",$_;
      next if $m<0;
      if ($r != 0) {
        $r = "($r)";
        $r = "2*$r" if $m>0;
        if (!defined($comp{$sort})) {
          $comp{$sort} .= "    $r*psi[ind($l,$m,0)]";
        } else {
          $comp{$sort} .= "\n    + $r*psi[ind($l,$m,0)]";
        }
      }
      if ($i != 0) {
        $i = -$i;
        $i = "($i)";
        $i = "2*$i" if $m>0;
        if (!defined($comp{$sort})) {
          $comp{$sort} .= "    $i*psi[ind($l,$m,1)]";
        } else {
          $comp{$sort} .= "\n    + $i*psi[ind($l,$m,1)]";
        }
      }
    }
    $comp{$sort} .= ";\n";
  }
  $comp{$sort} = "  a".join('',map("[$_]",@index))."=\n".$comp{$sort};

  my $i = $rank-1;
  while ($i>=0 && ++$index[$i]==3) {
    $index[$i] = 0;
    $i--;
  }
  last if $i<0;
}
foreach my $key (sort keys %comp) {
  print $comp{$key};
}

print "}\n";
