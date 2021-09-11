#/usr/bin/perl -w

use strict;

use Expect::Simple;

my $use_maxima = defined($ENV{"USE_MAXIMA"}) && $ENV{"USE_MAXIMA"};

my $rank = $ARGV[0];
die if $rank =~ /\D/ || $rank<=0; 

my %attr = (Cmd => "env TENSOR=1 perl expand-method.pl",
  Prompt => '> ',
  DisconnectCmd => "q"
);

my $methodtensor = Expect::Simple->new(\%attr);

sub comp {
  my $c = $_[0];
  $c =~ s/0/x/g;
  $c =~ s/1/y/g;
  $c =~ s/2/z/g;
  $methodtensor->send("c ".join "*",split"",$c);
  my $response = $methodtensor->before;
  $response =~ s/\r//g;
  $response =~ s/.*?\n//;
  chomp($response);
  return $response;
}

my @equ;
my @index;
for (my $i=0;$i<$rank;$i++) {
  $index[$i] = 0;
}
while (1) {
  my $equ = "a".join('',@index).($use_maxima?"=":"==");
  my $response = comp(join('',@index));
  foreach (split "\n", $response) {
    my ($l,$m,$r,$i) = split ";",$_;
    $m = "m".(-$m) if $m<0;
    $equ .= "+(($r)+I*($i))*psil${l}m$m"
  }
  @equ = (@equ,$equ);
  my $i = $rank-1;
  while ($i>=0 && $index[$i]==2) {
    $i--;
  }
  last if $i<0;
  $index[$i]++;
  while (++$i<$rank) {
    $index[$i] = $index[$i-1];
  }
}

undef($methodtensor);

my %var;
my $string = join ",",sort @equ;
while ($string =~ s/(psi\w+)//) {
  $var{$1}=1;
}

my $input;
my $return;
if ($use_maxima) {
  $input = "string(fullratsimp(solve([".join(",\n",@equ)."],\n[".join(",\n",keys %var)."])));\n";
  open(TMP,"> temp-file");
  print TMP "$input\n";
  close(TMP);
  $return = `maxima < temp-file`;
  $return =~ s/.*?\(\%o\d+\)(.*)/\1/s || die;
  $return =~ s/\s//g;
  $return =~ s/\r//g;
  $return =~ s/\\//g;
} else {
  $input = "InputForm[N[Expand[Solve[{".join(",\n",@equ)."},\n{".join(",\n",keys %var)."}]]]]\n";
  open(TMP,"> temp-file");
  print TMP "$input\n";
  close(TMP);
  $return = `math < temp-file`;
  unlink "temp-file";
  $return =~ s/.*?Out\[\d+\]\/\/InputForm\=\s+(.*)/\1/s || die;
  $return =~ s/\s//g;
  $return =~ s/\r//g;
}


print "/* Created by \"perl $0 $ARGV[0]\". */\n";
print "/* Adds to psi sqrt(4 pi) times moment tensor of rank $rank. */\n";
print "void reverse_tensor${rank}(double a", "[3]" x $rank,", COMPLEX *psi) {\n";

my @psi_equ;

if ($use_maxima) {
  while ($return =~ s/(psi.*?)\=(.*?)[\,\}]//) {
    my $lhs = $1;
    my $rhs = $2;
    next if $lhs =~ /mm/;
    $lhs =~ s/l(\d+)m(\d+)/\[ind\($1\,$2\)\]/ || die;
    $rhs =~ s/([\+\-])/\n                   $1/g;
    $rhs =~ s/\s+$//;
    $rhs =~ s/^\s+//;
    while ($rhs =~ s/(a\d*)(\d)/$1\[$2\]/) {}
    @psi_equ = (@psi_equ,"  $lhs += $rhs;\n");
  }
} else {
  while ($return =~ s/(psi.*?)\-\>(.*?)[\,\}]//) {
    my $lhs = $1;
    my $rhs = $2;
    next if $lhs =~ /mm/;
    $lhs =~ s/l(\d+)m(\d+)/\[ind\($1\,$2\)\]/ || die;
    $rhs =~ s/(a\d+)/$1\n                   /g;
    $rhs =~ s/\s+$//;
    while ($rhs =~ s/(a\d*)(\d)/$1\[$2\]/) {}
    @psi_equ = (@psi_equ,"  $lhs += $rhs;\n");
  }
}

print join'',sort @psi_equ;

print "}\n";
