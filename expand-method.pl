#/usr/bin/perl -w

# To use maxima instead of mathematica, set the environment variable USE_MAXIMA

use strict;

use Expect::Simple;

my $use_maxima = defined($ENV{"USE_MAXIMA"}) && $ENV{"USE_MAXIMA"};
my $tensor;

# Start up Mathematica or Maxima.
my %attr;
if ($use_maxima) {
  %attr = (Cmd => "maxima",
    Prompt => [-re => '\(\%i\d+\)'],
    DisconnectCmd => "quit();"
  );
} else {
  %attr = (Cmd => "math",
    Prompt => [-re => 'In\[\d+\]\:\=\s+'],
    DisconnectCmd => "Exit"
  );
}

my $mathcmd = Expect::Simple->new(\%attr);

# Convert expression to C.
sub cform {
  my $return;
  if ($use_maxima) {
    if ($tensor) {
      $mathcmd->send("$_[0], numer;");
      $return = $mathcmd->before;
      $return =~ s/.*?\(\%o\d+\)(.*)/\1/s || die;
      $return =~ s/\s//g;
      $return =~ s/\r//g;
      $return =~ s/\\//g;
    } else {
      $return = $_[0];
      $return =~ s/(\w+)\^(\d+)/pow\($1\,$2p\)/g;
      $return =~ s/(\d+)(?!p)/$1\./g;
      $return =~ s/(\d+)p/$1/g;
    }
  } else {
    if ($tensor) {
      $mathcmd->send("CForm[N[$_[0],20]]");
    } else {
      $mathcmd->send("CForm[Simplify[$_[0]]]");
    }
    $return = $mathcmd->before;
    $return =~ s/.*?Out\[\d+\]\/\/CForm\=\s+(.*)/\1/s || die;
    $return =~ s/\s//g;
    $return =~ s/\r//g;
    $return =~ s/Power/pow/g;
    $return =~ s/Sqrt/sqrt/g;
  }
  return $return;
}

# Convert expression to TeX.
sub texform {
  my $return;
  if ($use_maxima) {
    my $in = $_[0];
    $in =~ s/I/\%i/;
    $mathcmd->send("tex($in);");
    $return = $mathcmd->before;
    $return =~ s/.*\$\$(.*)\$\$.*/$1/s || die;
  } else {
    $mathcmd->send("TeXForm[Simplify[$_[0]]]");
    $return = $mathcmd->before;
    $return =~ s/.*?Out\[\d+\]\/\/TeXForm\=\s+(.*)/\1/s || die;
    $return =~ s/\s//g;
    $return =~ s/\r//g;
  }
  return $return;
}

# Simplify expression.
sub simplify {
  my $return;
  if ($use_maxima) {
    my $in = $_[0];
    $in =~ s/\[/\(/g;
    $in =~ s/\]/\)/g;
    $in =~ s/Sqrt/sqrt/g;
    $in =~ s/R/rationalize/g;
    $mathcmd->send("string(fullratsimp($in));");
    $return = $mathcmd->before;
    $return =~ s/.*?\(\%o\d+\)(.*)/\1/s || die;
    $return =~ s/\s//g;
    $return =~ s/\r//g;
    $return =~ s/\\//g;
  } else {
    my $in = $_[0];
    $in =~ s/R/Rationalize/g;
    $mathcmd->send("InputForm[Simplify[$in]]");
    $return = $mathcmd->before;
    $return =~ s/.*?Out\[\d+\]\/\/InputForm\=\s+(.*)/\1/s || die;
    $return =~ s/\s//g;
    $return =~ s/\r//g;
  }
  return $return;
}

# Simplify a whole computation.
sub simplifycomp {
  my @temp;
  my ($c_l,$c_m,$c_r,$c_i);
  my $first = 1;
  foreach my $term (sort split ";",$_[0]) {
    my ($l,$m,$r,$i) = split /\,/,$term;
    if (!$first && $c_l==$l && $c_m==$m) {
      $c_r .= "+($r)";
      $c_i .= "+($i)";
    } else {
      if (!$first) {
        push @temp, "$c_l,$c_m,$c_r,$c_i";
      }
      $c_l = $l;
      $c_m = $m;
      $c_r = "($r)";
      $c_i = "($i)";
      $first = 0;
    }
  }
  if (!$first) {
    push @temp, "$c_l,$c_m,$c_r,$c_i";
  }

  my @result;
  foreach (@temp) {
    my ($l,$m,$r,$i) = split /\,/,$_;
    $r = simplify($r);
    $i = simplify($i);
    if (($r ne '0' && $r ne '0.' && $r ne '0.0') || ($i ne '0' && $i ne '0.' && $i ne '0.0')) {
      push @result, "$l,$m,$r,$i";
    }
  }
  return join";",@result;
}

# Compute $_[0] composed with $_[1].
sub compose {
  my @result;
  foreach my $term1 (split";",$_[1]) {
    my ($l1,$m1,$r1,$i1) = split /\,/,$term1;
    foreach my $term0 (split";",$_[0]) {
      my ($l0,$m0,$r0,$i0) = split /\,/,$term0;
      $r0 =~ s/l/\(l\+\($l1\)\)/g;
      $r0 =~ s/m/\(m\+\($m1\)\)/g;
      $i0 =~ s/l/\(l\+\($l1\)\)/g;
      $i0 =~ s/m/\(m\+\($m1\)\)/g;
      push @result, ($l0+$l1).','.($m0+$m1).",($r0)*($r1)-($i0)*($i1),($r0)*($i1)+($r1)*($i0)";
    }
  }
  return simplifycomp(join";",@result);
}

my %symbol;
$symbol{"I"} = "0,0,0,1";
$symbol{"lm"} = "0,1,-Sqrt[l-m]*Sqrt[l+m+1],0";
$symbol{"lp"} = "0,-1,-Sqrt[l+m]*Sqrt[l-m+1],0";
$symbol{"lx"} = expr("0.5*(lp+lm)");
$symbol{"ly"} = expr("-0.5*I*(lp-lm)");
$symbol{"lz"} = "0,0,-m,0";
$symbol{"z"} = "-1,0,Sqrt[l+m]*Sqrt[l-m]/Sqrt[2*l-1]/Sqrt[2*l+1],0;".
               "1,0,Sqrt[l+m+1]*Sqrt[l-m+1]/Sqrt[2*l+1]/Sqrt[2*l+3],0";
$symbol{"x"} = expr("I*(z*ly-ly*z)");
$symbol{"y"} = expr("I*(lx*z-z*lx)");
$symbol{"dx"} = expr("I*(z*ly-y*lz)");
$symbol{"dy"} = expr("I*(x*lz-z*lx)");
$symbol{"dz"} = expr("I*(y*lx-x*ly)");
$symbol{"lap"} = "0,0,-l*(l+1),0";

sub expr {
  my @sub_expr;
  my $count = 0;

  my $line = $_[0];
# Evaluate expressions inside brackets.
  while ($line =~ s/\(([^\(]*?)\)/\@$count/) {
    return "Error: cannot use brackets in tensor mode" if $tensor;
    $sub_expr[$count] = expr($1);
    return $sub_expr[$count] if $sub_expr[$count] =~ /^Error\:/;
    $count++;
  }
# Evaluate numbers.
  while ($line =~s/(?<![\@\/\d\.])([\/\d\.]+)/\@$count/) {
    $sub_expr[$count] = simplifycomp "0,0,R[$1],0";
    $count++;
  }
# Evaluate symbols.
  while ($line =~s/(?<![\@\w])(\w+)/\@$count/) {
    if (defined $symbol{$1}) {
      $sub_expr[$count] = $symbol{$1};
      $count++;
    } else {
      return "Error: unknown symbol";
    }
  }
# Evaluate multiplications from right to left.
  $line = scalar reverse $line;
  while ($line =~s/(\d+)\@\*(\d+)\@/scalar reverse "\@$count"/e) {
    $sub_expr[$count] = compose($sub_expr[scalar reverse $2],$sub_expr[scalar reverse $1]);
    $count++;
  }
  $line = scalar reverse $line;
# Evaluate any unary subtraction at beginning of expression.
  if ($line =~ s/^\-\@(\d+)/\@$count/) {
    $sub_expr[$count] = compose("0,0,-1,0",$sub_expr[$1]);
    $count++;
  }
# Evaluate any unary addition at beginning of expression.
  if ($line =~ s/^\+\@(\d+)/\@$count/) {
    $sub_expr[$count] = $sub_expr[$1];
    $count++;
  }
# Evaluate addition and subtraction from left to right.
  while ($line =~s/\@(\d+)(\+|\-)\@(\d+)/\@$count/) {
    if ($2 eq '+') {
      $sub_expr[$count] = simplifycomp("$sub_expr[$1];$sub_expr[$3]");
    } else {
      $sub_expr[$count] = simplifycomp($sub_expr[$1].';'.compose("0,0,-1,0",$sub_expr[$3]));
    }
    $count++;
  }
# Return final answer.
  if ($line =~ /^\@(\d+)$/) {
    return $sub_expr[$1];
  } else {
    return "Error: syntax error";
  }
}

$tensor = defined($ENV{"TENSOR"}) && $ENV{"TENSOR"};

if ($tensor) {
  if ($use_maxima) {
    $mathcmd->send("l:0;m:0;");
  } else {
    $mathcmd->send("l=0;m=0");
  }
  my $junk = $mathcmd->before;
}

print "> ";
while (my $line=<>) {
  chomp $line;
  exit if $line eq 'q';

  my $expr;
  my $latex;
  my $cas;
  if ($line =~ /^(\w+)\s*\=\s*(.*)/) {
    $symbol{$1} = $expr = expr($2);
  } else {
    $latex = ($line =~ s/^t\s+//);
    $cas = ($line =~ s/^c\s+//);
    $expr = expr($line);
  }

  if ($expr =~ /^Error\:/) {
    print $expr;
  } else {
    my $first = 1;
    foreach my $c (split";",$expr) {
      my ($l,$m,$r,$i) = split ',',$c;
      if ($latex) {
        my $n = texform("$r+I*($i)");
        if (($n ne '0' && $n ne '0.' && $n ne '0.0')) {
          print "\n" if !$first;
          $first = 0;
          if ($tensor) {
            print "c_{$l,$m} = $n ";
          } else {
            $l = "+$l" if $l>0;
            $l = "" if $l==0;
            $m = "+$m" if $m>0;
            $m = "" if $m==0;
            print "c_{l,l$l}^{m,m$m} = $n ";
          }
        }
      } elsif ($cas) {
        print "\n" if !$first;
        $first = 0;
        print "$l;$m;$r;$i";
      } else {
        $r = cform($r);
        $i = cform($i);
        if (($r ne '0' && $r ne '0.' && $r ne '0.0') || ($i ne '0' && $i ne '0.' && $i ne '0.0')) {
          print "\n" if !$first;
          $first = 0;
          print "$l;$m;($r)+($i)*I";
        }
      }
    }
  }
  print "\n> ";
}
