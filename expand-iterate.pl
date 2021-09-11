#/usr/bin/perl -w

use strict;

use Expect::Simple;

sub plus {
  return $_[0] if $_[1]==0;
  return "$_[0]+$_[1]" if $_[1]>0;
  return "$_[0]$_[1]" if $_[1]<0;
}

my $filename = $ARGV[0];
open(IN,$filename);

my $main_program;
my @init_proc_assign;

my $nr_iterates = 0;
my @mc_count;

my %attr = (Cmd => "env TENSOR=0 perl expand-method.pl",
  Prompt => '> ',
  DisconnectCmd => "q"
);

my $methodcmd = Expect::Simple->new(\%attr);

while(my $line = <IN>) {
  my $text;
  if ($line=~/\@spherical_iterate\s*\{/) {
    while ($line = <IN>) {
      if ($line =~ /\}/) {
        last;
      } else {
        $text .= $line;
      }
    }
    $text =~s/^[\s^\n]*\n//;
    $text =~s/\s*$//;

    $mc_count[$nr_iterates] = 0;

    $main_program .= "  {\n";
    $main_program .= "    if (first_mc_$nr_iterates) initialize_mc_$nr_iterates();\n";
    $main_program .= "    /* Start-Threading */\n";
    $main_program .= "    int lstart=0, lend=max_order+2;\n";
    $main_program .= "    int l,m;\n";
    $main_program .= "    COMPLEX *mc;\n";
    $main_program .= "    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {\n";
    $main_program .= "      mc = mult_constant_$nr_iterates + mc_count_$nr_iterates*ind(l,m);\n";

    foreach $line (split "\n",$text) {
      $line =~ s/\@index/ind(l,m)/g;
      $line = "  $line";
      while ($line =~ s/^(.*?)\@method\((.*?)\,(.*?)\)//) {
        my $line_start = "$1(";
        my $var = $2;
        my $send_cmd = $3;
        while (1) {
          my $left = $send_cmd;
          my $right = $send_cmd;
          $left =~ s/[^\(]//g;
          $right =~ s/[^\)]//g;
          last if length($left) == length($right);
          $line =~ s/^(.*?)\)// || die;
          $send_cmd .= ")$1";
        }
        $methodcmd->send($send_cmd);
        my $response = $methodcmd->before;
        $response =~ s/\r//g;
        $response =~ s/.*?\n//;
        if ($response =~ /^Error\:/) {
          die "Error in \@method";
        }
        my @line_part;
        foreach my $comp(split "\n",$response) {
          my ($l,$m,$c) = split ';',$comp;
          $c =~ s/([lm])/\1\1/g;
          $init_proc_assign[$nr_iterates] .= "    if (abs(".plus("m",$m).")<=".plus("l",$l).")\n".
                                             "      mc[$mc_count[$nr_iterates]] = $c;\n".
                                             "    else\n".
                                             "      mc[$mc_count[$nr_iterates]] = 0;\n";
          push @line_part, "mc[$mc_count[$nr_iterates]]*index($var,".plus("l",$l).",".plus("m",$m).")";
          $mc_count[$nr_iterates]++;
        }
        if ($#line_part >= 0) {
          $main_program .= $line_start.join("+",@line_part).")";
        } else {
          $main_program .= "${line_start}0)";
        }
      }
      $main_program .= "$line\n";
    }
    $main_program .= "    }\n";
    $main_program .= "    /* Stop-Threading */\n";
    $main_program .= "  }\n";
    $nr_iterates++;
  } else {
    $main_program .= $line;
  }
}

print "/*\n";
print " * Created from $filename by $0.\n";
print " */\n";
print "\n";
print "\#include \"spherical.h\"\n";
print "\n";
print "extern int max_order;\n";

for (my $i=0;$i<$nr_iterates;$i++) {
  print "static int first_mc_$i = 1;\n";
  print "static COMPLEX *mult_constant_$i;\n";
  print "#define mc_count_$i $mc_count[$i]\n";
  print "static void initialize_mc_$i();\n";
  print "\n";
}

print "/* Declarations of threading go here. */\n\n";
print $main_program;

for (my $i=0;$i<$nr_iterates;$i++) {
  chomp $init_proc_assign[$i];
  print "\n";
  print "static void initialize_mc_$i() {\n";
  print "  int l,m;\n";
  print "  double ll,mm;\n";
  print "  COMPLEX *mc;\n";
  print "\n";
  print "  first_mc_$i = 0;\n";
  print "  mult_constant_$i = malloc(sizeof(COMPLEX)*length*mc_count_$i);\n";
  print "  for (l=0;l<=max_order;l+=2) for (m=0;m<=l;m++) {\n";
  print "    ll = l;\n";
  print "    mm = m;\n";
  print "    mc = mult_constant_$i + mc_count_$i*ind(l,m);\n";
  print "$init_proc_assign[$i]\n";
  print "  }\n";
  print "}\n";
}
