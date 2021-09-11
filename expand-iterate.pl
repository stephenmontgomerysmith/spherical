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

my %attr = (Cmd => "env TENSOR=0 OUTPUT_REAL=1 perl expand-method.pl",
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

    $main_program .= "  if (first_mc_$nr_iterates) initialize_mc_$nr_iterates(param);\n";
    $main_program .= "  {\n";
    $main_program .= "    /* Start-Threading */\n";
    $main_program .= "    int lstart=0, lend=param->max_order+2;\n";
    $main_program .= "    int l,m;\n";
    $main_program .= "    for (l=lstart;l<lend;l+=2) for (m=0;m<=l;m++) {\n";
    $main_program .= "      REAL *mc = mult_constant_$nr_iterates + mc_count_$nr_iterates*mc_ind(l,m);\n";
    my $main_program_r;
    my $main_program_i;
    foreach $line (split "\n",$text) {
      $line =~ s/\@index/ind(l,m,\@c)/g;
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
        my @line_part_r;
        my @line_part_i;
        foreach my $comp(split "\n",$response) {
          my ($l,$m,$r,$i) = split ';',$comp;
          $r =~ s/([lm])/\1\1/g;
          $i =~ s/([lm])/\1\1/g;
          if ($r ne '0') {
            $init_proc_assign[$nr_iterates] .= "    if (abs(".plus("m",$m).")<=".plus("l",$l).")\n".
                                               "      mc[$mc_count[$nr_iterates]] = $r;\n".
                                               "    else\n".
                                               "      mc[$mc_count[$nr_iterates]] = 0;\n";
            push @line_part_r, "mc[$mc_count[$nr_iterates]]*${var}[ind(".plus("l",$l).",".plus("m",$m).",0)]";
            push @line_part_i, "mc[$mc_count[$nr_iterates]]*${var}[ind(".plus("l",$l).",".plus("m",$m).",1)]";
            $mc_count[$nr_iterates]++;
          }
          if ($i ne '0') {
            $init_proc_assign[$nr_iterates] .= "    if (abs(".plus("m",$m).")<=".plus("l",$l).")\n".
                                               "      mc[$mc_count[$nr_iterates]] = $i;\n".
                                               "    else\n".
                                               "      mc[$mc_count[$nr_iterates]] = 0;\n";
            push @line_part_r, "(-mc[$mc_count[$nr_iterates]])*${var}[ind(".plus("l",$l).",".plus("m",$m).",1)]";
            push @line_part_i, "mc[$mc_count[$nr_iterates]]*${var}[ind(".plus("l",$l).",".plus("m",$m).",0)]";
            $mc_count[$nr_iterates]++;
          }
        }
        if ($#line_part_r >= 0) {
          $main_program_r .= $line_start.join("+",@line_part_r).")";
          $main_program_i .= $line_start.join("+",@line_part_i).")";
        } else {
          $main_program_r .= "${line_start}0)";
          $main_program_i .= "${line_start}0)";
        }
      }
      $main_program_r .= "$line\n";
      if ($line !~ /^\s*(int|REAL)/) {
        $main_program_i .= "$line\n";
      }
    }
    $main_program_r =~ s/\@c/0/g;
    $main_program_i =~ s/\@c/1/g;
    $main_program .= $main_program_r . $main_program_i;
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

for (my $i=0;$i<$nr_iterates;$i++) {
  print "static int first_mc_$i = 1;\n";
  print "static REAL *mult_constant_$i;\n";
  print "#define mc_count_$i $mc_count[$i]\n";
  print "static void initialize_mc_$i(param_list_t *param);\n";
  print "\n";
}

print "/* Declarations of threading go here. */\n\n";
print $main_program;

for (my $i=0;$i<$nr_iterates;$i++) {
  chomp $init_proc_assign[$i];
  print "\n";
  print "static void initialize_mc_$i(param_list_t *param) {\n";
  print "  int l,m;\n";
  print "  REAL ll,mm;\n";
  print "  REAL *mc;\n";
  print "\n";
  print "  first_mc_$i = 0;\n";
  print "  mult_constant_$i = (REAL*)malloc(sizeof(REAL)*mc_ind(param->max_order+2,0)*mc_count_$i);\n";
  print "  for (l=0;l<=param->max_order;l+=2) for (m=0;m<=l;m++) {\n";
  print "    ll = l;\n";
  print "    mm = m;\n";
  print "    mc = mult_constant_$i + mc_count_$i*mc_ind(l,m);\n";
  print "$init_proc_assign[$i]\n";
  print "  }\n";
  print "}\n";
}
