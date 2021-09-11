# This takes a parameters file as input, and replaces the gamma and w 
# variables with a shear moving in a random direction.  The output goes
# into a file with "-skew" placed before the ".txt" in its name.

# The purpose is to provide test scenarios in which gamma and w have entries
# that are all non-zero.


#!/usr/bin/perl

$param_file = $ARGV[0];
$new_param_file = $param_file;
$new_param_file =~ s/\.txt/\-skew\.txt/ || die
open(IN,$param_file) || die;
open(OUT,"> $new_param_file") || die;

$octave_commands = <<EOL;
format long e
Du = [0 1 0;0 0 0;0 0 0];
Q = orth(rand(3,3));
Du = Q*Du*Q';
gamma = Du+Du'
vort = Du-Du';
w = [vort(3,2) -vort(3,1) vort(2,1)]
EOL

$octave_commands =~ s/\:/\n/sg;
$numbers = `echo "$octave_commands" | perl -pe 's/\:/\n/g' | octave`;
$numbers =~ s/.*gamma\s\=//s;
$numbers =~ s/w\s\=//;
print $numbers;
$numbers =~ s/^\s*//;
@numbers = split '\s+',$numbers;

while (<IN>) {
  if (/gamma11/) {print OUT "gamma11=$numbers[0]\n"}
  elsif (/gamma12/) {print OUT "gamma12=$numbers[1]\n"}
  elsif (/gamma13/) {print OUT "gamma13=$numbers[2]\n"}
  elsif (/gamma22/) {print OUT "gamma22=$numbers[4]\n"}
  elsif (/gamma23/) {print OUT "gamma23=$numbers[5]\n"}
  elsif (/gamma33/) {print OUT "gamma33=$numbers[8]\n"}
  elsif (/w1/) {print OUT "w1=$numbers[9]\n"}
  elsif (/w2/) {print OUT "w2=$numbers[10]\n"}
  elsif (/w3/) {print OUT "w3=$numbers[11]\n"}
  else {print OUT $_}
}
