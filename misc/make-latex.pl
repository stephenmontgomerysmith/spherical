#!/usr/bin/perl

# Creates a latex file of coefficients for a method.  Example of usage:
# make-latex "x*dy+y*dx" > eg.tex ; latex eg

$text = `echo "t $ARGV[0]" | perl expand-method.pl`;
$text =~ s/\> //g;
$text =~ s/\s*$//;
$text =~ s/\=/\&\=/g;
$text =~ s/\n/\\\\\n/g;

$fulltext = << 'EOL';
% Created using perl $0 "$ARGV[0]"

\documentclass{amsart}

\begin{document}

\allowdisplaybreaks[4]

Output of {\tt expand-method.pl} applied to {\tt $ARGV[0]}.

\begin{align*}
$text
\end{align*}

\end{document}
EOL

$fulltext =~ s/\$ARGV\[0\]/$ARGV[0]/g;
$fulltext =~ s/\$text/$text/;
$fulltext =~ s/\$0/$0/;

print $fulltext;
