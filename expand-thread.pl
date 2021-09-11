#/usr/bin/perl -w

use strict;

my $main_text = "";
my $declarations;
$declarations .= "/* Converted into a threaded program by $0. */\n";
$declarations .= "\n";
$declarations .= "#include \"threads.h\"\n";
$declarations .= "\n";
$declarations .= "struct job_info_s {int lstart, lend; sem_t sem_start, sem_end;};\n";

my $thread_nr = 0;
my @thread_text;
my @thread_text_declarations;

my $arg;
my $inside_function=0;
my $inside_threading=0;
my $struct;
my $memcpy;
my $memcpy2;

while(my $line = <>) {
  if ($line =~ /^void\s+\S+\s*\((.*?)\).*$/) {
    $arg = $1;
    $arg = "  $arg;\n" if $arg =~ s/\,/\;/g;
    $inside_function = 1;
    $main_text .= $line;
  } elsif ($inside_function && $line =~ /^\}/) {
    $inside_function = 0;
    $main_text .= $line;
  } elsif ($line =~ /\/\* Start\-Threading \*\//) {
    $line = <>;
    $line = <>;
    $thread_text[$thread_nr] = "";
    $inside_threading = 1;
  } elsif ($line =~ /\/\* Stop\-Threading \*\//) {
    $inside_threading = 0;
    $struct = "";
    $memcpy = "";
    $memcpy2 = "";
    foreach my $a (split ';',$arg) {
      $a =~ s/^\s*(\S+)//;
      my $type = $1;
      if ($a =~ /\S/) {
        foreach my $b (split ',',$a) {
          $b =~ s/^\s*//;
          $b =~ s/\s*$//;
          my $c = $b;
          $c =~ s/\*//g;
          $c =~ s/\[.*?\]//g;
          if ($thread_text[$thread_nr] =~ /\b${c}\b/) {
            $struct .= "  $type $b;\n";
            $memcpy .= "    memcpy(&(thread_pass_args_$thread_nr.$c),&($c),sizeof($c));\n";
            $memcpy2 .= "    memcpy(&($c),&(thread_pass_args_$thread_nr.$c),sizeof($c));\n";
          }
        }
      }
    }
    $main_text .= "    int i;\n";
    $main_text .= "    /* Start threads the first time around. */\n";
    $main_text .= "    if (first_thread_$thread_nr) {\n";
    $main_text .= "      first_thread_$thread_nr = 0;\n";
    $main_text .= "      job_info_$thread_nr = (struct job_info_s*)malloc(param->nr_threads*sizeof(struct job_info_s));\n";
    $main_text .= "      for (i=0;i<param->nr_threads;i++) {\n";
    $main_text .= "        sem_init(job_info_${thread_nr}[i].sem_start);\n";
    $main_text .= "        sem_init(job_info_${thread_nr}[i].sem_end);\n";
    $main_text .= "        /* The ranges for l are assigned to each thread.  Note that the work\n";
    $main_text .= "           required to calculate for l between lstart and lend is\n";
    $main_text .= "           proportional to lend^2 - lstart^2. */\n";
    $main_text .= "        if (i==0)\n";
    $main_text .= "          job_info_${thread_nr}[i].lstart = 0;\n";
    $main_text .= "        else\n";
    $main_text .= "          job_info_${thread_nr}[i].lstart = job_info_${thread_nr}[i-1].lend;\n";
    $main_text .= "        if (i==param->nr_threads-1)\n";
    $main_text .= "          job_info_${thread_nr}[i].lend=param->max_order+1;\n";
    $main_text .= "        else\n";
    $main_text .= "          job_info_${thread_nr}[i].lend = (REAL)(param->max_order+1)*sqrt((i+1.)/(REAL)param->nr_threads);\n";
    $main_text .= "        job_info_${thread_nr}[i].lend += (job_info_${thread_nr}[i].lend)%2; /* round up to next even number. */\n";
    $main_text .= "        thread_create(thread_function_$thread_nr, job_info_$thread_nr+i);\n";
    $main_text .= "      }\n";
    $main_text .= "    }\n";
    $main_text .= "    /* Copy data to threads. */\n";
    $main_text .= $memcpy;
    $main_text .= "    /* Start the threads, and then wait for each thread to finish. */\n";
    $main_text .= "    for (i=0;i<param->nr_threads;i++)\n";
    $main_text .= "      sem_post(job_info_${thread_nr}[i].sem_start);\n";
    $main_text .= "    for (i=0;i<param->nr_threads;i++)\n";
    $main_text .= "      sem_wait(job_info_${thread_nr}[i].sem_end);\n";
    $declarations .= "static struct {\n$struct} thread_pass_args_$thread_nr;\n";
    $declarations .= "static thread_return_t thread_function_$thread_nr(void *arg);\n";
    $declarations .= "static int first_thread_$thread_nr = 1;\n";
    $declarations .= "static struct job_info_s *job_info_$thread_nr;\n";
    $struct =~ s/^/  /mg;
    $thread_text_declarations[$thread_nr] = "$struct$thread_text_declarations[$thread_nr]";
    $thread_text[$thread_nr] = "    /* Copy data to thread. */\n$memcpy2$thread_text[$thread_nr]";
    $thread_nr++;
  } else {
    if ($inside_threading) {
      $thread_text[$thread_nr] .= $line;
    } elsif ($inside_function) {
      $arg .= $line if $line =~ /^\s*(int|REAL|param\_list\_t)/;
      $main_text .= $line;
    } else {
      $main_text .= $line;
    }
  }
}

chomp $declarations;
$main_text =~ s/\/\* Declarations of threading go here\. \*\//$declarations/;
print $main_text;

for (my $i=0;$i<$thread_nr;$i++) {
  print "\n";
  print "static thread_return_t thread_function_$i(void *arg) {\n";
  print "  struct job_info_s* job_info = (struct job_info_s*)arg;\n";
  print "  while (1) {\n";
  print $thread_text_declarations[$i];
  print "    int lstart, lend;\n";
  print "    int l,m;\n";
  print "    /* Wait until thread is told to start. */\n";
  print "    sem_wait(job_info->sem_start);\n";
  print "    lstart = job_info->lstart;\n";
  print "    lend = job_info->lend;\n";
  print $thread_text[$i];
  print "    /* Broadcast that thread is finished. */\n";
  print "    sem_post(job_info->sem_end);\n";
  print "  }\n  return(0);\n}\n";
}
