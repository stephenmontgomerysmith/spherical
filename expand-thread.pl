#/usr/bin/perl -w

use strict;

my $main_text = "";
my $declarations;
$declarations .= "/* Converted into a threaded program by $0. */\n";
$declarations .= "#include <pthread.h>\n";
$declarations .= "extern int nr_threads;\n";
$declarations .= "static pthread_mutex_t job_m = PTHREAD_MUTEX_INITIALIZER;\n";
$declarations .= "static pthread_cond_t job_c = PTHREAD_COND_INITIALIZER;\n";
$declarations .= "struct job_info_s {int control, lstart, lend;};\n";

my $thread_nr = 0;
my @thread_text;

my $arg;
my $inside_function=0;
my $inside_threading=0;
my $struct;
my $memcpy;
my $memcpy2;

while(my $line = <>) {
  if ($line =~ /^void\s+\S+\s*\((.*?)\).*$/) {
    $arg = $1;
    $arg = "  $arg;\n" if $arg =~ s/\,/\;/;
    $inside_function = 1;
    $main_text .= $line;
  } elsif ($inside_function && $line =~ /^\}/) {
    $inside_function = 0;
    $main_text .= $line;
  } elsif ($line =~ /\/\* Start\-Threading \*\//) {
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
    $main_text .= "    pthread_t pid;\n";
    $main_text .= "    int i;\n";
    $main_text .= "    /* Start threads the first time around. */\n";
    $main_text .= "    if (first_thread_$thread_nr) {\n";
    $main_text .= "      first_thread_$thread_nr = 0;\n";
    $main_text .= "      job_info_$thread_nr = malloc(nr_threads*sizeof(struct job_info_s));\n";
    $main_text .= "      for (i=0;i<nr_threads;i++) {\n";
    $main_text .= "        job_info_${thread_nr}[i].control = 0;\n";
    $main_text .= "        /* The ranges for l are assigned to each thread.  Note that the work\n";
    $main_text .= "           required to calculate for l between lstart and lend is\n";
    $main_text .= "           proportional to lend^2 - lstart^2. */\n";
    $main_text .= "        if (i==0)\n";
    $main_text .= "          job_info_${thread_nr}[i].lstart = 0;\n";
    $main_text .= "        else\n";
    $main_text .= "          job_info_${thread_nr}[i].lstart = job_info_${thread_nr}[i-1].lend;\n";
    $main_text .= "        if (i==nr_threads-1)\n";
    $main_text .= "          job_info_${thread_nr}[i].lend=max_order+1;\n";
    $main_text .= "        else\n";
    $main_text .= "          job_info_${thread_nr}[i].lend = (double)(max_order+1)*sqrt((i+1.)/(double)nr_threads);\n";
    $main_text .= "        job_info_${thread_nr}[i].lend += (job_info_${thread_nr}[i].lend)%2; /* round up to next even number. */\n";
    $main_text .= "        pthread_create(&pid, NULL, thread_function_$thread_nr, job_info_$thread_nr+i);\n";
    $main_text .= "        pthread_detach(pid);\n";
    $main_text .= "      }\n";
    $main_text .= "    }\n";
    $main_text .= "    /* Copy data to threads. */\n";
    $main_text .= $memcpy;
    $main_text .= "    /* Start each thread, and then wait for each thread to finish. */\n";
    $main_text .= "    pthread_mutex_lock(&job_m);\n";
    $main_text .= "    for (i=0;i<nr_threads;i++)\n";
    $main_text .= "      job_info_${thread_nr}[i].control = 1;\n";
    $main_text .= "    pthread_mutex_unlock(&job_m);\n";
    $main_text .= "    pthread_cond_broadcast(&job_c);\n";
    $main_text .= "    pthread_mutex_lock(&job_m);\n";
    $main_text .= "    for (i=0;i<nr_threads;i++)\n";
    $main_text .= "      while (job_info_${thread_nr}[i].control != 0)\n";
    $main_text .= "        pthread_cond_wait(&job_c, &job_m);\n";
    $main_text .= "    pthread_mutex_unlock(&job_m);\n";
    $declarations .= "static struct {\n$struct} thread_pass_args_$thread_nr;\n";
    $declarations .= "static void *thread_function_$thread_nr(void *arg);\n";
    $declarations .= "static int first_thread_$thread_nr = 1;\n";
    $declarations .= "static struct job_info_s *job_info_$thread_nr;\n";
    $struct =~ s/^/  /mg;
    $thread_text[$thread_nr] = "$struct    /* Copy data to thread. */\n$memcpy2$thread_text[$thread_nr]";
    $thread_nr++;
  } else {
    if ($inside_threading) {
      $thread_text[$thread_nr] .= $line;
    } elsif ($inside_function) {
      $arg .= $line if $line =~ /^\s*(int|double|COMPLEX)/;
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
  print "\nstatic void * thread_function_$i(void *arg) {\n";
  print "  struct job_info_s* job_info = arg;\n";
  print "  while (1) {\n";
  print "    int lstart, lend;\n";
  print "    /* Wait until thread is told to start. */\n";
  print "    pthread_mutex_lock(&job_m);\n";
  print "    while (job_info->control != 1)\n";
  print "      pthread_cond_wait(&job_c, &job_m);\n";
  print "    lstart = job_info->lstart;\n";
  print "    lend = job_info->lend;\n";
  print "    pthread_mutex_unlock(&job_m);\n";
  print $thread_text[$i];
  print "    /* Broadcast that thread is finished. */\n";
  print "    pthread_mutex_lock(&job_m);\n";
  print "    job_info->control = 0;\n";
  print "    pthread_mutex_unlock(&job_m);\n";
  print "    pthread_cond_broadcast(&job_c);\n";
  print "  }\n  return(NULL);\n}\n";
}
