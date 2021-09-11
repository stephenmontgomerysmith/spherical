#include "spherical.h"

void get_parameters(int argc, const char **argv, param_list_t *param) {
  int i;
  int param_filename_set = 0;

  param->verbose = 0;
  param->verbose_print = 0;
  param->outfilename = (char*)malloc(strlen("s.out")+1);
  strcpy((char*)param->outfilename,"s.out");

  for (i=1;i<argc; i++) {
    if (strcmp(argv[i],"-p")==0 && i+1<argc) {
      i++;
      set_param_filename(argv[i]);
      param_filename_set = 1;
    }
    else if (strcmp(argv[i],"-o")==0 && i+1<argc) {
      i++;
      param->outfilename = argv[i];
    }
    else if (strcmp(argv[i],"-v")==0) {
      param->verbose = 1;
      if (i+1 < argc) {
        param->verbose_print = strtol(argv[i+1],NULL,10);
        if (param->verbose_print>=1) i++;
        if (param->verbose_print<=0) param->verbose_print = 0;
      }
    }
    else if (strncmp(argv[i],"-v",2)==0) {
      param->verbose = 1;
      param->verbose_print = strtol(argv[i]+2,NULL,10);
      if (param->verbose_print>=1) i++;
      if (param->verbose_print<=0) param->verbose_print = 0;
    }
    else {
      param_filename_set = 0;
      break;
    }
  }
  if (!param_filename_set) {
    printf("Usage: %s -p <parameters-file> [-o <output-file>] [-v [<number>]]\n"
           "  if -o is not set, <output-file> defaults to \"s.out\".\n"
           "  -v gives verbose output.\n"
           "  <number> is how often to print output to terminal.\n", argv[0]);
    exit(0);
  }
  set_param_verbose_level(param->verbose);

  if (getenv("NR_THREADS")==NULL)
    param->nr_threads = param_int("nr_threads");
  else {
    param->nr_threads = strtol(getenv("NR_THREADS"),NULL,10);
    if (param->verbose) printf("nr_threads=%d\n",param->nr_threads);
  }
  param_ignore("nr_threads");

  param->max_order = param_int("max_order");
  param->length = ind(param->max_order+2+PADDING,-PADDING,0);

  param->tstart = param_REAL("tstart");
  param->tend = param_REAL("tend");

  if (check_param("print_every"))
    param->print_every = param_int("print_every");
  else
    param->print_every = 1;

  param->ode_adams_bash_2 = param_bool("ode_adams_bash_2");
  param->ode_adams_bash_4 = param_bool("ode_adams_bash_4");
  param->ode_rk_4 = param_bool("ode_rk_4");
  param->ode_rkf_45 = param_bool("ode_rkf_45");
  if (check_param("h"))
    param->h = param_REAL("h");
  else
    param->h = 1e-3;
  if (check_param("tol"))
    param->tol = param_REAL("tol");
  else
    param->tol = 1e-3;

  param->do_vl = param_bool("do_vl");
  if (param->do_vl) {
    param->vl_fun = param_string("vl_fun");
  } else
    param->lambda = param_REAL("lambda");
  param->w[0] = param_REAL("w1");
  param->w[1] = param_REAL("w2");
  param->w[2] = param_REAL("w3");
  param->gamm[0*3+0] = param_REAL("gamma11");
  param->gamm[0*3+1] = param->gamm[1*3+0] = param_REAL("gamma12");
  param->gamm[0*3+2] = param->gamm[2*3+0] = param_REAL("gamma13");
  param->gamm[1*3+1] = param_REAL("gamma22");
  param->gamm[1*3+2] = param->gamm[2*3+1] = param_REAL("gamma23");
  param->gamm[2*3+2] = param_REAL("gamma33");
  param->do_koch = param_bool("do_koch");
  param->do_dd = param_bool("do_dd");
  param->do_vd = param_bool("do_vd");
  param->do_ard = param_bool("do_ard");
  i = 0;
  if (param->do_koch) i++;
  if (param->do_dd) i++;
  if (param->do_vd) i++;
  if (param->do_ard) i++;
  if (i>1) {
    fprintf(stderr,"Error: cannot use more than one type of diffusion.\n");
    exit(1);
  }
  if (param->do_koch) {
    param->C1 = param_REAL("C1");
    param->C2 = param_REAL("C2");
    param->C3 = param_REAL("C3");
    param->CI = 0;
  } else if (param->do_dd) {
    param->C1 = param_REAL("C1");
    param->C2 = param_REAL("C2");
    param->CI = 0;
  } else if (param->do_vd) {
    param->vd_fun = param_string("vd_fun");
    param->CI = 0;
  } else if (param->do_ard) {
    param->b1 = param_REAL("b1");
    param->b2 = param_REAL("b2");
    param->b3 = param_REAL("b3");
    param->b4 = param_REAL("b4");
    param->b5 = param_REAL("b5");
    param->CI = 0;
  } else
    param->CI = param_REAL("CI");
  param->do_rsc = param_bool("do_rsc");
  if (param->do_rsc) param->kappa = param_REAL("kappa");

  param->print_aij = param_bool("print_aij");
  param->print_aijkl = param_bool("print_aijkl");
  param->print_aijklmn = param_bool("print_aijklmn");
  param->print_daij = param_bool("print_daij");
  param->print_daijkl = param_bool("print_daijkl");
  param->print_daijklmn = param_bool("print_daijklmn");

  done_with_param();
}

static char zero[] = "0";
static int verbose=0;
static struct param_s {
  char *line;
  int used;
  struct param_s *next;
} *param = NULL;

void set_param_filename(const char *filename) {
  FILE *par_file;
  char s[1024];
  struct param_s *current=NULL;
  par_file = fopen(filename,"r");
  if (par_file==NULL) {
    perror("Error opening parameter file");
    exit(1);
  }
  while (fgets(s,sizeof(s)-1,par_file)!=NULL) {
    if (s[0] != '#') {
      if (current==NULL) {
        param = (struct param_s *)malloc(sizeof(struct param_s));
        current = param;
      } else {
        current->next = (struct param_s *)malloc(sizeof(struct param_s));
        current = current->next;
      }
      current->line = (char *)malloc(strlen(s)+1);
      strcpy(current->line,s);
      current->used = 0;
      current->next = NULL;
    }
  }
  fclose(par_file);
}

void set_param_verbose_level(int v) {
  verbose = v;
}

void done_with_param() {
  struct param_s *current=param, *next;
  int error=0;
  while (current!=NULL) {
    if (current->used==0) {
      fprintf(stderr,"Did not understand the parameter: %s",current->line);
      error = 1;
    }
    next = current->next;
    free(current);
    current = next;
  }
  if (error) exit(1);
}

static char *get_line(const char *p, int die_if_none) {
  struct param_s *current=param;
  while (current!=NULL) {
    if (strncmp(p,current->line,strlen(p))==0 && (current->line)[strlen(p)]=='=') {
      current->used = 1;
      return current->line+strlen(p)+1;
    }
    current = current->next;
  }
  if (die_if_none) {
    fprintf(stderr,"Parameter file is missing %s\n",p);
    exit(1);
  } else
    return zero;
}

int check_param(const char *p) {
  struct param_s *current=param;
  while (current!=NULL) {
    if (strncmp(p,current->line,strlen(p))==0 && (current->line)[strlen(p)]=='=') {
      return 1;
    }
    current = current->next;
  }
  return 0;
}


int param_bool(const char *p) {
  int r;
  char *p1,*p2;
  p1 = get_line(p,0);
  r = strtol(p1,&p2,10);
  if (p1==p2) {
    fprintf(stderr,"Syntax error with parameter %s.\n",p);
    exit(1);
  }
  if (verbose && r) printf("%s=%d\n",p,r);
  return r;
}

int param_int(const char *p) {
  int r;
  char *p1,*p2;
  p1 = get_line(p,1);
  r = strtol(p1,&p2,10);
  if (p1==p2) {
    fprintf(stderr,"Syntax error with parameter %s.\n",p);
    exit(1);
  }
  if (verbose) printf("%s=%d\n",p,r);
  return r;
}

REAL param_REAL(const char *p) {
  REAL r;
  char *p1,*p2;
  p1 = get_line(p,1);
  r = strtod(p1,&p2);
  if (p1==p2) {
    fprintf(stderr,"Syntax error with parameter %s.\n",p);
    exit(1);
  }
  if (verbose) printf("%s=%g\n",p,r);
  return r;
}

int param_choice(const char *p, ...) {
  va_list args;
  char *s,*t;
  int r;
  s = get_line(p,0);
  while (strlen(s)>0 && isspace(s[strlen(s)-1])) s[strlen(s)-1] = '\0';
  va_start(args,p);
  while (strcmp(t=va_arg(args,char*),"")!=0) {
    r = va_arg(args,int);
    if (strcmp(s,t)==0) {
      va_end(args);
      if (verbose) printf("%s=%s\n",p,s);
      return r;
    }
  }
  fprintf(stderr,"Parameter %s not set to acceptable value.\n",p);
  exit(1);
}

char *param_string(const char *p) {
  char *s, *ret_val;
  s = get_line(p,1);
  ret_val = (char *)malloc(strlen(s)+1);
  strcpy(ret_val,s);
  return ret_val;
}

void param_ignore(const char *p) {
  char *s;
  s = get_line(p,0);
}
