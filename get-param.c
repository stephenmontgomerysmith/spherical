#include "spherical.h"

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

int param_bool(const char *p) {
  int r;
  r = strtol(get_line(p,0),NULL,10);
  if (verbose && r) printf("%s=%d\n",p,r);
  return r;
}

int param_int(const char *p) {
  int r;
  r = strtol(get_line(p,1),NULL,10);
  if (verbose) printf("%s=%d\n",p,r);
  return r;
}

REAL param_REAL(const char *p) {
  REAL r;
  r = strtod(get_line(p,1),NULL);
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

void param_ignore(const char *p) {
  char *s;
  s = get_line(p,0);
}
