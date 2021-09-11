#include "spherical.h"

static REAL evaluate_string_recursive(char **s, char var, REAL x, int current_binding, int start_of_expression);

REAL evaluate_string(char *s, char var, REAL x) {
  return evaluate_string_recursive(&s, var, x, 0, 1);
}

/* binding
     m M: 10
     + -: 20
     * /: 30
     ^  : 40
*/
/* Only do operations with binding greater than current_binding */
static REAL evaluate_string_recursive(char **s, char var, REAL x, int current_binding, int start_of_expression) {
  REAL number1, number2;
  char *new_s;

//printf("s = %s b = %d\n",*s,current_binding);
    while (isspace(**s)) (*s)++;
    if (**s==0) {
      fprintf(stderr, "Syntax error\n");
      exit(1);
    }

    /* Get first number */
    number1 = strtod(*s,&new_s);
    if (new_s!=*s && **s!='-') {
      *s = new_s;
    }
    else if (**s==var) {
      number1 = x;
      (*s)++;
    }
    else if (**s=='-') { /* unary minus */
      if (!start_of_expression) {
        fprintf(stderr, "Syntax error\n");
        exit(1);
      }
      (*s)++;
      number1 = - evaluate_string_recursive(s, var, x, 20, 0);
    }
    else if (**s=='(') {
      (*s)++;
      number1 = evaluate_string_recursive(s, var, x, 0, 1);
      while (isspace(**s)) (*s)++;
      if (**s != ')') {
        fprintf(stderr, "Syntax error\n");
        exit(1);
      }
      (*s)++;
    }
    else {
      fprintf(stderr, "Syntax error\n");
      exit(1);
    }

  while (1) {
    while (isspace(**s)) (*s)++;
//printf("s2 = '%s' b = %d n = %g\n",*s,current_binding,number1);
    /* Get operator */
    if (**s==0 || **s==')') {
      return number1;
    }
    else if (**s=='M') {
      if (current_binding>=10) return number1;
      (*s)++;
      number2 = evaluate_string_recursive(s, var, x, 10, 0);
      if (number2>number1) number1 = number2;
    }
    else if (**s=='m') {
      if (current_binding>=10) return number1;
      (*s)++;
      number2 = evaluate_string_recursive(s, var, x, 10, 0);
      if (number2<number1) number1 = number2;
    }
    else if (**s=='+') {
      if (current_binding>=20) return number1;
      (*s)++;
      number2 = evaluate_string_recursive(s, var, x, 20, 0);
      number1 = number1 + number2;
    }
    else if (**s=='-') {
      if (current_binding>=20) return number1;
      (*s)++;
      number2 = evaluate_string_recursive(s, var, x, 20, 0);
      number1 = number1 - number2;
    }
    else if (**s=='*') {
      if (current_binding>=30) return number1;
      (*s)++;
      number2 = evaluate_string_recursive(s, var, x, 30, 0);
      number1 = number1 * number2;
    }
    else if (**s=='/') {
      if (current_binding>=30) return number1;
      (*s)++;
      number2 = evaluate_string_recursive(s, var, x, 30, 0);
      number1 = number1 / number2;
    }
    else if (**s=='^') {
      if (current_binding>=40) return number1;
      (*s)++;
      number2 = evaluate_string_recursive(s, var, x, 39, 0);
      number1 = pow(number1,number2);
    }
    else {
      fprintf(stderr,"Syntax error\n");
      exit(1);
    }
  }
}
