#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <kwset.h>


int main(int argc, char *argv[]) {
  
  int exit_code = 0;
  size_t blah = 0;

  size_t where = 10;
  size_t res = 0;
  struct kwsmatch kwsm;


  char *grpw;
  grpw = malloc(4);
  strcpy(grpw, "www");
  printf("%s\n", grpw);

  char *grpa;
  grpa = malloc(4);
  strcpy(grpa, "aaa");
  printf("%s\n", grpa);
  
  char *data = malloc(11);
  strcpy(data, "xxxaaaxbbb");
  printf("%s\n", data);

  kwset_t test = kwsalloc(NULL);
  kwsincr(test, grpw, 3);
  kwsincr(test, grpa, 3);
  kwsprep(test);
  
  res = kwsexec(test, data, where, &kwsm);
  printf("%d\n", (int)where);
  printf("%d\n", (int)res);
  printf("%d\n", (int)kwsm.index);
  printf("%d\n", (int)kwsm.offset[0]);
  printf("%d\n", (int)kwsm.size[0]);
  
  kwsfree (test);
  return exit_code;
  
}
