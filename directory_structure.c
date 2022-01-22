#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> // for mkdir, for example
#include <errno.h>
#include <string.h>
#include <fcntl.h>

// # This should be run first, to make sure the required packages are installed


char* makeNewDirectory(const char *dirname) {
  // Note that mkdir returns -1 when the error occurs and sets errno accordingly
  errno = 0;
  int success = mkdir(dirname); //S_IRWXU
  
  char* msg = malloc(150);
  if (success == -1) {
    if (errno == EEXIST) {
      strcpy(msg, "\nDirectory ");
      strcat(msg, dirname);
      strcat(msg, "' exists already");
      // printf("Directory %s exists already\n", dirname);
    }else if (errno == ENOENT) {
      strcpy(msg, "\nParent directory of '");
      strcat(msg, dirname);
      strcat(msg, "' does not exist");
    }else {
      strcpy(msg, "\nDid not make '");
      strcat(msg, dirname);
      strcat(msg, "' directory successfully");
    }
  }else {
    strcpy(msg, "\nDirectory ");
    strcat(msg, dirname);
    strcat(msg, "' created successfully");
  }
  
  return msg;
}

// S_IRWXU mode bits, which implies that the owner will have read/write/execute permissions on the directory
int main(void) {
  
  FILE *flogs;
  flogs = fopen("log_directory_structure.txt", "wb"); // move to logs/ at end
  
  char bar[96] = "||----------------------------------";
  strcat(bar, " Directory structure ");
  strcat(bar, "----------------------------------||");
  fputs(bar, flogs);
  
  fputs(makeNewDirectory("logs"), flogs);
  fputs(makeNewDirectory("inputs"), flogs);
  fputs(makeNewDirectory("inputs/processed"), flogs);

  fputs(makeNewDirectory("results"), flogs);

  fputs(makeNewDirectory("intermediate_data"), flogs);

  fputs(makeNewDirectory("report_specific"), flogs);
  fputs(makeNewDirectory("report_specific/epitables"), flogs);
  fputs(makeNewDirectory("report_specific/clustering"), flogs);
  fputs(makeNewDirectory("report_specific/heatmaps"), flogs);
  fputs(makeNewDirectory("report_specific/densityplots"), flogs);
  fputs(makeNewDirectory("report_specific/frequencies"), flogs);
  
  fclose(flogs);
  
  return 1;
}

// run using: $ mingw32-gcc.exe directory_structure.c -o directory_structure.exe
// option -Wall enables all compiler’s warning messages
// option -o is used to specify the output file name. 
//    If we do not use this option, then an output file with name a.out is generated.

// Refs: 
//  https://www.delftstack.com/howto/c/mkdir-in-c/
