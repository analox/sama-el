/* plotpoints.c

   This is a simple driver routine for generating plots of 3d data.  
   It takes two parameters, the first being the filename where the points
   are stored, and the second being the name of the destination file.
   The program will prompt the user for a title to add to the plot.
*/

#include <strings.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/wait.h>

/* Exec the named cmd as a child process, returning
 * two pipes to communicate with the process, and
 * the child's process ID */
int start_child(char *cmd, FILE **readpipe, FILE **writepipe) {
    int childpid, pipe1[2], pipe2[2];
    
    if ((pipe(pipe1) < 0) || (pipe(pipe2) < 0) ) {
       perror("pipe"); exit(-1);
    }

    if ((childpid = fork()) < 0) {
       perror("fork"); exit(-1);
    } else if (childpid > 0) {   /* Parent. */
       close(pipe1[0]); close(pipe2[1]);
       /* Write to child on pipe1[1], read from child on pipe2[0]. */
       *readpipe = fdopen(pipe2[0], "r");
       *writepipe = fdopen(pipe1[1], "w");
       setlinebuf(*writepipe);
       return childpid;

    } else { /* Child. */
       close(pipe1[1]); close(pipe2[0]);
       /* Read from parent on pipe1[0], write to parent on pipe2[1]. */
       dup2(pipe1[0],0); dup2(pipe2[1],1);
       close(pipe1[0]); close(pipe2[1]);

       if (execlp(cmd, cmd, NULL) < 0)
          perror("execlp");
       /* Never returns */
    }
    return childpid;
}


int main( int argc, char* argv[]) 
     //int argc; char* argv[];
{
  int childpid, i, e = 1;
  FILE *read_from, *write_to;
  char title[200];
  long p = 2;

  if( argc < 3 ){
    fprintf( stderr, "Wrong number of arguments.  Usage: plotpoints [-e] data_file_name output_file_name [dimension]\n");
  }

    read_from = NULL;
    write_to = NULL;

    childpid = start_child("gnuplot",&read_from,&write_to);
    fprintf(write_to, "set parametric\n");
    fprintf(write_to, "set data style lines\n");
    if( argv[1][0] == '-' ) {
      fprintf(write_to, "set terminal postscript color eps\n");
      e = 2;
    } else
      fprintf(write_to, "set terminal postscript color\n");
    fprintf(write_to, "set nokey\n");
    fprintf(write_to, "set output \"%s\"\n", argv[e+1]);
    printf( "Title of plot: ");
    fflush(stdout);
    fgets(title, 200, stdin);
    i = 0;
    while((title[i] != '\n')&&(i < 200)) i++;
    title[i] = '\0';
    fprintf(write_to, "set title \"%s\"\n", title);
    if( argc == e+3)
      p = strtol( argv[e+2], NULL, 10 );
    /*    printf("p = %ld\n", p);*/
    if( p == 2 ) {
      fprintf(write_to, "set contour\n");

      fprintf(write_to, "splot \"%s\"\n", argv[e]);
    }
    else if (p == 1)
      fprintf(write_to, "plot \"%s\"\n", argv[e]);
    fprintf(write_to, "quit\n");

    wait(NULL);
    fclose(read_from); fclose(write_to); 
    exit(0);
       
}







