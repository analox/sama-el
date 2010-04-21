/* Anthony Padula
kriggui.cc 

This file runs a gui interface for the krigifier.  It fork/execs the wish interpreter.

 *  The author of this software is Anthony D. Padula 
 *  Permission to use, copy, modify, and distribute this software  
 *  for any purpose without fee is hereby granted, provided that   
 *  this entire notice is included in all copies of any software   
 *  which is or includes a copy or modification of this software   
 *  and in all copies of the supporting documentation for such     
 *  software.  THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT    
 *  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, THE AUTHOR    
 *  OFFERS NO REPRESENTATION OR WARRANTY OF ANY KIND                   
 *  CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS    
 *  FITNESS FOR ANY PARTICULAR PURPOSE.   

*/

// filename should be "*.???" ie a dot with three letters after it in DOS fashion
char fname[20] = "./kriggui.out";

// The full pathname of your postscript viewer
char viewer[20] = "/usr/bin/X11/gv";

bool redisplay = false;

#include <cstdio>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include "krigify.h"
#include <signal.h>

using namespace std;

int guipid= 0, gvid = -1;

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
       //       setlinebuf(*writepipe);
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
 
void sighandler(int i) {
  if( guipid > 0) {
    kill( 9, guipid);
    wait(NULL);
  }
  if( gvid > 0 ) {
    kill( 9, gvid);
    wait(NULL);
  }
  exit(0);
}

int main( int argc, char* argv[], char* envp[] ) {

  double tau2, beta0, alpha, theta, sigma2, low, up ;
  ofstream outf;
  long seed, p, n, childpid = 0;
  int stream = 42;
  FILE *read_from, *write_to, *writetitle, *readcrap;
  Vector<double> beta1, x0, ain(1), tin(1), lower,upper;
  Matrix<double> beta2;
  char * args[5];
  int len, replots = 0;
  int pipe1[2], pipe2[2];
  randfunc* MyRandomFunction = NULL;
  char pstr[10];

  seed = -1;    

  if( argc >= 3 ) { 
    sscanf( argv[1], "%ld", &seed);
    sscanf( argv[2], "%d", &stream);
  }

  signal(9, sighandler);

  // Spawn a child for the gui
  guipid = start_child("wish",&read_from,&write_to);
  fprintf(write_to, "source kriggui.tcl\n");
  fflush(write_to);

  // Wait for parameters from the gui
  fscanf(read_from, "%ld", &p);
  
  while (p != -1) {
    //$p $n $tau2 $beta0 $alpha $theta $sigma2 $lower $upper
    fscanf(read_from, "%ld %lf %lf %lf %lf %lf %lf %lf %ld %d", 
	   &n, &tau2, &beta0, &alpha, &theta, &sigma2, &low, &up, &seed, &stream);
    //fprintf(stderr, "%ld %f %f %f %f %f %f %f", 
    //	   n, tau2, beta0, alpha, theta, sigma2, low, up);
    //if(childpid != 0)
    //  kill(childpid, 9);

    /*    if( seed != -1 ) */
      MyRandomFunction = new randfunc( seed, // seed for the random number generator
				       stream);       // stream for the rng
      /*  else
	  MyRandomFunction = new randfunc; */

    beta1.newsize(p);
    beta2.newsize(p,p);
    x0.newsize(p);
    lower.newsize(p);
    upper.newsize(p);
   
    beta1 = 0;
    beta2 = 0;
    x0 = 0;
    lower = low;
    upper = up;
    ain = alpha;
    tin = theta;

    MyRandomFunction->setvalues( p, n, beta0, beta1,
				 beta2, x0, 
				 ain, tin, 
				 sigma2, lower,
				 upper);

    MyRandomFunction->generateTrend(tau2, &beta0);
    
    fprintf(write_to, "set fname {Picture: creating picture}\n");
    fprintf(write_to, "displayfname\n");
    fflush(write_to);

    outf.open(fname);
    scanfunction(*MyRandomFunction,  // scan this function, 
		 outf,              // output results to this stream,
		 100);              // divide each axis into this many points
    // and evaluate on the resulting grid,
    outf.close();

    seed = MyRandomFunction->getSeed();
    stream = MyRandomFunction->getStream();
    //    fprintf(write_to, "set seedstring {Your Seed: %ld Your Stream: %d}\n", seed, stream);
    //fprintf(write_to, "displayseed\n");
    fprintf(write_to, "set FinalValue(0) %ld\n", seed);
    fprintf(write_to, "set FinalValue(1) %d\n", stream);
    fflush(write_to);

    seed = -1;
    delete MyRandomFunction;
    MyRandomFunction = NULL;

    // Convert data file to ps
    
    sprintf(pstr, "%ld", p);
    len = strlen(fname);
    args[3] = pstr;
    //    printf("%s\n", pstr);
    args[4] = NULL;
    new_array(args[0], 15);
    args[1] = fname;
    new_array(args[2],len+1);
    strcpy(args[0], "./plotpoints");
    strcpy(args[2], fname);
    args[2][len-3] = 'p';
    args[2][len-2] = 's';
    args[2][len-1] = '\0';
    //printf(" ps file is %s\n", args[2]);

    if ((pipe(pipe1) < 0) || (pipe(pipe2) < 0) ) {
      perror("pipe"); exit(-1);
    }
    
    if ((childpid = fork()) < 0) {
       perror("fork"); exit(-1);
    } else if (childpid > 0) {   /* Parent. */
       close(pipe1[0]); close(pipe2[1]);
       /* Write to child on pipe1[1], read from child on pipe2[0]. */
       readcrap = fdopen(pipe2[0], "r");
       writetitle = fdopen(pipe1[1], "w");
       //       setlinebuf(writetitle);
       fprintf(writetitle, "%ld-D Krigifier Function: n=%ld, alpha=%d, theta=%d, sigma^2=%d\n",
	       p, n, (int)alpha,(int) theta,(int)sigma2);
       fflush(writetitle);
       wait(NULL);
       close(pipe1[1]);
       close(pipe2[0]);
    } else { /* Child. */
       close(pipe1[1]); close(pipe2[0]);
       /* Read from parent on pipe1[0], write to parent on pipe2[1]. */
       dup2(pipe1[0],0); dup2(pipe2[1],1);
       close(pipe1[0]); close(pipe2[1]);
       if (execve(args[0], args, envp) < 0)
	 perror("execve plotter");
       /* Never returns */
    }
    args[1] = NULL;

    if( replots > 0) {   
      if( redisplay == true) {
	fprintf(write_to, "set fname {Picture: stored in file %s}\n", args[2]);
	kill(gvid, SIGHUP);
      } else {
	fprintf(write_to, 
		"set fname {Picture: stored in file %s - push redisplay on your viewer}\n",
		args[2]);
	}
      fprintf(write_to, "displayfname\n");
      fflush(write_to);
    }    
    else {
      fprintf(write_to, "set fname {Picture: stored in file %s}\n", args[2]);
      fprintf(write_to, "displayfname\n");
      fflush(write_to);
      /* Try to call a gv window to display */
      strcpy(args[0], viewer);
      args[1] = args[2];
      args[2] = NULL;

      if((gvid = fork()) < 0) {
	  perror("fork"); exit(-1);
      } else if (gvid == 0) { /* Child. */ 
	  if(execve(args[0], args, envp) < 0)
	    perror("execve viewer");
	} 
    }
    fprintf(write_to, "set plotting 0\n");
    fprintf(write_to, 
	    "set status {Status: Done plotting.  You may now change parameters and push Plot}\n");
    fprintf(write_to, "displaystatus\n");
    fprintf(write_to, "noise\n");
    fflush(write_to);
    replots++;

    if( args[0] != NULL)
      delete args[0];
    if( args[1] != NULL)
      delete args[1];
    if( args[2] != NULL)
      delete args[2];

    args[0] = NULL;
    args[1] = NULL;
    args[2] = NULL;

    fscanf(read_from, "%ld", &p);
  }
  if( gvid > 0) {
    kill(gvid, 9);
    // close(pipe3[1]);

    wait(NULL);
  }

  return 0;
}

