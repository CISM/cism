#include "cfortran.h"
#include "writestats.h"
#include "config.inc"

#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <sys/times.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>

FCALLSCSUB3(gc_writestats,GF_WRITESTATS,gf_writestats,STRING,STRING,DOUBLE)

#define CFG_LEN 35
#define BUFFER_LEN 400
#define PERM_FILE (S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)

void gc_writestats(const char *resname, const char *cfgname, double wallTime)
{
  struct tms runtime; 
  clock_t clck;
  double userTime,sysTime;
  time_t now;
  struct tm * timeAndDate;
  char dateStr[20];
  struct utsname unameData;
  char outBuffer[BUFFER_LEN+1];
  char hdrBuffer[BUFFER_LEN+1];
  int fd;
  int i,haveLock;
  struct flock fl = { F_WRLCK, SEEK_SET, 0,       0,     0 };
  off_t fileLength;


  /* get user and system time */
  clck = times(&runtime);
  userTime = ((double) runtime.tms_utime)/((double) sysconf(_SC_CLK_TCK));
  sysTime = ((double) runtime.tms_stime)/((double) sysconf(_SC_CLK_TCK));

  /* get the current data */
  now = time(NULL);
  timeAndDate = localtime(&now);
  snprintf(dateStr,20,"%4d-%02d-%02d_%02d:%02d",timeAndDate->tm_year+1900, timeAndDate->tm_mon+1, timeAndDate->tm_mday, timeAndDate->tm_hour, timeAndDate->tm_min);

  /* get host name and architecture */
  if ((uname(&unameData))!=0) {
    unameData.nodename[0]='\0';
    unameData.machine[0]='\0';
  }

  /* construct output line */
  snprintf(outBuffer,BUFFER_LEN,"%*s %9.2f %9.2f %8.2f %s %-10s %-6s %-10s \"%s\"\n",-CFG_LEN,cfgname, wallTime, userTime, sysTime, dateStr, \
	   unameData.nodename, unameData.machine, VERSION, GLIMMER_FCFLAGS);
  snprintf(hdrBuffer,BUFFER_LEN,"%*s %9s %9s %-8s %-16s %-10s %-6s %-10s %s\n",-CFG_LEN,"#cfg_file","wall_time","usr_time","sys_time","date","host","arch","version","FCFLAGS");

  /* open output file */
  if ((fd = open(resname, O_CREAT|O_WRONLY|O_SYNC,PERM_FILE)) == -1) {
    perror("opening result file");
    printf("%s\n",outBuffer);
    return;
  }

  /* get a lock on the file */
  i=0;
  while ((haveLock=fcntl(fd, F_SETLK, &fl))==-1 && i<100000)  i++;
  if (haveLock==-1) {
    close(fd);
    perror("getting lock");
    printf("%s\n",outBuffer);
    return;
  }

  /* go to the end of the file */
  fileLength = lseek(fd,0,SEEK_END);

  /* write data */
  if (fileLength == 0)
    write(fd,hdrBuffer,strlen(hdrBuffer));
  write(fd,outBuffer,strlen(outBuffer));

  /* release the lock */
  fl.l_type = F_UNLCK;
  if (fcntl(fd, F_SETLK, &fl) == -1) {
    perror("unlocking file");
    return;
  }
  
  close(fd);

  return;
}
