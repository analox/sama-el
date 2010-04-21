#include "utilities.h"
long long getTime()
{
   struct timeval tv;
   struct timezone tz;
   gettimeofday(&tv, &tz);

   long long result = 1000000LL * tv.tv_sec + tv.tv_usec;

   return result;
}
