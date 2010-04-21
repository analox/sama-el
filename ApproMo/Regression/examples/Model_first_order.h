#include "Array.h"
#include "nrutil.h"
#include "stdio.h"
#include "ludcmp.h"
#include "lubksb.h"
#include "sqr.h"
#include "math.h"

int Model_first_order(Array<double> InputData,
		      Array<double> TargetData,
		      int start,
		      int ende,
		      Array<double> &coefficient);
