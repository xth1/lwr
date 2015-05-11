#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <stddef.h>   /* for offsetof */
#include <unistd.h>
#include <string.h>

#include "mcstructs.h"
#include "uniinvdist.h"
#include "mcrdist.h"
#include "mcread_input.h"
#include "write_data.h"
#include "scenario.h"
#include "testrev.h"
#include "opt_scenario.h"
#include "unsigned.h"
#include "e_malloc.h"
#include "circ_align.h"
#include "texgraph.h"
#include "ext_function.h"
#include "countperms.h"
#include "unsignedhc.h"
#include "matrixmisc.h"

int main()
{
    printf("getFromGrimm\n");
    test_linkage();
}
