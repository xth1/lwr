/* mcmain.c
 *    Parse command line options and dispatch required functions.
 *
 * Copyright (C) 2001-2006 The Regents of the University of California
 * by Glenn Tesler
 *
 * Contains code from main.c in GRAPPA 1.02
 * Copyright (C) 2000-2001  The University of New Mexico and
 *                          The University of Texas at Austin
 * by David A. Bader, Bernard M.E. Moret, Tandy Warnow, Stacia K Wyman, Mi Yan
 *
 * See file COPYRIGHT for details.
 *****************************************************************************
 * This file is part of GRIMM.
 *
 * GRIMM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, Version 2,
 * dated June 1991, as published by the Free Software Foundation.
 *
 * GRIMM is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/* Last modified on Wed Aug 2, 2006, by Glenn Tesler
 */

/* Contains excerpts from GRAPPA 1.02: main.c */

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

/* put quotes around version number */
#define qstr1(s) #s
#define qstr(s) qstr1(s)
#define QVERS qstr(VERS)

FILE *outfile;

void do_fn_distmatrix(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
                      int circular,
                      int unsigned_dist,
                      int verbose,
                      struct genome_struct *genome_list);
void do_fn_distmatrix_v(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
                        int circular,
                        int unsigned_dist,
                        struct genome_struct *genome_list);



void print_usage(char *progname)
{
    fprintf(stderr,"\nGRIMM %s\n",QVERS);
    fprintf(stderr,"Copyright (C) 2001-2006 The Regents of the University of California\n");
    fprintf(stderr,"Contains code from GRAPPA 1.02 (C) 2000-2001 The University of New Mexico\n");
    fprintf(stderr,"                                         and The University of Texas at Austin\n");
    fprintf(stderr,"See file COPYRIGHT for full copyright and authorship details.\n");

    fprintf(stderr,
            "\nUsage: %s -f datafile [-o outfile] [other options]",
            progname);

    fprintf(stderr,
            "\n\n");

    fprintf(stderr,"   -f datafile: name of file containing all the genomes\n");
    fprintf(stderr,"   -o outfile: output file name.  Output goes to STDOUT unless -o used.\n");
    fprintf(stderr,"   -v: Verbose output.\n");
    fprintf(stderr,"       In pairwise mode, outputs all the graph parameters.\n");
    fprintf(stderr,"       In matrix mode, outputs matrices for all graph parameters.\n");
    fprintf(stderr,"       In special functions, outputs extra information.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"Genome type:\n");
    fprintf(stderr,"   -L: unichromosomal (directed) linear\n");
    fprintf(stderr,"   -C: unichromosomal circular\n");
    fprintf(stderr,"   Multichromosomal (undirected linear chromosomes) used unless -L/-C given\n");
    /* -L: "0" chromosomes -- linear -- the old unichromosomal case */
    /* -C: "0" chromosomes -- circular */
    /* -D2: pad it w/caps, but then do unichromosomal alg only on it;
       This has no mathematical usefulness whatsoever, it is just
       useful for debugging. */
    /* -M #: use user specified caps -- also just useful for debugging */
    fprintf(stderr,"\n");
    fprintf(stderr,"Genome selection:\n");
    fprintf(stderr,"   -g i,j: compare genomes i and j (counting from 1)\n");
    fprintf(stderr,"      Rearrangement scenario will go from genome i to genome j.\n");
    fprintf(stderr,"      In Hannenhalli-Pevzner papers, genome i is gamma and genome j is pi.\n");
    fprintf(stderr,"      For files with 3 or more genomes, unless -g is specified, a distance\n");
    fprintf(stderr,"      matrix is printed and options -d, -c, -z, -s are not available.\n");
    fprintf(stderr,"   -m: force matrix output even for 2 genomes\n");

    fprintf(stderr,"\n");

    fprintf(stderr,"Function: (default -d -c -s for multichromosomal, -d -s for unichromosomal)\n");
    fprintf(stderr,"   -d: distance\n");
    fprintf(stderr,"   -c: show capping\n");

    /* TODO: choose something better than -z */
    fprintf(stderr,"   -z: show capping and chromosome delimeters too\n");


    fprintf(stderr,"   -s: display an optimal scenario (see -S below)\n");
    fprintf(stderr,"   -u: distance of unsigned genomes; find signs giving optimal distance\n");
    fprintf(stderr,"   -U #: unsigned genomes approximation algorithm with # iterations\n");
    fprintf(stderr,"   -W file: optional file with weight matrix for -U\n");
    fprintf(stderr,"   -S #: produce/display optimal scenario in a certain format (#=1,...,7)\n");
    fprintf(stderr,"      #=1: operations may occur in any order\n");
    fprintf(stderr,"      #=2,...,7: operations are prioritized\n");
    fprintf(stderr,"        rev. short->long; flip chromos; fusion; transloc short->long; fission\n");
    fprintf(stderr,"      #=1,2,3,4: scenario respects initial capping\n");
    fprintf(stderr,"      #=5,6,7: scenario may not respect initial capping, so can avoid flips\n");
    fprintf(stderr,"      #=2: 1-line permutations with caps\n");
    fprintf(stderr,"      #=3: multiline permutations, show caps, / at end of chromo\n");
    fprintf(stderr,"      #=4: multiline permutations, hide caps, $ at end of chromo\n");
    fprintf(stderr,"      #=5,6,7: different greedy approaches to reversals first\n");
    fprintf(stderr,"        They allow recapping, so hide caps and show $ at end of each chromo\n");

    fprintf(stderr,"\nSpecial functions:\n");
    fprintf(stderr,"   -t#: test all 1-step reversals, list the ones reducing distance by 1,\n");
    fprintf(stderr,"       and give statistics\n");
    fprintf(stderr,"       #=1,2,3,4 specifies test method (debugging)\n");
    fprintf(stderr,"   -F n: tabulate statistics for given # genes n\n");
    fprintf(stderr,"   -T xxx: output graph in TeX form.  xxx is a combination of\n");
    fprintf(stderr,"      l: left-to-right format  c: component/cycle format  (only use one)\n");
    fprintf(stderr,"      b: bracket chromosomes\n");
    fprintf(stderr,"      s: signed labels\n");
    fprintf(stderr,"      u: unsigned labels 2x-1,2x  a: unsigned labels a,b  (only use one)\n");
    fprintf(stderr,"      r: list of raw components, cycles, paths, not TeX format\n");
    fprintf(stderr,"      1-6: stage in multichromosomal capping algorithm\n");
    fprintf(stderr,"      Some argument must be given, so use 1 if nothing else\n");
    fprintf(stderr,"\nDebugging (not for general use):\n");
    fprintf(stderr,"   -D1: delete caps from multichromosomal genomes, treat as linear\n");
    fprintf(stderr,"   -D2: put caps on multichromosomal genomes, then treat them as unichromosomal\n");
    fprintf(stderr,"   -M n: the data appears in unichromosomal format, but has been manually\n");
    fprintf(stderr,"       capped for n chromosomes\n");
    fprintf(stderr,"   -X arg: demo of extensions\n");

    fprintf(stderr,"\n");
    fprintf(stderr,"Format of datafile: each entry has the form\n");
    fprintf(stderr,"   >human\n");
    fprintf(stderr,"   # Chromosome 1\n");
    fprintf(stderr,"   5 -2 3 -1 $\n");
    fprintf(stderr,"   # Chromosome 2\n");
    fprintf(stderr,"   7 -4 6 $\n");
    fprintf(stderr,"   # Chromosome 3\n");
    fprintf(stderr,"   8\n");
    fprintf(stderr,"   ...\n");
    fprintf(stderr,"   >mouse\n");
    fprintf(stderr,"   ...\n");
    fprintf(stderr,"Each genome starts with a line \">genomename\".\n");
    fprintf(stderr,"The genes are separated by any whitespace, and may be on multiple lines.\n");
    fprintf(stderr,"The chromosomes are separated by \";\" or \"$\".  New lines between\n");
    fprintf(stderr,"chromosomes are recommended for neatness but are not required.\n");
    fprintf(stderr,"Comments begin with \"#\" and go to the end of the line.\n");
    fprintf(stderr,"\nWith the -U option only, strips of unknown orientation may be bracketed:\n");
    fprintf(stderr,"   5 [ -2 3 ] -1 $\n");
    fprintf(stderr,"means it's either -2 3 or -3 2 but the unbracketed genes are fixed.\n");
    fprintf(stderr,"If no [] with the -U option then unknown signs on all genes in genomes 2,3,...\n");
    fprintf(stderr,"Brackets {} () are reserved and currently ignored.\n");
    fprintf(stderr,"\n");
}


/****************************************************************************/
/*            Distance matrix for multiple genomes                          */
/****************************************************************************/


/* Stacia added this so that the ouput can be directly fed into other
 * programs such as the tds suite. It expects # of taxa on the first line
 * by itself, and then a full symmetric matrix with the taxa name as the
 * first thing on a line.
 */
#define GNAME_WIDTH 20
void print_full_distmatrix(int **distmatrix,int num_genomes,
                           struct genome_struct *genome_list)
{
    int i, j;     /* loop over matrix entries */
    int k;        /* loop to print genome name */
    int pad;

    int max_entry = 0;
    //int entry_width;
    for (i=0; i<num_genomes; i++) {
        for (j=0; j<i; j++) {
            if (max_entry > distmatrix[j][i]) {
                max_entry = distmatrix[j][i];
            }
        }
    }
    //entry_width = num_digits(max_entry);

    //fprintf(outfile,"\n   %d\n",num_genomes);

    for (i=0 ; i<num_genomes ; i++) {

        /* print genome name, truncating or padding if necessary */
        for (k=0, pad=FALSE; k<GNAME_WIDTH; k++) {
            if (!pad) {
                if (genome_list[i].gnamePtr[k] != '\0')
                    fputc(genome_list[i].gnamePtr[k], outfile);
                else
                    pad = TRUE;
            }
            if (pad) fputc(' ', outfile);
        }
        fputc('\t', outfile);
    }
    //fprintf(outfile,"\n");

    return;
}


/* print a matrix of the values of a particular field in the distance
 * structures */
void print_matrix_offset(char *matrix_name,
                         int field_offset,
                         graphstats_t **graphstats_matrix,
                         int num_genomes,
                         struct genome_struct *genome_list)
{
    int i, j;     /* loop over matrix entries */
    int k;        /* loop to print genome name */
    int pad;

    int max_entry = 0;
    int e;
    //int entry_width;
    for (i=0; i<num_genomes; i++) {
        for (j=0; j<i; j++) {
            e = *(int *)((char *)(&graphstats_matrix[i][j]) + field_offset);
            if (max_entry > e) {
                max_entry = e;
            }
        }
    }
    //entry_width = num_digits(max_entry);


    //fprintf(outfile,"\n");
    //fprintf(outfile,"%s Matrix:\n", matrix_name);
    //fprintf(outfile,"   %d\n",num_genomes);

    for (i=0 ; i<num_genomes ; i++) {

        /* print genome name, truncating or padding if necessary */
        for (k=0, pad=FALSE; k<GNAME_WIDTH; k++) {
            if (!pad) {
                if (genome_list[i].gnamePtr[k] != '\0')
                    fputc(genome_list[i].gnamePtr[k], outfile);
                else
                    pad = TRUE;
            }
            if (pad) fputc(' ', outfile);
        }
        fputc('\t', outfile);

        for (j=0 ; j<num_genomes ; j++) {
            //fprintf(outfile, "%*d ", entry_width,
//	      *(int *)((char *)(&graphstats_matrix[i][j]) + field_offset));
        }
        //fprintf(outfile,"\n");
    }
    //fprintf(outfile,"\n");

    return;
}




/*****************************************************************************/

/* circular option: the genomes have already been aligned so that
   linear signed inv dist = circular signed inv dist
   However, unsigned still needs to do more work */

void do_fn_distmatrix(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
                      int circular,
                      int unsigned_dist,
                      int verbose,
                      struct genome_struct *genome_list)
{
    int i;
    int **distmatrix;
    distmem_t distmem;
    int dist_error;

    if (verbose && !unsigned_dist) {
        do_fn_distmatrix_v(NUM_GENES, NUM_CHROMOSOMES, NUM_GENOMES,
                           circular, unsigned_dist,
                           genome_list);
        return;
    }


    //fprintf(outfile,"\nDistance Matrix:\n");

    /* allocate space for matrix */
    distmatrix = (int **) e_malloc((NUM_GENOMES)*sizeof(int *), "distmatrix");

    for (i=0; i < NUM_GENOMES; i++) {
        distmatrix[i] = (int *) e_malloc((NUM_GENOMES)*sizeof(int), "distmatrix");
    }

    /* allocate memory for distance computations */
    mcdist_allocmem(NUM_GENES,NUM_CHROMOSOMES,&distmem);


    /* compute the appropriate matrix */
    if (unsigned_dist) {
        set_unsigneddist_matrix(distmatrix, genome_list,
                                NUM_GENES, NUM_CHROMOSOMES, NUM_GENOMES,
                                circular,
                                &distmem,
                                &dist_error);
        if (dist_error) {
            //fprintf(outfile,"WARNING: These unsigned permutations are too complex;\n");
            //fprintf(outfile,"some answers are approximated.\n");
        }

    } else if (NUM_CHROMOSOMES == 0) {
        /* unichromosomal data */
        setinvmatrix(distmatrix,genome_list,
                     NUM_GENES,NUM_GENOMES,
                     &distmem,

                     /* CIRCULAR already converted to equiv linear problem */
                     FALSE);
    } else {
        /* multichromosomal data */
        setmcdistmatrix(distmatrix,genome_list,
                        NUM_GENES,NUM_CHROMOSOMES,NUM_GENOMES,
                        &distmem);
    }

    /* display the matrix */
    print_full_distmatrix(distmatrix,NUM_GENOMES,genome_list);

    /* free the memory for distance computations */
    mcdist_freemem(&distmem);

    /* free memory for the distance matrix */
    for (i=NUM_GENOMES-1; i >=0 ; i--)
        free(distmatrix[i]);
    free(distmatrix);
}


/* verbose matrix option: display matrices of all the various parameters */
/* TODO: unsigned_dist not supported at this time */

#define print_matrix_field(display,var) \
    print_matrix_offset(display, \
		    offsetof(graphstats_t, var), \
		    statsmatrix, \
		    NUM_GENOMES, \
		    genome_list)



void do_fn_distmatrix_v(int NUM_GENES, int NUM_CHROMOSOMES, int NUM_GENOMES,
                        int circular,
                        int unsigned_dist,
                        struct genome_struct *genome_list)
{
    int i,j;
    graphstats_t **statsmatrix;
    distmem_t distmem;


    /* allocate space for matrix */
    statsmatrix =
        (graphstats_t **) e_malloc((NUM_GENOMES)*sizeof(graphstats_t *),
                                   "statsmatrix");

    for (i=0; i < NUM_GENOMES; i++) {
        statsmatrix[i] = (graphstats_t *) e_malloc((NUM_GENOMES)*sizeof(graphstats_t),
                         "statsmatrix");
    }

    /* allocate memory for distance computations */
    mcdist_allocmem(NUM_GENES,NUM_CHROMOSOMES,&distmem);


    /* compute all pairwise graphs.
     * Some parameters aren't symmetric, so do all pairs (i,j) and (j,i).
     * TODO: compute diagonal ones simply
     */

    for (i = 0 ;  i < NUM_GENOMES ;  i++) {
        for (j = 0 ;  j < NUM_GENOMES ;  j++) {
            mcdist_noncircular(&genome_list[i], &genome_list[j],
                               NUM_GENES, NUM_CHROMOSOMES,
                               &distmem,
                               &statsmatrix[i][j]);
        }
    }


    /* display the matrices */
    if (NUM_CHROMOSOMES == 0) {
        print_matrix_field("Distance", d);
        print_matrix_field("Number of Breakpoints", br);
        print_matrix_field("Number of Cycles", c4);
        print_matrix_field("Number of Hurdles", h);
        print_matrix_field("Number of Fortresses", f);
    } else {
        print_matrix_field("Distance", d);
        print_matrix_field("Number of Black Edges", bl);
        print_matrix_field("Number of Cycles and Paths", cp);
        print_matrix_field("Number of Gamma-Gamma Paths", pgg);
        print_matrix_field("Number of Semi-knots", s);
#if 0
        print_matrix_field("Parameter rr", rr);
#endif
        print_matrix_field("Parameter r", r);
        print_matrix_field("Parameter fr", fr);
        print_matrix_field("Parameter gr", gr);

#if 0
        /* can't do this w/o forming capping/concat, which we didn't do */
        print_matrix_field("Number of bad bonds", badbonds);
#endif

        print_matrix_field("Number of internal breakpoints", bp_int);
        print_matrix_field("Number of external breakpoints", bp_ext);
    }



    /* free the memory for distance computations */
    mcdist_freemem(&distmem);

    /* free memory for the distance matrix */
    for (i=NUM_GENOMES-1; i >=0 ; i--)
        free(statsmatrix[i]);
    free(statsmatrix);
}
