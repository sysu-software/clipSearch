#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<getopt.h>
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
using namespace std;

#include "bioUtils.h"
#include "bedFile.h"
#include "faiFile.h"
#include "clipSearch.h"

char version[] = "clipSearch version 0.1";
void usage(void);

int main(int argc, char *argv[])
{
  char *outfile       = NULL;
  FILE *outfp         = NULL;
  FILE *genomefp      = NULL;
  FILE *faifp         = NULL;
  FILE *mirfp         = NULL;
  FILE *peakfp        = NULL;
  int showVersion     = 0;
  int showHelp        = 0;
  int i               = 0;
  int c               = 0;
  struct parameterInfo paraInfo;
  /* parse commmand line parameters */

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVbtPo:m:s:" ;

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "output" , required_argument , NULL, 'o' },
    { "max-mfe" , required_argument, NULL, 'm' },
    { "min-score" , required_argument, NULL, 's' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.verbose   = 0;
  paraInfo.maxMFE    = 0;
  paraInfo.minScore  = 0;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'o':
      outfile  = optarg;
      break;
    case 'm':
      paraInfo.maxMFE = atof(optarg);
      break;
    case 's':
      paraInfo.minScore = atoi(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  if (argc == optind || argc - 4 != optind) usage();
  genomefp = (FILE *) fopen(argv[optind], "r");
  if (genomefp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open genome %s\n", argv[optind]);
    usage();
  }
  faifp = (FILE *) fopen(argv[optind + 1], "r");
  if (faifp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open fai file %s\n", argv[optind + 1]);
    usage();
  }
  mirfp = (FILE *) fopen(argv[optind + 2], "r");
  if (mirfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open miR file %s\n", argv[optind + 2]);
    usage();
  }
  peakfp = (FILE *) fopen(argv[optind + 3], "r");
  if (peakfp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open peak file %s\n", argv[optind + 3]);
    usage();
  }


  if (outfile == NULL)
  {
    outfp = stdout;
  }
  else
  {
    outfp = (FILE *) fopen(outfile, "w");
    if (outfp == NULL)
    {
      fprintf(stderr, "ERROR: Can't open %s\n", outfile);
      usage();
    }
  }

  // help for version
  if (showVersion)
  {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }
  fprintf(stderr, "#clipSearch program start\n");
  scanMTI(&paraInfo, genomefp, faifp, outfp, mirfp, peakfp);
  fprintf(stderr, "#clipSearch program end\n");
  fclose(genomefp);
  fclose(faifp);
  fclose(mirfp);
  fclose(peakfp);
  fclose(outfp);
  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s", "Usage:  clipSearch [options] <genome file> <genome fai> <mir file, fasta> <peak, bed>\n\
peak format is bed6+\n\
[options]\n\
-v/--verbose                   : verbose information\n\
-V/--version                   : clipSearch version\n\
-h/--help                      : help informations\n\
-o/--output <string>           : output file\n\
-m/--max-mfe <double>          : maximum mfe for miR-target duplex [default < 0]\n\
-s/--min-score <int>           : minimum score for miR-target duplex [default > 0]\n\
");

  exit(1);
}
