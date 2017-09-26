/* clipSearch head file */

#ifndef clipSearch_HEAD_H
#define clipSearch_HEAD_H

#define GAP_OPEN -2.0
#define SEED_GAP_OPEN -5.0
#define GAP_CONT -2.0
#define MATCH 2.0
#define SEED_MATCH 5.0
#define MISMATCH -2.0
#define SEED_MISMATCH -5.0
#define GU_MATCH 1.0

#define H 3

typedef struct parameterInfo
{
  double maxMFE;
  int verbose;
  int minScore;
} parameterInfo;

typedef struct alignInfo
{
  char *mirSeq;
  char *tarSeq;
  char *pairStr;
} alignInfo;

void scanMTI(parameterInfo *paraInfo, FILE *gfp, FILE *faifp, FILE *outfp, FILE *mirfp, FILE *peakfp);

int searchTarget(parameterInfo *paraInfo,  FILE *outfp, FILE *gfp, faidx *fai, CBed6 *bed6, map<string, string> &mirHash, int targetNum);

const char *filterMTI(const char *seq, char *structure, struct parameterInfo *paraInfo);

const char *searchMirSeed(char *seq, char *mir, int start);

int encodeIntChar (char ch);

int RNApair(char bp1, char bp2);

void readMiRNAs(FILE *fp, map<string, string> &readHash);

double NeedlemanWunschAlign(struct parameterInfo *paraInfo, FILE *outfp, char *mirseq, char *seq, alignInfo *align);

double scorePair(char a, char b, int i, int seedLen);

double scoreGap(int i, int seedLen);

void freeFloatMatrix(double **matrix);

void freeIntMatrix(int **matrix);

int argmax(double m[], int len);

double max(double m[], int len);

void freeAlignInfo(alignInfo *align);

int filterPairs(char *pairStr);


#endif /* End clipSearch_HEAD_H */
