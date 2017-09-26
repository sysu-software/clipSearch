/* API for bed format */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<assert.h>
#include<math.h>
#include<time.h>
#include<limits.h>
extern "C" {
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
}
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

template<typename T>
string NumberToString(T Number)
{
  ostringstream ss;
  ss << Number;
  return ss.str();
}

void scanMTI(parameterInfo *paraInfo, FILE *gfp, FILE *faifp, FILE *outfp, FILE *mirfp, FILE *peakfp)
{
  int targetNum = 0;
  char *line = NULL;
  CBed6 *bed6 = NULL;
  faidxMap faiHash;
  map<string, string> mirHash;
  if(paraInfo->verbose) fprintf(stderr, "read fai file\n");
  readFai(faifp, faiHash);
  if(paraInfo->verbose) fprintf(stderr, "read mir file\n");
  readMiRNAs(mirfp, mirHash);
  while (line = getLine(peakfp))
  {
    if (feof(peakfp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    bed6 = parseBed6Line(line); // get bed information
    safeFree(line); // free memory

    string chrom(bed6->chrom);
    if (faiHash.find(chrom) == faiHash.end())
    {
      fprintf(stderr, "can't not find the chromosome %s, skip it.\n", bed6->chrom);
      freeBed6Item(bed6);
      continue;
    }
    faidx *fai = faiHash[chrom];
    targetNum = searchTarget(paraInfo, outfp, gfp, fai, bed6, mirHash, targetNum);
    freeBed6Item(bed6);
  }// loop-while
}

int searchTarget(parameterInfo *paraInfo,  FILE *outfp, FILE *gfp, faidx *fai, CBed6 *bed6, map<string, string> &mirHash, int targetNum)
{
  int seedLen = 7;
  int i = 0;
  int chromStart = 0;
  int chromEnd = 0;
  int seedStart = 0;
  int seedEnd = 0;
  char *seq = faidxFetchSeq(gfp, fai, bed6->chromStart, bed6->chromEnd, bed6->strand);
  int seqLen = strlen(seq);
  int trueSeedLen = seedLen;
  for (i = seedLen; i < seqLen - 1; i++)
  {
    for (map<string, string>::iterator curr = mirHash.begin(); curr != mirHash.end();
        curr++)
    {
      char *mirSeq = (char *)curr->second.c_str();
      const char *seedType = searchMirSeed(seq, mirSeq, i);
      if (strcmp(seedType, "non") != 0)
      {
        string foldSeq = mirSeq;
        int mirLen = strlen(mirSeq);
        double mfe = 0;
        if (bed6->strand == '+')
        {
          chromEnd = bed6->chromStart + i + 2;
          chromStart = chromEnd - mirLen;
          if (chromStart < 0 ) chromStart = 0;
        }
        else
        {
          chromStart = bed6->chromEnd - i - 2;
          chromEnd = chromStart + mirLen;
          if (chromStart < 0 ) chromStart = 0;
        }
        char *targetSeq = (char *)faidxFetchSeq(gfp, fai, chromStart, chromEnd, bed6->strand);
        foldSeq = foldSeq + "XXXXXXXXX" + targetSeq;
        const char *chimeraSeq = foldSeq.c_str();
        char *chimeraStruct = (char *)safeMalloc(strlen(chimeraSeq) + 1);
        mfe = fold(chimeraSeq, chimeraStruct);
        if (filterMTI(chimeraSeq, chimeraStruct, paraInfo) != NULL)
        {
          alignInfo *align = (alignInfo *)safeMalloc(sizeof(alignInfo));
          double alignScore = NeedlemanWunschAlign(paraInfo, outfp, mirSeq, targetSeq, align);
          if (mfe < paraInfo->maxMFE && alignScore > paraInfo->minScore && filterPairs(align->pairStr) != 0)
          {
            fprintf(outfp, ">%s\t%d\t%d\t%s:%s\t%.2f\t%c\t%s\t%.2f\t%.2f\n", bed6->chrom, chromStart, chromEnd,
                curr->first.c_str(), bed6->name, bed6->score, bed6->strand, seedType, mfe, alignScore);
            fprintf(outfp, "miRNA  3'-%s-5'\n          %s\ntarget 5'-%s-3'\n", align->mirSeq, align->pairStr, align->tarSeq);
          }
          freeAlignInfo(align);
          targetNum++;
        } // if filter mti
        free_arrays;
        safeFree(targetSeq);
        safeFree(chimeraStruct);
      }// if seed type
    } // for mirHash
  } // for seed
  return targetNum;
}

int filterPairs(char *pairStr)
{
  int pairLen = strlen(pairStr);
  int i = 0;
  int seedLen = 7;
  int pairNum = 1;
  int start = pairLen - seedLen;
  if (start < 0) start = 0;
  for(i = pairLen - seedLen; i < pairLen - 1; i++)
  {
    if (pairStr[i] == '.' || pairStr[i] == '-' || pairStr[i] == ':')
    {
      pairNum = 0;
      break;
    }
  }
  return pairNum;
}

const char *filterMTI(const char *seq, char *structure, struct parameterInfo *paraInfo)
{
  const char *seedType[6] = {"8mer", "7mer-m8", "7mer-A1", "6mer", "offset-6mer", "non-canonical"};
  int i = 0;
  int pairNum = 0;
  int guNum = 0;
  int seedLen = 8;
  short int *pairTable = make_pair_table(structure);
  int bulgeNum = 0;
  char c = '.';
  int j = 0;

  for(i = 1; i < seedLen && i < strlen(structure); i++)
  {
    j = i + 1;
    if(structure[i] == '(')
      pairNum++;
    if(structure[i] == ':')
      guNum++;
    if((structure[i] == '.' || structure[i] == ':' || structure[i] == ')') && i < seedLen - 1)
    {
      bulgeNum++;
    }
    if(structure[i] == '.' && pairTable[j - 1] > 0 && pairTable[j + 1] > 0 && (pairTable[j - 1] - pairTable[j + 1]) > 1 && i < seedLen - 2)
    {
      bulgeNum += pairTable[j - 1] - pairTable[j + 1] - 2;
    }
    if(pairTable[j] > 0 && pairTable[j + 1] > 0 && (pairTable[j] - pairTable[j + 1]) > 1 && i < seedLen - 2)
      bulgeNum += pairTable[j] - pairTable[j + 1] - 1;
    if(structure[i] == '(' && i == 1)
    {
      c = seq[pairTable[j]];
    }
  }
  safeFree(pairTable);
  if (pairNum == 7 && c == 'A' && bulgeNum == 0)
    return seedType[0];
  else if (pairNum == 7 && c != 'A' && bulgeNum == 0)
    return seedType[1];
  else if (pairNum == 6 && c == 'A' && bulgeNum == 0)
    return seedType[2];
  else if (pairNum == 6 && c != 'A' && bulgeNum == 0)
    return seedType[3];
  /*else if (pairNum == 6 && bulgeNum==1 && structure[1]=='.')
    return seedType[4];
  else if (pairNum >= 6 && totalPairs>=paraInfo->minPair && bulgeNum==1)
    return seedType[5];*/

  return NULL;
}

const char *searchMirSeed(char *seq, char *mir, int start)
{
  const char *seedType[5] = {"8mer", "7mer-m8", "7mer-A1", "6mer", "non"};
  int mirLen = strlen(mir);
  int seqLen = strlen(seq);
  int seedLen = 7;
  int i, j;
  int pairNum = 0;
  int tag = 0;
  char c;
  for(i = start, j = 1; i > 0 && j < seedLen + 1; i--, j++)
  {
    if (RNApair(seq[i], mir[j]) == 2)
    {
      pairNum += 1;
    }
    else
    {
      break;
    }
  }
  c = seq[start + 1];
  if (pairNum == 7 && c == 'A')
    return seedType[0];
  else if (pairNum == 7 && c != 'A')
    return seedType[1];
  else if (pairNum == 6 && c == 'A')
    return seedType[2];
  else if (pairNum == 6 && c != 'A')
    return seedType[3];
  else
    return seedType[4];
}

double NeedlemanWunschAlign(struct parameterInfo *paraInfo, FILE *outfp, char *mirseq, char *seq, alignInfo *align)
{
  double alignScore = 0;
  double **scoreMatrix;
  double tmp[3];
  int i, j, k;
  int seedLen = 8;
  int mirseqLen = strlen(mirseq);
  int seqLen = strlen(seq);
  int M = mirseqLen + 1;
  int N = seqLen + 1;
  char *newMirSeq = strClone(mirseq);
  reverseBytes(newMirSeq, mirseqLen);
  seedLen = mirseqLen - seedLen;
  char *seq1 = newMirSeq;
  char *seq2 = seq;
  char *qSeq = NULL;
  char *mSeq = NULL;
  char *tSeq = NULL;
  int allocSpace = 0;
  // allocate memory
  scoreMatrix = (double **)safeMalloc(sizeof(double *)*M);
  scoreMatrix[0] = (double *)safeMalloc(sizeof(double) * (M * N));
  for (i = 1; i < M; i++)
    scoreMatrix[i] = scoreMatrix[0] + N * i;

  // initialize scoreMatrix
  scoreMatrix[0][0] = 0;
  for (i = 1; i < N; i++)
    scoreMatrix[0][i] = scoreMatrix[0][i - 1] + scoreGap(i, seedLen);
  for (i = 1; i < M; i++)
    scoreMatrix[i][0] = scoreMatrix[i - 1][0] + scoreGap(i, seedLen);

  // matrix score
  j = 0;
  for (i = 1; i < M; i++)
  {
    for (j = 1; j < N; j++)
    {
      tmp[0] = scoreMatrix[i - 1][j - 1] + scorePair(seq1[i - 1], seq2[j - 1], i, seedLen); /* sequence index is i-1 and j-1*/
      tmp[1] = scoreMatrix[i - 1][j] + scoreGap(i, seedLen);
      tmp[2] = scoreMatrix[i][j - 1] + scoreGap(i, seedLen);
      scoreMatrix[i][j] = max(tmp, 3);
    }
  }

  // trace back
  allocSpace = M > N ? M : N;
  qSeq = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));
  mSeq = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));
  tSeq = (char *)safeMalloc(sizeof(char) * (allocSpace + 1));

  i = M - 1;
  j = N - 1;

  k = 0;
  for (;;)
  {
    if (i == 0 && j == 0) break;
    if (k >= allocSpace)
    {
      qSeq = (char *)xrealloc(qSeq, k + 2);
      mSeq = (char *)xrealloc(mSeq, k + 2);
      tSeq = (char *)xrealloc(tSeq, k + 2);
      allocSpace = k;
    }
    if (i >= 1 && j >= 1 && scoreMatrix[i][j] == scoreMatrix[i - 1][j - 1] + scorePair(seq1[i - 1], seq2[j - 1], i, seedLen))
    {
      qSeq[k] = seq1[i - 1];
      alignScore += scorePair(seq1[i - 1], seq2[j - 1], i, seedLen);
      if (RNApair(seq1[i - 1], seq2[j - 1]) == 2)
      {
        mSeq[k] = '|';
      }
      else if (RNApair(seq1[i - 1], seq2[j - 1]) == 1)
      {
        mSeq[k] = ':';
      }
      else
      {
        mSeq[k] = '.';
      }
      tSeq[k] = seq2[j - 1];
      i--;
      j--;
    }
    else if (i >= 1 && j >= 0 && scoreMatrix[i][j] == scoreMatrix[i - 1][j] + scoreGap(i, seedLen))
    {
      qSeq[k] = seq1[i - 1];
      mSeq[k] = '-';
      tSeq[k] = '-';
      alignScore += scoreGap(i, seedLen);
      i--;
    }
    else if (j >= 1 && i >= 0 && scoreMatrix[i][j] == scoreMatrix[i][j - 1] + scoreGap(i, seedLen))
    {
      qSeq[k] = '-';
      mSeq[k] = '-';
      tSeq[k] = seq2[j - 1];
      alignScore += scoreGap(i, seedLen);
      j--;
    }
    k++;
  }

  qSeq[k] = '\0';
  mSeq[k] = '\0';
  tSeq[k] = '\0';
  reverseBytes(qSeq, strlen(qSeq));
  reverseBytes(mSeq, strlen(mSeq));
  reverseBytes(tSeq, strlen(tSeq));

  //fprintf(outfp, "%s\n", qSeq);
  //fprintf(outfp, "%s\n", mSeq);
  //fprintf(outfp, "%s\n", tSeq);

  freeFloatMatrix(scoreMatrix);
  align->mirSeq = qSeq;
  align->tarSeq = tSeq;
  align->pairStr = mSeq;
  safeFree(newMirSeq);
  return alignScore;
}

void freeAlignInfo(alignInfo *align)
{
  safeFree(align->mirSeq);
  safeFree(align->tarSeq);
  safeFree(align->pairStr);
  safeFree(align);
}

int encodeIntChar (char ch)
{
  ch = toupper(ch);
  if (ch == 'A')
    return 0;
  else if (ch == 'C')
    return 1;
  else if (ch == 'G')
    return 2;
  else if (ch == 'T' || ch == 'U')
    return 3;
  else
    return 4;
}

int RNApair(char bp1, char bp2)
{
  int pairMatrix[5][5] =
  {
    /* A C G T N*/
    {0, 0, 0, 2, 0},
    {0, 0, 2, 0, 0},
    {0, 2, 0, 1, 0},
    {2, 0, 1, 0, 0},
    {0, 0, 0, 0, 0}
  };
  return pairMatrix[(int)(encodeIntChar(bp1))][(int)(encodeIntChar(bp2))];
}

void readMiRNAs(FILE *fp, map<string, string> &readHash)
{
  int fieldNum = 0;
  char **fields = NULL;
  char *line = NULL;
  char *seq = NULL;
  char *qualityName = NULL;
  char *quality = NULL;
  while (line = getLine(fp))
  {
    if (feof(fp) || line == NULL)
    {
      safeFree(line);
      break;
    }
    if (line[0] != '@' && line[0] != '>')
    {
      fprintf(stderr, "error read format: %c\n", line[0]);
      safeFree(line);
      continue;
    }
    fields = splitWhitespace(line + 1, &fieldNum);
    if (line[0] == '@')
    {
      seq = getLine(fp);
      qualityName = getLine(fp);
      quality = getLine(fp);
      string name(fields[0]);
      string sequence(seq);
      readHash[name] = sequence;
      safeFree(qualityName);
      safeFree(quality);
    }
    if (line[0] == '>')
    {
      seq = getLine(fp);
      string name(fields[0]);
      string sequence(seq);
      readHash[name] = sequence;
    }
    freeWords(fields, fieldNum);
    safeFree(line);
    safeFree(seq);
  } // gofile while loops
}

double scoreGap(int i, int seedLen)
{
  if (i > seedLen)
    return SEED_GAP_OPEN;
  else return GAP_OPEN;
  return GAP_OPEN;
}

double scorePair(char a, char b, int i, int seedLen)
/* score match and mismatch */
{
  int matchNum = RNApair(a, b);
  if (matchNum == 2)
  {
    if (i > seedLen)
      return SEED_MATCH;
    else return MATCH;
  }
  else if (matchNum == 1)
  {
    return GU_MATCH;
  }
  else
  {
    if (i > seedLen)
      return SEED_MISMATCH;
    else return MISMATCH;
  }
}

void freeFloatMatrix(double **matrix)
/* free float matrix */
{
  free(matrix[0]);
  free(matrix);
}

void freeIntMatrix(int **matrix)
/*free int matrix */
{
  free(matrix[0]);
  free(matrix);
}

double max(double m[], int len)
/* return maxinum value */
{
  int i = 0;
  double maxScore = m[0];

  for (i = 1; i < len; i++)
  {
    if (m[i] > maxScore)
    {
      maxScore = m[i];
    }
  }
  return maxScore;
}
