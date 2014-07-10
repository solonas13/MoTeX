/******************************************************************************/
/* SMILE v1.47 - Extraction of structured motifs common to several sequences  */
/* Copyright (C) 2004 L.Marsan (lama -AT- prism.uvsq.fr)                      */
/*                                                                            */
/* This program is free software; you can redistribute it and/or              */
/* modify it under the terms of the GNU General Public License                */
/* as published by the Free Software Foundation; either version 2             */
/* of the License, or (at your option) any later version.                     */
/*                                                                            */
/* This program is distributed in the hope that it will be useful,            */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with this program; if not, write to the Free Software                */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */
/******************************************************************************/

/* (C) by "coward" <coward -AT- ii.uib.no> from 1996-1999 */
/*   modified and extended by Laurent Marsan <lama@prism.uvsq.fr> 1999-2004.  */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <shufflet.h>



/******************************************************/
/* indexseq: convert it to a sequence of indices      */
/* input:    seq: sequence of upper-case letters      */
/*        seqlen: sequence length                     */
/* output:   seq: sequence of indices                 */ 

void indexseq(signed char seq[], int seqlen, int letter[128])
{
  int i;

  /* second pass: convert letters to indices */
  for (i = 0; i < seqlen; i++)
	{
	seq[i]=letter[(int)seq[i]];
	if (seq[i] == -1)
		{
		fprintf(stderr,"unknown character\n");
		exit(1);
		}
	}
}

/******************************************************/
/* kletcount: count k-lets in sequence of indices     */
/* input:    seq: sequence of indices (0,1,2,..,m-1)  */
/*        seqlen: sequence length                     */
/*             k: length of the k-lets                */
/*             m: alphabet size = highest index + 1   */
/* output: count: alphabetically ordered k-let counts */
/*                sufficient allocation assumed       */

void kletcount(signed char seq[], int seqlen, int k, int m, int count[])
{
  int i, hash = 0, chop;

  chop = pow(m,k-1); /* number of (k-1)-lets, used to chop off most significant
                        hash digit */
  for (hash = 0; hash < chop*m; count[hash++] = 0); /* initialize count */

  /* first k-1 letters: begin constructing hash, no k-lets to store yet */
  for (i = 0; i < k-1; i++)
    {
      assert(seq[i] >= 0 && seq[i] < m); /* if not: non-alphabetic index */
      hash = hash*m+seq[i];
    }   
  /* letter k-1 to end of sequence: update hash, store k-lets */
  for (i = k-1; i < seqlen; i++)
    {
      assert(seq[i] >= 0 && seq[i] < m); /* if not: non-alphabetic index */
      hash = (hash%chop)*m+seq[i]; /* chop%, shift*, add+ -> new hash */  
      count[hash]++;
    }
}

/******************************************************/
/* kletverify: check that k-let count in (shuffled)   */
/* sequence corresponds to stored count.              */
/* input:    seq: sequence of indices (0,1,2,..,m-1)  */
/*        seqlen: sequence length                     */
/*             k: length of the k-lets                */
/*             m: alphabet size = highest index + 1   */
/*        count0: reference k-let counts              */
/* output:count1: array for new k-let counts          */
/*                sufficient allocation assumed       */
/* return value: 1 of OK, 0 if mismatch               */

int kletverify(signed char seq[], int seqlen, int k, int m,
		int count0[], int count1[])
{
  int hash, nklets, errflag = 0;
  char klet[MAXORDER+1];

  nklets = pow(m,k);
  kletcount(seq,seqlen,k,m,count1);
  for (hash = 0; hash < nklets; hash++)
    if (count1[hash] != count0[hash])
      {
	fprintf(stderr,"Mismatch: %s: %i in orig, %i in shf\n",
	        hash2str(hash,k,m,klet),count0[hash],count1[hash]);
	errflag = 1;
      }
/*     else */
/*       { */
/*       printf("%c%c\t",alpha[hash/4],alpha[hash%4]); */
/*       printf("%d\t%f\n",count0[hash],(float)count0[hash]/seqlen); */
/*       } */
  return !errflag;
}		

/******************************************************/
/* edgecount: count outgoing edges from               */
/* each vertex, i.e. count (k-1)-lets                 */
/* input:      k: length of the k-lets                */
/*             m: alphabet size                       */
/*          last: hash key of last word (vertex)      */
/*         count: alphabetically ordered k-let counts */
/* output:  vdeg: alph. ordered (k-1)-let counts      */
/*                = in/out-degree of each vertex      */
/*                sufficient allocation assumed       */
/* return value: number of distinct (k-1)-lets        */
/*               (vertices)                           */
           
int edgecount(int k, int m, int last, int count[], int vdeg[])
{
  int i, j, nk1lets, sum, nver;

  nk1lets = pow(m,k-1);
  for (i = 0, nver = 0; i < nk1lets; i++)
    {
      for (j = m*i, sum = 0; j < m*(i+1); j++)
	sum += count[j];
      vdeg[i] = sum;
      if (sum > 0)
	nver++;
    }
  if (vdeg[last] == 0) /* count last vertex if not already counted */
    nver++;           
  return nver;
}      

/******************************************************/
/* kletoutput: output list of k-let counts            */
/* input:     fp: output file pointer                 */
/*             k: length of the k-lets                */
/*             m: alphabet size = highest index + 1   */
/*         count: alphabetically ordered k-let counts */

void kletoutput(FILE *fp, int k, int m, int count[])
{
  int hash, sum = 0, nklets;
  char klet[MAXORDER+1];

  nklets = pow(m,k);
  fprintf(fp,"\n");
  for (hash = 0; hash < nklets; hash++)
    {
      fprintf(fp,"%s %7i\n",hash2str(hash,k,m,klet),count[hash]);
      sum += count[hash];
    }
  fprintf(fp,"\nSum: %i\n",sum);
}

/******************************************************/
/* hash2str: k-let letter string from hash index      */
/* input:  hash: hash index                           */
/*            k: length of the k-let                  */
/*            m: alphabet size                        */
/* output: klet: k-let string (assumed big enough)    */
/* return value: klet corresponding to hash           */

char *hash2str(int hash, int k, int m, char klet[])
{
  int i;

  for (i = k-1; i >= 0; i--)
    {
      klet[i] = alpha[hash%m];
      hash /= m;
    }
  klet[k] = '\0';
  return klet;
}

/******************************************************/
/* ind2hash: hash index from k-let index string       */
/* input: klet: the k-let as string of indices        */
/*           k: length of the k-let                   */
/*           m: alphabet size                         */
/* return value: hash index corresponding to klet     */

int ind2hash(char klet[], int k, int m)
{
  int i, hash = 0;

  for (i = 0; i < k; i++)
    hash = hash*m+klet[i];
  return hash;
}

/******************************************************/
/* shuffle: make random shuffling given counts        */
/* cyclic counts assumed (k-lets wrapping end/start)  */
/* input:      m: alphabet size                       */
/*             k: length of the k-let                 */
/*          nver: number of DISTINCT (k-1)-lets       */
/*         count: alphabetically ordered k-let counts */
/*                CONTENT WILL BE DESTROYED!          */
/*          vdeg: alph. ordered (k-1)-let counts      */
/*                CONTENT WILL BE DESTROYED!          */
/*         first: hash key of first word (vertex)     */
/*          last: hash key of last word (vertex)      */
/*      lastedge: array to hold last edges (internal  */
/*                use). Required size: nver           */
/* output:   seq: random sequence (indices 0..m-1)    */
/*                where seq[0] is the k-th letter     */
/*                must be allocated to sequence length*/
/*                plus k                              */
/*                the first k-1 (invariable) letters  */
/*                are not assigned                    */
/* NB! For k=1, seq is assumed to contain a           */
/* permuatation of the sequence                       */

void shuffle(int m, int k, int nver, int count[], int vdeg[], int first, 
	     int last, int lastedge[], char seq[])
{
  int nk1lets;

  nk1lets = pow(m,k-1); /* number of (k-1)-let combinations */
/*   if (k == 1) */
/*     { */
/*       if (debugflag) */
/* 	fprintf(dlog,"simple letter shuffling\n"); */
/*       monoshuffle(vdeg[0],seq); */
/*     } */
/*   else */
/*     { */
      /* construct random inbound spanning tree with last as root */
      if (debugflag)
	fprintf(dlog,"constructing arborescence\n");
      arborescence(m,nk1lets,nver,count,vdeg,first,last,lastedge);
      /* find a random Euler trail from first, using the lastedge last */
      if (debugflag)
	fprintf(dlog,"\nmaking random trail\n");
      randomtrail(m,k,count,vdeg,first,lastedge,seq);
      if (debugflag)
	fprintf(dlog,"OK!\n");
/*     } */
}

/******************************************************/
/* monoshuffle: do simple 1-let permutation of seq    */
/* input: seqlen: sequence length                     */
/* input/output: seq: random sequence to shuffle      */

void monoshuffle(int seqlen, char seq[])
{
  int i, j;
  char temp;

  for (i = 0; i < seqlen; i++)
    {
      /* choose a random position j>i to swap */
      j = i + randomint(seqlen-i) - 1;
      temp = seq[i];
      seq[i] = seq[j];
      seq[j] = temp;
    }
}

/******************************************************/
/* arborescence: construct random inbound spanning    */
/* tree of the directed Eulerian graph having the     */
/* (k-1)-lets as vertices and the k-lets as edges,    */
/* starting at root.                                  */
/* NB! Eulerian graph <=> cyclic sequence             */
/*    <=> first (k-1)-let = last (k-1)-let            */ 
/* input:      m: alphabet size                       */
/*       nk1lets: number of (k-1)-let combinations    */
/*                = dimension of vdeg and branch      */
/*          nver: number of vertices                  */
/*                (distinct (k-1)-lets)               */
/*         count: alphabetically ordered k-let counts */
/*                (all the edge multiplicities)       */
/*          vdeg: alph. ordered (k-1)-let counts      */
/*                (# edges in/out from each vertex)   */
/*         first: index (in vdeg) of start vertex     */
/*          root: index (in vdeg) of root             */
/* output: branch: branch[i]=j represents a branch    */
/*                 from vertex i away from the root   */
/*                 corresponding to a following       */
/*                 letter index j (0<=j<m)            */
/*                 special value m for root           */
/*                 Must be allocated to nk1lets cells */


void arborescence(int m, int nk1lets, int nver, int count[], int vdeg[], 
		  int first, int root, int branch[])
{
  int vertex, remain, edgeno, edge, i;

  if (debugflag)
    fprintf(dlog,"m=%i, nver=%i, root=%i\n",m,nver,root);

  if (first != root) /* sequence not cyclic, i.e. graph not Eulerian? */
    vdeg[root]++;    /* add edge from root=last to first */

  for (i = 0; i < nk1lets; i++)
    branch[i] = -1;  /* meaning: vertex not yet hit */

  vertex = root;   /* start at the root */
  branch[root] = m; /* special value for root, also marking it as hit */
  remain = nver-1; /* number of vertices not yet hit (all except root) */
  if (debugflag)
    fprintf(dlog,"%i ",vertex);

  /* do random walk on graph until all vertices are hit */
  while (remain > 0)
    {
      /* choose a random edge from current vertex */
      edgeno = randomint(vdeg[vertex]); /* among vdeg[vertex] edges */
      if (debugflag)
	fprintf(dlog,"(%i/%i)",edgeno,vdeg[vertex]);
      if (vertex == first && first != root && edgeno == vdeg[vertex])
	{ /* special case: extra edge added from root to first */
	  if (debugflag)
	    fprintf(dlog,"xtra");
	  vertex = root;
	  continue;
	}
      edge = vertex;  /* first possible edge: same index as vertex */
      edgeno -= count[edge]; 
      while (edgeno > 0)   /* scan edges while counting down to zero */
	{
	  edge += nk1lets;
	  edgeno -= count[edge];
	}
      assert(edge < m*nk1lets); /* if fail: mismatch btw. count and vdeg! */
      vertex = edge/m;  /* first k-1 m-digits */
      if (debugflag)
	fprintf(dlog,"=%i ",vertex);
      if (branch[vertex] == -1) /* not yet visited? */
	{ 
	  branch[vertex] = edge%m;  /* store new branch (last m-digit) */
	  remain--;                     /* count it */
	  if (debugflag)
	    fprintf(dlog,"[%i]",remain);
	}
    }
  /* end of random walk */

  if (first != root)
    vdeg[root]--;   /* remove extra edge added */
}
 

/******************************************************/
/* randomtrail: construct random sequence by making   */
/* a random Euler trail at the k-let graph, at first  */
/* and using the edge in lastedge last for each vertex*/
/* input:      m: alphabet size                       */
/*             k: length of the k-let                 */
/*         count: alphabetically ordered k-let counts */
/*                (all the edge multiplicities)       */
/*                the branches in last must be removed*/
/*                CONTENT WILL BE DESTROYED!          */
/*          vdeg: alph. ordered (k-1)-let counts      */
/*                (# edges in/out from each vertex)   */
/*                CONTENT WILL BE DESTROYED!          */
/*         first: vertex to start walking             */
/*      lastedge: the last edge to use for each vertex*/
/*                indicating the next letter (0..m-1) */
/*                special value m for last vertex     */
/* output:   seq: random sequence (indices 0..m-1)    */
/*                where seq[0] is the k-th letter     */
/*                the first (k-1)-let (invariable and */
/*                equal to the last (k-1)-let) is     */
/*                NOT assigned                        */
/*                must be allocated to sequence length*/
/*                plus k                              */

void randomtrail(int m, int k, int count[], int vdeg[], int first,
		int lastedge[], char seq[])
{
  int nk1lets, vertex, nedges, letter, edge, edgeno, i;

  nk1lets = pow(m,k-1);

  assert(k >= 2);
  /* remove the last edges from the rest, stored in count */
  for (i = 0; i < nk1lets; i++)
    if (lastedge[i] >= 0 && lastedge[i] < m) 
      /* don't remove from unvisited vertices nor from last vertex */
      {
	count[m*i+lastedge[i]]--;
	vdeg[i]--;
      }

  vertex = first; /* start walking here */
  if (debugflag)
    fprintf(dlog,"%i ",vertex);
    
  /* construct random trail */
  for (i = k-1;; i++)  /* terminates by break statement */
    {
      nedges = vdeg[vertex]; /* number of remaining non-terminal edges */
      if (nedges == 0)
	{ /* last edge */
	  if (debugflag)
	    fprintf(dlog,"(*)");
	  letter = lastedge[vertex];
	  if (letter == m) /* special value for last vertex */
	    break;         /* trail complete */
	}
      else
	{ /* not last edge */
	  /* choose a random edge from current vertex */
	  edgeno = randomint(nedges);
	  if (debugflag)
	    fprintf(dlog,"(%i/%i)",edgeno,vdeg[vertex]);
	  edge = m*vertex;  /* first possible edge */
	  edgeno -= count[edge]; 
	  while (edgeno > 0)   /* scan edges while counting down to zero */
	    {
	      if (debugflag)
		fprintf(dlog,"%i,%i ",edge,count[edge]);
	      edge++;
	      edgeno -= count[edge];
	    }
	  letter = edge-m*vertex;
	  assert(letter < m); /* if fail: mismatch count/vdeg! */
	  count[edge]--;  /* this edge is now used */
	  vdeg[vertex]--; /* one less remaining */
	}	     
      seq[i] = letter;
      vertex = (m*vertex)%nk1lets+letter;
      if (debugflag)
	fprintf(dlog,"=%i ",vertex);
    }
}
