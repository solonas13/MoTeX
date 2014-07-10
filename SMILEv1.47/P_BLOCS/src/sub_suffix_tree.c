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

#include<global_fonctions.h>
#include<define.h>
#include<construction.h>
#include<unistd.h>
#include<libfasta.h>
#include<getopt.h>
 
#define SIZE 45


extern char *optarg;
extern int optind, opterr, optopt;

static struct option programme_options[] =
{
  {"seq", 1, NULL, 'f'},
  {"size", 1, NULL, 's'},
  {"print", 0, NULL, 'p'},
  {"alpha", 1, NULL, 'a'},
  {"help",0,NULL, 'h'},
  {"joker",0,NULL, 'j'},
  {"stat",0,NULL,'i'},
  {0, 0, 0, 0}
};

void printHelp(void)
{
  printf("option : --seq   -f <file sequence name>\n"
	 "         --size  -s <size of window>\n"
	 "         --print -p : print the tree\n"
	 "         --alpha -a <alphabet>\n"
	 "         --joker -j : whith joker\n"
	 "         --help  -h : print this help\n");
}
void readParametres(int argc,char **argv,int *windows_size,int *print,int *stat,char **fileName,int *joker,char **alphabet)
{
  int c;
  int index;
  
  *windows_size = DEFAULT_WINDOW_SIZE_OPTION;
  *print = 0;
  *joker = 0;
  *alphabet = NULL;
  *fileName = NULL;
  *stat=0;

  while(1)
    {
      c = getopt_long(argc,argv,"f:s:a:phji",
		      programme_options,
		      &index);
      switch(c)
	{
	case -1 :
	  if (*fileName==NULL)
	    {
	      fprintf(stderr,"option -f missing\n");
	      fprintf(stderr,"-h for more information\n");
	      exit(0);
	    }
	  if (*alphabet==NULL)
	    {
	      fprintf(stderr,"option -a missing\n");
	      fprintf(stderr,"-h for more information\n");
	      exit(0);
	    }
	  if (*windows_size<=3)
	    {
	      fprintf(stderr,"window size must be more then 3\n");
	      fprintf(stderr,"-h for more information\n");
	      exit(0);
	    }
	  return;
	case '?':
	  fprintf(stderr,"Unknow option %s \n",optarg);
	  fprintf(stderr,"-h for more information\n");
	  exit(0);
	  break;
	case 'a':
	  *alphabet =(char *)malloc(strlen(optarg)+1);
	  strcpy(*alphabet,optarg);
	  break;
	case 's':
	  *windows_size = atoi(optarg);
	  break;
	case 'i':
	  *stat = 1;
	  break;
	case 'p':
	  *print = 1;
	  break;
	case 'f':
	  *fileName = (char *)malloc(strlen(optarg)+1);
	  strcpy(*fileName,optarg);
	  break;
	case 'j':
	  *joker = 1;
	  break;
	case 'h' :
	  printHelp();
	  exit(0);
	  break;
	}
    };
}

void printParametres(int windows_size,int print,int stat,char *fileName,int joker,char *alphabet)
{
  printf("sequence file name : %s\n"
	 "alphabet           : %s\n"
	 "joker              : %s\n"
	 "windows size       : %d\n"
	 "print tree         : %s\n"
	 "statistiques       : %s\n",
	 fileName,alphabet,(joker)?"yes":"no",
	 windows_size,(print)?"yes":"no",
	 (stat)?"yes":"no");
}

int main(int argc, char **argv)
{
  Noeud *Arbre;
  FILE *fichier;
  FastaSequence	*seq;
  int ok;
  int i;
  int window_size = 4;
  int print;
  int stat;
  int joker;
  char *filename,*alphabet;


  // A  virer
  char *TEXT[5000];
  int indice_text = 0;
  //
  
  readParametres(argc,argv,&window_size,&print,&stat,&filename,&joker,&alphabet);
  printParametres(window_size,print,stat,filename,joker,alphabet);
  
  fichier = fopen(filename,"r");
  
  if (fichier==NULL)
    {
      fprintf(stderr,"Invalide sequence file name\n");
      exit(0);
    }
  
  do
    {
      seq = NewFastaSequence();
      ok = ReadFastaSequence(fichier,seq);
      if (ok)
	{
	  TEXT[indice_text] = (char *)malloc(sizeof(char )*seq->length+2);
	  memcpy(TEXT[indice_text],seq->seq,seq->length);
	  TEXT[indice_text][seq->length]='$';
	  TEXT[indice_text][seq->length+1]='\0';
	  indice_text++;
	}
    }
  while(ok);
  
  Init_All(alphabet,joker,indice_text);
  
  Arbre=Construction_Arbre(TEXT[0],window_size);
  
  //printf("Arbre de la premiere sequence %s :\n",TEXT[0]);
  //Print_Tree(Arbre,1,0);
  
  for(i=1;i<indice_text;i++)
    {
      fprintf(stdout,"\rajout de  %d ",i);
      AjouteSequence(Arbre,TEXT[i],window_size);
      //printf("Arbre apres ajout de %s : \n",TEXT[i]);
      //Print_Tree(Arbre,1,0);
    }
  //sleep(2);
  UpdateBit_TabForAllTree(Arbre);
  printf("\n");
  if ((print) || (stat))
    {
      printf("Arbre final\n");
      Print_Tree(Arbre,print,stat);
    }
  
  //sleep(20);
  //fprintf(stderr,"Liberation de l'arbre....");
  Free_Arbre(Arbre);
  //fprintf(stderr,"\tOK\n");
  //sleep(2);
  
  //fprintf(stderr,"Liberation des cellules de la liste....");
  Free_All_Liste_Cell();
  //fprintf(stderr,"\tOK\n");
  //sleep(2);
  
  //fprintf(stderr,"Liberation de la structure des listes de positions....");
  Free_ListePositions(Liste_positions_fin);
  //fprintf(stderr,"\tOK\n");
  //sleep(40);
  
  return 1;
}
