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

#include<bit_tab.h>


int NB_SEQUENCE=0;
int CHANGE_LIMITE=0;
int SIZE_STATIC_BIT_TAB=0;

void initBitTab(int nb_seq)
{
  NB_SEQUENCE=nb_seq;
  SIZE_STATIC_BIT_TAB=nb_seq/8 +1;
  CHANGE_LIMITE=SIZE_STATIC_BIT_TAB/2 -1;
  if (CHANGE_LIMITE<=0)
    CHANGE_LIMITE=1;
}

Bit_Tab *AllocBitTab(void)
{
  Bit_Tab *tmp = (Bit_Tab *)malloc(2);
  
#if DEBUG_JTREE
  nb_alloc_tab++;
#endif
  if (tmp==NULL)
    {
      fprintf(stderr,"No enougth space... \n Program aborded \n");
      exit(-2);
    }
  memset(tmp,0,2);
  tmp[0] = 0x80; /* --> tmp est en dynamic... */
  tmp[1] = 0;    /* --> pas d'element dans tmp. */
  return tmp;
}

void ReinitBitTab(Bit_Tab **bt)
{
  if ((*bt)[0]&0x80)
    {
      *bt = (Bit_Tab *)realloc(*bt,2);
      (*bt)[0] = 0x80;
      (*bt)[1] = 0;
    }
  else
    {
      *bt = (Bit_Tab *)memset(*bt,0,SIZE_STATIC_BIT_TAB);
    }
}

void CopyBitTab(Bit_Tab **dest,Bit_Tab *src)
{
  int nb_element;
  if (src[0] & 0x80)
    { /* Dynamic : */
      nb_element = ((src[0] & 0x7F) << 8) | src[1];
      (*dest) = (Bit_Tab *)realloc(*dest,2 + nb_element*sizeof(unsigned short int) );
      *dest = (Bit_Tab *)memcpy(*dest,src,2 + nb_element*sizeof(unsigned short int) );
    }
  else
    { /* Static : */
      free(*dest);
      *dest = AllocBitTabStatic();
      /*ReinitBitTab(dest); */
      *dest = (Bit_Tab *)memcpy(*dest,src,SIZE_STATIC_BIT_TAB);
    }
}

Bit_Tab *AllocBitTabStatic(void)
{
  Bit_Tab *tmp = (Bit_Tab *)malloc(SIZE_STATIC_BIT_TAB);
  if (tmp==NULL)
    {
      fprintf(stderr,"No enougth space... \n Program aborded \n");
      exit(-2);
    }
  tmp = (Bit_Tab *)memset(tmp,0,SIZE_STATIC_BIT_TAB);
  return tmp;
}

void addBitTabValue(Bit_Tab **tab,int value)
{
  if ((*tab)[0] & 0x80)
    addBitTabValueDynamic(tab,value);
  else
    addBitTabValueStatic(tab,value);
}

void fusionneBitTab(Bit_Tab **tab1,Bit_Tab *tab2) /* 1 <- 1 & 2 */
{
  if ((*tab1)[0]&0x80)
    fusionneBitTabDynamic(tab1,tab2);
  else
    fusionneBitTabStatic(tab1,tab2);
}

int  nbSequenceInBitTab(Bit_Tab *tab)
{
  if (tab[0] & 0x80)  
    return nbSequenceInBitTabDynamic(tab);
  return  nbSequenceInBitTabStatic(tab);
}

void printBitTab(Bit_Tab *tab)
{
  int i,value;
  int nb_element;
  unsigned char mask = 0;
  if (tab[0] & 0x80)
    {
      printf("d:");
      nb_element = nbSequenceInBitTabDynamic(tab);
      for(i=0;i<nb_element;i++)
	printf("%d,",((unsigned short int *)(tab + 2))[i]);
      return;
    }
  printf("s:");
  value = 0;
  for(mask = 0x40;mask != 0 ; mask>>=1,value++)
    if (tab[0] & mask)
      printf("%d,",value);
  for(i=1;i<SIZE_STATIC_BIT_TAB;i++)
    for(mask = 0x80;mask != 0 ; mask>>=1,value++)
      if (tab[i] & mask)
	printf("%d,",value);
}

/*---------------------------------------------------------*/

void convertBitTab(Bit_Tab **tab) 
{
  Bit_Tab *tmp = AllocBitTabStatic();
  int nb_elment = nbSequenceInBitTabDynamic(*tab);
  int i;

  for(i=0;i<nb_elment;i++)
    {
      addBitTabValueStatic(&tmp,((unsigned short int *)(*tab + 2))[i]);
    }
  free(*tab);
  *tab = tmp;
}

void addBitTabValueStatic(Bit_Tab **tab,int value)
{
  int position;
  int offset;
  
  position = (value+1) / 8;
  offset = (value+1) % 8;
  offset = (0x80) >> offset;
  (*tab)[position] = (*tab)[position] | offset;
}

void addBitTabValueDynamic(Bit_Tab **tab,int value)
{
  unsigned short int nb_element = 0,i;
  nb_element = (((*tab)[0] & 0x7F) << 8) | (*tab)[1];
  
  
  if (nb_element==CHANGE_LIMITE)
    {
      convertBitTab(tab);
      addBitTabValueStatic(tab,value);
      return;
    }
  for(i=0;(i<nb_element) && (((unsigned short int *)(*tab + 2))[i]!=value);i++);
  if ((i!=nb_element) && (((unsigned short int *)(*tab + 2))[i]==value))
    return;
  nb_element++;
  (*tab) = (Bit_Tab *) realloc(*tab,2 + nb_element*sizeof(unsigned short int) );
  ((unsigned short int *)(*tab + 2))[nb_element-1] = (unsigned short int )value;
  (*tab)[0] = 0x80 | (nb_element>> 8);
  (*tab)[1] = nb_element;
}


int  nbSequenceInBitTabStatic(Bit_Tab *tab)
{
  int i;
  unsigned char k;
  int nb_element=0;
  
  for(k=0x40;k!=0;k>>=1)
    if (tab[0] & k)
      nb_element++;
  
  for(i=1;i<SIZE_STATIC_BIT_TAB;i++)
    for(k=0x80;k!=0;k>>=1)
      if (tab[i] & k)
	nb_element++;
  return nb_element;
}

int  nbSequenceInBitTabDynamic(Bit_Tab *tab)
{
  return ((((tab)[0] & 0x7F) << 8) | (tab)[1]);
}

void fusionneBitTabStatic(Bit_Tab **tab1,Bit_Tab *tab2)
{
  int i=0;
  int nb_element;
  if ((tab2)[0] &0x80)
    {
      nb_element =  nbSequenceInBitTabDynamic(tab2);
      for(i=0;i<nb_element;i++)
	addBitTabValueStatic(tab1,((unsigned short int *)(tab2 + 2))[i]);
      return;
    }
  
  for(i=0;i<SIZE_STATIC_BIT_TAB;i++)
    (*tab1)[i] |= tab2[i];
}

void fusionneBitTabDynamic(Bit_Tab **tab1,Bit_Tab *tab2)
{
  int i,nb_element;
  int value;
  unsigned char mask;
  if (tab2[0] & 0x80)
    {
      nb_element =  nbSequenceInBitTabDynamic(tab2);
      for(i=0;i<nb_element;i++)
	addBitTabValue(tab1,((unsigned short int *)(tab2 + 2))[i]);
      return;
    }


  value = 0;  
  for(mask = 0x40;mask != 0 ; mask>>=1,value++)
    if (tab2[0] & mask)
      addBitTabValue(tab1,value);
  for(i=1;i<SIZE_STATIC_BIT_TAB;i++)
    for(mask = 0x80;mask != 0 ; mask>>=1,value++)
      if (tab2[i] & mask)
	addBitTabValue(tab1,value);
}


