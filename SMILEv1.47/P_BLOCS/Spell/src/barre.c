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

#include "barre.h"

/* Fonctions privees                                                          */
void    repChar(char, int);

/******************************************************************************/
/* repChar                                                                    */
/******************************************************************************/
void    repChar(char c, int nb)
{
int i;

for(i=0;i!=nb;i++)
    fprintf(stderr,"%c",c);
}

/******************************************************************************/
/* barre                                                                      */
/******************************************************************************/
void    barre(int n)
{
static int      pos=0;
static int      max=0;
static int      lastecrit=0;
static time_t   start;
int             nbpts;
int             sec, mn, hr;


if(n!=0)
    {
    start   = time(NULL);
    max     = n;
    lastecrit = 0;
    pos     = 1;
    fprintf(stderr,"  0%%");
#if MAXCOL != 0
    fprintf(stderr," [");
    repChar(' ',MAXCOL);
    fprintf(stderr,"]");
#endif
    return;
    }

if(pos>max)
    return;


if(pos==max)
    {
    fprintf(stderr,"\r100%%");
#if MAXCOL != 0
    fprintf(stderr," [");
    repChar(CHARBARRE,MAXCOL);
    fprintf(stderr,"]");
#endif
    fprintf(stderr," 00:00:00\n");
    pos++;
    return;
    }

sec = (int) difftime( time(NULL), start);
sec *= ((float) max-pos)/pos;
hr  = sec / 3600;
mn  = (sec - hr * 3600) / 60;
sec = sec - hr * 3600 - mn *60;

#if  MAXCOL != 0
nbpts   = pos*MAXCOL/max;
if(nbpts!=lastecrit)
    {
#endif
    fprintf(stderr,"\r%3d%%",pos*100/max);
#if  MAXCOL != 0
    fprintf(stderr," [");
    repChar(CHARBARRE,nbpts);
    repChar(' ',MAXCOL-nbpts);
    fprintf(stderr,"]");
#endif
    if(pos>3)
        fprintf(stderr, " %02d:%02d:%02d", hr,mn,sec);
#if  MAXCOL != 0
    lastecrit   = nbpts;
    }
#endif
pos++;
}

