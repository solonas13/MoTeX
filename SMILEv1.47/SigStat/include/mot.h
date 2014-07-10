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

#ifndef MOT_H
#define MOT_H

/******************************************************************************/
/* STRUCTURES                                                                 */
/******************************************************************************/
typedef struct struct_mot
    {
    char    *mot;
    char    *codes;
    float   quorum_reel;
    float   khi2;
    float   moyenne_shuffle;        /* sigstat */
    float   zscore;                 /* sigstat */
    float   sigma;                  /* sigstat */
    unsigned int    nbseq_vrai;
    unsigned int    nbseq_faux;     /* faux */
#if OCC
    unsigned int    nboccex_vrai;
    unsigned int    nboccex_faux;   /* faux */
    float   moyenne_shuffle_occ;     /* sigstat */
    float   khi2_occ;
    float   zscore_occ;             /* sigstat */
    float   sigma_occ;              /* sigstat */
    char    sign_occ;               /* faux */
#endif
    char    sign;                   /* faux */
    } Mot;


#endif

