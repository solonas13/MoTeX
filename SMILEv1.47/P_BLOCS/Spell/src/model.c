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

#include <model.h>

/******************************************************************************/
/* AfficheModel                                                               */
/******************************************************************************/
void AfficheModel(P_mod m)
{
int i;

for(i=0; i<m->lon; i++)
    printf("%d",m->name[i]);
printf("\n");
}

/******************************************************************************/
/* allocModel                                                                 */
/******************************************************************************/
/* Alloue la structure d'un modele                                            */
/******************************************************************************/
P_mod allocModel(void)
{
P_mod   model;

if ( !(model = (P_mod) malloc (sizeof(Mod)) ) )
    fatalError("allocModel: cannot allocate 'model'\n"); 

if ( !(model->name = (int *) calloc (GRAIN_SIZMOD, sizeof(int)) ) )
    fatalError("allocModel: cannot allocate 'model->name'\n");

model->taille  = GRAIN_SIZMOD;
model->lon     = 0;

return model;
}


/******************************************************************************/
/* changeModel                                                                */
/******************************************************************************/
/* Ajoute un symbole au modele                                                */
/******************************************************************************/
void changeModel(P_mod mod, int symbol)
{
if (mod->lon+2 >= mod->taille)
    {
    mod->taille += GRAIN_SIZMOD;
#if DEBUG_BASE
    printf("J'etends model a %d\n",mod->taille);
#endif
    mod->name   = (int *) realloc(mod->name, mod->taille * sizeof(int));
    if (!mod->name)
        fatalError("changeModel: cannot reallocate 'model->name'\n");
    }   

(mod->name)[mod->lon]   = symbol;
(mod->lon)++;
}


/******************************************************************************/
/* decrModel                                                                  */
/******************************************************************************/
void decrModel(P_mod model)
{
if(model->lon <= 0)
    return;

model->lon--;
}

