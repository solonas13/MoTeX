###############################################################################
# OPTIONS DU PROGRAMME ########################################################
#Affichage du nombre d'occurrences pour chaque modele trouve
NB_OCCS=1
#Affichage des positions des occurrences de chaque modele trouve
AFF_OCCS=0
###############################################################################

#Compilateur :
CC=gcc
OPT=-Wall -ansi -O3
#OPT=-g

#definition des repertoires des arbres:
INCL_DIR=include/
OBJ_DIR=obj/
SRC_DIR=src/
LIB_DIR=lib/
BIN_DIR=bin/
# definition de repertoires de spell :
SPELL_INCL_DIR=Spell/include/
SPELL_SRC_DIR=Spell/src/
SPELL_OBJ_DIR=Spell/obj/

# Options de compilation :
C_FLAG= $(OPT) -I$(INCL_DIR) -I$(SPELL_INCL_DIR)
O_FLAG= -c $(OPT) -I$(INCL_DIR) -I$(SPELL_INCL_DIR) -DOCC=$(NB_OCCS) -DAFF_OCC=$(AFF_OCCS)
LIB_FLAG= -lm
# Nom du programme , de l'archive
PROG_NAME=x-smile
ARCHIVE_NAME=x-smile.tar.gz
#Objet ...
OBJ_FILE=$(OBJ_DIR)*.o $(SPELL_OBJ_DIR)*.o

TREE_OBJ=$(OBJ_DIR)global_fonctions.o $(OBJ_DIR)construction.o $(OBJ_DIR)liste_pos.o $(OBJ_DIR)allocateurs.o $(OBJ_DIR)bit_tab.o $(OBJ_DIR)libfasta.o $(OBJ_DIR)global_variables.o  $(SPELL_OBJ_DIR)global.o $(SPELL_OBJ_DIR)criteres.o #$(OBJ_DIR)sub_suffix_tree.o

TREE_OBJ_DEB=$(OBJ_DIR)global_fonctions.o $(OBJ_DIR)construction.o $(OBJ_DIR)liste_pos2.o $(OBJ_DIR)allocateurs.o $(OBJ_DIR)bit_tab2.o $(OBJ_DIR)libfasta.o $(OBJ_DIR)global_variables.o  $(SPELL_OBJ_DIR)global.o $(SPELL_OBJ_DIR)criteres.o $(OBJ_DIR)sub_suffix_tree.o

SPELL_OBJ=$(SPELL_OBJ_DIR)global.o $(SPELL_OBJ_DIR)pile_occ.o $(SPELL_OBJ_DIR)spell.o $(SPELL_OBJ_DIR)occ.o $(SPELL_OBJ_DIR)criteres.o $(SPELL_OBJ_DIR)barre.o $(SPELL_OBJ_DIR)alphabet.o $(SPELL_OBJ_DIR)model.o

all:    $(BIN_DIR)$(PROG_NAME)
	echo ALL!!!

debug:  $(TREE_OBJ_DEB)
	$(CC) $(C_FLAG) $(OBJ_FILE) -o $(BIN_DIR)$(PROG_NAME) $(LIB_FLAG)

$(BIN_DIR)$(PROG_NAME): $(TREE_OBJ) $(SPELL_OBJ)
	$(CC) $(C_FLAG) $(OBJ_FILE) -o $(BIN_DIR)$(PROG_NAME) $(LIB_FLAG)

$(OBJ_DIR)sub_suffix_tree.o :  $(SRC_DIR)sub_suffix_tree.c
	$(CC) $(SRC_DIR)sub_suffix_tree.c -o $(OBJ_DIR)sub_suffix_tree.o $(O_FLAG) 

$(OBJ_DIR)global_fonctions.o : $(INCL_DIR)global_fonctions.h $(SRC_DIR)global_fonctions.c
	$(CC) $(SRC_DIR)global_fonctions.c -o $(OBJ_DIR)global_fonctions.o $(O_FLAG) 

$(OBJ_DIR)construction.o : $(INCL_DIR)construction.h $(SRC_DIR)construction.c
	$(CC) $(SRC_DIR)construction.c -o $(OBJ_DIR)construction.o $(O_FLAG) 

$(OBJ_DIR)liste_pos.o : $(INCL_DIR)liste_pos.h $(SRC_DIR)liste_pos.c
	$(CC) $(SRC_DIR)liste_pos.c -o $(OBJ_DIR)liste_pos.o $(O_FLAG) 

$(OBJ_DIR)liste_pos2.o : $(INCL_DIR)liste_pos.h $(SRC_DIR)liste_pos2.c
	$(CC) $(SRC_DIR)liste_pos2.c -o $(OBJ_DIR)liste_pos2.o $(O_FLAG) 

$(OBJ_DIR)allocateurs.o : $(INCL_DIR)allocateurs.h $(SRC_DIR)allocateurs.c
	$(CC) $(SRC_DIR)allocateurs.c -o $(OBJ_DIR)allocateurs.o $(O_FLAG) 

$(OBJ_DIR)bit_tab.o : $(INCL_DIR)bit_tab.h $(SRC_DIR)bit_tab.c
	$(CC) $(SRC_DIR)bit_tab.c -o $(OBJ_DIR)bit_tab.o $(O_FLAG) 

$(OBJ_DIR)bit_tab2.o : $(INCL_DIR)bit_tab.h $(SRC_DIR)bit_tab2.c
	$(CC) $(SRC_DIR)bit_tab2.c -o $(OBJ_DIR)bit_tab2.o $(O_FLAG) 

$(OBJ_DIR)libfasta.o : $(INCL_DIR)libfasta.h $(SRC_DIR)libfasta.c
	$(CC) $(SRC_DIR)libfasta.c -o $(OBJ_DIR)libfasta.o $(O_FLAG) 

$(OBJ_DIR)global_variables.o : $(SRC_DIR)global_variables.c
	$(CC) $(SRC_DIR)global_variables.c -o $(OBJ_DIR)global_variables.o $(O_FLAG) 

#Dependance pour le programme spell :

$(SPELL_OBJ_DIR)global.o : $(SPELL_INCL_DIR)global.h $(SPELL_SRC_DIR)global.c
	$(CC) $(SPELL_SRC_DIR)global.c -o $(SPELL_OBJ_DIR)global.o $(O_FLAG) 

$(SPELL_OBJ_DIR)criteres.o : $(SPELL_INCL_DIR)criteres.h $(SPELL_SRC_DIR)criteres.c $(SPELL_INCL_DIR)global.h
	$(CC) $(SPELL_SRC_DIR)criteres.c -o $(SPELL_OBJ_DIR)criteres.o $(O_FLAG) 

$(SPELL_OBJ_DIR)pile_occ.o : $(SPELL_INCL_DIR)pile_occ.h $(SPELL_SRC_DIR)pile_occ.c $(SPELL_INCL_DIR)global.h
	$(CC) $(SPELL_SRC_DIR)pile_occ.c -o $(SPELL_OBJ_DIR)pile_occ.o $(O_FLAG) 

$(SPELL_OBJ_DIR)spell.o : $(SPELL_INCL_DIR)spell.h $(SPELL_SRC_DIR)spell.c $(SPELL_INCL_DIR)global.h
	$(CC) $(SPELL_SRC_DIR)spell.c -o $(SPELL_OBJ_DIR)spell.o $(O_FLAG) 

$(SPELL_OBJ_DIR)occ.o : $(SPELL_INCL_DIR)occ.h $(SPELL_SRC_DIR)occ.c
	$(CC) $(SPELL_SRC_DIR)occ.c -o $(SPELL_OBJ_DIR)occ.o $(O_FLAG) 

$(SPELL_OBJ_DIR)barre.o : $(SPELL_INCL_DIR)barre.h $(SPELL_SRC_DIR)barre.c
	$(CC) $(SPELL_SRC_DIR)barre.c -o $(SPELL_OBJ_DIR)barre.o $(O_FLAG) 

$(SPELL_OBJ_DIR)alphabet.o : $(SPELL_INCL_DIR)alphabet.h $(SPELL_SRC_DIR)alphabet.c
	$(CC) $(SPELL_SRC_DIR)alphabet.c -o $(SPELL_OBJ_DIR)alphabet.o $(O_FLAG) 

$(SPELL_OBJ_DIR)model.o : $(SPELL_INCL_DIR)model.h $(SPELL_SRC_DIR)model.c
	$(CC) $(SPELL_SRC_DIR)model.c -o $(SPELL_OBJ_DIR)model.o $(O_FLAG) 



clean: clean_emacs clean_obj 
	echo CLEAN OK!!!

clean_emacs:
	find . -name "*~" -exec rm -f {} \;

clean_obj:
	rm -f $(OBJ_DIR)*o $(SPELL_OBJ_DIR)*o
	rm -f $(BIN_DIR)$(PROG_NAME)

clean_arch:
	rm ./$(ARCHIVE_NAME)

tgz: clean
	tar -zcvf $(ARCHIVE_NAME) ./*

depend:
	makedepend -o.o $(SPELL_INCL_DIR)*.h $(INCL_DIR)*.h $(SPELL_SRC_DIR)*.c $(SRC_DIR)*.c -I$(INCL_DIR) -I$(SPELL_INCL_DIR)

# DO NOT DELETE
