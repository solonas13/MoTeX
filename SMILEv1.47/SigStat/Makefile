# OPTIONS DU PROGRAMME GREP+ ##################################################
#Statistiques sur les occurrences exactes                        
OCC=1
###############################################################################

#Compilateur :
CC=gcc
#OPT=-g
OPT=-ansi -Wall -O3

# definition de repertoires de grep+ :
GREP_DIR=grep+/
GREP_INCL_DIR=$(GREP_DIR)Spell/include/
GREP_INCL_DIR2=$(GREP_DIR)include/
GREP_OBJ_DIR=$(GREP_DIR)Spell/obj/
GREP_OBJ_DIR2=$(GREP_DIR)obj/

#definition des repertoires de sigstat
SIG_INCL_DIR=include/
SIG_SRC_DIR=src/
SIG_OBJ_DIR=obj/
SIG_BIN_DIR=bin/

# Options de compilation :
C_FLAG= $(OPT) -I$(GREP_INCL_DIR) -I$(GREP_INCL_DIR2) -I$(SIG_INCL_DIR)
O_FLAG= -c $(OPT) -I$(GREP_INCL_DIR) -I$(GREP_INCL_DIR2) -I$(SIG_INCL_DIR) -DOCC=$(OCC)
LIB_FLAG= -lm
# Nom du programme, de l'archive
PROG_NAME=e-smile_shuffling
PROG_NAME2=e-smile_against
ARCHIVE_NAME=e-smile_shuffling.tar.gz
ARCHIVE_NAME2=e-smile_against.tar.gz

#Objets
SIG_OBJ=$(SIG_OBJ_DIR)euler.o $(SIG_OBJ_DIR)random.o  $(SIG_OBJ_DIR)shufflet.o $(SIG_OBJ_DIR)sigstat.o $(SIG_OBJ_DIR)seqio.o $(SIG_OBJ_DIR)io.o
FAUX_OBJ=$(SIG_OBJ_DIR)faux.o $(SIG_OBJ_DIR)io.o

TREE_OBJ=$(GREP_OBJ_DIR2)global_fonctions.o $(GREP_OBJ_DIR2)construction.o $(GREP_OBJ_DIR2)liste_pos.o $(GREP_OBJ_DIR2)allocateurs.o $(GREP_OBJ_DIR2)bit_tab.o $(GREP_OBJ_DIR2)libfasta.o $(GREP_OBJ_DIR2)global_variables.o

GREP_OBJ=$(GREP_OBJ_DIR)global.o $(GREP_OBJ_DIR)pile_occ.o $(GREP_OBJ_DIR)grep+.o $(GREP_OBJ_DIR)occ.o $(GREP_OBJ_DIR)criteres.o $(GREP_OBJ_DIR)barre.o $(GREP_OBJ_DIR)alphabet.o $(GREP_OBJ_DIR)model.o

OBJ_FILE=$(TREE_OBJ) $(GREP_OBJ) $(SIG_OBJ)
OBJ_FILE2=$(TREE_OBJ) $(GREP_OBJ) $(FAUX_OBJ)


# Regles

all:    grep+ $(SIG_BIN_DIR)$(PROG_NAME) $(SIG_BIN_DIR)$(PROG_NAME2)
	@echo OK.

grep+::
	@echo Making grep+...
	(cd $(GREP_DIR) && $(MAKE) obj)
	@echo Fin du Make de grep+.

$(SIG_BIN_DIR)$(PROG_NAME): $(SIG_OBJ)
	@echo Making sigstat...
	$(CC) $(C_FLAG) $(OBJ_FILE) -o $(SIG_BIN_DIR)$(PROG_NAME) $(LIB_FLAG)

$(SIG_BIN_DIR)$(PROG_NAME2): $(FAUX_OBJ)
	@echo Making faux...
	$(CC) $(C_FLAG) $(OBJ_FILE2) -o $(SIG_BIN_DIR)$(PROG_NAME2) $(LIB_FLAG)

#Sigstat
$(SIG_OBJ_DIR)euler.o : $(SIG_INCL_DIR)shufflet.h $(SIG_SRC_DIR)euler.c 
	$(CC) $(SIG_SRC_DIR)euler.c -o $(SIG_OBJ_DIR)euler.o $(O_FLAG)

$(SIG_OBJ_DIR)random.o : $(SIG_INCL_DIR)shufflet.h $(SIG_SRC_DIR)random.c 
	$(CC) $(SIG_SRC_DIR)random.c -o $(SIG_OBJ_DIR)random.o $(O_FLAG)

$(SIG_OBJ_DIR)shufflet.o : $(SIG_INCL_DIR)shufflet.h $(SIG_SRC_DIR)shufflet.c
	$(CC) $(SIG_SRC_DIR)shufflet.c -o $(SIG_OBJ_DIR)shufflet.o $(O_FLAG)

$(SIG_OBJ_DIR)seqio.o : $(SIG_INCL_DIR)shufflet.h $(SIG_SRC_DIR)seqio.c 
	$(CC) $(SIG_SRC_DIR)seqio.c -o $(SIG_OBJ_DIR)seqio.o $(O_FLAG)

$(SIG_OBJ_DIR)sigstat.o : $(SIG_INCL_DIR)mot.h $(SIG_INCL_DIR)shufflet.h $(SIG_SRC_DIR)sigstat.c
	$(CC) $(SIG_SRC_DIR)sigstat.c -o $(SIG_OBJ_DIR)sigstat.o $(O_FLAG)

$(SIG_OBJ_DIR)faux.o : $(SIG_INCL_DIR)mot.h $(SIG_SRC_DIR)faux.c
	$(CC) $(SIG_SRC_DIR)faux.c -o $(SIG_OBJ_DIR)faux.o $(O_FLAG)

$(SIG_OBJ_DIR)io.o : $(SIG_INCL_DIR)io.h $(SIG_SRC_DIR)io.c
	$(CC) $(SIG_SRC_DIR)io.c -o $(SIG_OBJ_DIR)io.o $(O_FLAG)



clean: clean_emacs clean_obj 
	(cd $(GREP_DIR) && $(MAKE) clean)
	@echo CLEAN OK!!!

clean_emacs:
	find . -name "*~" -exec rm -f {} \;

clean_obj:
	rm -f $(SIG_OBJ_DIR)*o
	rm -f $(SIG_BIN_DIR)$(PROG_NAME) $(SIG_BIN_DIR)$(PROG_NAME2)

clean_arch:
	rm ./$(ARCHIVE_NAME)

tgz: clean
	tar -zcvf $(ARCHIVE_NAME) ./*

