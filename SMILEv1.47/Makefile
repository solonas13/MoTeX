###############################################################################

# definition des repertoires
P_DIR=P_BLOCS/
PD_DIR=P_BLOCS+DELTA/
S_DIR=SigStat/
C_DIR=Converter/

all:    clean
	@echo "**** Making P_BLOCS..."
	(cd $(P_DIR) && $(MAKE) all)
	@echo "==> P_BLOCS OK."
	@echo "**** Making P_BLOCS+DELTA..."
	(cd $(PD_DIR) && $(MAKE) all)
	@echo "==> P_BLOCS+DELTA OK."
	@echo "**** Making SigStat..."
	(cd $(S_DIR) && $(MAKE) all)
	@echo "==> SigStat OK."

clean:
	@echo Cleaning...
	(cd $(P_DIR) && $(MAKE) clean)
	(cd $(PD_DIR) && $(MAKE) clean)
	(cd $(S_DIR) && $(MAKE) clean)
	@echo Clean OK.

