# Makefile

SRC_DIR := src
BIN_DIR := bin

all := ProcessAlignment 

#-------------------------------------

ProcessAlignment:
		@cd $(SRC_DIR); make
		@ln -s ../$(SRC_DIR)/$(all) $(BIN_DIR)

#-------------------------------------

.PHONY:
		clean

clean:
		cd $(SRC_DIR); make clean
		cd $(BIN_DIR); rm -f $(all)
