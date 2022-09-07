APP_NAME := lulesh

APP_ENTRY := $(APP_NAME).cpp
APP_ENTRY_MPI := $(APP_ENTRY)

OP2_LIBS_WITH_HDF5 := true

include ${OP2_INSTALL_PATH}/makefiles/common.mk
include ${OP2_INSTALL_PATH}/makefiles/c_app.mk