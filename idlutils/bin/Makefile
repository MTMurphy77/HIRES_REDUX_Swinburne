###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

SHELL = /bin/sh
#
all :;
#
# Install things in their proper places in $(IDLUTILS_DIR)

EXECS1 = Makefile build_links evilmake eviltarball idlutils_version install_ups_table memtrack pstogif sxpar
EXECS2 = balkanize drangle harmonize map poly2poly polyid ransack snap unify
EXECS3 = mrb_balkanize mrb_snap mrb_unify weight pixelize

EXECS = $(EXECS1) $(EXECS2) $(EXECS3)

install :
	@ for f in $(EXECS); do \
		(/bin/cp $$f $(IDLUTILS_DIR)/bin/ ); \
	done


clean :
	- /bin/rm -f *~ core
