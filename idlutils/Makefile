###############################################################################
# Sloan Digital Sky Survey (SDSS)
# IDL support code for products: idlmapper, idlspec2d
#
# S. Burles & D. Schlegel
###############################################################################

#
# IDL support utilities for spectro2d and the fibermapper
#
SHELL = /bin/sh
#
.c.o :
	$(CC) -c $(CCCHK) $(CFLAGS) $*.c
#
CFLAGS  = $(SDSS_CFLAGS) -DCHECK_LEAKS -I../include

SUBDIRS = bin data doc goddard include lib pro src ups

all :	include/export.h
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) all ); \
	done

#
# Install things in their proper places in $(IDLUTILS_DIR)
#
install :
	@echo "You should be sure to have updated before doing this."
	@echo ""
	@if [ "$(IDLUTILS_DIR)" = "" ]; then \
		echo You have not specified a destination directory >&2; \
		exit 1; \
	fi 
	@if [ -e $(IDLUTILS_DIR) ]; then \
		echo The destination directory already exists >&2; \
		exit 1; \
	fi 
	@echo ""
	@echo "You will be installing in \$$IDLUTILS_DIR=$$IDLUTILS_DIR"
	@echo "I'll give you 5 seconds to think about it"
	@sleep 5
	@echo ""
	@ rm -rf $(IDLUTILS_DIR)
	@ mkdir $(IDLUTILS_DIR)
	@ for f in $(SUBDIRS); do \
		(mkdir $(IDLUTILS_DIR)/$$f; cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) install ); \
	done
	- cp Makefile $(IDLUTILS_DIR)
	- cp RELEASE_NOTES $(IDLUTILS_DIR)

clean :
	- /bin/rm -f *~ core
	@ for f in $(SUBDIRS); do \
		(cd $$f ; echo In $$f; $(MAKE) $(MFLAGS) clean ); \
	done

# Some versions of IDL have non-backward compatible export.h files. So use whatever is current.
include/export.h:
	@if test -z "$$IDL_DIR"; then \
		echo "IDL_DIR environment variable is not set -- it must point to the top of the IDL product directory, "; \
		echo "   e.g. /usr/local/itt/idl"; \
		exit 1; \
	else \
		echo "Linking in IDL export.h, using IDL_DIR=$$IDL_DIR"; \
		test -r "$$IDL_DIR/external/export.h" || (echo "no valid $$IDL_DIR/external/export.h"; exit 1); \
		( cd include; rm export.h; ln -s $$IDL_DIR/external/export.h ); \
	fi
