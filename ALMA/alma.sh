#
# This is 'alma.sh' - Last modified by GS 6/5/2007 
# 
ECHO " +++++++++++++++++++++++++++++++++++++++++++++++"
ECHO "                    A L M A                     "                 
ECHO "     Version 1.0 - last modified June 2007      "                   
ECHO "              Copyright (C) 2007                " 
ECHO "  Giorgio Spada -Email: giorgio.spada@gmail.com "
ECHO " +++++++++++++++++++++++++++++++++++++++++++++++"
#
#
#
g95 alma.f90 -ffree-form -c -o alma.o
g95  alma.o FMSAVE.o FMZM90.o FM.o -o alma.exe
./alma.exe 
#
#
ECHO ""
ECHO " +++++++++++++++++++++++++++++++++++++++++++++++"
ECHO "                    A L M A                     "                 
ECHO "     Version 1.0 - last modified June 2007      "                   
ECHO "              Copyright (C) 2007                " 
ECHO "  Giorgio Spada -Email: giorgio.spada@gmail.com "
ECHO " +++++++++++++++++++++++++++++++++++++++++++++++"
#
#
