# BEHR
Bayesian Estimation of Hardness Ratios

BEHR is a standalone command-line C program designed to quickly estimate the hardness ratios and their uncertainties for astrophysical sources. It is especially useful in the Poisson regime of low counts, and computes the proper uncertainty regardless of whether the source is detected in both passbands or not.

The code

The latest version is ver.12-12-2013

The code is a combination of Fortran and C-based programs. The theory behind it is described in detail in the BEHR paper (Park et al., 2006, ApJ, 652, 610; this is the citation to use if you used this code in your analysis). It has been successfully compiled and run on SunOS 5.8, various flavors of macOS and Linux, using g77, gfortran, and gcc. To install, first download the tar file that contains all the necessary files to a suitable directory (e.g., /path/to/BEHR/installation/directory/) and then type

	cd /path/to/BEHR/installation/directory/
	tar xvf BEHR.tar
	cd BEHR
	make clean
	make BEHR
which will create an executable named BEHR. It may be necessary to delete the "-arch i386" flags on newer systems in the Makefile. Also, make sure that the Makefile is using the right gcc and gfortran. The executable can be run from the command line. If run without arguments, it will print a short usage syntax and exit. If run as
	prompt> ./BEHR help
it will print a more detailed description of the parameters before exiting.

NOTES ON COMPILATION ERRORS: People that try to compile with gfortran may encounter some compilation errors on macOS. Here is a brief rundown of the kind of things to try to fix them.
-- The gfortran compiler listed in the Makefile in the tar file is v4.3 of the macports port. Change it to whichever one you have installed. A good repository of the latest versions is at hpc.sourceforge.net. (Be sure to update gcc too.)
-- The Makefile targets the i386 architecture. The -arch i386 flags to gcc can be omitted from newer versions of macOS, or explicitly included with -arch x86_64.
-- If you get errors of the sort
pquad.f:1.72:
      subroutine pquad (f, a, b, absacc, relacc, minrul, maxrul,ruleno, 
Error: Invalid character in name at (1)
This is due to what appears to be a different interpretation of the continuation character. To fix the errors, edit the file pquad.f and move all the &'s from the first column to the sixth column, i.e., from
      subroutine pquad (f, a, b, absacc, relacc, minrul, maxrul,ruleno,
&     quads, nevals, ercode)
to
      subroutine pquad (f, a, b, absacc, relacc, minrul,maxrul,ruleno,
     &quads, nevals, ercode)
-- If you have macports installed, the gcc assembler in /opt/local/bin might interfere with the compilation. To get around it, do
prompt> set newpath = `echo $PATH | sed 's,/opt/local/bin:,,g'`
prompt> set path = ( /usr/bin $newpath )
prompt> make clean ; make BEHR
