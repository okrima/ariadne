#[=======================================================================[:
# FindGnuplot
#-----------
#
# this module looks for gnuplot
#
#
#
# Once done this will define
#
#
#
#  GNUPLOT_FOUND - system has Gnuplot
#  GNUPLOT_EXECUTABLE - the Gnuplot executable
#  GNUPLOT_VERSION_STRING - the version of Gnuplot found (since CMake 2.8.8)
#
#
#
#  GNUPLOT_VERSION_STRING will not work for old versions like 3.7.1.
#]=======================================================================]

find_program(GNUPLOT_EXEC gnuplot)

if(NOT GNUPLOT_EXEC)
    message(FATAL_ERROR "[Gnuplot] gnuplot not found")
else()
  set(GNUPLOT_FOUND TRUE)
endif()
mark_as_advanced(GNUPLOT_EXEC)
