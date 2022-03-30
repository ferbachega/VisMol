#!/bin/bash
# . Bash environment variables and paths to be added to a user's ".bash_profile" file.
# . The root of the program.
VISMOL_HOME=/home/igorchem/VisMol ; export VISMOL_HOME
# . Additional paths.
VISMOL_GLCORE=$VISMOL_HOME/glCore                                    ; export VISMOL_GLCORE
VISMOL_GLWIDGET=$VISMOL_HOME/glWidget                                ; export VISMOL_GLWIDGET
