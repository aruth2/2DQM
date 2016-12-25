# GOMonteCarlo
This software package simulates the reduction of graphene oxide and calculates the optical response of the generated sp2 domains

![alt tag](https://raw.githubusercontent.com/aruth2/2DQM/master/GUI.png)

Prerequisites:

    GTK+3
    gnuplot
    Clapack
    Blas


Installation Instructions:

Linux:
compile with this command

    gcc 2DQM.c -o 2DQM `pkg-config --cflags --libs gtk+-3.0` -lm -lpthread -llapack -lblas

To Run:

    ./2DQM


Mac:
Use this script to install gnuplot and GTK+3:
    
    ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
    brew install gtk+3
    install gnuplot
    BASEDIR="$( dirname "$0" )"
    cd "$BASEDIR"

Then compile

    gcc 2DQM.c -o 2DQM `pkg-config --cflags --libs gtk+-3.0` -lm -lpthread -llapack -lblas

And Run
    
    ./2DQM

Windows:
There are some issues with using GTK+ in Windows, but I have compiled GTK+ code in Windows before.
If you would like to run this program on a Windows machine please contact me.
