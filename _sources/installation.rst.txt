.. _installation:

############
Installation
############

..  .. contents::
        :local:
        :depth: 2

.. note::

    Want to try ecoevolity without installing it?
    :ref:`Try out our Docker image <docker_install>`.
    If you run into troubles below building and installing |eco|, the Docker
    image provides an alternative.


.. _prerequisites:

*************
Prerequisites
*************

Compiling |eco| requires |cmake|_ and a new-ish C++ compiler (one that
supports the C++11 standard).

We also strongly recommend using |git|_ to acquire the source code.
If you use |git|_, the |eco| tools will report the version of
|eco| you are using much more precisely, which will make your
work more reproducible.

While not required, the |pyco|_ Python package can be useful for summarizing
the output of |eco|.
|Pyco| will be used in tutorials.
:ref:`See below <pycoevolity_install>`
for how to install |pyco|.

Linux
=====

If you are on a Linux machine, use your package manager to update/install g++,
|cmake|, and |git|.
On Ubuntu (or other Debian-based Linux distribution) you can simply use::

    $ sudo apt-get install cmake g++ git

Mac
===

If you are on a Mac, you can install Xcode command line tools (that'll get you
a C++ compiler and |git|) and download and install |cmake| from
https://cmake.org/download/.
Download the ``.dmg`` file for the latest version for Mac.
Once installed, open CMake, and from the "Tools" menu click 
"How to Install For Command Line Use." This should pop up a window that has a
line like::

    sudo "/Applications/CMake.app/Contents/bin/cmake-gui" --install

The line might be slightly different on your Mac.
Copy that line from the CMake pop up window (not the line above), and
paste it into the command line (Terminal).

Windows
=======

We have not compiled |eco| in Windows.
It's possible to do, but I suspect it will take a fair bit of tweaking of the
|cmake|_ configuration.
If you have a Windows machine, and you are not experienced with compiling C/C++ code
in Windows, you have some options:

#.  If you have Windows 10,
    `you can install Ubuntu <https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0>`_.
    Once installed, then you can follow the instructions above for Linux.

#.  Try running our Docker image
    (:ref:`see below <docker_install>`).


*****************************
Getting the |eco| source code
*****************************

If you have |git|_ installed, download the |eco| repository with the ``clone``
command::

    $ git clone https://github.com/phyletica/ecoevolity.git

If you prefer not to use |git|_, you can download an archive of the 
`source code here <https://github.com/phyletica/ecoevolity/archive/master.zip>`_.
Once downloaded, you'll need to unzip the archive.


***********************
Basic build and install 
***********************

Next, move into the downloaded directory and run the build script::

    $ cd ecoevolity
    $ ./build.sh

If the build was successful, the |eco| executables should now be in the
``./build/bin`` directory, and you should be able to run::

    $ ./build/bin/ecoevolity -h

and see the |eco| help menu, the beginning of which should looks something
like::

    ======================================================================
                                  Ecoevolity
                      Estimating evolutionary coevality
            Version 0.2.0 (dev d4a4d48: 2018-05-23T15:22:40-05:00)
    ======================================================================
    
    Usage: ecoevolity [OPTIONS] YAML-CONFIG-FILE
    
    Ecoevolity: Estimating evolutionary coevality
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit

The executables in ``./build/bin`` are ready to use, but you'll probably want
to put them in your PATH (a list of directories that your shell looks in to
find the commands you type on the command line). You can do this via::

    $ sudo cp ./build/bin/* /usr/local/bin

If this worked, you're good to go; you can try ``ecoeovlity -h`` to be sure.

If it didn't work, you probably don't have admin privileges.
If so, you can create a bin folder in your home folder and put the tools
there::

    $ mkdir -p "${HOME}/bin"
    $ cp ./build/bin/* "${HOME}/bin"

Then, you can add this directory to your PATH (if it's not already there; you
can check with ``echo $PATH``)::

    $ export PATH="${PATH}:${HOME}/bin"

Note, this update to PATH is only for your current terminal window.  If you
want this to be permanent (work for all future terminal windows), add ``export
PATH="${PATH}:${HOME}/bin"`` to your ``.bashrc`` or ``.bash_profile`` file in
your home directory.


********************
Install during build
********************

If you want to build and install in one go, you just need to specify where you
want the installation to go, for example::

    $ sudo ./build.sh --prefix /usr/local


*****************************
Building the threaded version
*****************************

If you want to install a version of |eco| that performs the likelihood
calculations across multiple threads, you just need to add the ``--threads``
flag::

    $ ./build.sh --threads

In my opinion, you're usually better off running multiple independent chains
rather than multithreading, but the option is there.


.. _pycoevolity_install:

**********************
Installing pycoevolity
**********************

|Pyco|_ is a Python package for summarizing the output of |eco|.
It should work with Python 2 or 3.
If you have
`Python <https://www.python.org/>`_
and 
`pip <https://pypi.org/project/pip/>`_
installed, you can install |Pyco| via::

    $ pip install git+git://github.com/phyletica/pycoevolity.git

If this isn't working, try the
`manual installation instructions here <https://github.com/phyletica/pycoevolity>`_.
Also, |pyco| uses the
`R <https://www.r-project.org/>`_
packages
`ggplot2 <http://ggplot2.tidyverse.org/>`_
and
`ggridges <https://github.com/clauswilke/ggridges>`_
for creating some plots.
So, if you want plotting by pycoevolity to be fully functional,
and you don't already have
`R <https://www.r-project.org/>`_
installed, you'll need to install it.
Once
`R <https://www.r-project.org/>`_
is in place, you can install the packages from the
`R <https://www.r-project.org/>`_
prompt using:::

    install.packages(c("ggplot2", "ggridges"))


.. _docker_install:

***********************************
Using ecoevolity without installing
***********************************

Docker provides a nice way of sharing lightweight containers that act like a
virtual machine.
We have created a Docker container with |eco| built in.
To get started, you first need to 
`install Docker <https://www.docker.com/community-edition>`_.
If you're on a Mac or Windows machine, might need to sign up for a free Docker
account to download it.
Once Docker is installed and running pull down our Docker image::

    $ docker pull phyletica/ecoevolity-docker

.. note::

    Depending on your system and how Docker is configured, you may need to use
    ``sudo`` to run Docker commands. If you received a "permission denied"
    message when you ran the command above, try::
    
        $ sudo docker pull phyletica/ecoevolity-docker

This download could take several minutes depending on your internet connection.
After it completes, run and enter the docker container::

    $ docker run -it phyletica/ecoevolity-docker bash

.. note::

    Again, you might need to prefix this command with ``sudo``.

That's it, you are now in a virtual container with 
a fully functioning |eco| ecosystem
(|eco| and |pyco| are installed, along with example data).
Try typing::

    $ ecoevolity -h

This should display the |eco| help menu.
Next, ``cd`` into the example data directory::

    $ cd ecoevolity-example-data
    $ ls

There you will find an |eco| configuration file and nexus-formatted data files.
Go ahead and run an |eco| analysis::

    $ ecoevolity --relax-missing-sites --relax-triallelic-sites --ignore-data ecoevolity-config.yml

To exit the container, simply type::

    $ exit

Docker will keep the |eco| image around, so you can always jump
back in anytime via::

    $ docker run -it phyletica/ecoevolity-docker bash

However, any files you created on your last visit will be gone.
So, if you want to analyze *your* data and keep the results around, ``cd``
to the directory where you want to run |eco|, then jump into
the Docker container using::

    $ docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash

Then, once inside, type::

    $ cd portal
    $ ls

You should see the files that were in the directory on *your* computer.
Now you can run |eco| on data in this directory, and all output files will be
on your computer when you exit the container (magic!).
