.. _installation:

************
Installation
************

.. .. note::
.. 
..     Want to try ecoevolity without installing it?
..     :ref:`Try out our Docker image <docker_install>`.
..     If you run into trouble below building and installing |eco|, the Docker
..     image provides an alternative.

.. note::

    While not required, the |pyco|_ Python package can be useful for
    summarizing the output of |eco|, and is used in some of the tutorials.
    :ref:`See below <pycoevolity_install>`
    for how to install |pyco|.

There are three options for installing |eco|: Downloading precompiled |eco|
programs, compiling from the source code, or using the Docker image.
Use the outline below to help navigate the instructions for your preferred
installation option.

.. contents::
    :local:
    :depth: 3


Downloading precompiled programs
================================

Download and extract precompiled programs 
-----------------------------------------

Linux
^^^^^

Use the following link that matches your CPU architecture to download |eco|'s
precompiled programs:

-   |eco_latest_linux_amd_binaries|_
-   |eco_latest_linux_arm_binaries|_

To determine your CPU architecture, you can run the
``uname -m``
or
``dpkg --print-architecture``
command.

The download is a gzipped archive.
To extract it, use one of these two commands, depending on which architecture
you downloaded:

.. code-block:: shell
    :substitutions:

    # For AMD arch download
    tar xzf ecovolity-|eco_latest_tag|-amd64.tar.gz
    # For ARM arch download
    tar xzf ecovolity-|eco_latest_tag|-arm64.tar.gz

This should result in a directory (folder) named
:substitution-code:`ecoevolity-|eco_latest_tag|`

Try running the ``ecoevolity`` program in this directory to get the help menu:

.. code-block:: shell
    :substitutions:

    ./ecoevolity-|eco_latest_tag|/bin/ecoevolity -h

If it works, the beginning of help menu output should look something like::

    ======================================================================
                                  Ecoevolity
                      Estimating evolutionary coevality
           Version 1.1.0 (HEAD d64cc24: 2026-01-01T14:52:21-06:00)
    ======================================================================
    
    Usage: ecoevolity [OPTIONS] YAML-CONFIG-FILE
    
    Ecoevolity: Estimating evolutionary coevality
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit

Didn't work? Try the troubleshooting tips below.

.. admonition:: **Troubleshooting missing linked libraries**

    If you get an error when trying to run ``ecoevolity`` about a missing library,
    you might need to install the ``build-essential`` package.
    For example::
    
        sudo apt-get install build-essential
    
    This should install the missing libraries.

Mac OS X
^^^^^^^^

Use the following link that matches your Mac's hardware to download |eco|'s
precompiled programs:

-   |eco_latest_mac_intel_binaries|_
-   |eco_latest_mac_arm_binaries|_

The download is a gzipped archive. To extract it, use one of these two commands
depending on which architecture you downloaded:

.. code-block:: shell
    :substitutions:

    # For Intel download
    tar xzf ecovolity-|eco_latest_tag|-mac-intel64.tar.gz
    # For Apple silicon download
    tar xzf ecovolity-|eco_latest_tag|-mac-arm64.tar.gz

This should result in a directory (folder) named
:substitution-code:`ecoevolity-|eco_latest_tag|`

Try running the ``ecoevolity`` program in this directory to get the help menu:

.. code-block:: shell
    :substitutions:

    ./ecoevolity-|eco_latest_tag|/bin/ecoevolity -h

If it works, the beginning of the help menu output should look something like::

    ======================================================================
                                  Ecoevolity
                      Estimating evolutionary coevality
           Version 1.1.0 (HEAD d64cc24: 2026-01-01T14:52:21-06:00)
    ======================================================================
    
    Usage: ecoevolity [OPTIONS] YAML-CONFIG-FILE
    
    Ecoevolity: Estimating evolutionary coevality
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit

Didn't work? Try the troubleshooting tips below.

.. admonition:: **Troubleshooting Mac OS security**

    When you try running the ``ecoeovlity`` program, there's a good chance that
    Mac OS will state that it can't be opened because it cannot check it for
    malicious software (or can't identify the developer).
    Click "OK" and then go into the "Security & Privacy" section of Sytem
    Preferences and click "Allow Anyway" where it states that "ecoevolity" was
    blocked.
    The next time you try to run
    :substitution-code:`./ecoevolity-|eco_latest_tag|/bin/ecoevolity -h`,
    the OS will complain again, but now there should be an "Open" button to
    click.
    After clicking "Open" ecoevolity should run and you should see the help menu.
    When you run ``ecoevolity`` in the future, the OS shouldn't complain anymore.
    
    You will have to repeat the above for all the programs in 
    :substitution-code:`./ecoevolity-|eco_latest_tag|/bin`
    you would like to use.

.. admonition:: **Troubleshooting missing linked libraries**

    When trying to run the binary, if you get an error about a library not being
    found, you might need to install Xcode command line tools.
    For example::
    
        xcode-select --install
    
    The Xcode installation should include the missing libraries.


Install the precompiled programs 
--------------------------------

If running ``ecoevolity`` to view the help menu above worked, the programs in
:substitution-code:`./ecoevolity-|eco_latest_tag|/bin`
are ready to use, but you'll probably
want to put them in your PATH (a list of directories that your shell searches
to find the commands you type on the command line).

Installing locally (no admin privileges required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install the programs locally within your home folder, you can use:

.. code-block:: shell
    :substitutions:

    mkdir -p ~/.local/bin
    cp ./ecoevolity-|eco_latest_tag|/bin/* ~/.local/bin

The ``~/.local/bin`` directory is likely not in your PATH
variable, but you can check using
``echo $PATH``, which will show a ":"-separated list of directories in your
PATH.
To add the ``~/.local/bin`` directory to your PATH, you can use::

    export PATH="${PATH}:${HOME}/.local/bin"

Now, running ``ecoevolity -h`` on the command line should work.

.. note::

    This update to PATH is only for your current terminal window.
    If you want this to be permanent (work for all future terminal windows),
    add the line:: 
    
        export PATH="${PATH}:${HOME}/.local/bin"
    
    to your shell's startup file. E.g.,:
    
    -   ``~/.bashrc`` file in your home directory if you are using the bash shell
    -   ``~/.zshrc`` file in your home directory if you are using the zsh shell

Installing globally (admin privileges required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have admin privileges and want to install the programs to a location
that is almost certainly in your PATH, you can use:

.. code-block:: shell
    :substitutions:

    sudo cp ./ecoevolity-|eco_latest_tag|/bin/* /usr/local/bin

You can check things are working by running ``ecoevolity -h``, which should
show a help menu.


Compiling ecoevolity from source code
=====================================

.. _prerequisites:

Installing prerequisites
------------------------

Compiling |eco| requires |cmake|_ and a new-ish C++ compiler (one that
supports the C++11 standard).

.. We also strongly recommend using |git|_ to acquire the source code.
   If you use |git|_, the |eco| tools will report the version of
   |eco| you are using much more precisely, which will make your
   work more reproducible.

Linux
^^^^^

If you are on a Linux machine, use your package manager to install the required
software for compiling |eco|.
On Ubuntu (or other Debian-based Linux distribution) you can use::

    sudo apt-get install build-essential cmake git

Mac
^^^

If you are on a Mac, you can install Xcode command line tools using the following command::

    xcode-select --install

This will get you a C++ compiler and |git|.
You also need to install |cmake|.
If you have |homebrew| installed you can install |cmake| via::

    brew install cmake

If you don't have |homebrew|, you can download and install |cmake| from
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
^^^^^^^

We have not compiled |eco| in Windows.
It's possible to do, but I suspect it will take a fair bit of tweaking of the
|cmake|_ configuration.
If you have a Windows machine, and you are not experienced with compiling C/C++ code
in Windows, you have some options:

#.  You can use
    `Windows Subsystem for Linux to install Ubuntu <https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0>`_.
    Once installed, then you can follow the instructions above for Linux.
#.  Try running our Docker image (:ref:`see below <docker_install>`).


Getting the |eco| source code
-----------------------------

You can use the following link to download the latest release (|eco_latest_tag|) of the
|eco| source code:

-   |eco_latest_source_zip|_.

Once downloaded, you'll need to unzip the archive.

Alternatively, if you have |git|_ installed, you can download the |eco|
repository with the following command::

    git clone https://github.com/phyletica/ecoevolity.git


Build programs from source code 
-------------------------------

Next, move into the downloaded directory. If you downloaded the release
archive, use:

.. code-block:: shell
    :substitutions:

    cd ecoevolity-|eco_latest_tag_short|

If you cloned the git repository, use:

.. code-block:: shell

    cd ecoevolity

Then, run the build script::

    ./build.sh --threads

.. note:: 

    The ``--threads`` option builds a version of |eco| that can perform
    likelihood calculations across multiple threads.   
    If you don't plan to use multithreading, you can remove this option
    (i.e., just run ``./build.sh``).
    In my opinion, you're usually better off running multiple independent MCMC
    chains rather than multithreading, but if you have plenty of CPUs to spare,
    you can do both.

If the build was successful, the |eco| programs should now be in the
``./build/bin`` directory, and you should be able to run::

    ./build/bin/ecoevolity -h

and see the |eco| help menu. The beginning of the help menu should look
something like::

    ======================================================================
                                  Ecoevolity
                      Estimating evolutionary coevality
           Version 1.1.0 (HEAD d64cc24: 2026-01-01T14:52:21-06:00)
    ======================================================================
    
    Usage: ecoevolity [OPTIONS] YAML-CONFIG-FILE
    
    Ecoevolity: Estimating evolutionary coevality
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit


Install the programs 
--------------------

If running ``ecoevolity`` to get the help menu above worked, the programs in
``./build/bin`` are ready to use, but you'll probably want to install them in
your PATH (a list of directories that your shell looks in to find the commands
you type on the command line) for easier use.

Installing locally (no admin privileges required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install the programs locally within your home folder, use::

    mkdir -p ~/.local/bin
    cp ./build/bin/* ~/.local/bin"

Then, try opening a new terminal terminal window and run ``ecoevolity -h``.
If this doesn't work, the ``~/.local/bin`` directory is likely not in your PATH
variable.
You can check what directories are in your PATH using ``echo $PATH``.
To add the ``~/.local/bin`` directory to your
PATH, you can use::

    export PATH="${PATH}:${HOME}/.local/bin"

Now, ``ecoevolity -h`` should work.

.. note::

    This update to PATH is only for your current terminal window.
    If you want this to be permanent (work for all future terminal windows),
    add the line:: 
    
        export PATH="${PATH}:${HOME}/.local/bin"
    
    to your shell's startup file. E.g.,:
    
    -   ``~/.bashrc`` file in your home directory if you are using the bash shell
    -   ``~/.zshrc`` file in your home directory if you are using the zsh shell


Installing globally (admin privileges required)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have admin privileges and want to install the programs to a location
that is almost certainly in your PATH, you can use::

    sudo cp ./build/bin/* /usr/local/bin

You can check things are working by running ``ecoevolity -h``, which should
show the help menu for ``ecoevolity``.


Build and install with one command 
----------------------------------

If you want to build and install in one go, you just need to specify where you
want the installation to go using the ``--prefix`` argument.
If you want to install the programs locally in your home directory (this does
not require admin privileges), use::

    ./build.sh --threads --prefix ~/.local


If you want to install the programs "globally" for all users of the computer
system (this requires admin privileges), use::

    sudo ./build.sh --threads --prefix /usr/local


.. _docker_install:

Using ecoevolity via Docker image
=================================

Docker provides a nice way of sharing lightweight containers that act like a
virtual machine.
We have created a Docker container with |eco| built in.
To get started, you first need to 
`install Docker <https://docs.docker.com/install/>`_.
To do so, go to `<https://docs.docker.com/install>`_ and scroll down and click on
the download button for your platform.
Once Docker is installed and running, pull down our Docker image::

    docker pull phyletica/ecoevolity-docker

.. note::

    Depending on your system and how Docker is configured, you may need to use
    ``sudo`` to run Docker commands. If you received a "permission denied"
    message when you ran the command above, try::
    
        sudo docker pull phyletica/ecoevolity-docker

This download could take several minutes depending on your internet connection.
After it completes, run and enter the docker container::

    docker run -it phyletica/ecoevolity-docker bash

.. note::

    Again, you might need to prefix this command with ``sudo``.

That's it, you are now in a virtual container with 
a fully functioning |eco| ecosystem
(|eco| and |pyco| are installed, along with example data).
Try typing::

    ecoevolity -h

This should display the |eco| help menu.

.. Next, ``cd`` into the example data directory::

..     cd ecoevolity-example-data
..     ls

.. There you will find an |eco| configuration file and nexus-formatted data files.
.. Go ahead and run an |eco| analysis::

..     ecoevolity --relax-missing-sites --relax-triallelic-sites --ignore-data ecoevolity-config.yml

To exit the container, simply type::

    exit

Docker will keep the |eco| image around, so you can always jump
back in anytime via::

    docker run -it phyletica/ecoevolity-docker bash

However, any files you created on your last visit will be gone.
So, if you want to analyze *your* data and keep the results around, ``cd``
to the directory where you want to run |eco|, then jump into
the Docker container using::

    docker run -v "$(pwd)":/portal -it phyletica/ecoevolity-docker bash

Then, once inside, type::

    cd portal
    ls

You should see the files that were in the directory on *your* computer.
Now you can run |eco| on data in this directory, and all output files will be
on your computer when you exit the container (magic!).


Using ecoevolity via Singularity image
======================================

We've also created a Singularity image you can download here:

-   |eco_latest_singularity|_

Once downloaded, you can jump on the command line within the container via:

.. code-block:: shell
    :substitutions:

    singularity shell ecoevolity-|eco_latest_tag|-linux64-intel-singularity.sif

|Eco| is preinstalled in the container, so you have access to all of it's
programs; e.g., try ``ecoevolity -h``.
Singularity mounts your home folder and current working directory into the
container by default, so you have access to the files on *your* computer from
inside the container.


Getting some example data
=========================

We have created |git| repositories containing example data and configuration
files.
You can download example files for ``ecoevolity`` using::

    git clone https://github.com/phyletica/ecoevolity-example-data.git

You can download example files for ``phycoeval`` using::

    git clone https://github.com/joaks1/phycoeval-example-data.git

If you prefer not to use |git|_, you can download an archive of:

-   `The ecoevolity example files here <https://github.com/phyletica/ecoevolity-example-data/archive/master.zip>`_.
-   `The phycoeval example files here <https://github.com/joaks1/phycoeval-example-data/archive/master.zip>`_.

.. note::

    If you are using the
    :ref:`Docker image <docker_install>`,
    these example data are included in the container.
    But, for the tutorials, it will be helpful to follow the instructions above
    to get a copy on your computer, outside of the container.


.. _pycoevolity_install:

Installing pycoevolity
======================

|Pyco|_ is a Python package for summarizing the output of |eco|.
It should work with Python 2 or 3.
If you have
`Python <https://www.python.org/>`_
and 
`pip <https://pypi.org/project/pip/>`_
installed, you can install |Pyco| via::

    pip install git+https://github.com/phyletica/pycoevolity.git@master

If this isn't working, try the
`manual installation instructions here <https://github.com/phyletica/pycoevolity>`_.

Also, |pyco| the option of using the
`R <https://www.r-project.org/>`_
packages
`ggplot2 <http://ggplot2.tidyverse.org/>`_
and
`ggridges <https://github.com/clauswilke/ggridges>`_
for creating some plots.
So, if you want this optional plotting functionality,
and you don't already have
`R <https://www.r-project.org/>`_
installed, you'll need to install it.
Once
`R <https://www.r-project.org/>`_
is in place, you can install the packages from the
`R <https://www.r-project.org/>`_
prompt using::

    install.packages(c("ggplot2", "ggridges"))
