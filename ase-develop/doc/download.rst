.. _download_and_install:

=========================
Installation requirements
=========================

The following packages are required for basic ASE functionality:

1) Python_ 2.4 - 2.7.
2) NumPy_.

.. _Python: http://www.python.org
.. _NumPy: http://www.scipy.org/NumPy

It is highly recommended (but not required) to install also these two:

3) matplotlib_.
4) pygtk_.

Matplotlib is needed for :mod:`writing png and eps files <io>`, and
both packages are needed for ASE's simple GUI (called **ag**, see :mod:`gui`).
Some of these packages may already be installed on your system.

.. _matplotlib: http://matplotlib.sourceforge.net
.. _pygtk: http://www.pygtk.org


Specific information for different operating systems 
is provided at :ref:`installation`.

.. _download:

========
Download
========

.. highlight:: bash

.. _latest_stable_release:

Latest stable release
=====================

The latest stable release can be obtained from ``svn`` or as a ``tarball``.

.. note::

   The recommended installation path is :envvar:`$HOME`.

When using svn please set the following variable:

- bash::

   export ASE_TAGS=https://svn.fysik.dtu.dk/projects/ase/tags/

- csh/tcsh::

   setenv ASE_TAGS https://svn.fysik.dtu.dk/projects/ase/tags/

======= =========== ============================================ =============================
Release Date        Retrieve as svn checkout                     Retrieve as tarball
======= =========== ============================================ =============================
 3.6.0_ Feb 24 2012 ``svn co -r 2515 $ASE_TAGS/3.6.0 ase-3.6.0`` python-ase-3.6.0.2515.tar.gz_
 3.5.1_ May 24 2011 ``svn co -r 2175 $ASE_TAGS/3.5.1 ase-3.5.1`` python-ase-3.5.1.2175.tar.gz_
 3.4.1_ Aug 11 2010 ``svn co -r 1765 $ASE_TAGS/3.4.1 ase-3.4.1`` python-ase-3.4.1.1765.tar.gz_
 3.4.0_ Apr 23 2010 ``svn co -r 1574 $ASE_TAGS/3.4.0 ase-3.4.0`` python-ase-3.4.0.1574.tar.gz_
 3.3.1_ Jan 20 2010 ``svn co -r 1390 $ASE_TAGS/3.3.1 ase-3.3.1`` python-ase-3.3.1.1390.tar.gz_
 3.2.0_ Sep 4 2009  ``svn co -r 1121 $ASE_TAGS/3.2.0 ase-3.2.0`` python-ase-3.2.0.1121.tar.gz_
 3.1.0_ Mar 27 2009 ``svn co -r 846 $ASE_TAGS/3.1.0 ase-3.1.0``  python-ase-3.1.0.846.tar.gz_
 3.0.0_ Nov 13 2008 ``svn co -r 657 $ASE_TAGS/3.0.0 ase-3.0.0``  python-ase-3.0.0.657.tar.gz_
======= =========== ============================================ =============================

.. _3.6.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.6.0

.. _python-ase-3.6.0.2515.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.6.0.2515.tar.gz

.. _3.5.1:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.5.1

.. _python-ase-3.5.1.2175.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.5.1.2175.tar.gz

.. _3.4.1:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.4.1

.. _python-ase-3.4.1.1765.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.4.1.1765.tar.gz

.. _3.4.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.4.0

.. _python-ase-3.4.0.1574.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.4.0.1574.tar.gz

.. _3.3.1:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.3.1

.. _python-ase-3.3.1.1390.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.3.1.1390.tar.gz

.. _3.2.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.2.0

.. _python-ase-3.2.0.1121.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.2.0.1121.tar.gz

.. _3.1.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.1.0

.. _python-ase-3.1.0.846.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.1.0.846.tar.gz

.. _3.0.0:
    https://trac.fysik.dtu.dk/projects/ase/browser/tags/3.0.0

.. _python-ase-3.0.0.657.tar.gz:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.0.0.657.tar.gz



.. _latest_development_release:

Latest development release
==========================

The latest revision can be obtained like this::

  $ svn checkout https://svn.fysik.dtu.dk/projects/ase/trunk ase

or from the daily snapshot: `<python-ase-snapshot.tar.gz>`_.

.. note::

   The recommended checkout path is :envvar:`$HOME`.

.. _installation:

============
Installation
============

After performing the installation do not forget to :ref:`running_tests`!

.. _pm_installation:

Installation with package manager on Linux
==========================================

Install the binaries with the software package manager of your Linux distribution.

This is **the preferred** way to install on a Linux system.

If you prefer to install from sources follow :ref:`manual_installation`.

The currently supported systems include (issue the commands below **as root**):

- RHEL/CentOS 6::

    yum install wget
    cd /etc/yum.repos.d/
    wget http://download.opensuse.org/repositories/home:/dtufys/CentOS_CentOS-6/home:dtufys.repo
    yum install python-ase
    yum install python-matplotlib # optionally

- Fedora 17::

    yum install wget
    cd /etc/yum.repos.d/
    wget http://download.opensuse.org/repositories/home:/dtufys/Fedora_17/home:dtufys.repo
    yum install python-ase
    yum install python-matplotlib # optionally

- openSUSE 12.2::

    zypper ar -f http://download.opensuse.org/repositories/home:/dtufys/openSUSE_12.2/home:dtufys.repo
    yast -i python-ase
    yast -i python-matplotlib # optionally

- Debian 6.0::

    sudo bash -c 'echo "deb http://widehat.opensuse.org/repositories/home:/dtufys/Debian_6.0 /" > /etc/apt/sources.list.d/home_dtufys.sources.list'
    wget http://widehat.opensuse.org/repositories/home:/dtufys/Debian_6.0/Release.key && sudo apt-key add Release.key && rm Release.key
    sudo apt-get update
    sudo apt-get install python-ase
    sudo apt-get install python-matplotlib # optionally

- Ubuntu 12.04::

    sudo bash -c 'echo "deb http://widehat.opensuse.org/repositories/home:/dtufys/xUbuntu_12.04 /" > /etc/apt/sources.list.d/home_dtufys.sources.list'
    wget http://widehat.opensuse.org/repositories/home:/dtufys/xUbuntu_12.04/Release.key && sudo apt-key add Release.key && rm Release.key
    sudo apt-get update
    sudo apt-get install python-ase
    sudo apt-get install python-matplotlib # optionally

  .. note::

    Alternative packages for ubuntu are provided at
    `Ubuntu package <https://wiki.fysik.dtu.dk/gpaw/install/Linux/Ubuntu_ppa.html#ubuntupackage>`_.

For the full list of supported distributions check
https://build.opensuse.org/package/show?package=python-ase&project=home%3Adtufys

.. note::

   Explore the repositories - more software packages are available!

.. note::

   If you prefer to install manually, proceed to :ref:`manual_installation`, or
   alternatively, manually unpack the RPMS, e.g.::

     # download the packages + dependencies (you can do that also manually!)
     $ yumdownloader --resolve python-ase
     # unpack into the current directory
     $ find . -name "*.rpm" | xargs -t -I file sh -c "rpm2cpio file | cpio -idm"
     # modify profile.d environment scripts
     $ find . -name "*.*sh" | xargs -t -I file sh -c 'sed -i "s#PA=/usr#PA=$PWD/usr#" file'
     # modify environment modules scripts
     $ find . -name "*.modules" | xargs -t -I file sh -c 'sed -i "s# /usr# $PWD/usr#" file'
     # make scripts executable
     $ find . -name "*.*sh" | xargs -t -I file sh -c "chmod u+x file"
     # source the scripts (example for bash)
     $ for f in `find . -name "*.sh"`; do source $f; done
     # verify the desired installation location is used
     $ python -c "import ase; print ase.__file__"

   This method works for all the RPM packages from the repository (like gpaw),
   however not for the external, distribution provided packages,
   which may require manually creating the environment scripts.


OSX
===

For Apple users, the MacPorts_ Project provides a straight-forward
route to obtain all necessary requirements. Unfortunately, MacPorts
does not install the `gtk` bindings to matplotlib_ by default, which
are required to open the GUI. To get all the ASE prerequisites for
python 2.7 in one single command anyway, install MacPorts and then run::

  $ port install py27-matplotlib +gtk2

Use the `sudo` command if you have root access and if you require 
a system-wide install. Once finished, please follow :ref:`manual_installation`.

.. _MacPorts: http://www.macports.org/

Windows
=======

.. note::

   ASE is not yet fully functional on Windows!
   https://trac.fysik.dtu.dk/projects/ase/ticket/62

On Windows the following packages need to installed.
On the command prompt:

.. note:: installation assumes the python TARGETDIR `C:\\Python27` ,
          leave also the default `C:\\Program Files\\pythonxy` .

-  pythonxy_. Download the `exe` installer and install with::

     Python(x,y)-2.7.2.2.exe /Log="%TMP%\pythonxy_install.log" /S

.. note::

   Open Task Manager and control when the process in finished.

- pygtk_win32_. Download the `msi` **pygtk-all-in-one** installer.
  Specify the correct TARGETDIR and install::

     pygtk-all-in-one-2.24.2.win32-py2.7.msi TARGETDIR="%HOMEDRIVE%\Python27" ALLUSERS=1 /l*vx "%TMP%\pygtk_install.log" /passive

.. note::

   If performing clicking-installation make sure that the default
   python Windows TARGETDIR is selected.

- Download the python-ase-win32.msi_ installer and install with::

     python-ase-X.X.X.win32.msi /l*vx "%TMP%\python-ase_install.log" /passive

.. note::

   You can build the `msi` ASE package on Windows with::

      python setup.py bdist_msi

   The `msi` package will be created under the `dist` directory.

.. _pythonxy: http://code.google.com/p/pythonxy
.. _pygtk_win32: http://ftp.gnome.org/pub/GNOME/binaries/win32/pygtk/2.24/

.. _python-ase-win32.msi:
    https://wiki.fysik.dtu.dk/ase-files/python-ase-3.6.1.2627.win32.msi

.. _manual_installation:

Manual installation
===================

After the :ref:`download` of ASE source create the link
to the requested version, e.g.:

- if retrieved from ``svn``::

   $ cd $HOME
   $ ln -s ase-3.5.1 ase
    
- if retrieved as ``tarball``::

   $ cd $HOME
   $ tar zxf python-ase-3.5.1.2175.tar.gz
   $ ln -s python-ase-3.5.1.2175 ase

It is sufficient to
put the directory :file:`$HOME/ase` in your :envvar:`PYTHONPATH`
environment variable, and the directory :file:`$HOME/ase/tools` in
your :envvar:`PATH` environment variable.  Do this permanently in
your :file:`~/.bashrc` file::

  export PYTHONPATH=$HOME/ase:$PYTHONPATH
  export PATH=$HOME/ase/tools:$PATH

or your :file:`~/.cshrc` file::

  setenv PYTHONPATH ${HOME}/ase:${PYTHONPATH}
  setenv PATH ${HOME}/ase/tools:${PATH}

Instead of :envvar:`HOME`, you may use any other directory.

.. index:: test

Optional, **NOT** recommended way of installing ASE system-wide is::

  $ cd ase
  $ sudo python setup.py install

This is one of the best ways to ruin a Linux system.

.. _running_tests:

Run the tests
=============

Make sure that everything works by running the :mod:`test
suite <test>`.  This will create many files, so run the tests in a new
directory (preferably using bash)::

  $ bash
  $ mkdir /tmp/testase.$$; cd /tmp/testase.*
  $ testase.py 2>&1 | tee testase.log

.. note:: 

   In the development version of ASE, and in future stable versions,
   the test script is just named :file:`testase`.

.. note::

   The last test :trac:`ase/test/COCu111.py` requires closing
   the graphics windows to terminate the whole test-suite.

If any of the tests fail,
then please send us :file:`testase.log` (see :ref:`bugs`).

.. note::

   If matplotlib_ or pygtk_ is not installed, one of the tests will
   fail - avoid this with::

     $ testase.py --no-display

Video tutorial
==============

In the video: :ref:`overview` of the features of ASE,
followed by a :ref:`manual_installation` of ASE on a Linux system.

.. note::

   Use "Right Click -> Play" to play.

.. raw:: html

        <p></p>
        <object width="800" height="600">
        <embed src="https://wiki.fysik.dtu.dk/ase-files/oi_en_800x600.swf"
        type="application/x-shockwave-flash"
        allowFullScreen="false"
        allowscriptaccess="never"
        loop="false"
        play="false"
        width="800" height="600">
        <p></p>
        Video not playing? Download avi <a href="https://wiki.fysik.dtu.dk/ase-files/oi_en.avi">file</a> instead.
        </embed></object>
        <p></p>
