"README"
#Xrecon - External Reconstruction

This is Xrecon, a collection of magnetic resonance image reconstruction
routines tuned to read data acquired on Varian MRI Systems.
It has been developed to work externally and independently of the Varian data
acquisition software VnmrJ.

The routines provided here link against the GNU Scientific Library (GSL) and
the Fastest Fourier Transform in the West (FFTW) library of subroutines.
The GSL is free software published under the GNU General Public License (GPL).
The FFTW is free software published under the GNU General Public License (GPL).
As such it is a requirement that the routines provided here are also published
under the GPL.

To publish under the GPL two elements must be added to each source file:
1. A copyright notice with the names of all contributors.
2. A statement of copying permission saying the program is distributed under
   the terms of the GNU General Public License.

A copy of the license itself should also be included.

Xrecon is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Xrecon is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

##Travis-CI
[![Build Status](https://travis-ci.org/timburrow/xrecon.svg?branch=master)](https://travis-ci.org/timburrow/xrecon)

##Installation

Please consult the INSTALL file in this distribution for instructions.


##Files

COPYING         The GPL
INSTALL         Installation instructions
README          This file
Makefile        Xrecon makefile
data.h          Stripped down Agilent data file handler
Xrecon.c        Main External Reconstruction
Xrecon.h        Global includes and defines

###recon1D
recon1D.c       1D recon
default1D.c     Default 1D recon
profile1D.c     1D profile recon

###recon2D
recon2D.c       2D recon
default2D.c     Default 2D recon
asltest2D.c     2D ASL test mode
proj2D.c        2D projection recon

###recon3D
recon3D.c       3D recon
default3D.c     Default 3D recon
multiblock3D.c  Default multiblock 3D recon

###reconEPI
reconEPI.c      EPI recon
defaultEPI.c    Default EPI recon
dprocEPI.c      EPI data processing
prescanEPI.c    EPI prescan recon

##reconCSI
reconCSI.c      CSI recon
reconCSI2D.c    2D CSI recon
reconCSI3D.c    3D CSI recon
dprocCSI.c      CSI data processing

###common1D
dproc1D.c       1D data processing
dread1D.c       1D data read
dutils1D.c      1D data utilities
write1D.c       1D phasefile and datafaile writing

###common2D
dmask2D.c       2D data masking
dproc2D.c       2D data processing
dread2D.c       2D data read
dutils2D.c      2D data utilities
noise2D.c       2D routines for noise matrix measurements
fdfwrite2D.c    2D fdf writing
rawwrite2D.c    2D raw binary writing
tifwrite2D.c    2D TIFF writing

###common3D
dproc3D.c       3D data processing
fdfwrite3D.c    3D fdf writing
rawIO3D.c       3D raw binary writing

###common
dhead.c         Data header routines
dproc.c         Data processing
dread.c         Data read
dutils.c        Data utilities
options.c       Input options
pars.c          Parameter read routines
utils.c         General utilities

###nifti
niftiwrite.c    NIFTI-1/Analyze7.5 writing
nifti1.h        NIFTI-1 header


##Contacts

Paul Kinchesh,
Agilent Technologies, 6 Mead Road, Oxford Industrial Park, Oxford OX5 1QU, UK
~~paul.kinchesh@agilent.com~~


##Bug Reports

~~Paul Kinchesh, paul.kinchesh@agilent.com~~


##Contributors

Paul Kinchesh, Magnetic Resonance Systems, Agilent Technologies.
Martyn Klassen, Robarts Research Institute, University of Western Ontario.
Alexandr Khrapichev, Department of Physiology, Oxford University.
Margaret Kritzer, Magnetic Resonance Systems, Agilent Technologies.
Lana Kaiser, Magnetic Resonance Systems, Agilent Technologies.


