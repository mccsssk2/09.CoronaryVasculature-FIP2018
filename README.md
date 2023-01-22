# Repository.  

Human coronary vasculature. 2018.

# Background.  
 
This code generates a 3D coronary vasculature aligned with potentially imaged cardiac geometry, and
provided blood flow distribution from perfusion CT imaging. This modular code is prepared for application
in imaging data interpretation.  

The article related to this code is:  
Kharche et al. https://www.frontiersin.org/articles/10.3389/fphys.2018.00511/full  

# Dependencies.

This serial code requires GNU C.

# Install.

None. The driver function is srkVasculature.c .

# Sources and data description.

## Data.  

The data that drove this study was clinical resolution perfusion CT imaging of the heart, primarily the left ventricle.
A MATLAB code to stack the registered CT slices into a 3D volume is provided. Whereas fractal dimension (FD, box dimension)
is an established measure of heterogeniety, this code provides the means to estimate FD in the imaging.

## Source.  

The source code uses following algorithms and methods to simulate steady state blood flow in the human heart:  
* The vasculature is treated as a binary tree where each node may give rise to up to 2 daughter vessels.  
* The vasculature is distributed in a 3D volume representing the heart ventricles (see generator codes and example VTKs). To do so,
individual nodes are assigned coordinates that keep them inside the volume while distance from nearest neighbours is maximized.  
* Steady state blood flow and metrics such as vessel resistance and pressure are calculated using  Poiseuille's law.
* Output is written to ASCII files and to VTK files to assist with post processing.  
* The means to do non-interactive postprocessing using VTK and custom codes is provided.  

# Use.

The driver function in the source is srkVasculature.c. The user can replace the morphometry and connectivty to use
specific experimental data. A geometry is generated, in this case idealized cardiac ventricles. Compiling is using
the provided makefiles. Since the optimization is serial and takes a long time, multiple instances are generated
to allow averaging over an ensemble. The program outputs the optimial structure with data for blood flow, pressure,
resistance, and other metrics. Code for postprocessing is provided. The end result is an indicative blood flow
distribution under imposed conditions, e.g. stenosis of LAD.  

# Maintainer.

These codes are maintained by SR Kharche. To find wider applicability, it can be further developed and provided as
precompiled binaries for users.

# Acknowledements.

This project was generously funded by the Kidney Unit in Lawson, London Ontario, Canada, Compute Canada, and Western University.

# Licence.

BSD 3-Clause License

Copyright (c) 2023, Sanjay R. Kharche, Ph.D.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Dr. Sanjay R. Kharche, Ph.D. (App. Maths).  
January 23rd, 2023.  
Email: Sanjay.Kharche@lhsc.on.ca , skharche@uwo.ca , srk.london@yahoo.com .  
Phone: +15198780685 (mobile).  
Website: https://kcru.lawsonresearch.ca/research/srk/index.html  

