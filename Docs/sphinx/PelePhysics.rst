
 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

 .. _PelePhysics:

PelePhysics
===========
`PelePhysics` is a repository of physics databases and implementation code for use with the `Pele` suite of of codes. There are several instances of shared physics between PeleC and PeleLM that are stored in PelePhysics. Specifically, equations of state and thermodynamic property evaluation, transport coefficients evaluation, and chemistry rate computation.

FUEGO
-----
The `PelePhysics` repository contains the FUEGO source code generation tool in order to support the inclusion
of chemical models specified in the CHEMKIN-II format.  Typically, such models include ASCII files specifying :


* chemistry (elements, chemical species, Arrhenius reaction parameters) in CHEMKIN-II format
* thermodynamics (NASA polynomial expressions)
* transport (input data for the CHEMKIN function `tranfit` to compute transport coefficients)
* ...optionally, there may be multiple source files to evaluate the species production rates using various reduced/QSS approximations

Note that the CHEMKIN-II specification assumes a mixture of ideal gases.  If the `Pele` codes are built with `Eos_dir = Fuego`, the variable `Chemistry_Model` must be set to one of the models existing in the `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models` folder (named after the folder containing the model files there).

The repository provides models for inert air, and reacting hydrogen, methane, ethane, dodecane and di-methyl ether, and more can be added, using the procedure below.  Note that although all provided models can be regenerated following these procedures, we include the resulting source files for convenience.  As a result, if the generation procedure is modified, the complete set of models will need to be updated manually.  There is currently no procedure in place to keep the models up to date with the input files in the `PelePhysics` repository.

Model Generation Procedures
~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Ensure that the environment variable `PELE_PHYSICS_HOME` is set to the root folder of your local `PelePhysics` clone containing this file.
#. Set up some FUEGO variables by doing:

   .. code-block:: bash

      . ${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/setup.sh # for bash
      source ${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/setup.csh # for csh


#. Check that all is setup correctly by reproducing one of the existing .c files, for example for the LiDryer mechanism:

   .. code-block:: bash

      cd ${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/LiDryer  
      ./make-LiDryer.sh

      
   This script uses the CHEMKIN-II thermodynamics database files expected by the usual CHEMKIN programs as well as a chemistry input file modified to include a transport section (see below) and generates a c source file.  If this file "almost" matches the one in git, you're all set... By "almost" we mean minus reordering in the egtransetEPS, egtransetSIG, egtransetDIP, egtransetPOL, egtransetZROT and egtransetNLIN functions. If not, and there are no error messages, it might be that the git version needs to be updated.  Contact `MSDay@lbl.gov` to verify the correct actions to take.

   .. note::

      The parser requires that the line endings be 'unix-style', i.e., LF only; if you have input files with 'dos-style' (CR LF) endings you may need to run a utility such as dos2unix or something like :code:`sed -i 's/^M$//' chem.inp` to clean them up. 

#. To build and use a new CHEMKIN-II based model, make a new model folder in `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models`,  say `XXX`, and copy your CHEMKIN input files there.  Additionally copy `make-LiDryer.sh` from the `LiDryer` model folder there as a template, rename it to `make-XXX.sh`, and edit the filenames at the top; this file contains the following:

   .. code-block:: bash

      CHEMINP=LiDryer.mec
      THERMINP=LiDryer.therm
      FINALFILE=LiDryer.cpp

      FMC=${PELE_PHYSICS_HOME}/Support/Fuego/Pythia/products/bin/fmc.py
      HEADERDIR=${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/header

      CHEMLK=chem.asc
      LOG=chem.log
      CHEMC=chem.c

      ${FUEGO_PYTHON} ${FMC} -mechanism=${CHEMINP} -thermo=${THERMINP}  -name=${CHEMC}
      echo Compiling ${FINALFILE}...
      cat ${CHEMC} \
      ${HEADERDIR}/header.start\
      ${HEADERDIR}/header.mec   ${CHEMINP}\
      ${HEADERDIR}/header.therm ${THERMINP}\
      ${HEADERDIR}/header.end > ${FINALFILE}
      rm -f ${CHEMC} ${CHEMLK} ${LOG}


   In addition, you must modify the chemistry input file to include a transport section. To do so, call the script located under `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models` labeled `script_trans.py` with your chem.inp and trans.dat as arguments (in that order). The script should return a text file labeled `TRANFILE_APPEND.txt`. Open this file, copy all the lines and paste them in your chemistry input file, say right above the REACTION section. Next, run the `make-XXX.sh` script file, and if all goes well, the software will generate a `chem.c` file that gets concatenated with a few others, puts everything into a result file, `YYY.c`, and cleans up its mess.

#. Add a `Make.package` text file in that model folder that will be included by the AMReX make system.  In most case, this file will contain a single line, `cEXE_sources+=YYY.c` (see the other models for examples if there are auxiliary files in languages other than c to include in the build of this model).

#. Finally, edit the `GNUmakefile` where you want to use this (in, e.g., `PeleC/Exec`) so that `Chemistry_Model` is `XXX`.  In `PeleC/Exec/Make.PeleC`, the model is expected to be in the folder `${PELE_PHYSICS_HOME}/Support/Fuego/Mechanism/Models/$(Chemistry_Model)`, and it is expected that the folder contains a `Make.package` file to include, so make sure things are where they need to be.


Equation of State
-----------------

PeleC allows the user to use different equation of state (eos) as the constitutive equation and close the compressible Navier-Stokes system of equations. All the routines needed to fully define an eos are implemented through PelePhysics module. Available models include:

* An ideal gas mixture model (similar to the CHEMKIN-II approach)
* A simple `GammaLaw` model
* Cubic models such as `Soave-Redlich-Kwong`; `Peng-Robinson` support was started but is currently stalled.

Examples of eos implementation can be seen in `PelePhysics/Eos`. The following sections will fully describe the implementation of Soave-Redlich-Kwong, a non-ideal cubic eos, for a general mixture of species. Integration with the FUEGO, for a chemical mechanism described in a chemkin format, will also be highlighted. For an advanced user interested in implementing a new eos this chapter should provide a good starting point.

Soave-Redlich-Kwong (SRK)
~~~~~~~~~~~~~~~~~~~~~~~~~

The cubic model is built on top of the ideal gas models, and is selected by specifying its name as the `Eos_dir` during the build (the `Chemistry_Model` must also be specified).  Any additional parameters (e.g., attractions, repulsions, critical states) are either included in the underlying FUEGO database used to generate the source file model implementation, or else are inferred from the input model data.


SRK EOS as a function of Pressure (p), Temperature(T), and :math:`\tau` (specific volume) is given by

.. math::
   p = R T \sum \frac{Y_k}{W_k} \frac{1}{\tau - b_m} - \frac{a_m}{\tau(\tau + b_m)}

where :math:`Y_k` are species mass fractions, :math:`R` is the universal gas constant, and
:math:`b_m` and :math:`a_m` are mixture repulsion and attraction terms, respectively.

Mixing rules
############

For a mixture of species, the following mixing rules are used to compute :math:`b_m` and :math:`a_m`.

.. math::
   a_m = \sum_{ij} Y_i Y_j \alpha_i \alpha_j \;\;\;  b_m = \sum_k Y_k b_k

where :math:`b_i` and :math:`a_i` for each species is defined using critical pressure and temperature.

.. math::
   a_i(T) = 0.42748 \frac{\left(R T_{c,i} \right)^2}{W_i^2 p_{c,i}} \bar{a}_i \left(T/T_{c,i}\right) \;\;\;
   b_i = 0.08664 \frac{R T_{c,i}}{W_i p_{c,i}}  

where

.. math::
   \bar{a}_i (T/T_{c,i}) = \left(1 + \mathcal{A} \left[ f\left( \omega_i \right) \left(1-\sqrt{T/T_{c,i}} \right ) \right] \right)^2

where :math:`\omega_i` are the accentric factors and

.. math::
   f\left( \omega_i \right) = 0.48508 + 1.5517 \omega_i - 0.151613 \omega_{i}^2

For chemically unstable species such as radicals, critical temperatures and pressures are not available.  
For species where critical properties are not available, we use the Lennard-Jones potential for that species to construct attractive and repulsive coefficients.

.. math::
   T_{c,i} = 1.316 \frac{\epsilon_i}{k_b} \;\;\;  a_i(T_{c,i}) = 5.55 \frac{\epsilon_i \sigma_i^3}{m_i^2} \;\;\;
   \mathrm{and} \;\;\; b_i = 0.855 \frac{\sigma_i^3}{m_i} 

where :math:`\sigma_i`, :math:`\epsilon_i` are the Lennard-Jones potential molecular diameter and well-depth, respectively,
:math:`m_i` the molecular mass, and :math:`k_b` is Boltzmann's constant.

In terms of implementation, a routine called `MixingRuleAmBm` can be found in the SRK eos implementation. The following code block shows the subroutine which receives species mass fractions and temperature as input. The outputs of this routine are :math:`b_m` and :math:`a_m` .

.. code-block:: fortran
		
   do i = 1, nspecies
     Tr = T*oneOverTc(i)
     amloc(i) =  (1.0d0 + Fomega(i)*(1.0d0-sqrt(Tr))) *sqrtAsti(i)

     bm = bm + massFrac(i)*Bi(i)

   enddo
   do j = 1, nspecies
      do i = 1, nspecies
        
         am = am + massFrac(i)*massFrac(j)*amloc(i)*amloc(j)
   
      end do
   end do

Thermodynamic Properties
########################

Most of the thermodynamic properties can be calculated from the equation of state and involve derivatives of various thermodynamic quantities and of EOS parameters. In the following, some of these thermodynamic properties for SRK and the corresponding routines are presented. 

Specific heat 
"""""""""""""
For computing mixture specific heat at constant volume and pressure, the ideal gas contribution and the departure from the ideal gas are computed. Specific heat at constant volume can be computed using the following

.. math::
   c_v = \left( \frac{\partial e_m}{\partial T}\right)_{\tau,Y}

For SRK EOS, the formula for :math:`c_v` reduces to

.. math::
   c_v = c_v^{id} - T \frac{\partial^2 a_m}{\partial T^2} \frac{1}{b_m} ln ( 1 + \frac{b_m}{\tau})

where :math:`c_v^{id}` is the specific heat at constant volume. Mixture specific heat at constant volume is implemented through the routine `SRK_EOS_GetMixtureCv`

.. code-block:: fortran

   subroutine SRK_EOS_GetMixtureCv(state)
   implicit none
   type (eos_t), intent(inout) :: state
   real(amrex_real) :: tau, K1

   state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

   call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

   tau = 1.0d0/state%rho

   ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
   call Calc_dAmdT(state%T,state%massFrac,state%am,state%dAmdT)
  
   ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
   call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)

   ! Ideal gas specific heat at constant volume
   call ckcvbs(state%T, state % massfrac, iwrk, rwrk, state % cv)

   ! Real gas specific heat at constant volume
   state%cv = state%cv + state%T*state%d2AmdT2* (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
  
   end subroutine SRK_EOS_GetMixtureCv
		
Specific heat at constant pressure is given by

.. math::
   
   c_p = \left( \frac{\partial h_m}{\partial T}\right)_{p,Y}   \;\; \\
   c_p =  \frac{\partial h_m}{\partial T} - \frac {\frac{\partial h}{\partial \tau}} {\frac{\partial p}{\partial \tau}} \frac{\partial p}{\partial T}

where all the derivatives in the above expression for SRK EOS are given by

.. math::

   \frac{\partial p}{\partial T} = \sum Y_k / W_k  \frac{R}{\tau-b_m} - \frac{\partial a_m}{\partial T} \frac{1}{\tau(\tau +b_m)} \\
   \frac{\partial p}{\partial \tau} = -\sum Y_k / W_k  \frac{R T}{(\tau-b_m)^2} + \frac{a_m (2 \tau + b_m)}{[\tau(\tau +b_m)]^2} \\
   \frac{\partial h_m}{\partial \tau} = -\left(T \frac{\partial a_m}{\partial T}  - a_m \right) \frac{1}{\tau(\tau+b_m)} + \frac{a_m}{(\tau+b_m)^2} -\sum Y_k / W_k  \frac{R T b_m}{(\tau-b_m)^2}  \\
   \frac{\partial h_m}{\partial T} = c_p^{id} +T \frac{\partial^2 a_m}{\partial T^2} \frac{1}{b_m} ln ( 1 + \frac{b_m}{\tau}) - \frac{\partial a_m}{\partial T} \frac{1}{\tau+b_m} +\sum Y_k / W_k  \frac{R b_m}{\tau-b_m} 

.. code-block:: fortran

    subroutine SRK_EOS_GetMixtureCp(state)
    implicit none
    type (eos_t), intent(inout) :: state
    real(amrex_real) :: tau, K1
    real(amrex_real) :var: : Cpig
    real(amrex_real) :: eosT1Denom, eosT2Denom, eosT3Denom 
    real(amrex_real) :: InvEosT1Denom,InvEosT2Denom,InvEosT3Denom
    real(amrex_real) :: dhmdT,dhmdtau
    real(amrex_real) :: Rm

    state % wbar = 1.d0 / sum(state % massfrac(:) * inv_mwt(:))

    call MixingRuleAmBm(state%T,state%massFrac,state%am,state%bm)

    tau = 1.0d0/state%rho
  
    ! Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
    call Calc_dAmdT(state%T,state%massFrac,state%dAmdT)
  
    ! Second Derivative of the EOS AM w.r.t Temperature - needed for calculating enthalpy, Cp, Cv and internal energy
    call Calc_d2AmdT2(state%T,state%massFrac,state%d2AmdT2)
  
    K1 = (1.0d0/state%bm)*log(1.0d0+state%bm/tau)
    
    eosT1Denom = tau-state%bm
    eosT2Denom = tau*(tau+state%bm)
    eosT3Denom = tau+state%bm

    InvEosT1Denom = 1.0d0/eosT1Denom
    InvEosT2Denom = 1.0d0/eosT2Denom
    InvEosT3Denom = 1.0d0/eosT3Denom

    Rm = (Ru/state%wbar)
  
    ! Derivative of Pressure w.r.t to Temperature
    state%dPdT = Rm*InvEosT1Denom - state%dAmdT*InvEosT2Denom

    ! Derivative of Pressure w.r.t to tau (specific volume)
    state%dpdtau = -Rm*state%T*InvEosT1Denom*InvEosT1Denom + state%am*(2.0*tau+state%bm)*InvEosT2Denom*InvEosT2Denom

    ! Ideal gas specific heat at constant pressure
    call ckcpbs(state % T, state % massfrac, iwrk, rwrk,Cpig)
  
    ! Derivative of enthalpy w.r.t to Temperature
    dhmdT = Cpig + state%T*state%d2AmdT2*K1 - state%dAmdT*InvEosT3Denom + Rm*state%bm*InvEosT1Denom
  
    ! Derivative of enthalpy w.r.t to tau (specific volume)
    dhmdtau = -(state%T*state%dAmdT - state%am)*InvEosT2Denom + state%am*InvEosT3Denom*InvEosT3Denom - &
       Rm*state%T*state%bm*InvEosT1Denom*InvEosT1Denom

    ! Real gas specific heat at constant pressure
    state%cp = dhmdT - (dhmdtau/state%dpdtau)*state%dPdT
  
    end subroutine SRK_EOS_GetMixtureCp

Internal energy and Enthalpy 
""""""""""""""""""""""""""""

Similarly mixture internal energy for SRK EOS is given by

.. math::
   e_m = \sum_k Y_k e_k^{id} + \left( T  \frac{\partial a_m}{\partial T}  - a_m \right)\frac{1}{b_m} ln \left( 1 + \frac{b_m}{\tau}\right)

and mixture enthalpy :math:`h_m = e + p \tau`

.. math::
   h_m = \sum_k Y_k h_k^{id} + \left ( T \frac{\partial a_m}{\partial T} - a_m \right) \frac{1}{b_m} \ln \left( 1 + \frac{b_m}{\tau}\right) + R T \sum \frac{Y_k}{W_k} \frac{b_m}{\tau -b_m} - \frac{a_m}{\tau + b_m}

and the implementation can be found in the routine `SRK_EOS_GetMixture_H`.

Speed of Sound
""""""""""""""

The sound speed for SRK EOS is given by

.. math::

   a^2 = -\frac{c_p}{c_v} \tau^2  \frac{\partial p}{\partial \tau}

Species enthalpy
""""""""""""""""

For computation of kinetics and transport fluxes we will also need the species partial enthalpies and the chemical potential.  The species enthalpies for SRK EOS are given by

.. math::
   
   h_k = \frac{\partial h_m}{\partial Y_k } - \frac {\frac{\partial h}{\partial \tau}} {\frac{\partial p}{\partial \tau}} \frac{\partial p}{\partial Y_k}

where

.. math::
   \frac{\partial h_m}{\partial Y_k } &=  h_k^{id} + (T \frac{\partial^2 a_m}{\partial T \partial Y_k}  - \frac{\partial a_m }{\partial Y_k}) \frac{1}{b_m} \ln\left(1+ \frac{b_m}{\tau}\right) \\&-\left(T \frac{\partial a_m}{\partial T}  - a_m \right) \left[ \frac{1}{b_m^2} \ln\left(1+ \frac{b_m}{\tau}\right) - \frac{1}{b_m(\tau+b_m)} \right ] \frac{\partial b_m}{\partial Y_k} \nonumber \\&+ \frac{a_m}{(\tau+b_m)^2}  \frac{\partial b_m}{\partial Y_k} - \frac{1}{\tau+b_m}  \frac{\partial a_m}{\partial Y_k} + 1 / W_k  \frac{R T b_m}{\tau-b_m}\\&+\sum_i \frac{Y_i}{W_i} R T \left( \frac{1}{\tau -b_m} + \frac{b_m}{(\tau-b_m)^2} \right) \frac{ \partial b_m}{\partial Y_k} 

.. math::
   
   \frac{\partial p}{\partial Y_k} &= R T \frac{1}{W_k} \frac{1}{\tau - b_m} - \frac{\partial a_m}{\partial Y_k} \frac{1}{\tau(\tau + b_m)} \\&+\left(R T \sum \frac{Y_i}{W_i} \frac{1}{(\tau - b_m)^2} + \frac{a_m}{\tau(\tau + b_m)^2} \right ) \frac{\partial b_m}{\partial Y_k}

Chemical potential 
""""""""""""""""""
The chemical potentials are the derivative of the free energy with respect to composition.  Here the free energy `f`` is given by

.. math::
   f &= \sum_i Y_i (e_i^{id} - T s_i^{id,*}) +  \sum_i \frac{Y_i R T}{W_i} ln (\frac{Y_i R T}{W_i \tau p^{st}})  \nonumber \\ &+ \sum_i \frac{Y_i R T}{W_i} ln (\frac{\tau}{\tau-b_m}) -  a_m \frac{1}{b_m}ln (1+ \frac{b_m}{\tau})  \nonumber \\ &= \sum_i Y_i (e_i^{id} - T s_i^{id,*}) +  \sum_i \frac{Y_i R T}{W_i} ln (\frac{Y_i R T}{W_i (\tau-b_m) p^{st}} )- a_m \frac{1}{b_m} ln (1+ \frac{b_m}{\tau})  \nonumber 

Then

.. math::
   
   \mu_k &= \frac{\partial f}{\partial Y_k} = e_k^{id} - T s_k^{id,*}  + \frac{RT}{W_k} ln (\frac{Y_k R T}{W_k (\tau-b_m) p^{st}}) + \frac{RT}{W_k} +  \frac{RT}{\bar{W}} \frac{1}{\tau-b_m} \frac {\partial b_m}{\partial Y_k} \nonumber \\
   &- \frac{1}{b_m} ln(1 + \frac{b_m}{\tau}) \frac{\partial a_m}{\partial Y_k}+ \frac{a_m}{b_m^2} ln(1 + \frac{b_m}{\tau}) \frac{\partial b_m}{\partial Y_k}- \frac{a_m}{b_m} \frac{1}{\tau+b_m} \frac{\partial b_m}{\partial Y_k}

Other primitive variable derivatives
""""""""""""""""""""""""""""""""""""

The Godunov (FV) algorithm also needs some derivatives to express source terms in terms of primitive variables. In particular one needs

.. math::

   \left . \frac{\partial p}{\partial \rho} \right|_{e,Y} =-\tau^2 \left( \frac{\partial p}{\partial \tau}- \frac {\frac{\partial e}{\partial \tau}} {\frac{\partial e}{\partial T}} \frac{\partial p}{\partial T} \right )

and

.. math::

   \left . \frac{\partial p}{\partial e} \right|_{\rho,Y} = \frac{1}{c_v} \frac{\partial p}{\partial T}

All of the terms needed to evaluate this quantity are known except for

.. math::

   \frac{\partial e}{\partial \tau} = \frac{1}{\tau ( \tau + b_m)} \left( a_m - T  \frac{\partial a_m}{\partial T}  \right) \;\; .




Transport
---------
.. note:: Placeholder, to be written


Chemistry
---------

.. note:: Placeholder, to be written


Usage
-----

This section contains an evolving set of usage notes for PelePhysics.


1. In the FUEGO version of the chemistry ODE solver (`Fuego/actual_reactor.F90`)  the algorithm will attempt to reuse the Jacobian from the previous cell; this can improve performance significantly for many problems. However, this can cause slight differences in the solution when the cells are visited in a different order. While not physically significant these differences can make it difficult to debug potential problems when trying to verify correct operation of code changes that impact the order in which cells are visited such as loop reordering or using tiling. To assist in this process the re-use can be switched off by setting the following flag as part of the extern namelist in the relevant probin file: ::

     &extern
       new_Jacobian_each_cell = 1
     /
     
  When this namelist variable is set PelePhysics will treat each cell as a brand new problem and therefore be independent of the order cells are visited.
