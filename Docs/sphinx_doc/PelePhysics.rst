
 .. role:: cpp(code)
    :language: c++
 
 .. role:: fortran(code)
    :language: fortran

 .. _BCs:

PelePhysics
-----------

Equation of State
^^^^^^^^^^^^^^^^^
PeleC allows the user to use different equation of state (eos) as the constitutive equation and close the compressible Navier-Stokes system of equations. All the routines needed to fully define an eos are implemented through PelePhysics module. Examples of eos implementation can be seen in `PelePhysics/Eos`. The following sections will fully describe the implementation of Soave-Redlich-Kwong, a non-ideal cubic eos, for a general mixture of species. Integration with the Fuego, for a chemical mechanism described in a chemkin format, will also be highlighted. For an advanced user interested in implementing a new eos this chapter should provide a good starting point.

Soave-Redlich-Kwong (SRK)
^^^^^^^^^^^^^^^^^^^^^^^^^
SRK EOS as a function of Pressure (p), Temperature(T), and :math:`\tau` (specific volume) is given by

.. math::
   p = R T \sum \frac{Y_k}{W_k} \frac{1}{\tau - b_m} - \frac{a_m}{\tau(\tau + b_m)}

where :math:`Y_k` are species mass fractions, :math:`R` is the universal gas constant, and
:math:`b_m` and :math:`a_m` are mixture repulsion and attraction terms, respectively.

Mixing rules
""""""""""""
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
For species where critical properties are not available, we use the Lennard-Jones potential for that species to construct attractive and replusive coefficients.

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

Derivatives of EOS
""""""""""""""""""



Thermodynamic Properties
""""""""""""""""""""""""

Specific heat 
"""""""""""""
For computing mixture specific heat at constant volume and pressure, the ideal gas contribution and the departure from the ideal gas are computed. Specific heat at constant volume can be computed using the following

.. math::
   c_v = \left( \frac{\partial e_m}{\partial T}\right)_{\tau,Y}

For SRK EOS, the formula for :math:`c_v` reduces to

.. math::
   c_v = c_v^{id} - T \frac{\partial^2 a_m}{\partial T^2} \frac{1}{b_m} ln ( 1 + \frac{b_m}{\tau})

where :math:`c_v^{id}` is the specific heat at constant volume. 
   
	

