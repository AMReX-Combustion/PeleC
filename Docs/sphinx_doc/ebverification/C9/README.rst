C9. Acoustic wave in cylindrical channel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This is the case of an acoustic wave propagating in a cylindrical
channel. The geometry is a :math:`x`-direction aligned circular
channel with periodic boundary conditions. The acoustic wave
propagates along the :math:`x` direction.

The density pulse at the center of the channel is defined as:

.. math::
   \rho(x) = \rho_0 + \alpha \exp(-(x/\sigma)^2)

The background pressure is set to :math:`100000 erg/cm^3`, the
background density is set to :math:`\rho_0 = 0.0014 g/cm^3`,
:math:`\alpha=10^{-4} g/cm^3`, :math:`\sigma=10cm`, and the background
velocity is set to 0. The length of the channel is 100cm and the
radius is 5cm. The simulations are performed until the acoustic pulse
comes back to the center of the channel, :math:`t=0.01s`.

.. image:: ./ebverification/C9/pulse.png
   :height: 200pt


Density profiles in the centerline at t=0.01s
#############################################

.. image:: ./ebverification/C9/rho.png
   :height: 300pt

:math:`L_2` error norm of velocity
##################################

The :math:`L_2` error norm for a quantity :math:`s` is defined as

.. math::
   e_s = \sqrt{ \frac{\int_{-L/2}^{L} (s^h-s^*)^2 \mathrm{d}r}{n_x}}

where :math:`s^h` is the numerical solution, :math:`s^*` is the exact
solution, and :math:`n_x` is the number of cells in the
:math:`x`-direction.

.. image:: ./ebverification/C9/error.png
   :height: 300pt
