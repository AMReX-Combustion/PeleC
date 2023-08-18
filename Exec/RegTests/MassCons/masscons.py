# ========================================================================
#
# Imports
#
# ========================================================================
import os
import numpy.testing as npt
import pandas as pd
import unittest


# ========================================================================
#
# Test definitions
#
# ========================================================================
class ConsTestCase(unittest.TestCase):
    """Tests for conservation in Pele."""

    def test_conservation(self):
        """Are mass and energy conserved?"""

        # Load the data
        fdir = os.path.abspath(".")
        fname = os.path.join(fdir, "datlog")
        df = pd.read_csv(fname, delim_whitespace=True)
        npt.assert_allclose(df.mass, df.mass[0], rtol=1e-13)
        npt.assert_allclose(df.rho_E, df.rho_E[0], rtol=1e-13)


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":
    unittest.main()
