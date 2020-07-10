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
class MMSTestCase(unittest.TestCase):
    """Tests for symmetry in Pele."""

    def test_symmetry(self):
        """Is the MMS error symmetric with symmetric IC (u,v,w)?"""

        # Load the data
        fdir = os.path.abspath(".")
        fname = os.path.join(fdir, "mmslog")
        df = pd.read_csv(fname, delim_whitespace=True)
        npt.assert_allclose(df.u_mms_err, df.v_mms_err, rtol=1e-13)
        npt.assert_allclose(df.u_mms_err, df.w_mms_err, rtol=1e-13)
        npt.assert_allclose(df.rhou_residual, df.rhov_residual, rtol=1e-13)
        npt.assert_allclose(df.rhou_residual, df.rhow_residual, rtol=1e-13)


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":
    unittest.main()
