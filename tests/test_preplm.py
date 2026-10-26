import math
import unittest

from cmaq_port.preplm import preplm
from cmaq_port.plume import real32


def _expected_wind(u: float, v: float) -> float:
    """Replicate the source's binary32 wind-speed arithmetic."""
    u = real32(u)
    v = real32(v)
    return real32(math.sqrt(real32(real32(u * u) + real32(v * v))))


def _predict_layer(value: float, zf) -> int:
    """First 1-based layer L whose top face zf[L-1] (0-based) >= value."""
    for layer in range(1, len(zf) + 1):
        if zf[layer - 1] >= value:
            return layer
    return len(zf)


def _base_kwargs() -> dict:
    """A fully valid, physically reasonable non-fire configuration."""
    emlays = 6
    zf = (100.0, 200.0, 300.0, 400.0, 500.0, 600.0)
    zh = (50.0, 150.0, 250.0, 350.0, 450.0, 550.0)
    ta = tuple(real32(289.5 - 0.01 * (z - 50.0)) for z in zh)
    uw = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    vw = (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    qv = (0.005,) * emlays
    ddzf = (0.01,) * emlays
    pres = (0.0, 1000.0, 900.0, 800.0, 700.0, 600.0, 500.0)
    return dict(
        fireflg=False,
        emlays=emlays,
        hmix=450.0,
        hts=250.0,
        psfc=1000.0,
        ts=290.0,
        ddzf=ddzf,
        qv=qv,
        ta=ta,
        uw=uw,
        vw=vw,
        zh=zh,
        zf=zf,
        pres=pres,
    )


class PreplmWindSpeedTests(unittest.TestCase):
    def test_wind_speed_per_layer(self) -> None:
        kwargs = _base_kwargs()
        uw = (3.0, 6.0, 5.0, 8.0, 0.0, 9.0)
        vw = (4.0, 8.0, 12.0, 15.0, 7.0, 12.0)
        kwargs["uw"] = uw
        kwargs["vw"] = vw
        result = preplm(**kwargs)
        self.assertEqual(len(result.wind_speed), kwargs["emlays"])
        for layer in range(kwargs["emlays"]):
            self.assertAlmostEqual(
                result.wind_speed[layer],
                _expected_wind(uw[layer], vw[layer]),
                places=5,
            )


class PreplmFireBranchTests(unittest.TestCase):
    def test_fire_bypasses_interpolation(self) -> None:
        kwargs = _base_kwargs()
        kwargs["fireflg"] = True
        # ts is far outside the TA profile; if POLY were used, the stack
        # temperature would land near the ~285-290 profile, not at ts.
        kwargs["ts"] = 400.0
        result = preplm(**kwargs)
        self.assertEqual(result.stack_temperature, real32(400.0))
        self.assertEqual(result.stack_wind, result.wind_speed[0])

    def test_fire_allows_fewer_layers_than_poly(self) -> None:
        # emlays=2 is below the non-fire POLY minimum of 4, proving the fire
        # path never requires the interpolation points.
        emlays = 2
        zf = (100.0, 200.0)
        zh = (50.0, 150.0)
        result = preplm(
            fireflg=True,
            emlays=emlays,
            hmix=150.0,
            hts=80.0,
            psfc=1000.0,
            ts=295.0,
            ddzf=(0.01, 0.01),
            qv=(0.005, 0.005),
            ta=(288.0, 287.0),
            uw=(3.0, 4.0),
            vw=(4.0, 3.0),
            zh=zh,
            zf=zf,
            pres=(0.0, 1000.0, 900.0),
        )
        self.assertEqual(result.stack_temperature, real32(295.0))
        self.assertEqual(result.stack_wind, result.wind_speed[0])


class PreplmLayerIndexTests(unittest.TestCase):
    def test_lstk_lpbl_configurations(self) -> None:
        zf = _base_kwargs()["zf"]
        for hts, hmix in [(250.0, 450.0), (150.0, 550.0), (350.0, 250.0)]:
            with self.subTest(hts=hts, hmix=hmix):
                kwargs = _base_kwargs()
                kwargs["hts"] = hts
                kwargs["hmix"] = hmix
                result = preplm(**kwargs)
                self.assertEqual(result.lstk, _predict_layer(hts, zf))
                self.assertEqual(result.lpbl, _predict_layer(hmix, zf))


class PreplmGradientTests(unittest.TestCase):
    def test_first_gradient_overridden_by_second(self) -> None:
        result = preplm(**_base_kwargs())
        self.assertEqual(result.theta_gradient[0], result.theta_gradient[1])

    def test_gradient_length_matches_layers(self) -> None:
        kwargs = _base_kwargs()
        result = preplm(**kwargs)
        self.assertEqual(len(result.theta_gradient), kwargs["emlays"])


class PreplmGuardTests(unittest.TestCase):
    def test_emlays_below_two(self) -> None:
        kwargs = _base_kwargs()
        kwargs["emlays"] = 1
        with self.assertRaises(ValueError):
            preplm(**kwargs)

    def test_non_fire_emlays_below_four(self) -> None:
        emlays = 3
        with self.assertRaises(ValueError):
            preplm(
                fireflg=False,
                emlays=emlays,
                hmix=250.0,
                hts=150.0,
                psfc=1000.0,
                ts=290.0,
                ddzf=(0.01,) * emlays,
                qv=(0.005,) * emlays,
                ta=(288.0, 287.0, 286.0),
                uw=(3.0, 4.0, 5.0),
                vw=(4.0, 3.0, 2.0),
                zh=(50.0, 150.0, 250.0),
                zf=(100.0, 200.0, 300.0),
                pres=(0.0, 1000.0, 900.0, 800.0),
            )

    def test_profile_wrong_length(self) -> None:
        kwargs = _base_kwargs()
        kwargs["ta"] = kwargs["ta"][:-1]
        with self.assertRaises(ValueError):
            preplm(**kwargs)

    def test_pres_wrong_length(self) -> None:
        kwargs = _base_kwargs()
        kwargs["pres"] = kwargs["pres"][:-1]
        with self.assertRaises(ValueError):
            preplm(**kwargs)

    def test_fireflg_not_bool(self) -> None:
        kwargs = _base_kwargs()
        kwargs["fireflg"] = 1
        with self.assertRaises(ValueError):
            preplm(**kwargs)


class PreplmRegressionTests(unittest.TestCase):
    def test_stack_values_finite_and_in_range(self) -> None:
        kwargs = _base_kwargs()
        result = preplm(**kwargs)

        self.assertTrue(math.isfinite(result.stack_temperature))
        self.assertTrue(math.isfinite(result.stack_wind))

        ta = kwargs["ta"]
        tol = 1e-3
        self.assertGreaterEqual(result.stack_temperature, min(ta) - tol)
        self.assertLessEqual(result.stack_temperature, max(ta) + tol)

        wspd = result.wind_speed
        self.assertGreaterEqual(result.stack_wind, min(wspd) - tol)
        self.assertLessEqual(result.stack_wind, max(wspd) + tol)


if __name__ == "__main__":
    unittest.main()
