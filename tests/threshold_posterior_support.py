"""Small threshold posteriors by independent integration, not Gibbs updates.

These references deliberately construct scalar normal densities directly and
never use NWKIT's conditional-parameter or transition routines.
"""

import math

import numpy as np
from scipy.integrate import quad, quad_vec
from scipy.special import ndtr

_SQRT_2PI = math.sqrt(2 * math.pi)


def _normal(x, variance=1.0):
    return math.exp(-x * x / (2 * variance)) / (_SQRT_2PI * math.sqrt(variance))


def density_moments(density):
    """Normalized low-state probability, mean and variance for threshold zero."""
    mass = quad(density, -np.inf, np.inf, epsabs=1e-11)[0]
    mean = quad(lambda x: x * density(x), -np.inf, np.inf, epsabs=1e-11)[0] / mass
    second = quad(lambda x: x * x * density(x), -np.inf, np.inf, epsabs=1e-11)[0] / mass
    low = quad(density, -np.inf, 0, epsabs=1e-11)[0] / mass
    return {
        "probabilities": np.array([low, 1 - low]),
        "mean": mean,
        "variance": second - mean * mean,
    }


def binary_root_reference():
    """(A:1,B:2)R, A low and B high, unit root variance."""
    return density_moments(lambda r: _normal(r) * ndtr(-r) * ndtr(r / math.sqrt(2)))


def binary_internal_reference(a=1.0, b=0.5, c=2.0, stem=0.7):
    """((A:a,B:b)I:stem,C:c)R; low/high/high. Integrate root analytically.

    I has prior variance 1+stem. Conditional on I=u, root is normal with
    mean u/(1+stem) and variance stem/(1+stem); thus C can be integrated too.
    """
    variance = 1 + stem
    return density_moments(
        lambda u: (
            _normal(u, variance)
            * ndtr(-u / math.sqrt(a))
            * ndtr(u / math.sqrt(b))
            * ndtr((u / variance) / math.sqrt(c + stem / variance))
        )
    )


def ordinal_star_reference(*, ambiguous=False, tolerance=1e-8):
    """Three-state star, lengths 1,1,1 (plus ambiguous A|C tip of length 2).

    Root N(0,1), free threshold t>0 with flat density, observations low, middle,
    high. Return root probabilities/moments and posterior mean threshold.
    Category integration splits at t=r, avoiding gridded indicator error.
    """

    def integrand(t, r):
        zero = ndtr(-r)
        upper = ndtr(t - r)
        value = zero * (upper - zero) * ndtr(r - t)
        if ambiguous:
            value *= ndtr(-r / math.sqrt(2)) + ndtr((r - t) / math.sqrt(2))
        return np.array([value, t * value])

    def root_integrand(r):
        mass, threshold_moment = quad_vec(
            lambda t: integrand(t, r),
            0,
            np.inf,
            epsabs=tolerance / 100,
            epsrel=tolerance,
        )[0]
        high_mass = (
            0.0
            if r <= 0
            else quad(lambda t: integrand(t, r)[0], 0, r, epsabs=tolerance / 100)[0]
        )
        return _normal(r) * np.array(
            [
                mass,
                r * mass,
                r * r * mass,
                mass if r <= 0 else 0,
                high_mass,
                threshold_moment,
            ]
        )

    result = sum(
        quad_vec(root_integrand, lo, hi, epsabs=tolerance / 100, epsrel=tolerance)[0]
        for lo, hi in ((-np.inf, 0), (0, np.inf))
    )
    mass, first, second, low, high, threshold_moment = result
    mean = first / mass
    return {
        "probabilities": np.array([low, mass - low - high, high]) / mass,
        "mean": mean,
        "variance": second / mass - mean * mean,
        "threshold_mean": threshold_moment / mass,
    }


def fixed_ordinal_root_reference():
    """(A:1,B:1,C:1)R, observed low/middle/high, fixed thresholds 0,1."""

    def density(r):
        return _normal(r) * ndtr(-r) * (ndtr(1 - r) - ndtr(-r)) * ndtr(r - 1)

    result = density_moments(density)
    mass = quad(density, -np.inf, np.inf, epsabs=1e-11)[0]
    result["probabilities"] = np.array(
        [
            quad(density, lo, hi, epsabs=1e-11)[0] / mass
            for lo, hi in ((-np.inf, 0), (0, 1), (1, np.inf))
        ]
    )
    return result
