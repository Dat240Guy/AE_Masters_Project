import numpy as np

import numpy as np

def kirsch_polar_stress(r, theta, a, sigma_inf):
    """
    Kirsch solution for an infinite plate with a circular hole of radius a,
    remote uniaxial tension sigma_inf applied in the x-direction.

    Parameters
    ----------
    r : array_like
        Radial coordinate (r >= a).
    theta : array_like
        Angle [rad], measured from +x axis (direction of applied tension).
    a : float
        Hole radius.
    sigma_inf : float
        Far-field tensile stress in x-direction.

    Returns
    -------
    sigma_rr, sigma_tt, tau_rt : ndarray
        Radial, hoop (circumferential), and shear stresses in polar coords.
    """
    r = np.asarray(r, dtype=float)
    theta = np.asarray(theta, dtype=float)

    if np.any(r < a):
        raise ValueError("Kirsch solution is defined only for r >= a (outside the hole).")

    a2_r2 = (a**2) / (r**2)
    a4_r4 = (a**4) / (r**4)

    cos2t = np.cos(2.0 * theta)
    sin2t = np.sin(2.0 * theta)

    sigma_rr = 0.5 * sigma_inf * (1 - a2_r2) \
               + 0.5 * sigma_inf * (1 + 3*a4_r4 - 4*a2_r2) * cos2t

    sigma_tt = 0.5 * sigma_inf * (1 + a2_r2) \
               - 0.5 * sigma_inf * (1 + 3*a4_r4) * cos2t

    tau_rt  = -0.5 * sigma_inf * (1 - 3*a4_r4 + 2*a2_r2) * sin2t

    return sigma_rr, sigma_tt, tau_rt



def kirsch_cartesian_stress(x, y, a, sigma_inf):
    """
    Kirsch solution: stresses around a circular hole in an infinite plate,
    returned in Cartesian components (sigma_x, sigma_y, tau_xy).
    
    Parameters
    ----------
    x, y : array_like
        Cartesian coordinates.
    a : float
        Hole radius.
    sigma_inf : float
        Remote uniform tensile stress in x-direction.
    
    Returns
    -------
    sigma_x, sigma_y, tau_xy : ndarray
        Cartesian stress components.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    # Get polar stresses from Kirsch solution
    sigma_rr, sigma_tt, tau_rt = kirsch_polar_stress(r, theta, a, sigma_inf)

    # Transform polar → Cartesian
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    cos2t = np.cos(2.0 * theta)
    sin2t = np.sin(2.0 * theta)

    # Standard plane stress transformation from (r,theta) to (x,y)
    # (assuming r direction aligned with theta as usual polar coordinates)
    sigma_x = (sigma_rr + sigma_tt) / 2.0 + (sigma_rr - sigma_tt) / 2.0 * cos2t - tau_rt * sin2t
    sigma_y = (sigma_rr + sigma_tt) / 2.0 - (sigma_rr - sigma_tt) / 2.0 * cos2t + tau_rt * sin2t
    tau_xy = (sigma_rr - sigma_tt) / 2.0 * sin2t + tau_rt * cos2t

    # Mask points inside the hole (r < a) as NaN (solution not defined there)
    inside = r < a
    if np.any(inside):
        sigma_x = np.where(inside, np.nan, sigma_x)
        sigma_y = np.where(inside, np.nan, sigma_y)
        tau_xy = np.where(inside, np.nan, tau_xy)

    return sigma_x, sigma_y, tau_xy


if __name__ == "__main__":
    a = 1.0          # hole radius
    sigma_inf = 1000.0  # remote tension

    # Sample a ring outside the hole
    r = 1.0 * a
    theta = np.linspace(0, 2*np.pi, 361)

    srr, stt, trt = kirsch_polar_stress(r, theta, a, sigma_inf)

    # Check hoop stress at the boundary (r = a)
    # Expect: sigma_theta(theta=pi/2) = 3*sigma_inf
    print("sigma_theta at theta=0 deg:", stt[theta == 0])
    print("sigma_theta at theta=90 deg:", stt[theta == np.pi/2])
