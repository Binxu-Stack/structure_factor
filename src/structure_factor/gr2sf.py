## from https://gist.githubusercontent.com/by256/b747e0bb9693c913249e83d30ace9dc2/raw/2f9ccf21b2b6e6a761d26e7824e20b2a4366244f/structure_factor.py
def compute_structure_factor(g_r, radii, rho, qmax=1.0, n=512):
    """
    Compute structure factor S(q) from a radial distribution function g(r).
    The calculation of S(q) can be further simplified by moving some of 
    the terms outside the sum. I have left it in its current form so that 
    it matches S(q) expressions found commonly in books and in literature.
    
    Parameters
    ----------
    
    g_r : np.ndarray
        Radial distribution function, g(r).
    radii : np.ndarray
        Independent variable of g(r).
    rho : float
        Average number density of particles.
    qmax : float
        Maximum value of momentum transfer (the independent variable of S(q)).
    n : int
        Number of points in S(q).
        
    Returns
    -------
    S_q : np.ndarray
        Structure factor
    Q : np.ndarray
        Momentum transfer (the independent variable of S(q)).
    
    """
    n_r = len(g_r)
    Q = np.linspace(0.0, qmax, n)
    S_q = np.zeros_like(Q)
    
    dr = radii[1] - radii[0]
    h_r = g_r - 1
    
    for q_idx, q in enumerate(Q):
        if q_idx == 0:
            S_q[q_idx] = np.sum([4*np.pi*(r**2)*h_r[r_idx]*dr for r_idx, r in enumerate(radii)])
        else:
            S_q[q_idx] = np.sum([4*np.pi*(r**2)*h_r[r_idx]*dr*np.sin(q*r)/(q*r) for r_idx, r in enumerate(radii)])
    
    S_q = 1 + rho * S_q / n_r
    
    return S_q, Q
