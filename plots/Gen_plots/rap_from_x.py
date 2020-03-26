import numpy as np

sqrt_s = 13000.
def rap_from_xf(xf, m, pt):
    mt = (m**2 + pt**2)**(0.5)
    rap = np.arcsinh(sqrt_s * xf / (2. * mt))
    return rap

def xf_from_rap(rap, m, pt):
    mt = (m**2 + pt**2)**(0.5)
    xf = 2. * mt / sqrt_s * np.sinh(rap)
    return xf




print(rap_from_xf(0.24, 500.,0.))
print(xf_from_rap(1.4, 170., 0.))
