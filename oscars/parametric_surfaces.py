from math import sin, cos, sqrt, pi

class PSRectangle:
    """A Parametric surface - rectangle"""

    # This shape specific parameters
    L = 1
    W = 1
    

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = -L/2.
    ustop  = +L/2.
    vstart = -W/2.
    vstop  = +W/2.
    nu = 10
    nv = 10
 
    def __init__ (self, L=1, W=1, ustart=-L/2., ustop=L/2., nu=10, vstart=-W/2., vstop=W/2., nv=10):
        self.L = L
        self.W = W
        self.ustart = -L/2.
        self.ustop = L/2.
        self.nu = nu
        self.vstart = -W/2.
        self.vstop = W/2.
        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        x = u
        y = v
        z = 0
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        xn = 0
        yn = 0
        zn = 1
        
        return [xn, yn, zn]




class PSTorus:
    """A Parametric surface - torus"""

    # This shape specific parameters
    R = 1
    r = 1
    

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = 0
    ustop  = 2 * pi
    vstart = 0
    vstop  = 2 * pi
    nu = 10
    nv = 10
 
    def __init__ (self, R=1, r=1, ustart=0, ustop=2*pi, nu=10, vstart=0, vstop=2*pi, nv=10):
        self.R = R
        self.r = r
        self.ustart = ustart
        self.ustop = ustop
        self.nu = nu
        self.vstart = vstart
        self.vstop = vstop
        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        x = (self.R + self.r * cos(u)) * cos(v)
        y = (self.R + self.r * cos(u)) * sin(v)
        z = self.r * sin(u)
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        xn = -self.r * cos(u) * (self.R + self.r * cos(u)) * cos(v)
        yn = -self.r * cos(u) * (self.R + self.r * cos(u)) * sin(v)
        zn = -self.r * (self.R + self.r * cos(u)) * sin(u)
        
        mag = sqrt(xn*xn + yn*yn + zn*zn)
        
        return [xn / mag, yn / mag, zn / mag]




class PSCylinder:
    """A Parametric surface - cylinder with no top or bottom"""

    # This shape specific parameters
    R = 1
    L = 1

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = -L/2.
    ustop  = +L/2.
    vstart = 0
    vstop  = 2 * pi
    nu = 10
    nv = 10
 
    def __init__ (self, R=1, L=1, ustart=None, ustop=None, nu=10, vstart=None, vstop=None, nv=10):
        self.R = R
        self.L = L
        self.ustart = -L/2.
        self.ustop  = +L/2.
        self.nu = nu
        if vstart is None:
            self.vstart = 0
        else:
            self.vstart = vstart
        if vstop is None:
            self.vstop = 2.*pi
        else:
            self.vstop = vstop

        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        x = self.R * cos(v)
        y = self.R * sin(v)
        z = u
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        xn = -cos(v)
        yn = -sin(v)
        zn = 0
        
        mag = sqrt(xn*xn + yn*yn + zn*zn)
        
        return [xn / mag, yn / mag, zn / mag]


class PSSphere:
    """A Parametric surface - sphere"""

    # This shape specific parameters
    R = 1

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = 0.
    ustop  = 2. * pi
    vstart = -pi/2.
    vstop  = pi/2.
    nu = 10
    nv = 10
 
    def __init__ (self, R=1, ustart=None, ustop=None, nu=10, vstart=None, vstop=None, nv=10):
        self.R = R
        self.nu = nu
        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        x = self.R * cos(u) * cos(v)
        y = self.R * sin(u) * cos(v)
        z = self.R * sin(v)
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        xn = -self.R * self.R * cos(u) * cos(v) * cos(v)
        yn = -self.R * self.R * sin(u) * cos(v) * cos(v)
        zn = -self.R * self.R * cos(v) * sin(v)
        
        mag = sqrt(xn*xn + yn*yn + zn*zn)
        
        return [xn / mag, yn / mag, zn / mag]


class PSDisk:
    """A Parametric surface - Disk"""

    # This shape specific parameters
    r0 = 0
    r1 = 1
    

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = 0
    ustop  = 1
    vstart = 0
    vstop  = 2 * pi
    nu = 10
    nv = 10
 
    def __init__ (self, r0=1, r1=1, nu=10, vstart=0, vstop=2*pi, nv=10):
        self.r0 = r0
        self.r1 = r1
        self.ustart = r0
        self.ustop = r1
        self.nu = nu
        self.vstart = vstart
        self.vstop = vstop
        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        x = u * cos(v)
        y = u * sin(v)
        z = 0
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        return [0, 0, 1]




class PSCone:
    """A Parametric surface - cylinder with no top or bottom"""

    # This shape specific parameters
    r0 = 0
    r1 = 1
    L = 1

    # Required for all PS shapes
    # Start, stop, and number of points for the u and v parameters
    # All PSShapes must have these defined
    ustart = -L/2.
    ustop  = +L/2.
    vstart = 0
    vstop  = 2 * pi
    nu = 10
    nv = 10
 
    def __init__ (self, r0=0, r1=1, L=1, ustart=None, ustop=None, nu=10, vstart=None, vstop=None, nv=10):
        self.r0 = r0
        self.r1 = r1
        self.L = L
        self.ustart = -L/2.
        self.ustop  = +L/2.
        self.nu = nu
        if vstart is None:
            self.vstart = 0
        else:
            self.vstart = vstart
        if vstop is None:
            self.vstop = 2.*pi
        else:
            self.vstop = vstop

        self.nv = nv
       
    
    def position (self, u, v):
        """Return the position in 3D at this u and v"""

        R = self.r0 + (self.r1 - self.r0)*(u + self.L/2) / self.L
        x = R * cos(v)
        y = R * sin(v)
        z = u
        
        return [x, y, z]
    
    
    
    def normal (self, u, v):
        """Return a unit normal in 3D at this u and v position"""

        R = self.r0 + (self.r1 - self.r0)*(u + self.L/2) / self.L
        xn = -R * cos(v)
        yn = -R * sin(v)
        zn = (self.r1 - self.r0) / self.L * cos(v) * R * cos(v) + (self.r1 - self.r0) / self.L * sin(v) * R * sin(v)
        
        mag = sqrt(xn*xn + yn*yn + zn*zn)
        
        return [xn / mag, yn / mag, zn / mag]


