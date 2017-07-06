def setup_beamline(facility, beamline, osr=None, oth=None):
    """
    Setup a known beamline.  Load magnetic field data tables if known, setup
    beam parameters, load summary tables.

    Parameters
    ----------
    facility : str
        Name of the facility.  See below for currently supported.

    beamline : str
        Name of beamline.  See below for currently supported.

    osr : oscars.sr object
        oscars.sr object to load

    oth : oscars.th object
        oscars.th object to load

    Currently Supported - facility beamline
        NSLSII LIX
        NSLSII HXN

    Returns
    -------
    None
    """

    if facility == 'NSLSII':
        if beamline == 'LIX':
            pass
        elif beamline == 'HXN':
            pass
        else:
          raise  ValueError('beamline not recognized: ' + str(beamline))
    else:
      raise ValueError('facility not recognized: ' + str(facility))

    return
