"""
Module: plots/common.py
Description: Checks shared by the generic plots before they draw: absorption needs a single origin for the
             observer distances (a wind farm has none), and the PSD and the 1/3 octave bands need a frequency
             bin width (frequencies that are not uniformly spaced, e.g. the tones of frequency-domain sources,
             have none).

Author: Raul Sanz Ramirez (raul.sanz.ramirez@upm.es / raul.sanz.ramirez@gmail.com)
Institution: Universidad Politecnica de Madrid - ETSIAE
Date: 10/2026
"""


def no_absorption(case      : dict = None,     # [-] Case
                  absorption: bool = True,     # [-] Absorption requested by the call
                  name      : str  = None):    # [-] Name of the calling plot, for the message
    """
    Absorption needs the distance travelled by the sound, measured from a single origin: without one (e.g. a
    wind farm, several sources) it is switched off with a message.

    Returns
    -------
    absorption : bool
    """

    if absorption and case["origin"] is None:
        print(f"{name}(): no single origin for the distances of {case['Name']}: absorption = False")
        return False

    return absorption

def psd_not_applicable(cases: list = None,     # [-] Cases drawn by the plot
                       mode : str  = None,     # [-] Mode of the plot
                       name : str  = None):    # [-] Name of the calling plot, for the message
    """
    The PSD divides by the frequency bin width: a case without one (frequencies not uniformly spaced) cannot
    be drawn in mode 'PSD'.

    Returns
    -------
    skip : bool
        True (with a message) if the plot is not applicable.
    """

    if str(mode).upper() != "PSD":
        return False

    tonal = [case for case in cases if case["bin_width"] is None]
    if not tonal:
        return False

    print(f"{name}(): mode 'PSD' needs uniformly spaced frequencies (a bin width) and those of {tonal[0]['Name']} "
          f"are not: not applicable, skipping")
    return True

def octave_not_applicable(cases : list = None,     # [-] Cases drawn by the plot
                          octave: bool = True,     # [-] 1/3 octave bands requested by the call
                          name  : str  = None):    # [-] Name of the calling plot, for the message
    """
    The 1/3 octave bands of frequencies without a bin width (not uniformly spaced) would be wrong: the bands
    holding none of them take the level of the nearest one (AcousticUtils.to_third_octave).

    Returns
    -------
    skip : bool
        True (with a message) if the plot is not applicable.
    """

    if not octave:
        return False

    tonal = [case for case in cases if case["bin_width"] is None]
    if not tonal:
        return False

    print(f"{name}(): 1/3 octave bands need uniformly spaced frequencies (a bin width) and those of "
          f"{tonal[0]['Name']} are not: not applicable, skipping")
    return True
